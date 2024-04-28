package gestalt.evolution.substitutionmodel;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import cern.jet.random.AbstractDiscreteDistribution;
import cern.jet.random.NegativeBinomial;
import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;
import gestalt.evolution.alignment.*;
import org.apache.commons.math3.util.Pair;
import org.jblas.DoubleMatrix;

import java.util.*;

import static java.lang.Math.max;
import static org.jblas.DoubleMatrix.concatHorizontally;
import static org.jblas.MatrixFunctions.*;

@Description("Mutation model of GESTALT barcode evolution, as described in Feng et Al, 2021. Stores model parameters and processes them for calculation of the tree likelihood. ")
@Citation(value =
        "Feng J et Al. (2021) Estimation of cell lineage trees a\n" +
                "  by maximum-likelihood phylogenetics \n" +
                "  15.1:343-362.", DOI = "10.1214/20-AOAS1400", year = 2021, firstAuthorSurname = "feng")

public class gestaltGeneral extends SubstitutionModel.Base {


    //metadata inputs
    public Input<String> uneditedBarcodeInput = new Input<>("barcodeSequence",
            "The exact sequence of the unedited barcode, split in targets and spacers.");
    public Input<Integer> cutSiteInput = new Input<>("cutSite",
            "How far from the end of the target in bp the cut site is.");
    public Input<IntegerParameter> crucialPosInput = new Input<>("crucialPos",
            "Position range left and right of the cut that can deactivate the target.", new IntegerParameter("6 6"));
    public Input<Integer> maxSumStepsInput = new Input<>("maxSumSteps",
            "Max total number of steps from the parsimony states on each branch allowed to be taken");
    public Input<Integer> maxExtraStepsInput = new Input<>("maxExtraSteps",
            "Max steps past parsimony state taken");

    //target cut parameters input
    final public Input<List<RealParameter>> cutRatesInput = new Input<>("cutRates",
            "Rates at which each target is cut in the barcode",
            new ArrayList<>());
    final public Input<RealParameter> doubleCutWeightInput = new Input<>("doubleCutWeight",
            "Rate at which editing regions are lost (weigth)", new RealParameter("0.3"));

    final public Input<List<RealParameter>> longTrimFactorsInput = new Input<>("longTrimScalingFactors",
            "Scaling factor to get the rate at which long deletions left and right of the cut sites", new ArrayList<>());

    //indel parameters input
    //probability of zero deletions
    final public Input<List<RealParameter>> trimZeroProbsInput = new Input<>("trimZeroProbs",
            "Probablilties of a length zero deletion for indel types: focal-left, intertarget-left, focal right, intertarget-right", new ArrayList<>());
    final public Input<RealParameter> insertZeroProbInput = new Input<>("insertZeroProb",
            "Length zero insertion probability", new RealParameter("0.5"));
    //Log(means) for the distributions on the short deletions lengths left and right of the cut
    final public Input<List<RealParameter>> trimShortParamsInput = new Input<>("trimShortParams",
            "Log(Mean) x2 for the short deletion length distribution, left and right of the cut (poisson)", new ArrayList<>());
    //Log(means) for the distributions on the short deletions lengths left and right of the cut
    final public Input<List<RealParameter>> trimLongParamsInput = new Input<>("trimLongParams",
            "Log(Mean) x2 for the long deletion length distribution, left and right of the cut (poisson)", new ArrayList<>());
    final public Input<List<RealParameter>> insertParamsInput = new Input<>("insertParams",
            "Log(Mean) for the insert length distribution (poisson)", new ArrayList<>());

    //penalization parameters input
    public Input<Double> branchPenParamInput = new Input<>("branchPenParam",
            "Penalty parameter for branch lengths", 0.01);
    public Input<Double> targetLamPenParamInput = new Input<>("cutRatePenParam",
            "Penalty parameter for target cut rates", 0.01);

    //processed metadata
    static Integer numTargets;
    public BarcodeMeta metaData;

    //target cut type parameters
    public List<RealParameter> cutRates;
    public RealParameter doubleCutWeight;
    public List<RealParameter> longTrimFactors;

    //indel parameters
    public List<RealParameter> trimZeroProbs;
    public List<RealParameter> trimShortParams;
    public List<RealParameter> trimLongParams;
    public RealParameter insertZeroProb;
    public List<RealParameter> insertParams;

    //reshaped trimzero probabilities
    public DoubleMatrix trimZeroProbsDict = DoubleMatrix.zeros(2, 2);

    //choice of distribution for indel lengths
    boolean usePoisson = true;

    //penalty parameters
    public Double cutRatesPenalty = 0.0;
    public Double branchLensPenalty = 0.0;

    //indel distributions
    public Poisson insertDist;
    public List<AbstractDiscreteDistribution> delShortDist;
    public List<AbstractDiscreteDistribution> delLongDist;

    //parameterised event probabilities
    public Hashtable<TargetStatus, Double> hazardAwayDict;
    static Hashtable<IndelSet.Singleton, Integer> singletonIndexDict = new Hashtable<>();

    //conditional probabilities of the trims, this needs to be updated for each move involving parameter changes
    static DoubleMatrix singletonCondProb;

    //dictionary for storing hazards between target statuses -- assuming all moves are possible
    public static Hashtable<TargetStatus, Hashtable<TargetStatus, List<IndelSet.TargetTract>>> targStatTransitionsDict;

    //  hazards for all the target tracts need to be updated for each move involving parameter changes
    public DoubleMatrix targetTractHazards;
    // precalculcated hazards dictionary for all the target tracts to speed up overall calculation
    public Hashtable<IndelSet.TargetTract, Integer> targetTractDict;

    //this is to store hazards to potentially avoid recalculation
    public Hashtable<TargetStatus, Hashtable<TargetStatus, DoubleMatrix>> targStatTransitionHazardsDict;

    @Override
    public void initAndValidate() {

        //metadata
        String[] barcodeSplit = uneditedBarcodeInput.get().split(" ");
        int[] crucial = new int[2];
        crucial[0] = crucialPosInput.get().getValue(0);
        crucial[1] = crucialPosInput.get().getValue(1);
        metaData = new BarcodeMeta(Arrays.asList(barcodeSplit), cutSiteInput.get(), crucial, maxSumStepsInput.get(), maxExtraStepsInput.get());
        //toDo check why this needs to be in both metadata and here
        numTargets = metaData.nTargets;

        //target cut parameters
        cutRates = cutRatesInput.get();
        doubleCutWeight = doubleCutWeightInput.get();
        longTrimFactors = longTrimFactorsInput.get();

        //parameters controlling indel lengths
        insertZeroProb = insertZeroProbInput.get();
        trimZeroProbs = trimZeroProbsInput.get();
        trimShortParams = trimShortParamsInput.get();
        trimLongParams = trimLongParamsInput.get();
        insertParams = insertParamsInput.get();

        //penalization parameters
        branchLensPenalty = branchPenParamInput.get();
        cutRatesPenalty = targetLamPenParamInput.get();

        //simple checks:
        for (int i = 0; i < cutRates.size(); i++) {

            if (cutRates.get(i).getValue() < 0) {
                throw new RuntimeException("All cut rates must be positive!");
            }
        }
        if (doubleCutWeight.getValue() < 0) {
            throw new RuntimeException("Double cut weight must be positive!");
        }

        for (int i = 0; i < longTrimFactors.size(); i++) {
            if (longTrimFactors.get(i).getValue() < 0) {
                throw new RuntimeException("long trim factors must be positive!");
            }
        }
        for (int i = 0; i < trimZeroProbs.size(); i++) {
            if ((trimZeroProbs.get(i).getValue() >= 1) || (trimZeroProbs.get(i).getValue() <= 0)) {
                throw new RuntimeException("The probability of deleting no sequence must be >=0 and <=1!");
            }

        }
        if ((insertZeroProb.getValue() >= 1) || (insertZeroProb.getValue() <= 0)) {
            throw new RuntimeException("The probability of inserting no sequence must be >=0 and <=1!");
        }

        if (metaData.nTargets <= 0) {
            throw new RuntimeException("Empty barcodes are impossible, number of targets must be positive!");
        }

        //initialise reparameterised event probabilities

        //this only needs to happen once
        targStatTransitionsDict = TargetStatus.getAllTransitions(numTargets);
//        Pair<DoubleMatrix, Hashtable<IndelSet.TargetTract, Integer>> doubleMatrixHashtablePair = createAllTargetTractHazards();
//
//        targetTractHazards = doubleMatrixHashtablePair.getFirst();
//        targetTractDict = doubleMatrixHashtablePair.getSecond();
//
//        targStatTransitionHazardsDict = new Hashtable<>();
//
//        for (TargetStatus stat : targStatTransitionsDict.keySet()) {
//            Hashtable<TargetStatus, DoubleMatrix> empty = new Hashtable<>();
//            targStatTransitionHazardsDict.put(stat, empty);
//        }
        //createTrimInsertDistributions(10);
    }


    public static DoubleMatrix reorderLikelihoods(DoubleMatrix orderedDownProbs, TransitionWrap newWrapper, TransitionWrap oldWrapper) {
        List<Pair<Integer, Double>> indexVals = new ArrayList<>();
        for (TargetStatus newstat : newWrapper.transStatuses) {
            if (oldWrapper.statusMap.containsKey(newstat)) {
                indexVals.add(new Pair(newWrapper.statusMap.get(newstat), orderedDownProbs.get(oldWrapper.statusMap.get(newstat))));
            }
        }
        DoubleMatrix downProbs = new DoubleMatrix(newWrapper.numStatuses + 1, 1);
        for (int i = 0; i < indexVals.size(); i++) {
            downProbs.put(indexVals.get(i).getFirst(), indexVals.get(i).getSecond());
        }
        return downProbs;
    }

    /**
     * Calculates the rates at which all possible target tracts are introduced in the barcode, and maps each target tract rate to its index.
     */
    public Pair<DoubleMatrix, Hashtable<IndelSet.TargetTract, Integer>> createAllTargetTractHazards() {
        TargetStatus targetStatusAllActive = new TargetStatus();

        //this is done on a fully active target status, which means it is not tree dependent, but only dependent on mutation parameters
        List<IndelSet.TargetTract> allTargetTracts = targetStatusAllActive.getPossibleTargetTracts(new ArrayList<>(), metaData.nTargets);
        Hashtable<IndelSet.TargetTract, Integer> ttDict = new Hashtable<>();
        int size = allTargetTracts.size();
        for (int i = 0; i < size; i++) {
            ttDict.put(allTargetTracts.get(i), i);
        }
        int nTTs = allTargetTracts.size();
        DoubleMatrix minTargets = new DoubleMatrix(nTTs);
        DoubleMatrix maxTargets = new DoubleMatrix(nTTs);
        DoubleMatrix longLeftStatuses = new DoubleMatrix(nTTs);
        DoubleMatrix longRightStatuses = new DoubleMatrix(nTTs);
        size = allTargetTracts.size();
        for (int i = 0; i < size; i++) {
            IndelSet.TargetTract currentTT = allTargetTracts.get(i);
            minTargets.put(i, currentTT.getminTarg());
            maxTargets.put(i, currentTT.getmaxTarg());
            longLeftStatuses.put(i, currentTT.isLeftLong() ? 1.0 : 0.0);
            longRightStatuses.put(i, currentTT.isRightLong() ? 1.0 : 0.0);

        }
        DoubleMatrix allHazards = createHazardTargetTract(minTargets, maxTargets, longLeftStatuses, longRightStatuses);
        return new Pair(allHazards, ttDict);

    }


    /**
     * calculates the rates at which target tracts are introduced in the barcode given the targets cut, and the length of the deletion (short/long) at each side of the cut
     */
    public DoubleMatrix createHazardTargetTract(DoubleMatrix minTargets, DoubleMatrix maxTargets, DoubleMatrix longLeftStatuses, DoubleMatrix longRightStatuses) {


        DoubleMatrix logLeftTrimFactor = new DoubleMatrix(longLeftStatuses.length);
        DoubleMatrix logRightTrimFactor = new DoubleMatrix(longRightStatuses.length);
        int length = longLeftStatuses.length;

        //we operate in log space to apply (multiplicative) factors to the base cut rates

        //if there are deletions, add corresponding factors to the right and left
        for (int i = 0; i < length; i++) {
            if (longLeftStatuses.get(i) == 1.0) {
                logLeftTrimFactor.put(i, log(longTrimFactors.get(0).getValue(0)));
            } else {
                logLeftTrimFactor.put(i, 0.0);
            }
            if (longRightStatuses.get(i) == 1.0) {
                logRightTrimFactor.put(i, log(longTrimFactors.get(0).getValue(1)));
            } else {
                logRightTrimFactor.put(i, 0.0);
            }
        }


        DoubleMatrix gatheredTargetLamMin = new DoubleMatrix(minTargets.length);
        DoubleMatrix gatheredTargetLamMax = new DoubleMatrix(minTargets.length);
        length = minTargets.length;

        //collect the base cut rates for each target tract, where the cut occurs left and right
        for (int i = 0; i < length; i++) {
            gatheredTargetLamMin.put(i, cutRates.get(0).getValue((int) minTargets.get(i)));
            gatheredTargetLamMax.put(i, cutRates.get(0).getValue((int) maxTargets.get(i)));

        }


        DoubleMatrix gatherMinandMax = DoubleMatrix.zeros(minTargets.length);

        //combine the left and right cut rate: rate left cut + rate right cut
        gatherMinandMax.addi(gatheredTargetLamMax);
        gatherMinandMax.addi(gatheredTargetLamMin);

        //add the factor for the deletion types to the left and right
        DoubleMatrix logFocalLambdaPart = DoubleMatrix.zeros(minTargets.length);
        logFocalLambdaPart.addi(logLeftTrimFactor);
        logFocalLambdaPart.addi(logRightTrimFactor);
        logFocalLambdaPart.addi(logi(gatheredTargetLamMin));

        //if there is a double cut, add the cut weight
        DoubleMatrix logDoubleLambdaPart = DoubleMatrix.zeros(minTargets.length);
        logDoubleLambdaPart.addi(logLeftTrimFactor);
        logDoubleLambdaPart.addi(logRightTrimFactor);
        logDoubleLambdaPart.addi(logi(gatherMinandMax));
        logDoubleLambdaPart.addi(log(doubleCutWeight.getValue()));

        //hazard is (cut_rate_left + cut_rate_right)*double_cut_factor*long_cut_factors)
        DoubleMatrix hazard = new DoubleMatrix(minTargets.length);
        length = minTargets.length;
        for (int i = 0; i < length; i++) {
            if (minTargets.get(i) == maxTargets.get(i)) {
                hazard.put(i, logFocalLambdaPart.get(i));
            } else {
                hazard.put(i, logDoubleLambdaPart.get(i));
            }

        }

        //exit log scale
        expi(hazard);
        return hazard;

    }


    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {


    }

    public DoubleMatrix getTransitionProbabilities(Node node, TransitionWrap wrap, double branchLength, double rate) {


        //Create the probability matrix exp(Qt)
        final double branchTime = branchLength * rate;
        DoubleMatrix rateM = this.createRateMatrix(wrap);
        DoubleMatrix ptma = (expm(rateM.muli(branchTime)));
        return ptma;


    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Nucleotide;
    }


    /**
     * Creates the instantaneous transition matrix of the aggregated Markov chain
     */
    public DoubleMatrix createRateMatrix(TransitionWrap wrapper) {

        DoubleMatrix targTractLeft = createMarginalTransitionMatrixLeft(wrapper);
        DoubleMatrix trimProbsLeft = createTrimInstantProbMatrixLeft(wrapper);

        DoubleMatrix qMatrix = targTractLeft.muli(trimProbsLeft);
        DoubleMatrix lastColumn = qMatrix.rowSums();
        lastColumn.muli(-1);
        qMatrix = concatHorizontally(qMatrix, lastColumn);

        return qMatrix;

    }

    /**
     * Creates the matrix for rates for the introduction of target tracts into barcodes with specific satuses.
     */
    public DoubleMatrix createMarginalTransitionMatrixLeft(TransitionWrap wrapper) {
        //first component in the factorization
        //                of the instantaneous transition rates, the part that deals with the target status
        //                and target tracts. The last row corresponds to the sink state.
        //                We omit the last column since its parent function will handle its creation
        Set specialTts = new HashSet();
        for (IndelSet sgwc : wrapper.ancStates.getSingletonWCs()) {
            specialTts.add(sgwc.getTargetTract());
        }

        Set possibleStates = new HashSet();
        possibleStates.addAll(wrapper.transStatuses);

        List<Pair> singleTtSparseIndices = new ArrayList<>();
        List<Integer> singleTtGatherIndices = new ArrayList<>();
        List<Pair<Integer, Integer>> sparseIndices = new ArrayList<>();
        List<Double> sparseValues = new ArrayList<>();


        for (TargetStatus startState : wrapper.transStatuses) {
            Integer startKey = wrapper.statusMap.get(startState);

            Double hazardAway = hazardAwayDict.get(startState);
            sparseIndices.add(new Pair(startKey, startKey));
            sparseValues.add(hazardAway * (-1));
            Set<TargetStatus> allEndStates = targStatTransitionsDict.get(startState).keySet();
            //should be the intersection
            List<TargetStatus> possibleEndStates = new ArrayList(allEndStates);
            possibleEndStates.retainAll(possibleStates);
            for (TargetStatus endState : possibleEndStates) {
                Integer endKey = wrapper.statusMap.get(endState);
                List<IndelSet.TargetTract> targetTractsTransition = targStatTransitionsDict.get(startState).get(endState);

                List<IndelSet.TargetTract> matchingTTs = IndelSet.TargetTract.intersectList(new ArrayList<>(specialTts), targetTractsTransition);

                if (matchingTTs != null && matchingTTs.size() != 0) {
                    assert matchingTTs.size() == 1;
                    IndelSet.TargetTract matchingTT = matchingTTs.get(0);

                    singleTtGatherIndices.add(targetTractDict.get(matchingTT));
                    singleTtSparseIndices.add(new Pair(startKey, endKey));
                } else {
                    Double hazard = 0.0;
                    if (targStatTransitionHazardsDict.get(startState).containsKey(endState)) {
                        hazard = targStatTransitionHazardsDict.get(startState).get(endState).get(0);
                    } else {
                        List<Integer> hazard_idxs = new ArrayList<>();
                        for (IndelSet.TargetTract tt : targetTractsTransition) {
                            hazard_idxs.add(targetTractDict.get(tt));
                        }
                        for (int i : hazard_idxs) {
                            hazard = hazard + targetTractHazards.get(i);
                        }
                        Hashtable<TargetStatus, DoubleMatrix> toAppend = targStatTransitionHazardsDict.get(startState);
                        DoubleMatrix entry = new DoubleMatrix(1);
                        entry.put(0, hazard);
                        toAppend.put(endState, entry);
                        targStatTransitionHazardsDict.put(startState, toAppend);

                    }
                    sparseIndices.add(new Pair(startKey, endKey));
                    sparseValues.add(hazard);


                }
            }
        }

        Integer matrixLength = wrapper.numStatuses + 1;

        DoubleMatrix qSingleTtMatrix = DoubleMatrix.zeros(matrixLength, matrixLength - 1);
        if (singleTtGatherIndices.size() != 0) {

            DoubleMatrix singleTtSparseVals = new DoubleMatrix(singleTtGatherIndices.size());
            for (int i = 0; i < singleTtGatherIndices.size(); i++) {
                singleTtSparseVals.put(i, targetTractHazards.get(singleTtGatherIndices.get(i)));
            }
            for (int i = 0; i < singleTtSparseIndices.size(); i++) {
                Pair<Integer, Integer> indices = singleTtSparseIndices.get(i);
                qSingleTtMatrix.put(indices.getFirst(), indices.getSecond(), singleTtSparseVals.get(i));
            }

        }


        DoubleMatrix qAllTtMatrix = DoubleMatrix.zeros(matrixLength, matrixLength - 1);

        for (int i = 0; i < sparseIndices.size(); i++) {

            Pair<Integer, Integer> indices = sparseIndices.get(i);
            qAllTtMatrix.put(indices.getFirst(), indices.getSecond(), sparseValues.get(i));
        }
        DoubleMatrix qmatrixLeft = qAllTtMatrix;

        qmatrixLeft.addi(qSingleTtMatrix);
        return qmatrixLeft;

    }

    /**
     * Creates the matrix of conditional probabilities of each trim.
     * The entry in (i,j) is the trim probability for transitioning from target status i
     * to target status j.
     */

    public static DoubleMatrix createTrimInstantProbMatrixLeft(TransitionWrap childWrap) {
        //The second component when factorizing the instantaneous transition rates
        //                between the meta-states. This is the matrix of conditional probabilities of each trim.
        //                So the entry in (i,j) is the trim probability for transitioning from target status i
        //                to target status j. There is a trim prob associated with it because this must correspond to
        //                the introduction of a singleton in the AncState.
        //                Last row corresponds to sink state.
        //                We omit the last column corresponding to sink state. Its parent function will handle
        //                the creation.
        //
        //                note: this is used to generate the matrix for a specific branch in the tree
        //                If the transition is not possible, we fill in with trim prob 1.0 since it doesnt matter

        List<IndelSet.Singleton> childSingletons = childWrap.ancStates.getSingletons();

        List<IndelSet.Singleton> targetToSingleton = new ArrayList<>();
        for (int i = 0; i < numTargets; i++) {
            targetToSingleton.add(null);
        }
        for (IndelSet.Singleton sg : childSingletons) {
            targetToSingleton.set(sg.getminTarg(), sg);
            targetToSingleton.set(sg.getmaxTarg(), sg);
        }

        Set<TargetStatus> possibleStates = new HashSet<>();
        possibleStates.addAll(childWrap.transStatuses);
        List<Pair<Integer, Integer>> sparseIndices = new ArrayList<>();
        List<Double> sparseVals = new ArrayList<>();
        for (TargetStatus startTargetStatus : possibleStates) {
            Set<TargetStatus> allEndStates = new HashSet<>();
            allEndStates.addAll(targStatTransitionsDict.get(startTargetStatus).keySet());
            Set<TargetStatus> possibleEndState = new HashSet<>(allEndStates);
            possibleEndState.retainAll(possibleStates);

            for (TargetStatus endTargetStatus : possibleEndState) {
                Set<Integer> newDeacTargs = endTargetStatus.minus(startTargetStatus);
                Set<IndelSet.Singleton> singletonSet = new HashSet<>();

                for (Integer deactTarg : newDeacTargs) {
                    if (targetToSingleton.size() != 0 && targetToSingleton.get(deactTarg) != null) {
                        singletonSet.add(targetToSingleton.get(deactTarg));
                    }
                }
                Integer startKey = childWrap.statusMap.get(startTargetStatus);
                Integer endKey = childWrap.statusMap.get(endTargetStatus);
                if (singletonSet.size() == 1) {
                    //CHECK THAT
                    sparseIndices.add(new Pair(startKey, endKey));
                    IndelSet.Singleton sg = (IndelSet.Singleton) singletonSet.toArray()[0];
                    Double trimProbVal = singletonCondProb.get(singletonIndexDict.get(sg));
                    sparseVals.add(trimProbVal);


                } else if (singletonSet.size() > 1) {
                    sparseIndices.add(new Pair(startKey, endKey));
                    sparseVals.add(0.0);
                }

            }


        }
        int outputLength = childWrap.numStatuses + 1;
        int[] outputShape = new int[]{outputLength, outputLength - 1};
        if (sparseVals.size() != 0) {
            DoubleMatrix qMatrixLeft = DoubleMatrix.ones(outputShape[0], outputShape[1]);
            for (int i = 0; i < sparseIndices.size(); i++) {
                qMatrixLeft.put(sparseIndices.get(i).getFirst(), sparseIndices.get(i).getSecond(), sparseVals.get(i));

            }
            return qMatrixLeft;
        } else {
            return DoubleMatrix.ones(outputShape[0], outputShape[1]);
        }

    }

    /**
     * returns a matrix of hazards of transitioning away from specific target statuses
     */
    public DoubleMatrix createHazardAwayTargetStatuses(List<TargetStatus> targetStatuses) {

        //For each target status create a corresponding vector containing the associated cut rate or 0 if it is already
        //inactive
        DoubleMatrix activeTargetHazards = DoubleMatrix.zeros(targetStatuses.size(), numTargets);
        DoubleMatrix activeMasks = DoubleMatrix.zeros(targetStatuses.size(), numTargets);
        for (int i = 0; i < targetStatuses.size(); ++i) {
            Integer[] binaryrep = targetStatuses.get(i).getBinaryStatus(numTargets);
            for (int j = 0; j < numTargets; j++) {
                activeMasks.put(i, j, 1 - binaryrep[j]);
                activeTargetHazards.put(i, j, activeMasks.get(i, j) * cutRates.get(0).getValue(j));
            }
        }

        if (numTargets == 1) {
            return activeTargetHazards.getColumn(0);
        }


        int[] indicesInner = new int[activeTargetHazards.getColumns() - 2];
        for (int i = 1; i < activeTargetHazards.getColumns() - 1; ++i) {
            indicesInner[i - 1] = i;
        }

        DoubleMatrix focalHazards = ((activeTargetHazards.getColumn(0)).muli(1 + longTrimFactors.get(0).getValue(1))).addi(((activeTargetHazards.getColumns(indicesInner)).rowSums()).muli((1 + longTrimFactors.get(0).getValue(1)) * (1 + longTrimFactors.get(0).getValue(0)))).addi((activeTargetHazards.getColumn(activeTargetHazards.getColumns() - 1)).muli(1 + longTrimFactors.get(0).getValue(0)));
        DoubleMatrix middleHazards = (activeTargetHazards.getColumns(indicesInner).rowSums());
        DoubleMatrix numMiddle = (activeMasks.getColumns(indicesInner).rowSums());
        DoubleMatrix middleDoubleCutHazards = middleHazards.muli((numMiddle.subi(1)).muli((1 + longTrimFactors.get(0).getValue(0)) * (1 + longTrimFactors.get(0).getValue(1))));
        DoubleMatrix numAfterStarts = (activeMasks.getColumns(indicesInner).rowSums());
        DoubleMatrix startToOtherHazard = activeMasks.getColumn(0).muli(1 + longTrimFactors.get(0).getValue(1)).muli((numAfterStarts.muli(activeTargetHazards.getColumn(0))).addi((activeTargetHazards.getColumns(indicesInner).rowSums())));
        DoubleMatrix numBeforeEnd = (activeMasks.getColumns(indicesInner).rowSums());
        DoubleMatrix endToOtherHazard = (activeMasks.getColumn(activeMasks.getColumns() - 1)).muli(1 + longTrimFactors.get(0).getValue(0)).muli((numBeforeEnd).muli((activeTargetHazards.getColumn(activeTargetHazards.getColumns() - 1))).addi((activeTargetHazards.getColumns(indicesInner).rowSums())));
        DoubleMatrix startToEndHazards = (activeMasks.getColumn(0)).muli(activeMasks.getColumn(activeMasks.getColumns() - 1)).muli((activeTargetHazards.getColumn(0).addi(activeTargetHazards.getColumn(activeTargetHazards.getColumns() - 1))));


        DoubleMatrix hazardAwayNodes = focalHazards.addi((middleDoubleCutHazards.addi(startToOtherHazard).addi(endToOtherHazard).addi(startToEndHazards)).muli(doubleCutWeight.getValue()));

        return hazardAwayNodes;


    }

    /**
     * get all possible singletons from the computed set of ancestral states  for a specific tree and leaf sequences
     */

    public static List<IndelSet.Singleton> getAllSingletons(beast.base.evolution.tree.TreeInterface tree, Hashtable<Integer, AncStates> statesDict) {
        Set<IndelSet.Singleton> singletonSet = new HashSet();
        for (Node node : tree.getExternalNodes()) {
            if (node.isLeaf()) {
                AncStates leaf_state = statesDict.get(node.getNr());
                for (IndelSet sgwc : leaf_state.getSingletonWCs()) {
                    IndelSet.Singleton sing = sgwc.getSingleton();
                    singletonSet.add(sing);
                }
            }
        }

        List<IndelSet.Singleton> singletonList = new ArrayList<>();
        singletonList.addAll(singletonSet);
        return singletonList;


    }


    /**
     * get all possible conditional probabilities of the trims (all at once to speed the computation)
     */

    public void initSingletonProbs(List<IndelSet.Singleton> singletons) {


        for (int i = 0; i < singletons.size(); i++) {
            singletonIndexDict.put(singletons.get(i), i);
        }
        createTrimInsertDistributions();
        singletonCondProb = expi(createLogIndelProbs(singletons));

    }

    /**
     * Creates the basic trim + insert helper distributions
     */

    public void createTrimInsertDistributions() {

        DoubleMatrix trimLongParamReshaped = new DoubleMatrix(2, 1);
        //left log(mean)
        trimLongParamReshaped.put(0, trimLongParams.get(0).getValue(0));
        //right log(mean)
        trimLongParamReshaped.put(1, trimLongParams.get(0).getValue(1));

        DoubleMatrix trimShortParamsReshaped = new DoubleMatrix(2, 1);
        //left log(mean)
        trimShortParamsReshaped.put(0, trimShortParams.get(0).getValue(0));
        //right log(mean)
        trimShortParamsReshaped.put(1, trimShortParams.get(0).getValue(1));


        delShortDist = makeDelDist(trimShortParamsReshaped, 2, usePoisson);
        delLongDist = makeDelDist(trimLongParamReshaped, 2, usePoisson);
        if (usePoisson) {
            insertDist = new Poisson(exp(insertParams.get(0).getValue(0)), RandomEngine.makeDefault());
        }

        //in original ML implementation, the option of having Negative Binomial distributions for insert lengths was available
       /*else {

            insertDist = new NegativeBinomial(insertParams[1].intValue(),exp(insertParams[0]),RandomEngine.makeDefault());
        }*/

    }


    /**
     * Creates a basic trim helper distribution
     */

    public List<AbstractDiscreteDistribution> makeDelDist(DoubleMatrix params, int nTrimTypes, boolean usePoisson) {
        List<AbstractDiscreteDistribution> delDistList = new ArrayList<>();
        if (usePoisson) {
            for (int i = 0; i < nTrimTypes; i++) {
                delDistList.add(new Poisson(exp(params.get(i)), RandomEngine.makeDefault()));
            }
        } else {
            for (int i = 0; i < nTrimTypes; i++) {
                for (int j = 0; j < numTargets; j++) {
                    delDistList.add(new NegativeBinomial((int) params.get(i, 1), exp(params.get(i, 0)), RandomEngine.makeDefault()));
                }
            }

        }
        return delDistList;

    }



    /**
     * Create the conditional probability of indels
     */
    public DoubleMatrix createLogIndelProbs(List<IndelSet.Singleton> singletons) {

        trimZeroProbsDict.put(0, 0, trimZeroProbs.get(0).getValue(0));
        trimZeroProbsDict.put(0, 1, trimZeroProbs.get(0).getValue(1));
        trimZeroProbsDict.put(1, 0, trimZeroProbs.get(0).getValue(2));
        trimZeroProbsDict.put(1, 1, trimZeroProbs.get(0).getValue(3));

        if (singletons.size() == 0) {
            return new DoubleMatrix();
        } else {
            // assemble individual probabilities:
            DoubleMatrix leftDelProb = createLeftDelProbs(singletons);
            DoubleMatrix rightDelProb = createRightDelProbs(singletons);
            DoubleMatrix insertProb = createInsertProbs(singletons);

            //combine everything
            DoubleMatrix allLogProbs = (logi(leftDelProb)).addi(logi(rightDelProb)).addi(logi(insertProb));

            //initialise normalization
            Double logShortFocalNormalization = log(1 - (trimZeroProbsDict.get(0, 0) * trimZeroProbsDict.get(1, 0) * insertZeroProb.getValue()));

            boolean[] isLongIndel = new boolean[singletons.size()];
            for (int i = 0; i < singletons.size(); i++) {
                isLongIndel[i] = (singletons.get(i).isLeftLong() || singletons.get(i).isRightLong() || singletons.get(i).isIntertarget());
            }
            DoubleMatrix fullLogIndelProbs = new DoubleMatrix(singletons.size());
            for (int i = 0; i < singletons.size(); i++) {
                if (isLongIndel[i]) {
                    fullLogIndelProbs.put(i, allLogProbs.get(i));
                } else {
                    fullLogIndelProbs.put(i, allLogProbs.get(i) - logShortFocalNormalization);
                }

            }
            return fullLogIndelProbs;

        }


    }

    /**
     * Create insert probabilities
     */
    public DoubleMatrix createInsertProbs(List<IndelSet.Singleton> singletons) {

        List<Integer> insertLengths = new ArrayList<>();
        for (IndelSet.Singleton sg : singletons) {
            insertLengths.add(sg.getInsertLength());
        }
        /*DoubleMatrix insertSeqProb = new DoubleMatrix();
        for (int length:insertLengths) {
            insertSeqProb.add(pow(4,length));
        }*/
        DoubleMatrix insertLenProb = new DoubleMatrix(insertLengths.size());
        for (int i = 0; i < insertLengths.size(); i++) {
            insertLenProb.put(i, insertDist.pdf(max((insertLengths.get(i) - 1), 0)));
        }
        DoubleMatrix Probs = new DoubleMatrix(insertLengths.size());
        for (int i = 0; i < insertLengths.size(); i++) {
            if (insertLengths.get(i) == 0) {
                Probs.put(i, insertZeroProb.getValue());
            } else {
                Probs.put(i, (1 - insertZeroProb.getValue()) * insertLenProb.get(i));
            }
        }

        return Probs;
    }


    /**
     * Create deletion probabilities for the left side of the cut
     */
    public DoubleMatrix createLeftDelProbs(List<IndelSet.Singleton> singletons) {
        List<Integer> minTargets = new ArrayList<>();
        DoubleMatrix isLeftLongs = new DoubleMatrix(singletons.size());
        DoubleMatrix isIntertargets = new DoubleMatrix(singletons.size());
        DoubleMatrix startpos = new DoubleMatrix(singletons.size());
        int size = singletons.size();
        for (int i = 0; i < size; i++) {
            IndelSet.Singleton sg = singletons.get(i);
            minTargets.add(sg.getminTarg());
            isLeftLongs.put(i, sg.isLeftLong() ? 1.0 : 0.0);
            isIntertargets.put(i, sg.isIntertarget() ? 1.0 : 0.0);
            startpos.put(i, sg.getStartPos());
        }

        DoubleMatrix minTargetSites = new DoubleMatrix(minTargets.size());
        DoubleMatrix leftTrimLongMin = new DoubleMatrix(minTargets.size());
        DoubleMatrix leftTrimLongMax = new DoubleMatrix(minTargets.size());
        size = minTargets.size();
        for (int i = 0; i < size; i++) {
            minTargetSites.put(i, metaData.absCutSites.get(minTargets.get(i)));
            leftTrimLongMin.put(i, metaData.leftLongTrimMin.get(minTargets.get(i)));
            leftTrimLongMax.put(i, metaData.leftMaxTrim.get(minTargets.get(i)));
        }
        DoubleMatrix leftTrimLens = minTargetSites.subi(startpos);

        return createDelProbs(leftTrimLens, isLeftLongs, isIntertargets, leftTrimLongMin, leftTrimLongMax, false);


    }


    /**
     * Create deletion probabilities for the right side of the cut
     */
    public DoubleMatrix createRightDelProbs(List<IndelSet.Singleton> singletons) {
        List<Integer> maxTargets = new ArrayList<>();
        DoubleMatrix isRightLongs = new DoubleMatrix(singletons.size());
        DoubleMatrix isIntertargets = new DoubleMatrix(singletons.size());
        DoubleMatrix endpos = new DoubleMatrix(singletons.size());
        int size = singletons.size();
        for (int i = 0; i < size; i++) {
            IndelSet.Singleton sg = singletons.get(i);
            maxTargets.add(sg.getmaxTarg());
            isRightLongs.put(i, sg.isRightLong() ? 1.0 : 0.0);
            isIntertargets.put(i, sg.isIntertarget() ? 1.0 : 0.0);
            endpos.put(i, sg.getEndPos());
        }
        DoubleMatrix maxTargetSites = new DoubleMatrix(singletons.size());
        DoubleMatrix rightTrimLongMin = new DoubleMatrix(singletons.size());
        DoubleMatrix rightTrimLongMax = new DoubleMatrix(singletons.size());
        size = maxTargets.size();
        for (int i = 0; i < size; i++) {
            maxTargetSites.put(i, metaData.absCutSites.get(maxTargets.get(i)));
            rightTrimLongMin.put(i, metaData.rightLongTrimMin.get(maxTargets.get(i)));
            rightTrimLongMax.put(i, metaData.rightMaxTrim.get(maxTargets.get(i)));
        }
        //todo check that the lengths add up correctly
        DoubleMatrix rightTrimLen = endpos.subi(maxTargetSites);
        return createDelProbs(rightTrimLen, isRightLongs, isIntertargets, rightTrimLongMin, rightTrimLongMax, true);
    }


    /**
     * calculate the log conditional probability of the deletions found in
     * each of the singletons
     */
    public DoubleMatrix createDelProbs(DoubleMatrix trimLen, DoubleMatrix isLongs, DoubleMatrix isIntertargets, DoubleMatrix trimMins, DoubleMatrix trimMaxs, boolean isRight) {
        int esRight = (isRight ? 1 : 0);

        DoubleMatrix trimZeroProbs = new DoubleMatrix(isIntertargets.length);
        int length = isIntertargets.length;
        for (int i = 0; i < length; i++) {
            if (isIntertargets.get(i) == 1) {
                trimZeroProbs.put(i, trimZeroProbsDict.get(esRight, 1));
            } else {
                trimZeroProbs.put(i, trimZeroProbsDict.get(esRight, 0));
            }
        }

        Poisson delShortD = (Poisson) delShortDist.get(esRight);
        DoubleMatrix shortNonZeroProb = new DoubleMatrix(trimZeroProbs.length);
        length = trimZeroProbs.length;
        for (int i = 0; i < length; i++) {
            int pdfInput = (int) max(trimLen.get(i) - 1, 0);
            int cdfInput = (int) trimMaxs.get(i) - 1;
            shortNonZeroProb.put(i, (1 - trimZeroProbs.get(i)) * delShortD.pdf(pdfInput) / (delShortD.cdf(cdfInput)));
        }

        Poisson delLongD = (Poisson) delLongDist.get(esRight);
        DoubleMatrix LongProb = new DoubleMatrix(trimZeroProbs.length);
        length = trimZeroProbs.length;
        for (int i = 0; i < length; i++) {
            int pdfInput = (int) max(trimLen.get(i) - trimMins.get(i), 0);
            int cdfInput = (int) max(trimMaxs.get(i) - trimMins.get(i), 0);
            LongProb.put(i, delLongD.pdf(pdfInput) / (delLongD.cdf(cdfInput)));
        }

        DoubleMatrix Probs = new DoubleMatrix(trimLen.length);
        length = trimLen.length;
        for (int i = 0; i < length; i++) {
            if (trimLen.get(i) == 0) {
                Probs.put(i, trimZeroProbs.get(i));
            } else {
                if (isLongs.get(i) == 1) {
                    Probs.put(i, LongProb.get(i));
                } else {
                    Probs.put(i, shortNonZeroProb.get(i));
                }
            }
        }

        return Probs;

    }

    /**
     * Returns:
     * Dictionary mapping all possible TargetStatus to matrix for the hazard away,
     */
    public Hashtable<TargetStatus, Double> createHazardAwayDict() {
        Set<TargetStatus> targetStatuses = targStatTransitionsDict.keySet();

        List<TargetStatus> targetStatusesList = new ArrayList<>();
        targetStatusesList.addAll(targetStatuses);
        DoubleMatrix hazardAwayNodes = createHazardAwayTargetStatuses(targetStatusesList);

        Hashtable<TargetStatus, Double> hazardAwayDict = new Hashtable<>();
        int size = targetStatusesList.size();
        for (int i = 0; i < size; i++) {
            hazardAwayDict.put(targetStatusesList.get(i), hazardAwayNodes.get(i));
        }
        return hazardAwayDict;

    }

}


