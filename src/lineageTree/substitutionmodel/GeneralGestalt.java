package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import cern.jet.random.AbstractDiscreteDistribution;
import cern.jet.random.NegativeBinomial;
import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;
import org.jblas.DoubleMatrix;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.stream.Stream;
import java.util.*;
import org.apache.commons.math3.util.Pair;
import org.jblas.MatrixFunctions;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static org.jblas.DoubleMatrix.*;
import static org.jblas.MatrixFunctions.*;


public class GeneralGestalt extends SubstitutionModel.Base {

        final public Input<List<RealParameter>> cutRatesInput = new Input<>("cutRates",
                "Rates at which each target is cut in the barcode",
                new ArrayList<>());
        final public Input<List<RealParameter>> longTrimScalingInput = new Input<>("longTrimScaling",
                "Scaling factor for long deletions left and right of the cut sites", new ArrayList<>());
        final public Input<RealParameter> doubleCutWeightInput = new Input<>("doubleCutWeight",
            "Rate, at which editing regions are lost", Input.Validate.REQUIRED);

       public Input<Double> branchPenParamInput = new Input<>("branchPenParam",
               "Penalty parameter for branch lengths",0.01);
       public Input<Double> targetLamPenParamInput = new Input<>("cutRatePenParam",
            "Penalty parameter for target cut rates",0.01);

        public Input<String> uneditedBarcodeInput = new Input<>("barcodeSequence",
                "The exact sequence of the unedited barcode, split in targets and spacers.");
        public Input<Integer> cutSiteInput = new Input<>("cutSite",
                "How far from the end of the target in bp the cut site is.");
        public Input<IntegerParameter> crucialPosInput = new Input<>("crucialPos",
                "Position range left and right of the cut that can deactivate the tarfet.");
        public Input<Integer> maxSumStepsInput = new Input<>("maxSumSteps",
            "Position range left and right of the cut that can deactivate the tarfet.");
        public Input<Integer> maxExtraStepsInput = new Input<>("maxExtraSteps",
            "Position range left and right of the cut that can deactivate the tarfet.");


        /**
         * flag to indicate matrix is up to date *
         */
        protected boolean updateMatrix = false;

      /*  double editingHeight;
        double editingDuration;*/
        double[] frequencies;
        double[][] rateMatrix;

        //Model parameters
        Integer numTargets;
        public Double[] cutRates = {1.09459452, 1.03371947, 1.02624685};

        public Double[] longTrimFactors = {0.04,0.04};

        public Double[] trimZeroProbs = {0.5,0.5,0.5,0.5};
        public Double[] trim_short_params = {3.0,3.0};
        public Double[] trim_long_params = {3.0,3.0};
        Double insertZeroProb = 0.5;
         Double[] insertParams = {2.0};
        //toy example clt_calc: Double[] insertParams = {2.0};
        //tune_topology: Double[] insertParams = {1.0};
        public Double doubleCutWeight = 0.03;
        //toy example clt_calc: public Double doubleCutWeight = 0.3;
        //tune_topology: public Double doubleCutWeight = 0.1;
        public Double cutRatesPenalty = 0.0;
        public Double branchLensPenalty = 0.0;
        boolean usePoisson= true;

        //Computed/rearranged model parameters
        DoubleMatrix trim_short_params_reshaped;
        DoubleMatrix trim_long_params_reshaped;
        public Hashtable<TargetStatus, Double> hazardAwayDict;
        Hashtable<IndelSet.Singleton,Integer> singletonIndexDict = new Hashtable<>();
        DoubleMatrix singletonLogCondProb;
        DoubleMatrix singletonCondProb;
        DoubleMatrix trimZeroProbsDict = DoubleMatrix.zeros(2,2);
        Poisson insertDist;
        List<AbstractDiscreteDistribution> delShortDist;
        List<AbstractDiscreteDistribution> delLongDist;
        Hashtable<TargetStatus, Hashtable<TargetStatus, List<IndelSet.TargetTract>>> targStatTransitionsDict;
        public DoubleMatrix targetTractHazards;
        public Hashtable<IndelSet.TargetTract,Integer> targetTractDict;
        Hashtable<TargetStatus,Hashtable<TargetStatus,DoubleMatrix>> targ_stat_transition_hazards_dict;



        public BarcodeMeta metaData;


    @Override
    public void initAndValidate() {
        //TO DO: remove loss and scarRates to test the rest!
        String[] barcodeSplit = uneditedBarcodeInput.get().split(" ");
        int[] crucial = new int[2];
        crucial[0] = crucialPosInput.get().getValue(0);
        crucial[1] = crucialPosInput.get().getValue(1);
        metaData = new BarcodeMeta(Arrays.asList(barcodeSplit),cutSiteInput.get(),crucial,maxSumStepsInput.get(),maxExtraStepsInput.get());

        numTargets = metaData.nTargets;
        //nr of states changes
        nrOfStates = cutRatesInput.get().get(0).getDimension() + 2;
        rateMatrix = new double[nrOfStates][nrOfStates];
        cutRates = cutRatesInput.get().get(0).getValues();
        //longTrimFactors = longTrimScalingInput.get().get(0).getValues();
        doubleCutWeight = doubleCutWeightInput.get().getValue();
        branchLensPenalty = branchPenParamInput.get();
        cutRatesPenalty =  targetLamPenParamInput.get();


        // assert positive rates
        double sumScarringRates = 0;

        // check cut rates
        for (int i=0; i<cutRates.length; i++){

            if (cutRates[i] <= 0) {
                throw new RuntimeException("All cut rates must be positive!");
            }

        }

        // check cut rates
        for (int i=0; i<longTrimFactors.length; i++){

            if (longTrimFactors[i] > 1 || longTrimFactors[i] < 0) {
                throw new RuntimeException("long trim factors are assumed to be less than 1");
            }

        }

        if (doubleCutWeight <= 0) {
            throw new RuntimeException("Double cut weight must be positive!");
        }

        if (metaData.nTargets <= 0) {
            throw new RuntimeException("Empty barcodes are impossible, number of targets must be positive!");
        }




        // center root frequency on unedited state
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;

      /*  editingHeight = editingHeightInput.get();
        editingDuration = editingDurationInput.get();
*/

        //TODO check that the indices match

        trimZeroProbsDict.put(0,0,trimZeroProbs[0]);
        trimZeroProbsDict.put(0,1,trimZeroProbs[1]);
        trimZeroProbsDict.put(1,0,trimZeroProbs[2]);
        trimZeroProbsDict.put(1,1,trimZeroProbs[3]);


        //CREATING THE targStatTransitionsDict is empty!!! have to redo that
        targStatTransitionsDict = TargetStatus.getAllTransitions(numTargets);
        Pair<DoubleMatrix,Hashtable<IndelSet.TargetTract,Integer>> doubleMatrixHashtablePair = createAllTargetTractHazards();

        targetTractHazards = doubleMatrixHashtablePair.getFirst();
        targetTractDict = doubleMatrixHashtablePair.getSecond();
        targ_stat_transition_hazards_dict = new Hashtable<>();

        for(TargetStatus stat:targStatTransitionsDict.keySet()) {
            Hashtable<TargetStatus,DoubleMatrix> empty = new Hashtable<>();
            targ_stat_transition_hazards_dict.put(stat,empty);
        }

        hazardAwayDict = createHazardAwayDict();


        trim_long_params_reshaped = new DoubleMatrix(2,1);
        trim_long_params_reshaped.put(0,trim_long_params[0]);
        trim_long_params_reshaped.put(1,trim_long_params[1]);

        trim_short_params_reshaped = new DoubleMatrix(2,1);
        trim_short_params_reshaped.put(0,trim_short_params[0]);
        trim_short_params_reshaped.put(1,trim_short_params[1]);











    }

    public Double createTopologyLogLikelihood(beast.evolution.tree.TreeInterface tree, Hashtable<Integer, TransitionWrap> transitionWrappers) {
        Integer rootNodeID = tree.getRoot().getNr();
        Hashtable<Integer, DoubleMatrix> transMats = new Hashtable<>();
        //TODO check that
        //Hashtable<Integer, DoubleMatrix> Ddiags = new Hashtable<>();
        Hashtable<Integer, DoubleMatrix> ptMats = new Hashtable<>();
        Hashtable<Integer, DoubleMatrix> down_probs_dict= new Hashtable<>();
        Hashtable<Integer, DoubleMatrix> LProb = new Hashtable<>();
        Hashtable<Integer, Double> logScalingTerms = new Hashtable<>();
        for (Node node:tree.listNodesPostOrder(null,null)) {
            if(node.isLeaf()) {
                Log.info.println("NODE IS LEAF");

                TransitionWrap nodeWrap = transitionWrappers.get(node.getNr());
                DoubleMatrix probArray = DoubleMatrix.zeros(nodeWrap.numStatuses+1,1);
                Integer observedKey = nodeWrap.statusMap.get(nodeWrap.leafState);
                probArray.put(observedKey,1.0);
                LProb.put(node.getNr(),probArray);
            }
            else {
                Log.info.println("NODE IS NOT LEAF");

                TransitionWrap nodeWrap = transitionWrappers.get(node.getNr());
                DoubleMatrix logLprobNode = initializeLowerLogProb(nodeWrap,node);
                Log.info.println("initial logLprobNode"+ logLprobNode);

                DoubleMatrix hasPosProb = DoubleMatrix.ones(1);
                for (Node child: node.getChildren()) {

                    Log.info.println("CHILD NODE");

                    TransitionWrap childNodeWrap = transitionWrappers.get(child.getNr());
                    //TODO check the trim probs init
                    //trans_mats[child.node_id], trim_probs[child.node_id] = self._create_transition_matrix(
                    //                        child_wrapper)
                    DoubleMatrix childMat = createRateMatrix(childNodeWrap);
                    Log.info.println("childMat " + childMat);

                    transMats.put(child.getNr(),childMat);
                   /* DoubleMatrix expInput = new DoubleMatrix((childMat.toArray2()));*/
                    childMat.muli(node.getHeight() - child.getHeight());

                    DoubleMatrix ptMat = expm(childMat);
                    Log.info.println("ptMat " + ptMat);
                    //TODO check AFTER THIS

                    ptMats.put(child.getNr(),ptMat);
                    //TO DO check matrix multiplication IS CORRECT, JUST CAREFUL WITH DECIMALS !

                    DoubleMatrix chOrderedDownProbs = ptMat.mmul(LProb.get(child.getNr()));

                    DoubleMatrix downProbs = new DoubleMatrix();
                    if (! node.isRoot()) {
                        nodeWrap = transitionWrappers.get(node.getNr());
                        downProbs = reorderLikelihoods(chOrderedDownProbs,nodeWrap,childNodeWrap);

                    //TODO CHECK BEFORE THIS
                    }
                    //ROOT BLOCK!!!!
                    else {
                        Log.info.println("ROOT block !");
                        Integer chId = childNodeWrap.statusMap.get(new TargetStatus());
                        double downProb = chOrderedDownProbs.get(chId);
                        downProb = max(downProb,0.0);
                        downProbs = new DoubleMatrix(1,1,downProb);
                        down_probs_dict.put(chId,downProbs);
                        Log.info.println("en of ROOT block !");
                        //CORRECT UNTIL NOW!!!!
                    }

                    /*if (child.isLeaf()) {
                    //TODO change the abundance weight
                    Double leafAbundanceWeight = 1.0;
                    }*/
                    Double leafAbundanceWeight = 1.0;
                    hasPosProb = downProbs.ge(0);
                    //TODO LOSS OF DOWN PROBS DIMENSION!!!!

                    //ISSUE WITH THE ADDI ??
                    Log.info.println("logLprobNode before" + logLprobNode );
                    Log.info.println("downProbs before" + downProbs );
                    Log.info.println("hasPosProb" + hasPosProb );
                    Double dummy = downProbs.add((hasPosProb.neg()).add(1).mul( 1e-30 )).get(0);
                    NumberFormat formatter = new DecimalFormat();

                    Log.info.println("downProbs.addi((hasPosProb.neg()).add(1).mul( 1e-30 )) non zero!" + (formatter.format(dummy)));
                    Log.info.println("is it  > 0 ?" + (dummy > 0) );
                    Log.info.println("is it  = 0 ?" + (dummy == 0) );
                    Log.info.println("log of that" + log(dummy) );
                    Log.info.println("log of that" + log(dummy) );
                    //proctection against states with zero probability.

                    logLprobNode = logLprobNode.addi(logi(downProbs.addi((hasPosProb.neg()).add(1).mul( 1e-30 ))))  ;
                    Log.info.println("logLprobNode after" + logLprobNode);
                    Log.info.println("END CHILD NODE ITER");

                }
                Log.info.println("logLprobNode" + logLprobNode);
                Double logScalingTerm = logLprobNode.max();
                LProb.put(node.getNr(),expi((logLprobNode.subi(logScalingTerm)).muli(hasPosProb)));
                logScalingTerms.put(node.getNr(), logScalingTerm);
            }



        }
        Log.info.println("LProb" + LProb.get(0).get(1) *1000000000);
        Collection<Double> logScalingTermsAll = logScalingTerms.values();
        Log.info.println("logScalingTermsAll" + logScalingTermsAll);

        List<Double> terms = new ArrayList<>(logScalingTermsAll);
        Double fullScaling = 0.0;
        for(int i = 0; i < terms.size();++i) {
            fullScaling = fullScaling + terms.get(i);
        }
        DoubleMatrix logLikAlleles = log(LProb.get(rootNodeID)).addi(fullScaling);
        Log.info.println("logLikAlleles" + logLikAlleles);
        Log.info.println("fullScaling" + fullScaling);
        /*Lprob = Lprob;
        down_probs_dict = down_probs_dict;
        pt_matrix = pt_matrix;
       trans_mats = trans_mats;
        trim_probs = trim_probs;*/
        return logLikAlleles.get(0);

    }


    DoubleMatrix reorderLikelihoods(DoubleMatrix orderedDownProbs, TransitionWrap newWrapper, TransitionWrap oldWrapper) {
      List<Pair<Integer,Double>> indexVals = new ArrayList<>();
        for (TargetStatus newstat: newWrapper.transStatuses ) {
            if (oldWrapper.statusMap.containsKey(newstat)) {
                indexVals.add(new Pair(newWrapper.statusMap.get(newstat),orderedDownProbs.get(oldWrapper.statusMap.get(newstat))));
            }
        }
        DoubleMatrix downProbs = new DoubleMatrix(newWrapper.numStatuses + 1,1);
        for (int i=0; i<indexVals.size();i++) {
            downProbs.put(indexVals.get(i).getFirst(),indexVals.get(i).getSecond());
        }
      return downProbs;
    }



    public DoubleMatrix initializeLowerLogProb(TransitionWrap nodeWrap, Node node) {
        //HERE WE ARE WORKING WITH BIFURCATIONS, THERE ARE NO UNRESOLVED MULTIFURACTIONS, as in the original code.
        DoubleMatrix zero = DoubleMatrix.zeros(1);
        return zero;
    }

    public Pair<DoubleMatrix,Hashtable<IndelSet.TargetTract,Integer>> createAllTargetTractHazards() {
        TargetStatus targetStatusAllActive = new TargetStatus();
        List<IndelSet.TargetTract> allTargetTracts = targetStatusAllActive.getPossibleTargetTracts(new ArrayList<>(),metaData.nTargets);
        Hashtable<IndelSet.TargetTract,Integer> ttDict = new Hashtable<>();
        for(int i = 0;i<allTargetTracts.size();i++) {
            ttDict.put(allTargetTracts.get(i),i);
        }
        int nTTs = allTargetTracts.size();
        DoubleMatrix minTargets = new DoubleMatrix(nTTs);
        DoubleMatrix maxTargets = new DoubleMatrix(nTTs);
        DoubleMatrix longLeftStatuses = new DoubleMatrix(nTTs);
        DoubleMatrix longRightStatuses = new DoubleMatrix(nTTs);
        for(int i=0;i<allTargetTracts.size();i++ ) {
            IndelSet.TargetTract currentTT = allTargetTracts.get(i);
            minTargets.put(i,currentTT.getminTarg());
            maxTargets.put(i,currentTT.getmaxTarg());
            longLeftStatuses.put(i,currentTT.isLeftLong() ? 1.0 : 0.0 );
            longRightStatuses.put(i,currentTT.isRightLong() ? 1.0 : 0.0 );

        }
        DoubleMatrix allHazards = createHazardTargetTract(minTargets, maxTargets, longLeftStatuses, longRightStatuses);
        return new Pair(allHazards,ttDict);


    }





    public DoubleMatrix createHazardTargetTract(DoubleMatrix minTargets, DoubleMatrix maxTargets,DoubleMatrix longLeftStatuses, DoubleMatrix longRightStatuses ) {

        //TODO CAREFUL BECAUSE DOUBLEMATRIX FUNCTIONS ARE IN PLACE FUNCTIONS!!
        DoubleMatrix logLeftTrimFactor = new DoubleMatrix(longLeftStatuses.length);
        DoubleMatrix logRightTrimFactor = new DoubleMatrix(longRightStatuses.length);
        for(int i=0;i< longLeftStatuses.length;i++ ){
            if(longLeftStatuses.get(i) == 1.0) {
                logLeftTrimFactor.put(i,log(longTrimFactors[0])); }
            else { logLeftTrimFactor.put(i,0.0); }
            if(longRightStatuses.get(i) == 1.0) {
                logRightTrimFactor.put(i,log(longTrimFactors[1])); }
            else { logRightTrimFactor.put(i,0.0); }
        }



        DoubleMatrix gatheredTargetLamMin = new DoubleMatrix(minTargets.length);
        DoubleMatrix gatheredTargetLamMax = new DoubleMatrix(minTargets.length);
        for(int i=0;i<minTargets.length;i++) {
            gatheredTargetLamMin.put(i,cutRates[(int) minTargets.get(i)]);
            gatheredTargetLamMax.put(i,cutRates[(int) maxTargets.get(i)]);

        }


        DoubleMatrix gatherMinandMax = DoubleMatrix.zeros(minTargets.length);
        gatherMinandMax.addi(gatheredTargetLamMax);
        gatherMinandMax.addi(gatheredTargetLamMin);
        DoubleMatrix logFocalLambdaPart = DoubleMatrix.zeros(minTargets.length);
        logFocalLambdaPart.addi(logLeftTrimFactor);
        logFocalLambdaPart.addi(logRightTrimFactor);
        logFocalLambdaPart.addi(logi(gatheredTargetLamMin));
        DoubleMatrix logDoubleLambdaPart = DoubleMatrix.zeros(minTargets.length);

        logDoubleLambdaPart.addi(logLeftTrimFactor);
        logDoubleLambdaPart.addi(logRightTrimFactor);
        logDoubleLambdaPart.addi(logi(gatherMinandMax));
        logDoubleLambdaPart.addi(log(doubleCutWeight));
        //HERE SHOULD GET SOME POSITIVE VALUES!!

       /* (log_left_trim_factor
                + tf.log(tf.gather(self.target_lams, min_target) + tf.gather(self.target_lams, max_target))
                + log_right_trim_factor
                + tf.log(self.double_cut_weight))*/



        DoubleMatrix hazard = new DoubleMatrix(minTargets.length);

        for(int i=0;i<minTargets.length;i++) {
            if(minTargets.get(i) == maxTargets.get(i)) {
                hazard.put(i,logFocalLambdaPart.get(i));
            }
            else {
                hazard.put(i,logDoubleLambdaPart.get(i));
            }

        }

        expi(hazard);
        return hazard;

    }



    /*
    Calculate transition probability matrix for loss without editing
    The transition probability matrix for loss with editing differs only in the first row.
     */
    public void getLossProbabilities(double[] matrix, double expOfDeltaLoss){

        // fill diagonal and final column
        for (int i=0; i<nrOfStates; i++){
            for (int j=0; j<nrOfStates; j++){

                if ( i==j ){
                    matrix[i*nrOfStates + j] = expOfDeltaLoss;
                }else if(j == nrOfStates-1){
                    matrix[i*nrOfStates + j] = 1 - expOfDeltaLoss;
                }else{
                    matrix[i*nrOfStates + j] = 0;
                }
            }
        }
        // set final diagonal element to 1
        matrix[nrOfStates * nrOfStates - 1] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {

       /* double delta = startTime - endTime;
        double expOfDeltaLoss = Math.exp(-delta * lossRate);

        // calculate transition probabilities for loss process
        getLossProbabilities(matrix, expOfDeltaLoss);

        // for loss & editing, add the editing transition probabilities
        if ( (endTime >= (editingHeight - editingDuration)) & (endTime < editingHeight)){

            Stream<Double> scarSum = Stream.of(scarRates);
            Double scarRateSum = scarSum.reduce(0.0, (subtotal, element) -> subtotal + element);
            //.sum();

            // fill first row
            matrix[0] = Math.exp(-delta * (lossRate + scarRateSum));
            for (int i=0; i<nrOfStates-2; i++){
                matrix[i+1] = (scarRates[i] * expOfDeltaLoss - scarRates[i] * Math.exp(-delta * (lossRate + scarRateSum))) / scarRateSum;
            }
            matrix[nrOfStates-1] = 1 - expOfDeltaLoss;
        }*/
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {return null;}

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    /*public double getEditingHeight(){return editingHeight;}

    public double getEditingDuration(){return editingDuration;}*/

    public double[][] getRateMatrix(){return rateMatrix;}

    public DoubleMatrix createRateMatrix(TransitionWrap wrapper) {

        Integer matrxLen = wrapper.numStatuses + 1;
        DoubleMatrix targTractLeft = createMarginalTransitionMatrixLeft(wrapper);
        DoubleMatrix trimProbsLeft = createTrimInstantProbMatrixLeft(wrapper);

        //TODO CHECK IF TENSORS ARE FINITE


        DoubleMatrix qMatrix = targTractLeft.muli(trimProbsLeft);
        DoubleMatrix lastColumn = qMatrix.rowSums();
        lastColumn.muli(-1);
        qMatrix = concatHorizontally(qMatrix,lastColumn);

        return qMatrix;

    }


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
                    if (targ_stat_transition_hazards_dict.get(startState).containsKey(endState)) {
                        hazard = targ_stat_transition_hazards_dict.get(startState).get(endState).get(0);
                    } else {
                        List<Integer> hazard_idxs = new ArrayList<>();
                        for (IndelSet.TargetTract tt : targetTractsTransition) {
                            hazard_idxs.add(targetTractDict.get(tt));
                        }
                        for (int i : hazard_idxs) {
                            hazard = hazard + targetTractHazards.get(i);
                        }
                        Hashtable<TargetStatus, DoubleMatrix> toAppend = targ_stat_transition_hazards_dict.get(startState);
                        DoubleMatrix entry = new DoubleMatrix(1);
                        entry.put(0,hazard);
                        toAppend.put(endState, entry);
                        targ_stat_transition_hazards_dict.put(startState, toAppend);

                    }
                    sparseIndices.add(new Pair(startKey, endKey));
                    sparseValues.add(hazard);


                }
            }
        }

        Integer matrixLength = wrapper.numStatuses + 1;

        DoubleMatrix qSingleTtMatrix = DoubleMatrix.zeros(matrixLength,matrixLength-1);
        if (singleTtGatherIndices.size() != 0) {

            DoubleMatrix singleTtSparseVals = new DoubleMatrix(singleTtGatherIndices.size());
            for(int i=0;i<singleTtGatherIndices.size();i++) {
               singleTtSparseVals.put(i,targetTractHazards.get(singleTtGatherIndices.get(i)));
            }
            for (int i=0;i<singleTtSparseIndices.size();i++) {
                Pair<Integer,Integer> indices = singleTtSparseIndices.get(i);
                qSingleTtMatrix.put(indices.getFirst(),indices.getSecond(),singleTtSparseVals.get(i));
            }

        }


        DoubleMatrix qAllTtMatrix = DoubleMatrix.zeros(matrixLength,matrixLength-1);

        for (int i=0;i<sparseIndices.size();i++) {

            Pair<Integer, Integer> indices = sparseIndices.get(i);
            qAllTtMatrix.put(indices.getFirst(), indices.getSecond(), sparseValues.get(i));
        }
        DoubleMatrix qmatrixLeft = qAllTtMatrix;

        qmatrixLeft.addi(qSingleTtMatrix);
        return qmatrixLeft;

    }

    public DoubleMatrix createTrimInstantProbMatrixLeft(TransitionWrap childWrap){
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
        for(int i = 0;i<numTargets;i++) {
            targetToSingleton.add(null);
        }
            for (IndelSet.Singleton sg: childSingletons) {
                targetToSingleton.set(sg.getminTarg(),sg);
                targetToSingleton.set(sg.getmaxTarg(),sg);


            }

            Set<TargetStatus> possibleStates = new HashSet<>();
            possibleStates.addAll(childWrap.transStatuses);
            List<Pair<Integer,Integer>> sparseIndices = new ArrayList<>();
            List<Double> sparseVals = new ArrayList<>();
            for (TargetStatus startTargetStatus : possibleStates) {
                Set<TargetStatus> allEndStates = new HashSet<>();
                allEndStates.addAll(targStatTransitionsDict.get(startTargetStatus).keySet());
                Set<TargetStatus> possibleEndState = new HashSet<>(allEndStates);
                possibleEndState.retainAll(possibleStates);

                for(TargetStatus endTargetStatus:possibleEndState) {
                    Set<Integer> newDeacTargs = endTargetStatus.minus(startTargetStatus);
                    Set<IndelSet.Singleton> singletonSet = new HashSet<>();

                    for (Integer deactTarg: newDeacTargs) {
                        if(targetToSingleton.size() != 0 && targetToSingleton.get(deactTarg) != null) {
                            singletonSet.add(targetToSingleton.get(deactTarg));
                        }
                    }
                    Integer startKey = childWrap.statusMap.get(startTargetStatus);
                    Integer endKey = childWrap.statusMap.get(endTargetStatus);
                    if( singletonSet.size() ==1) {
                        //CHECK THAT
                        sparseIndices.add(new Pair(startKey,endKey));
                        IndelSet.Singleton sg = (IndelSet.Singleton) singletonSet.toArray()[0];
                        Double trimProbVal = singletonCondProb.get(singletonIndexDict.get(sg));
                        sparseVals.add(trimProbVal -1 );
                        Log.info.println(" trimProbVal" + trimProbVal);
                        Log.info.println(" trimProbVal -1" + (trimProbVal -1));
                        //CHECK THAT

                    }
                    else if (singletonSet.size()>1) {
                        sparseIndices.add(new Pair(startKey,endKey));
                        sparseVals.add(-1.0);
                    }

                }





            }
        int outputLength = childWrap.numStatuses+1;
        int[] outputShape = new int[] {outputLength,outputLength -1};
        if(sparseVals.size() !=0) {
        DoubleMatrix zeroes = DoubleMatrix.zeros(outputShape[0],outputShape[1]);
        DoubleMatrix tobeReturned = DoubleMatrix.ones(outputShape[0],outputShape[1]);
        for(int i=0;i<sparseIndices.size();i++) {
            Log.info.println(" tobereturned sparsevals" + (sparseVals.get(i)));
            Log.info.println("filling tobereturned" + (1.0 + sparseVals.get(i)));
            tobeReturned.put(sparseIndices.get(i).getFirst(),sparseIndices.get(i).getSecond(),1.0 + sparseVals.get(i));

        }
        tobeReturned.max(zeroes);
        return tobeReturned;
        }
        else{
            return DoubleMatrix.ones(outputShape[0],outputShape[1]);
        }

    }


    public DoubleMatrix createHazardAwayTargetStatuses(List<TargetStatus> targetStatuses) {

        //For each target status create a corresponding vector containing the associated cut rate or 0 if it is already
        //inactive

        DoubleMatrix activeTargetHazards = DoubleMatrix.zeros(targetStatuses.size(),numTargets);
        DoubleMatrix activeMasks = DoubleMatrix.zeros(targetStatuses.size(),numTargets);
        for(int i=0;i<targetStatuses.size();++i) {
            Integer[] binaryrep = targetStatuses.get(i).getBinaryStatus(numTargets);
            for(int j=0; j<numTargets; j++) {
                activeMasks.put(i,j, 1 - binaryrep[j]);
                activeTargetHazards.put(i,j,activeMasks.get(i,j) * cutRates[j]);
            }
        }


        // TODO CHECK THAT if the barcode only contains a single target, just return the cut rates as is:
        if(numTargets == 1) {
            return activeTargetHazards.getColumn(0);
        }


        int[] indicesInner = new int[activeTargetHazards.getColumns()-2];
        for (int i=1;i <activeTargetHazards.getColumns()-1;++i) {
            indicesInner[i-1] = i;
        }


        DoubleMatrix focalHazards = ((activeTargetHazards.getColumn(0)).muli(1+ longTrimFactors[1])).addi(((activeTargetHazards.getColumns(indicesInner)).rowSums()).muli((1+ longTrimFactors[1])*(1+longTrimFactors[0]))).addi((activeTargetHazards.getColumn(activeTargetHazards.getColumns()-1)).muli(1+ longTrimFactors[0]));
        DoubleMatrix middleHazards = (activeTargetHazards.getColumns(indicesInner).rowSums()) ;
        DoubleMatrix numMiddle = (activeMasks.getColumns(indicesInner).rowSums()) ;
        DoubleMatrix middleDoubleCutHazards = middleHazards.muli((numMiddle.subi(1)).muli((1 + longTrimFactors[0])*(1 + longTrimFactors[1])));
        //(1 + self.trim_long_factor[0]) * (1 + self.trim_long_factor[1]) * (num_in_middle - 1) * middle_hazards

        //CORRECT UNTIL HERE

        DoubleMatrix numAfterStarts = (activeMasks.getColumns(indicesInner).rowSums());
        DoubleMatrix startToOtherHazard = activeMasks.getColumn(0).muli(1 + longTrimFactors[1]).muli((numAfterStarts.muli(activeTargetHazards.getColumn(0))).addi((activeTargetHazards.getColumns(indicesInner).rowSums())));
       /* start_to_other_hazard = active_masks[:,0] * (1 + self.trim_long_factor[1]) * (
                num_after_start * active_targ_hazards[:, 0] + tf.reduce_sum(active_targ_hazards[:, 1:-1], axis=1))*/

        DoubleMatrix numBeforeEnd = (activeMasks.getColumns(indicesInner).rowSums());

        DoubleMatrix endToOtherHazard = (activeMasks.getColumn(activeMasks.getColumns()-1)).muli(1 + longTrimFactors[0]).muli((numBeforeEnd).muli((activeTargetHazards.getColumn(activeTargetHazards.getColumns()-1))).addi((activeTargetHazards.getColumns(indicesInner).rowSums())));

       /* active_masks[:,-1] * (1 + self.trim_long_factor[0]) * (
                num_before_end * active_targ_hazards[:, -1] + tf.reduce_sum(active_targ_hazards[:, 1:-1], axis=1))*/
        //correct!
        DoubleMatrix startToEndHazards = (activeMasks.getColumn(0)).muli(activeMasks.getColumn(activeMasks.getColumns()-1)).muli((activeTargetHazards.getColumn(0).addi(activeTargetHazards.getColumn(activeTargetHazards.getColumns()-1))));
        /*active_masks[:,0] * active_masks[:,-1] * (
                active_targ_hazards[:, 0] + active_targ_hazards[:, -1])*/

        DoubleMatrix hazardAwayNodes = focalHazards.addi((middleDoubleCutHazards.addi(startToOtherHazard).addi(endToOtherHazard).addi(startToEndHazards)).muli(doubleCutWeight));

        return hazardAwayNodes;



    }

    public static List<IndelSet.Singleton> getAllSingletons(beast.evolution.tree.TreeInterface tree,Hashtable<Integer, AncStates> statesDict ) {
        Set<IndelSet.Singleton> singletonSet = new HashSet();
        for (Node node: tree.getExternalNodes()) {
            if(node.isLeaf()) {
                AncStates leaf_state = statesDict.get(node.getNr());
                for(IndelSet sgwc :leaf_state.getSingletonWCs()) {
                    IndelSet.Singleton sing = sgwc.getSingleton();
                    singletonSet.add(sing);
                }
            }
        }

        List<IndelSet.Singleton> singletonList = new ArrayList<>();
        singletonList.addAll(singletonSet);
        return singletonList;



    }

    /*def get_all_singletons(topology: CellLineageTree):
    singletons = set()
        for leaf in topology:
            for leaf_anc_state in leaf.anc_state_list:
            for singleton_wc in leaf_anc_state.indel_set_list:
    sg = singleton_wc.get_singleton()
            singletons.add(sg)
            return singletons*/


    public void initSingletonProbs(List<IndelSet.Singleton> singletons) {


        for (int i=0; i<singletons.size();i++) {
            singletonIndexDict.put(singletons.get(i),i);
        }
        createTrimInsertDistributions(singletons.size());
        singletonLogCondProb = createLogIndelProbs(singletons);
        Log.info.println("singletonLogCondProb" + singletonLogCondProb);
        singletonCondProb = expi(singletonLogCondProb);
        Log.info.println("singletonCondProb" + singletonCondProb);


    }

    public void createTrimInsertDistributions(int numSingletons) {

        delShortDist = makeDelDist(trim_short_params_reshaped, 2, usePoisson);
        delLongDist = makeDelDist(trim_long_params_reshaped, 2, true);
        if (usePoisson) {

             insertDist = new Poisson(exp(insertParams[0]),RandomEngine.makeDefault());
       }
       /*else {

            insertDist = new NegativeBinomial(insertParams[1].intValue(),exp(insertParams[0]),RandomEngine.makeDefault());
        }*/

    }

    public List<AbstractDiscreteDistribution> makeDelDist(DoubleMatrix params, int nTrimTypes, boolean usePoisson) {
        List<AbstractDiscreteDistribution> delDistList = new ArrayList<>();
        if (usePoisson) {
            for(int i=0;i < nTrimTypes;i++) {
                delDistList.add(new Poisson(exp(params.get(i)),RandomEngine.makeDefault()));
            }
        }
        else {
            for(int i=0;i < nTrimTypes;i++) {
                for(int j=0;j< numTargets;j++) {
                    delDistList.add(new NegativeBinomial((int) params.get(i,1),exp(params.get(i,0)), RandomEngine.makeDefault()));
                }
            }

        }
        return delDistList;

    }

    public DoubleMatrix createLogIndelProbs(List<IndelSet.Singleton> singletons) {
        if(singletons.size() == 0) {
            return new DoubleMatrix();
        }
        else {
            // assemble individual probabilities:
            DoubleMatrix leftDelProb = createLeftDelProbs(singletons);
            DoubleMatrix rightDelProb = createRightDelProbs(singletons);
            DoubleMatrix insertProb = createInsertProbs(singletons);
            Log.info.println("insertProb" + insertProb);
            Log.info.println("rightDelProb" + rightDelProb);
            Log.info.println("leftDelProb" + leftDelProb);
            //combine everything
            DoubleMatrix allLogProbs = (logi(leftDelProb)).addi(logi(rightDelProb)).addi(logi(insertProb));

            Log.info.println("allLogProbs" + allLogProbs);

            //SOMETHING TO INITIALIZE
            Double logShortFocalNormalization = log(1-(trimZeroProbsDict.get(0,0)*trimZeroProbsDict.get(1,0)*insertZeroProb));
            Log.info.println("logShortFocalNormalization" + logShortFocalNormalization);

            boolean[] isLongIndel = new boolean[singletons.size()];
            for (int i=0; i<singletons.size();i++) {
                isLongIndel[i] =  (singletons.get(i).isLeftLong() || singletons.get(i).isRightLong() || singletons.get(i).isIntertarget());
            }
            DoubleMatrix fullLogIndelProbs = new DoubleMatrix(singletons.size());
            for (int i=0;i<singletons.size();i++ ) {
                if (isLongIndel[i]) {
                    fullLogIndelProbs.put(i,allLogProbs.get(i));
                }
                else {
                    fullLogIndelProbs.put(i,allLogProbs.get(i) - logShortFocalNormalization);
                }

            }
           return fullLogIndelProbs;

        }


    }


    public DoubleMatrix createInsertProbs(List<IndelSet.Singleton> singletons) {

        List<Integer> insertLengths = new ArrayList<>();
        for(IndelSet.Singleton sg:singletons) {
            insertLengths.add(sg.getInsertLength());
        }
        /*DoubleMatrix insertSeqProb = new DoubleMatrix();
        for (int length:insertLengths) {
            insertSeqProb.add(pow(4,length));
        }*/
        DoubleMatrix insertLenProb = new DoubleMatrix(insertLengths.size());
        for(int i=0; i<insertLengths.size();i++) {
            insertLenProb.put(i,insertDist.pdf(max((insertLengths.get(i)-1),0)));
        }
        DoubleMatrix Return =new DoubleMatrix(insertLengths.size());
        for(int i=0;i<insertLengths.size();i++) {
            if(insertLengths.get(i) == 0) {
                Return.put(i,insertZeroProb);
            }
            else { Return.put(i,(1-insertZeroProb)*insertLenProb.get(i));}
        }

        return Return;
    }



    public DoubleMatrix createLeftDelProbs(List<IndelSet.Singleton> singletons) {
        List<Integer> minTargets = new ArrayList<>();
        DoubleMatrix isLeftLongs = new DoubleMatrix(singletons.size());
        DoubleMatrix isIntertargets = new DoubleMatrix(singletons.size());
        DoubleMatrix startpos = new DoubleMatrix(singletons.size());
        for(int i=0;i<singletons.size();i++) {
            IndelSet.Singleton sg = singletons.get(i);
            minTargets.add(sg.getminTarg());
            isLeftLongs.put(i,sg.isLeftLong() ? 1.0 : 0.0);
            isIntertargets.put(i,sg.isIntertarget() ? 1.0 : 0.0);
            startpos.put(i,sg.getStartPos());
        }

        DoubleMatrix minTargetSites = new DoubleMatrix(minTargets.size());
        DoubleMatrix leftTrimLongMin = new DoubleMatrix(minTargets.size());
        DoubleMatrix leftTrimLongMax = new DoubleMatrix(minTargets.size());
        for(int i=0;i<minTargets.size();i++) {
            minTargetSites.put(i,metaData.absCutSites.get(minTargets.get(i)));
            leftTrimLongMin.put(i,metaData.leftLongTrimMin.get(minTargets.get(i)));
            leftTrimLongMax.put(i,metaData.leftMaxTrim.get(minTargets.get(i)));
        }
        DoubleMatrix leftTrimLens = minTargetSites.subi(startpos);



        return createDelProbs(leftTrimLens,isLeftLongs,isIntertargets,leftTrimLongMin,leftTrimLongMax,false);


    }



    public DoubleMatrix createRightDelProbs(List<IndelSet.Singleton> singletons) {
        List<Integer> maxTargets = new ArrayList<>();
        DoubleMatrix isRightLongs = new DoubleMatrix(singletons.size());
        DoubleMatrix isIntertargets = new DoubleMatrix(singletons.size());
        DoubleMatrix endpos = new DoubleMatrix(singletons.size());
        for(int i=0;i<singletons.size();i++) {
            IndelSet.Singleton sg = singletons.get(i);
            maxTargets.add(sg.getmaxTarg());
            isRightLongs.put(i,sg.isRightLong() ? 1.0 : 0.0);
            isIntertargets.put(i,sg.isIntertarget() ? 1.0 : 0.0);
            endpos.put(i,sg.getEndPos());
        }
        DoubleMatrix maxTargetSites = new DoubleMatrix(singletons.size());
        DoubleMatrix rightTrimLongMin = new DoubleMatrix(singletons.size());
        DoubleMatrix rightTrimLongMax = new DoubleMatrix(singletons.size());
        for(int i=0;i<maxTargets.size();i++) {
            maxTargetSites.put(i,metaData.absCutSites.get(maxTargets.get(i)));
            rightTrimLongMin.put(i,metaData.rightLongTrimMin.get(maxTargets.get(i)));
            rightTrimLongMax.put(i,metaData.rightMaxTrim.get(maxTargets.get(i)));
        }
        //CAREFUL HERE!
        DoubleMatrix rightTrimLen = endpos.subi(maxTargetSites);


       /* max_target_sites = tf.constant([self.bcode_meta.abs_cut_sites[mt] for mt in max_targets], dtype=tf.float64)
        right_trim_len = del_ends - max_target_sites

        right_trim_long_min = tf.constant([self.bcode_meta.right_long_trim_min[mt] for mt in max_targets], dtype=tf.float64)
        right_trim_long_max = tf.constant([self.bcode_meta.right_max_trim[mt] for mt in max_targets], dtype=tf.float64)*/

        return createDelProbs(rightTrimLen,isRightLongs,isIntertargets,rightTrimLongMin,rightTrimLongMax,true);


    }




    public DoubleMatrix createDelProbs(DoubleMatrix trimLen,DoubleMatrix isLongs,DoubleMatrix isIntertargets,DoubleMatrix trimMins,DoubleMatrix trimMaxs,boolean isRight) {
    int esRight = (isRight ? 1 : 0);

    DoubleMatrix trimZeroProbs = new DoubleMatrix(isIntertargets.length);
    for(int i=0;i<isIntertargets.length;i++ ) {
        if (isIntertargets.get(i) == 1) {
            trimZeroProbs.put(i, trimZeroProbsDict.get(esRight, 1));
        } else {
            trimZeroProbs.put(i, trimZeroProbsDict.get(esRight, 0));
        }
    }

/*

        del_short_dist = self.del_short_dist[is_right]
        short_nonzero_prob = (1 - trim_zero_prob) * del_short_dist.prob(tf.maximum(trim_len - 1, 0))/del_short_dist.cdf(trim_maxs - 1)
*/

        Poisson delShortD = (Poisson) delShortDist.get(0);
        DoubleMatrix shortNonZeroProb = new DoubleMatrix(trimZeroProbs.length);
        for(int i=0; i<trimZeroProbs.length;i++) {
            int pdfInput = (int) max(trimLen.get(i)-1,0);
            int cdfInput = (int) trimMaxs.get(i)-1;
            shortNonZeroProb.put(i,(1- trimZeroProbs.get(i))* delShortD.pdf(pdfInput)/(delShortD.cdf(cdfInput)));
        }
        //(1- trimZeroProbs.get(i))* delShortD.pdf(pdfInput)/(delShortD.cdf(cdfInput))


        Poisson delLongD = (Poisson) delLongDist.get(esRight);
        DoubleMatrix LongProb = new DoubleMatrix(trimZeroProbs.length);
        for(int i=0; i<trimZeroProbs.length;i++) {
            int pdfInput = (int) max(trimLen.get(i) - trimMins.get(i), 0);
            int cdfInput = (int) max(trimMaxs.get(i) - trimMins.get(i), 0);
            LongProb.put(i,delLongD.pdf(pdfInput)/(delLongD.cdf(cdfInput)));
        }

        DoubleMatrix Return = new DoubleMatrix(trimLen.length);
        for(int i=0;i<trimLen.length;i++) {
            if(trimLen.get(i) == 0) {
                Return.put(i,trimZeroProbs.get(i));
            }
            else {
                if(isLongs.get(i) == 1) {
                    Return.put(i,LongProb.get(i));
                }
                else {
                    Return.put(i,shortNonZeroProb.get(i));
                }
            }
        }

    return Return;

    }




    public Hashtable<TargetStatus,Double> createHazardAwayDict() {
        Set<TargetStatus> targetStatuses = targStatTransitionsDict.keySet();

        List<TargetStatus> targetStatusesList = new ArrayList<>();
        targetStatusesList.addAll(targetStatuses);
        DoubleMatrix hazardAwayNodes = createHazardAwayTargetStatuses(targetStatusesList);

        Hashtable<TargetStatus, Double> hazardAwayDict = new Hashtable<>();
        for(int i=0;i<targetStatusesList.size();i++) {
            hazardAwayDict.put(targetStatusesList.get(i),hazardAwayNodes.get(i));
        }
        return hazardAwayDict;

    }

    public Double createTopologyLogLikBarcode(Hashtable<Integer,TransitionWrap> transitionWrapper) {






        return 0.0;
    }

}


