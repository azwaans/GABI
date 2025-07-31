package gestalt.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.State;
import gestalt.evolution.alignment.*;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.apache.commons.math3.util.Pair;
import org.jblas.DoubleMatrix;

import java.io.IOException;
import java.util.*;

import static org.jblas.MatrixFunctions.logi;


@Description("Generic tree likelihood for an alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model")
// Use this as base class to define any non-standard TreeLikelihood.
// Override Distribution.calculatLogP() to make this class functional.
// TODO: This could contain a generic traverse() method that takes dirty trees in account.

public class gestaltTreeLikelihood extends GenericTreeLikelihood {

    protected gestaltGeneral substitutionModel;

    protected SiteModel.Base m_siteModel;

    protected BranchRateModel.Base branchRateModel;

    protected double[] m_branchLengths;

    protected double[] storedBranchLengths;

    protected List<IndelSet.Singleton> singletonList;

    protected Hashtable<Integer, String> processedAlignment;


    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    protected gestaltLikelihoodCore likelihoodCore;


    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    public void initAndValidate() {
        // sanity check: site model should be an instance of the base site model class
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        int nodeCount = treeInput.get().getNodeCount();

        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());

        substitutionModel = (gestaltGeneral) m_siteModel.substModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        likelihoodCore = new gestaltLikelihoodCore();

        //this only depends on parameters that govern indel types.
        //this is NOT tree dependent it accounts for ALL TRANSITIONS
        Pair<DoubleMatrix, Hashtable<IndelSet.TargetTract, Integer>> doubleMatrixHashtablePair = substitutionModel.createAllTargetTractHazards();
        substitutionModel.targetTractHazards = doubleMatrixHashtablePair.getFirst();
        substitutionModel.targetTractDict = doubleMatrixHashtablePair.getSecond();

        //update target status rates: this depends only on indel type parameters
        //this is also not tree dependent, because it accounts FOR ALL TRANSITIONS
        substitutionModel.hazardAwayDict = substitutionModel.createHazardAwayDict();
        initCore();

        //if the input format is GSM, process the sequences
        if (substitutionModel.GSMformat) {
            try {
                processedAlignment = processGSMAlignment(treeInput.get(), dataInput.get());
            } catch (IOException e) {
                e.printStackTrace();
            }

        }

    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.init(nodeCount);
        hasDirt = Tree.IS_FILTHY;

    }

    public double calculateLogP() {
        final TreeInterface tree = treeInput.get();
        //recording time
        //long start1 = System.nanoTime();

        //if the tree topology has changed, update everything
        //todo make this a traversal/part of the traversal
        if (likelihoodCore.transitionWraps == null || isSomeNodeFilthy()) {
            //this.initAndValidate();
            hasDirt = Tree.IS_FILTHY;
            //finds the set of likely ancestral states at each internal node based on leaf sequences
            Hashtable<Integer, AncStates> statesDict = new Hashtable<>();
            if (substitutionModel.GSMformat) {
                statesDict = TransitionWrap.createStatesDictGSM(tree, processedAlignment, substitutionModel.metaData.posSites, substitutionModel.metaData.nTargets);
            } else {
                statesDict = TransitionWrap.createStatesDict(tree, dataInput.get(), substitutionModel.metaData.posSites, substitutionModel.metaData.nTargets);

            }

            //extract all single indel events from the possible ancestral states to compute conditional probabilities
            singletonList = gestaltGeneral.getAllSingletons(tree, statesDict);

            //each possible state at each node create a wrap of metadata
            likelihoodCore.transitionWraps = likelihoodCore.createTransitionWraps(tree, substitutionModel.metaData, statesDict);

        }
        //update combined parameters of the substitution model.
        //todo enable storing and restoring state for the following
        substitutionModel.initSingletonProbs(singletonList);

        //update target status rates: this depends only on mutation parameters
        substitutionModel.hazardAwayDict = substitutionModel.createHazardAwayDict();

        //update target tract rates: this depends only on mutation parameters
        Pair<DoubleMatrix, Hashtable<IndelSet.TargetTract, Integer>> doubleMatrixHashtablePair = substitutionModel.createAllTargetTractHazards();

        substitutionModel.targetTractHazards = doubleMatrixHashtablePair.getFirst();
        substitutionModel.targetTractDict = doubleMatrixHashtablePair.getSecond();

        //empty targStatTransitionHazardsDict
        substitutionModel.targStatTransitionHazardsDict = new Hashtable<>();

        for (TargetStatus stat : gestaltGeneral.targStatTransitionsDict.keySet()) {
            Hashtable<TargetStatus, DoubleMatrix> empty = new Hashtable<>();
            substitutionModel.targStatTransitionHazardsDict.put(stat, empty);
        }
        //long end1 = System.nanoTime();
        //System.out.println("Elapsed Time in seconds: "+ (end1-start1)*0.000000001);

        try {
            //populate partial likelihoods with postorder traversal
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN) {
                calcLogP();
            }

        } catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }
        return logP;
    }


    /* Assumes there IS a branch rate model as opposed to traverse() */
    protected int traverse(final Node node) {

        int update = hasDirt;

        if (node != null) {
            update = (node.isDirty() | hasDirt);
            final int nodeIndex = node.getNr();

            final double branchRate = branchRateModel.getRateForBranch(node);
            final double branchTime = node.getLength() * branchRate;

            // First update the transition probability matrix(ices) for this branch
            if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
                m_branchLengths[nodeIndex] = branchTime;
                likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
                TransitionWrap nodeWrap = likelihoodCore.transitionWraps.get(node.getNr() + 1);
                DoubleMatrix ptMat = substitutionModel.getTransitionProbabilities(node.getParent(), nodeWrap, node.getLength(), branchRate);
                likelihoodCore.setNodeMatrix(nodeIndex, ptMat);

                update |= Tree.IS_DIRTY;
            }

            // Only update the partial likelihoods if the node is a leaf
            if (node.isLeaf()) {
                likelihoodCore.setLeafPartials(node);
            } else {

                // Traverse down the two child nodes
                final Node child1 = node.getLeft(); //Two children
                final int update1 = traverse(child1);
                final Node child2 = node.getRight();
                final int update2 = traverse(child2);

                // If either child node was updated then update this node too
                if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                    update |= (update1 | update2);
                    likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                    DoubleMatrix nodeLogPartials = likelihoodCore.calculatePartials(node);
                    likelihoodCore.setNodePartials(node.getNr(), nodeLogPartials);


                }
            }
        }
        return update;
    } // traverseWithBRM


    void calcLogP() {

        final TreeInterface tree = treeInput.get();

        // take care of scaling terms
        Double fullScaling = 0.0;
        for (Node node : tree.getInternalNodes()) {
            fullScaling = fullScaling + likelihoodCore.logScalingTerms.get(likelihoodCore.makeCachingIndexScalingTerm(node.getNr()));
        }

        //get the likelihood from the root partial
        int rootIndex = tree.getRoot().getNr();
        DoubleMatrix logLikAlleles = logi(likelihoodCore.getNodePartials(rootIndex)).addi(fullScaling);
        //todo does this means that we may condition on unedited allele at the root (rather than at the origin?)
        logP = logLikAlleles.get(0);


    }


    public String processSequence(String GSMSequence) throws IOException {
        //create a list of target sequences
        String[] targetStrings = GSMSequence.split(",");

        List<GestaltEvent> ObservedAlignedSeq = processObservedSeqFormat7B(targetStrings, substitutionModel.mergingThreshold, substitutionModel.paddingLength, substitutionModel.minPosition);
        if (!ObservedAlignedSeq.isEmpty()) {
            String EventFormatSeq = concatenateObservedAlignedSeqString(ObservedAlignedSeq);
            return EventFormatSeq;
        } else {
            return "None,";
        }


    }

    public String concatenateObservedAlignedSeqString(List<GestaltEvent> AlignedEventSeq) {
        String concatenated = "";
        for (GestaltEvent i : AlignedEventSeq) {
            concatenated = concatenated + i.getString() + ",";
        }
        return concatenated;
    }

    public List<GestaltEvent> processObservedSeqFormat7B(String[] targetStrings, int mergeThresh, int padLength, int minPos) throws IOException {

        //Find list of targets for each event
        Hashtable<String, Pair<Integer, Integer>> evtTargetDict = new Hashtable<>();

        int targetIdx = 0;
        for (String targetString : targetStrings) {
            String[] targetEventStrings = targetString.split("&");
            for (String eventString : targetEventStrings) {
                if (!evtTargetDict.containsKey(eventString)) {
                    Pair<Integer, Integer> targets = new Pair(targetIdx, targetIdx);
                    evtTargetDict.put(eventString, targets);
                } else {
                    Pair<Integer, Integer> minTargMaxTarg = evtTargetDict.get(eventString);
                    Integer minTarg = java.lang.Math.min(minTargMaxTarg.getFirst(), targetIdx);
                    Integer maxTarg = java.lang.Math.max(minTargMaxTarg.getSecond(), targetIdx);
                    evtTargetDict.put(eventString, new Pair(minTarg, maxTarg));

                }
            }
            targetIdx++;
        }

        List<GestaltEvent> events = new ArrayList<>();
        for (String eventString : evtTargetDict.keySet()) {
            Pair<Integer, Integer> minTargMaxTarg = evtTargetDict.get(eventString);

            //handle unedited/unkwown
            if (!(eventString.equals("NONE") || eventString.equals("UNKNOWN"))) {
                GestaltEvent eventFormat7B = processEventFormat7B(eventString, minTargMaxTarg.getFirst(), minTargMaxTarg.getSecond(), minPos - padLength);
                events.add(eventFormat7B);

            }


        }
        //sort by startPos
        if (!events.isEmpty()) {
            Collections.sort(events);
            List<GestaltEvent> cleanedEvents = new ArrayList(events);
            int index = 0;
            for (GestaltEvent evt : events) {
                int newStartTarget = evt.getMinTarg();
                int newEndTarget = evt.getMaxTarg();
                if (evt.getStartPos() > substitutionModel.metaData.absCutSites.get(evt.getMinTarg())) {
                    //start position if after the min target. Proposal is to shift the min cut target
                    newStartTarget = evt.getMinTarg() + 1;

                }
                if (evt.getDelEnd() < substitutionModel.metaData.absCutSites.get(evt.getMaxTarg())) {
                    //end position if before the max target. Proposal is to shift the max cut target
                    newEndTarget = evt.getMaxTarg() - 1;

                }

                if (newStartTarget > newEndTarget) {
                    //If we shifted such that the targets don't make sense, we need to start over.
                    // This is probably a focal deletion that is not aligned with a cut site.
                    // Determine the new target
                    int newTarget = 0;
                    if (newStartTarget > substitutionModel.metaData.nTargets - 1) {
                        newTarget = newEndTarget;
                    } else if (newEndTarget < 0) {
                        newTarget = 0;
                    } else {
                        int cutSite1 = (int) substitutionModel.metaData.absCutSites.get(newStartTarget);
                        int cutSite2 = (int) substitutionModel.metaData.absCutSites.get(newEndTarget);
                        int minDist1 = java.lang.Math.min(java.lang.Math.abs(evt.getStartPos() - cutSite1), java.lang.Math.abs(evt.getDelEnd() - cutSite1));
                        int minDist2 = java.lang.Math.min(java.lang.Math.abs(evt.getStartPos() - cutSite2), java.lang.Math.abs(evt.getDelEnd() - cutSite2));
                        if (minDist1 < minDist2) {
                            newTarget = newStartTarget;
                        } else {
                            newTarget = newEndTarget;
                        }
                    }
                    newStartTarget = newTarget;
                    newEndTarget = newTarget;
                    int newStartPos = java.lang.Math.min(evt.getStartPos(), (int) substitutionModel.metaData.absCutSites.get(newStartTarget));
                    int newDelLen = java.lang.Math.max(evt.getDelEnd() - newStartPos, (int) substitutionModel.metaData.absCutSites.get(newStartTarget) - newStartPos);
                    //dummy insert sequence
                    String dummy = "a".repeat(newDelLen - evt.getDelLen());
                    String newInsertStr = evt.getInsSeq() + dummy;
                    GestaltEvent event = new GestaltEvent(newStartPos + "_" + newDelLen + "_" + newStartTarget + "_" + newEndTarget + "_" + newInsertStr);
                    //to do check if we need to add at certain index
                    cleanedEvents.set(index, event);


                } else {
                    GestaltEvent event = new GestaltEvent(evt.getStartPos() + "_" + evt.getDelLen() + "_" + newStartTarget + "_" + newEndTarget + "_" + evt.getInsSeq());
                    cleanedEvents.set(index, event);
                }
                ++index;


            }

            // Deal with multiple/compound events affecting the same targets.
            // We will merge these clashing events into a single event.
            List<GestaltEvent> nonClashingEvents = new ArrayList<>(cleanedEvents.subList(0, 1));
            for (GestaltEvent evt : cleanedEvents.subList(1, cleanedEvents.size())) {

                GestaltEvent prevEvt = nonClashingEvents.get(nonClashingEvents.size() - 1);
                int minDeactTarg = evt.getMinMaxDeactTargets(substitutionModel.metaData.posSites, (substitutionModel.metaData.nTargets)).getFirst();
                int maxDeactTarg = prevEvt.getMinMaxDeactTargets(substitutionModel.metaData.posSites, (substitutionModel.metaData.nTargets)).getSecond();

                boolean doMergeCauseClose = ((maxDeactTarg == minDeactTarg) && ((evt.getStartPos() - prevEvt.getDelEnd()) <= mergeThresh));
                boolean doMergeCauseImpossible = prevEvt.getMaxTarg() >= evt.getMinTarg();
                if (doMergeCauseClose || doMergeCauseImpossible) {

                    int newOmitStrLen = java.lang.Math.max(0, (evt.getStartPos() - prevEvt.getStartPos() >= prevEvt.getDelLen()) ? 1 : 0);

                    String dummy = "a".repeat(newOmitStrLen);
                    GestaltEvent newEvent = new GestaltEvent(prevEvt.getStartPos() + "_" + ((evt.getStartPos() - prevEvt.getStartPos()) + evt.getDelLen()) + "_" + prevEvt.getMinTarg() + "_" + java.lang.Math.max(prevEvt.getMaxTarg(), evt.getMaxTarg()) + "_" + prevEvt.getInsSeq() + evt.getInsSeq() + dummy);
                    nonClashingEvents.set(nonClashingEvents.size() - 1, newEvent);


                } else {
                    nonClashingEvents.add(evt);
                }
                //if not empty
                // Make sure the right trim length for the right-most target is not too long
                if (!(nonClashingEvents.size() == 0)) {
                    GestaltEvent lastEvt = nonClashingEvents.get(nonClashingEvents.size() - 1);
                    if (lastEvt.getMaxTarg() == substitutionModel.metaData.nTargets - 1) {
                        if (lastEvt.getDelEnd() > substitutionModel.metaData.origLength) {

                            nonClashingEvents.set(nonClashingEvents.size() - 1, new GestaltEvent(lastEvt.getStartPos() + "_" + (substitutionModel.metaData.origLength - lastEvt.getStartPos()) + "_" + lastEvt.getMinTarg() + "_" + lastEvt.getMaxTarg() + "_" + lastEvt.getInsSeq()));
                        }
                    }

                }

            }

            return nonClashingEvents;
        } else {
            return events;
        }

    }


    public static GestaltEvent processEventFormat7B(String eventString, Integer minTarget, Integer maxTarget, Integer minPos) throws IOException {
        //example "140D+141"
        String[] eventSplit = eventString.split("\\+");
        // ["140D", "141"]
        String eventTypeString = eventSplit[0].substring(eventSplit[0].length() - 1);
        // "D"
        Integer eventPos = Integer.parseInt(eventSplit[1]) - minPos;
        if (eventTypeString.equals("D")) {
            Integer delLen = Integer.parseInt(eventSplit[0].substring(0, eventSplit[0].length() - 1));
            //todo check final "_
            String constructedEventString = eventPos + "_" + delLen + "_" + minTarget + "_" + maxTarget + "_";
            return new GestaltEvent(constructedEventString);
        } else if (eventTypeString.equals("I")) {
            String constructedEventString = eventPos + "_" + 0 + "_" + minTarget + "_" + maxTarget + "_" + eventSplit[2];
            return new GestaltEvent(constructedEventString);
        } else {
            throw new IOException("\nUnexpected event while processing sequence!");

        }

    }

    /**
     * in original implementation: annotate_ancestral_states
     * finds all possible Ancestral States at internal nodes
     */
    public Hashtable<Integer, String> processGSMAlignment(beast.base.evolution.tree.TreeInterface tree,
                                                          Alignment alinmt) throws IOException {
        Hashtable<Integer, String> processedAlignmentDict = new Hashtable<>();
        Log.info.println("Processing GSM sequences: \n");
        for (Node node : tree.listNodesPostOrder(null, null)) {

            if (node.isLeaf()) {
                String leafSeqGSMformatWithTaxon = alinmt.sequenceInput.get().get(node.getNr()).toString();
                String[] taxonAndLeafSeq = leafSeqGSMformatWithTaxon.split(":");
                String leafSeqGSMformat = taxonAndLeafSeq[1];
                String leafSeqGABIformat = processSequence(leafSeqGSMformat);
                String leafSeqWithTaxon = taxonAndLeafSeq[0] + ":" + leafSeqGABIformat;
                Log.info.println("Taxon " + taxonAndLeafSeq[0]);
                Log.info.println("Input sequence: " + leafSeqGSMformat);
                Log.info.println("Processed sequence: " + leafSeqGABIformat);

                processedAlignmentDict.put(node.getNr(), leafSeqWithTaxon);
            }

        }
        return processedAlignmentDict;
    }

//penalization factors from Feng et al. Not used

//    protected double penalization() {
//
//        for(int i =0;  i< treeInput.get().getNodeCount(); i++) {
//            Node node = treeInput.get().getNode(i);
//            final double branchRate = branchRateModel.getRateForBranch(node);
//            final double branchTime = node.getLength() * branchRate;
//            m_branchLengths[i] = branchTime;
//
//            if (branchTime < 0.0) {
//                    throw new RuntimeException("Negative branch length: " + branchTime);
//            }
//
//        }
//
//        //Log.info.println("branch lengths"+ Arrays.toString(m_branchLengths));
//        DoubleMatrix branchLengths = new DoubleMatrix(m_branchLengths);
//        DoubleMatrix cutRates = new DoubleMatrix(Arrays.asList(substitutionModel.cutRates));
//        DoubleMatrix logbranchLenghts = logi(branchLengths);
//        DoubleMatrix logcutRates = logi(cutRates);
//        Double branchLengPen = substitutionModel.branchLensPenalty*(powi(logbranchLenghts.subi(logbranchLenghts.mean()),2)).sum();
//        Double cutRatesPen = substitutionModel.cutRatesPenalty*(powi(logcutRates.subi(logcutRates.mean()),2)).sum();
//
//        return branchLengPen + cutRatesPen;
//
//    }

    @Override
    public void store() {
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }


    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    public boolean isSomeNodeFilthy() {
        Node[] nodes = treeInput.get().getNodesAsArray();
        boolean isFilthy = false;
        for (Node node : nodes) {
            if (node.isDirty() == Tree.IS_FILTHY) {
                isFilthy = true;
            }
        }

        return isFilthy;
    }


}
