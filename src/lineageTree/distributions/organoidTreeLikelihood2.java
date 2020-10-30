package lineageTree.distributions;

import beast.core.Input;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.*;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.GeneralScarringLoss;

import java.util.Arrays;

public class organoidTreeLikelihood2 extends GenericTreeLikelihood {

    final public Input<Frequencies> rootFrequenciesInput =
            new Input<>("rootFrequencies", "prior state frequencies at root, optional", Input.Validate.OPTIONAL);

    /**
     * calculation engine *
     */
    protected LikelihoodCore likelihoodCore;
    public LikelihoodCore getLikelihoodCore() {return likelihoodCore;}

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected GeneralScarringLoss substitutionModel;
    protected double scarringStart;
    protected double scarringStop;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials and helper node partials *
     */
    protected double[] m_fRootPartials;
    protected double[][] helperNodePartials;
    protected double[] nodePartials;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    protected int matrixSize;
    /**
     * flag to indicate ascertainment correction should be applied *
     */
    protected boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */

    // number of categories across site -> number of transition rate matrices
    protected int nrOfMatrices;
    // in the alignment matrix a pattern is a unique column of alignment entries; columns == sites with the same pattern
    // have the same likelihood. thus we save computation by not recomputing the site likelihood for identical patterns
    protected int nrOfPatterns;
    // number of states in substitution rate matrix
    protected int nrOfStates;



    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        //TODO check type conversion
        substitutionModel = (GeneralScarringLoss) m_siteModel.substModelInput.get();
        scarringStart = substitutionModel.getScarringHeight();
        scarringStop = scarringStart - substitutionModel.getScarringDuration();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        nrOfStates = dataInput.get().getMaxStateCount();
        nrOfPatterns = dataInput.get().getPatternCount();
        nrOfMatrices = m_siteModel.getCategoryCount();

        if (nrOfStates == 4) {
            likelihoodCore = new BeerLikelihoodCore4();
        } else {
            likelihoodCore = new BeerLikelihoodCore(nrOfStates);
        }

        Node rootNode = treeInput.get().getRoot();

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        initCore();

        patternLogLikelihoods = new double[nrOfPatterns];
        m_fRootPartials = new double[nrOfPatterns * nrOfStates];
        helperNodePartials = new double[4][nrOfPatterns * nrOfStates];
        nodePartials = new double[nrOfPatterns * nrOfStates];
        matrixSize = (nrOfStates + 1) * (nrOfStates + 1);
        probabilities = new double[(nrOfStates + 1) * (nrOfStates + 1)];
        Arrays.fill(probabilities, 1.0);



    }

    protected void initCore() {
        // currently do not handle ambiguities
        boolean m_useAmbiguities = false;
        boolean m_useTipLikelihoods = false;

        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities || m_useTipLikelihoods) {
            //setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
            //placeholder
            int m = 0;
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /* Assumes there IS a branch rate model as opposed to traverse()*/
    protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    // check whether the branch crosses windows for different rate matriced
                    // if yes then use different rate matrices to calculate the Partials
                    // if branch crosses scarring window
                    boolean parentBeforeChildrenAfterScarringHeight = (node.getHeight() >= scarringStart ) &
                            ((child1.getHeight() < scarringStart) | (child2.getHeight() < scarringStart));
                    boolean parentBeforeChildrenAfterScarringStop = (node.getHeight() >= scarringStop) &
                            ((child1.getHeight() < scarringStop)
                                    | (child2.getHeight() < scarringStop));
                    if (parentBeforeChildrenAfterScarringHeight | parentBeforeChildrenAfterScarringStop) {

                        //calculate partials
                        nodePartials = calculatePartialsForCrossBranches(nodePartials, node, child1, child2, parentBeforeChildrenAfterScarringHeight, parentBeforeChildrenAfterScarringStop);

                        //set partials at node
                        likelihoodCore.setNodePartials(nodeIndex, nodePartials);
                    }
                    // else use default calculation engine
                    }else {

                        likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                    }
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                   /* if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : constantPattern) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }*/

                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }
                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }

            }
        return update;
        }

    double[] calculatePartialsForCrossBranches(double[] nodePartials, Node parent, Node child1, Node child2, boolean bool1, boolean bool2) {

        int child1Index = child1.getNr();
        int child2Index = child2.getNr();

        final double branchRate1 = branchRateModel.getRateForBranch(child1);
        final double branchRate2 = branchRateModel.getRateForBranch(child2);

        // TODO remove all this varying rates across sites stuff!
        double jointBranchRate1 = m_siteModel.getRateForCategory(0, child1) * branchRate1;
        double jointBranchRate2 = m_siteModel.getRateForCategory(0, child2) * branchRate1;

        double[] probs1 = new double[(nrOfStates +1) * (nrOfStates + 1)];
        double[] probs2 = new double[(nrOfStates + 1) * (nrOfStates +1)];


        // if parent is above scarring window
        if (bool1){
            if (child1.getHeight() < scarringStop){

                if ( child2.getHeight() < scarringStop){

                    if (child1.isLeaf()){

                        int[] states1 = new int[nrOfPatterns];
                        likelihoodCore.getNodeStates(child1Index, states1);

                        //calculate transition probability for all intermediate regions until parent
                        // from leaf to first helper node
                        substitutionModel.getTransitionProbabilities(null, scarringStop, child1.getHeight(), jointBranchRate1, probs1);
                        helperNodePartials[0] = calculateStatesPruning(states1, probs1, helperNodePartials[0]);

                        // from first to second helper node
                        substitutionModel.getTransitionProbabilities(null, scarringStart, scarringStop, jointBranchRate1, probs1);
                        helperNodePartials[1] =  calculatePartialsPruning(helperNodePartials[0], probs1, helperNodePartials[1]);

                        substitutionModel.getTransitionProbabilities(null, parent.getHeight(), scarringStart, jointBranchRate1, probs1);

                        if(child2.isLeaf()){
                            int[] states2 = new int[nrOfPatterns];
                            likelihoodCore.getNodeStates(child2Index, states2);

                            //calculate transition probability for all intermediate regions until parent
                            //from leaf to first helper node
                            substitutionModel.getTransitionProbabilities(null, scarringStop, child2.getHeight(), jointBranchRate2, probs2);
                            helperNodePartials[2] = calculateStatesPruning(states2, probs2, helperNodePartials[2]);

                            // from first to second helper node
                            substitutionModel.getTransitionProbabilities(null, scarringStart, scarringStop, jointBranchRate2, probs2);
                            helperNodePartials[3] = calculatePartialsPruning(helperNodePartials[2], probs2, helperNodePartials[3]);

                            //from second helper node to parent
                            substitutionModel.getTransitionProbabilities(null, parent.getHeight(), scarringStart, jointBranchRate2, probs2);
                            nodePartials = calculatePartialsPartialsPruning(helperNodePartials[1], probs1, helperNodePartials[3], probs2, nodePartials);

                        }else{
                            double [] partials2 = new double[nrOfPatterns];
                            likelihoodCore.getNodePartials(child2Index, partials2);
                            helperNodePartials[0] = calculateStatesPartialsPruning(states1, probs1, partials2, probs2, helperNodePartials[0]);
                        }
                    }else{
                        double [] partials1 = new double[nrOfPatterns];
                        likelihoodCore.getNodePartials(child1Index, partials1);
                        if(child2.isLeaf()){
                            int[] states2 = new int[nrOfPatterns];
                            likelihoodCore.getNodeStates(child2Index, states2);
                            helperNodePartials[0] = calculateStatesPartialsPruning(states2, probs2, partials1, probs1, helperNodePartials[0]);
                        }else{
                            double [] partials2 = new double[nrOfPatterns];
                            likelihoodCore.getNodePartials(child2Index, partials2);
                            helperNodePartials[0] = calculatePartialsPartialsPruning(partials1, probs1, partials2, probs2, helperNodePartials[0]);
                        }
                    }

                    }
                }


                // call special function for child 1

                // call beerlikelihood function for child 2
            }else if(child2.getHeight() < substitutionModel.getScarringHeight()){
                // call special function for child 2

                // call beerlikelihood function for child 1
            }

        return nodePartials;
    }



        /**
         * Calculates partial likelihoods at a node when both children have states.
         * From BeerLikelihoodCore
         */
    protected double [] calculateStatesStatesPruning(int[] stateIndex1, double[] matrices1,
                                                int[] stateIndex2, double[] matrices2,
                                                double[] partials3) {
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int state1 = stateIndex1[k];
                int state2 = stateIndex2[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates && state2 < nrOfStates) {

                    for (int i = 0; i < nrOfStates; i++) {

                        partials3[v] = matrices1[w + state1] * matrices2[w + state2];

                        v++;
                        w += nrOfStates;
                    }

                } else if (state1 < nrOfStates) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < nrOfStates; i++) {

                        partials3[v] = matrices1[w + state1];

                        v++;
                        w += nrOfStates;
                    }
                } else if (state2 < nrOfStates) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < nrOfStates; i++) {

                        partials3[v] = matrices2[w + state2];

                        v++;
                        w += nrOfStates;
                    }
                } else {
                    // both children have a gap or unknown state so set partials to 1

                    for (int j = 0; j < nrOfStates; j++) {
                        partials3[v] = 1.0;
                        v++;
                    }
                }
            }
        }
        return partials3;
    }

    protected double [] calculateStatesPruning(int[] stateIndex1, double[] matrices1,
                                                     double[] partials3) {
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int state1 = stateIndex1[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates) {

                    for (int i = 0; i < nrOfStates; i++) {

                        partials3[v] = matrices1[w + state1];

                        v++;
                        w += nrOfStates;
                    }

                } else {
                    // single child has a gap or unknown state so set partials to 1
                    for (int j = 0; j < nrOfStates; j++) {
                        partials3[v] = 1.0;
                        v++;
                    }
                }
            }
        }
        return partials3;
    }


    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     * From BeerLikelihoodCore
     */
    protected double[] calculateStatesPartialsPruning(int[] stateIndex1, double[] matrices1,
                                                  double[] partials2, double[] matrices2,
                                                  double[] partials3) {

        double sum, tmp;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {

                int state1 = stateIndex1[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates) {


                    for (int i = 0; i < nrOfStates; i++) {

                        tmp = matrices1[w + state1];

                        sum = 0.0;
                        for (int j = 0; j < nrOfStates; j++) {
                            sum += matrices2[w] * partials2[v + j];
                            w++;
                        }

                        partials3[u] = tmp * sum;
                        u++;
                    }

                    v += nrOfStates;
                } else {
                    // Child 1 has a gap or unknown state so don't use it

                    for (int i = 0; i < nrOfStates; i++) {

                        sum = 0.0;
                        for (int j = 0; j < nrOfStates; j++) {
                            sum += matrices2[w] * partials2[v + j];
                            w++;
                        }

                        partials3[u] = sum;
                        u++;
                    }

                    v += nrOfStates;
                }
            }
        }
        return partials3;
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected double [] calculatePartialsPartialsPruning(double[] partials1, double[] matrices1,
                                                    double[] partials2, double[] matrices2,
                                                    double[] partials3) {
        double sum1, sum2;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < nrOfStates; i++) {

                    sum1 = sum2 = 0.0;

                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[w] * partials1[v + j];
                        sum2 += matrices2[w] * partials2[v + j];

                        w++;
                    }

                    partials3[u] = sum1 * sum2;
                    u++;
                }
                v += nrOfStates;
            }
        }
        return partials3;
    }

    // calculate partials for single child nodes
    protected double [] calculatePartialsPruning(double[] partials1, double[] matrices1,
                                                         double[] partials3) {
        double sum1;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < nrOfStates; i++) {

                    sum1 =  0.0;

                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += matrices1[w] * partials1[v + j];

                        w++;
                    }

                    partials3[u] = sum1;
                    u++;
                }
                v += nrOfStates;
            }
        }
        return partials3;
    }


    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
    }
    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {

        final TreeInterface tree = treeInput.get();

        try {
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
                calcLogP();
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
            //TODO deleted possibility to turn off scaling
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 ) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }
    }





}