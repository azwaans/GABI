package gestalt.evolution.likelihood;


import beast.base.core.Log;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.tree.Node;
import gestalt.evolution.alignment.*;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.jblas.DoubleMatrix;

import java.util.*;

import static gestalt.evolution.alignment.TransitionWrap.getCloseTransitionWrap;
import static java.lang.Math.max;
import static org.jblas.DoubleMatrix.ones;
import static org.jblas.DoubleMatrix.zeros;
import static org.jblas.MatrixFunctions.expi;
import static org.jblas.MatrixFunctions.logi;

/**
 * standard likelihood core, uses no caching *
 */
public class gestaltLikelihoodCore extends LikelihoodCore {
    protected int nrOfNodes;

    protected Hashtable<Integer, TransitionWrap> transitionWraps;
    protected Hashtable<Integer, TransitionWrap> storedtransitionWraps;

    protected Hashtable<Integer, DoubleMatrix> ptMat;
    protected Hashtable<Integer, DoubleMatrix> partials;
    protected Hashtable<Integer, Double> logScalingTerms;

    protected int currentWrapIndex;
    protected int storedWrapIndex;

    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;

    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected int[] currentScalingTermsIndex;
    protected int[] storedScalingTermsIndex;




    public gestaltLikelihoodCore() {
    } // c'tor

    public void initCore(int nodeCount) {

        nrOfNodes = nodeCount;

        partials = new Hashtable<>();
        ptMat = new Hashtable<>();
        logScalingTerms = new Hashtable<>();

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        currentScalingTermsIndex = new int[nodeCount];
        storedScalingTermsIndex = new int[nodeCount];


    }


    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    @Override
    public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
//        int v = 0;
//        for (int k = 0; k < nrOfPatterns; k++) {
//
//            double sum = 0.0;
//            for (int i = 0; i < nrOfStates; i++) {
//
//                sum += frequencies[i] * partials[v];
//                v++;
//            }
//            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
//        }
    }

    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

    }


    /**
     * initializes partial likelihood hashmaps.
     *
     * @param nodeCount the number of nodes in the tree
     */
    public void init(int nodeCount) {

        nrOfNodes = nodeCount;

        this.partials = new Hashtable<>();
        this.ptMat = new Hashtable<>();
        this.logScalingTerms = new Hashtable<>();
        this.transitionWraps = null;

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

        currentScalingTermsIndex = new int[nodeCount];
        storedScalingTermsIndex = new int[nodeCount];

        currentWrapIndex = 0;
        storedWrapIndex = 0;

    }

    @Override
    public void initialize(int nodeCount, int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities) {

    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
    public void finalize() throws Throwable {
        nrOfNodes = 0;

        currentPartialsIndex = null;
        storedPartialsIndex = null;

        currentMatrixIndex = null;
        storedMatrixIndex = null;

        currentWrapIndex = 0;
        storedWrapIndex = 0;

        currentScalingTermsIndex = null;
        storedScalingTermsIndex = null;

        partials = null;
        ptMat = null;

        logScalingTerms = null;
    }

    /**
     * Sets partials for a node
     */
    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {

//        if (this.partials[0][nodeIndex] == null) {
//            createNodePartials(nodeIndex);
//        }
//        if (partials.length < partialsSize) {
//            int k = 0;
//            for (int i = 0; i < nrOfMatrices; i++) {
//                System.arraycopy(partials, 0, this.partials[0][nodeIndex], k, partials.length);
//                k += partials.length;
//            }
//        } else {
        //       System.arraycopy(partials, 0, this.partials[0][nodeIndex], 0, partials.length);
        //}
    }

    /**
     * set leaf partials in likelihood core *
     */
    protected void setLeafPartials(Node node) {
        TransitionWrap nodeWrap = transitionWraps.get(node.getNr() +1);
        DoubleMatrix leafPartials = zeros(nodeWrap.numStatuses + 1, 1);
        Integer observedStateKey = nodeWrap.statusMap.get(nodeWrap.leafState);
        leafPartials.put(observedStateKey, 1.0);
        partials.put(makeCachingIndexPartials(node.getNr()) , leafPartials);
    }

    public void setNodePartials(int nodeIndex, DoubleMatrix partial) {
        partials.put(makeCachingIndexPartials(nodeIndex) , partial);
    }

    public void setScalingTerm(int nodeIndex, Double ScalingTerm) {
        logScalingTerms.put(makeCachingIndexScalingTerm(nodeIndex) , ScalingTerm);
    }
    public void setNodeMatrix(int nodeIndex, DoubleMatrix Mat) {
        ptMat.put(makeCachingIndexMatrix(nodeIndex), Mat);
    }

    //todo find a more elegant way of unambiguously caching all values
    int makeCachingIndexMatrix(int nodeIndex) {
        int node = nodeIndex + 1;
        String forHashing = node + "" +  currentMatrixIndex[nodeIndex] + ""+ node;
        return forHashing.hashCode();

    }

    int makeCachingIndexScalingTerm(int nodeIndex) {
        int node = nodeIndex + 1;
        String forHashing = node + "" +  currentScalingTermsIndex[nodeIndex] + ""+ node;
        return forHashing.hashCode();

    }

    int makeCachingIndexPartials(int nodeIndex) {
        int node = nodeIndex + 1;
        String forHashing = node + "" +  currentPartialsIndex[nodeIndex] + ""+ node;
        return forHashing.hashCode();

    }


    /**
     * Gets probability matrix for a node
     */
    @Override
    public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
    }

    @Override
    public void setUseScaling(double scale) {

    }

    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
    }

    public void setLogScalingTermsForUpdate(int nodeIndex) {
        currentScalingTermsIndex[nodeIndex] = 1 - currentScalingTermsIndex[nodeIndex];
    }


    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    @Override
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {


    }


    public DoubleMatrix calculatePartials(Node node) {


        DoubleMatrix nodeLogPartials = DoubleMatrix.zeros(1);
        DoubleMatrix hasPosProb = ones(1);

        //here, combine child 1 + child 2 partial
        for (Node child : node.getChildren()) {

            int childNum = child.getNr();
            TransitionWrap childNodeWrap = transitionWraps.get((childNum + 1));
            DoubleMatrix Mat = ptMat.get(makeCachingIndexMatrix(childNum));

            //Get the probability for the data descended from the child node, assuming that the node
            //has a particular target tract repr.
            //These down probs are ordered according to the child node's numbering of the TTs states
            DoubleMatrix chOrderedDownProbs = Mat.mmul(partials.get(makeCachingIndexPartials(childNum)));
            DoubleMatrix downProbs = new DoubleMatrix();

            if (!node.isRoot()) {

                // Reorder summands according to node's numbering of tract_repr states
                downProbs = gestaltGeneral.reorderLikelihoods(chOrderedDownProbs, transitionWraps.get(node.getNr() + 1), childNodeWrap);

            }

            if (node.isRoot()) {
                //For the root node, we just want the probability where the root node is unmodified
                //No need to reorder

                // the index is null:
                Integer childUnmodifiedIndex = childNodeWrap.statusMap.get(new TargetStatus());
                double downProb = chOrderedDownProbs.get(childUnmodifiedIndex);
                downProb = max(downProb, 0.0);
                downProbs = new DoubleMatrix(1, 1, downProb);
            }
            //todo: implement node abundances
            //if (child.isLeaf()) {
            //Double leafAbundanceWeight = 1.0;
            hasPosProb = downProbs.ge(0);

            //protection against states with zero probability.
            nodeLogPartials = nodeLogPartials.addi(logi(downProbs.addi((hasPosProb.neg()).add(1).mul(1e-30))));


        }
        Double logScalingTerm = nodeLogPartials.max();
        setLogScalingTermsForUpdate(node.getNr());
        setScalingTerm(node.getNr(), logScalingTerm);

        return expi((nodeLogPartials.subi(logScalingTerm)).muli(hasPosProb));

    }

    /**
     * Creates a transition wrap for all nodes in the tree
     */

    public Hashtable<Integer, TransitionWrap> createTransitionWraps(beast.base.evolution.tree.TreeInterface tree,
                                      BarcodeMeta metaData,
                                      Hashtable<Integer,
                                              AncStates> statesDict) {

        Hashtable<Integer, TransitionWrap> wrap = new Hashtable<>();

        for (Node node : tree.listNodesPostOrder(null, null)) {
            //creating an empty transition wrap, avoiding to use NULL as we need empty entries( not null)
            List<List<IndelSet.TargetTract>> initNull = new ArrayList<>();
            List<IndelSet.TargetTract> innerinitNull = new ArrayList<>();
            initNull.add(innerinitNull);
            TransitionWrap temp = new TransitionWrap(initNull, statesDict.get(node.getNr()), node.isLeaf());
            wrap.put(node.getNr() + 1, temp);
        }
        //dictionary of all singleton states (not node assigned)
        Hashtable<Integer, List<IndelSet.Singleton>> parsimDict = new Hashtable<>();
        for (Integer key : statesDict.keySet()) {
            parsimDict.put(key, statesDict.get(key).getSingletons());
        }

        //traverse the tree to fill up the dictionary of wraps
        List<Node> preorderList = Arrays.asList(tree.listNodesPostOrder(null, null));
        for (int reverseIt = preorderList.size() - 1; reverseIt >= 0; reverseIt--) {
            Node parentNode = preorderList.get(reverseIt);

            List<List<IndelSet.TargetTract>> parentNodeTuples = wrap.get((parentNode.getNr() + 1)).targetTractsTuples;
            List<List<IndelSet.TargetTract>> filteredtargetTractsTupless = parentNodeTuples;

            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap filterWrap = getCloseTransitionWrap(parentNode, childNode, statesDict, parsimDict, parentNodeTuples, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);
                filteredtargetTractsTupless = IndelSet.TargetTract.intersect(filteredtargetTractsTupless, filterWrap.targetTractsTuples);
                //removing duplicates
                Set noDup = new LinkedHashSet();
                noDup.addAll(filteredtargetTractsTupless);
                filteredtargetTractsTupless.clear();
                filteredtargetTractsTupless.addAll(noDup);

            }

            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap finalWrap = getCloseTransitionWrap(parentNode, childNode, statesDict, parsimDict, filteredtargetTractsTupless, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);
                List<TargetStatus> cleanTTUPLES = finalWrap.transStatuses;
                Collections.reverse(cleanTTUPLES);
                finalWrap.transStatuses = cleanTTUPLES;
                wrap.put((childNode.getNr() + 1) , finalWrap);

            }
        }
        return wrap;

    }


    /**
     * Store current state
     */
    @Override
    public void restore() {
//        // Rather than copying the stored stuff back, just swap the pointers...
        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int[] tmp2 = currentPartialsIndex;
        currentPartialsIndex = storedPartialsIndex;
        storedPartialsIndex = tmp2;

        int[] tmp3 = currentScalingTermsIndex;
        currentScalingTermsIndex = storedScalingTermsIndex;
        storedScalingTermsIndex = tmp3;

        transitionWraps = new Hashtable<>(storedtransitionWraps);

    }

    @Override
    public void unstore() {
        Log.info.println("UNSTORE");
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
        System.arraycopy(storedScalingTermsIndex, 0, currentScalingTermsIndex, 0, nrOfNodes);

    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
        System.arraycopy(currentScalingTermsIndex, 0, storedScalingTermsIndex, 0, nrOfNodes);
        storedtransitionWraps = new Hashtable<>(transitionWraps);
    }

    public DoubleMatrix getNodePartials(int nodeIndex) {
        return partials.get(makeCachingIndexPartials(nodeIndex));
    }

    @Override
    public boolean getUseScaling() {
        return true;
    }

    @Override
    public double getLogScalingFactor(int patternIndex_) {
        return 0;
    }

    @Override
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
    }


    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
    }

    /**
     * Allocates states for a node
     */
    public void createNodeStates(int nodeIndex) {
    }

    /**
     * Sets states for a node
     */
    @Override
    public void setNodeStates(int nodeIndex, int[] states) {
    }

    @Override
    public void getNodeStates(int nodeIndex, int[] states) {

    }

    /**
     * Allocates partials for a node
     */
    @Override
    public void createNodePartials(int nodeIndex) {

    }

    @Override
    public void getNodePartials(int nodeIndex, double[] partialsOut) {
    }


} // class gestaltCore
