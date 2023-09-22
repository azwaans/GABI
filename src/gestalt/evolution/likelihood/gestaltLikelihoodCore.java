package gestalt.evolution.likelihood;


import beast.base.core.Log;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.tree.Node;
import gestalt.evolution.alignment.TargetStatus;
import gestalt.evolution.alignment.TransitionWrap;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.jblas.DoubleMatrix;

import java.util.Hashtable;

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

    protected Hashtable<Integer, DoubleMatrix> ptMat;
    protected Hashtable<Integer, DoubleMatrix> partials;
    protected Hashtable<Integer, Double> logScalingTerms;


    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;

    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;


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

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

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

        partials = null;
        ptMat = null;

        logScalingTerms = null;
    }

    /**
     * Sets partials for a node
     */
    @Override
    public void setNodePartials(int nodeIndex, double[] partials) {

    }

    /**
     * set leaf partials in likelihood core *
     */
    protected void setLeafPartials(Node node) {


        TransitionWrap nodeWrap = transitionWraps.get(node.getNr());
        DoubleMatrix leafPartials = zeros(nodeWrap.numStatuses + 1, 1);
        Integer observedStateKey = nodeWrap.statusMap.get(nodeWrap.leafState);
        leafPartials.put(observedStateKey, 1.0);
        setNodePartials(node.getNr(), leafPartials);
    }


    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
    }


    public void setNodePartials(int nodeIndex, DoubleMatrix partial) {
        partials.put(nodeIndex + currentPartialsIndex[nodeIndex] * nodeIndex, partial);
    }


    public void setNodeMatrix(int nodeIndex, DoubleMatrix Mat) {
        ptMat.put(nodeIndex + currentPartialsIndex[nodeIndex] * nodeIndex, Mat);
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
            TransitionWrap childNodeWrap = transitionWraps.get(childNum);
            DoubleMatrix Mat = ptMat.get(childNum + currentMatrixIndex[childNum] * childNum);
            //Get the probability for the data descended from the child node, assuming that the node
            //has a particular target tract repr.
            //These down probs are ordered according to the child node's numbering of the TTs states
            DoubleMatrix chOrderedDownProbs = Mat.mmul(partials.get(childNum + currentPartialsIndex[childNum] * childNum));
            DoubleMatrix downProbs = new DoubleMatrix();

            if (!node.isRoot()) {

                // Reorder summands according to node's numbering of tract_repr states
                downProbs = gestaltGeneral.reorderLikelihoods(chOrderedDownProbs, transitionWraps.get(node.getNr()), childNodeWrap);

            }

            if (node.isRoot()) {
                //For the root node, we just want the probability where the root node is unmodified
                //No need to reorder

                // the index is null:
                Integer childUnmodifiedIndex = childNodeWrap.statusMap.get(new TargetStatus());
                //Log.info.println(chOrderedDownProbs);
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
        logScalingTerms.put(node.getNr(), logScalingTerm);
        return expi((nodeLogPartials.subi(logScalingTerm)).muli(hasPosProb));

    }


    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        int[] tmp1 = currentMatrixIndex;
        currentMatrixIndex = storedMatrixIndex;
        storedMatrixIndex = tmp1;

        int[] tmp2 = currentPartialsIndex;
        currentPartialsIndex = storedPartialsIndex;
        storedPartialsIndex = tmp2;

    }

    @Override
    public void unstore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);

    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
    }

    public DoubleMatrix getNodePartials(int nodeIndex) {
        return partials.get(nodeIndex + currentPartialsIndex[nodeIndex] * nodeIndex);
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
