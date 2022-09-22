package lineageTree.distributions;


import beast.evolution.alignment.TargetStatus;
import beast.evolution.alignment.TransitionWrap;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.tree.Node;
import lineageTree.substitutionmodel.GeneralGestalt;
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
public class GapmlLikelihoodCore extends LikelihoodCore {
    protected int nrOfStates;
    protected int nrOfNodes;

    protected Hashtable<Integer, TransitionWrap> transitionWrappers;

    protected Hashtable<Integer, DoubleMatrix> ptMat;
    protected Hashtable<Integer, DoubleMatrix> partials;
    protected Hashtable<Integer, Double> logScalingTerms ;


    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;

    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;



    public GapmlLikelihoodCore() {
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
     * @param nodeCount           the number of nodes in the tree
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
    //CHANGE BLOCK
    protected void setLeafPartials(Node node) {
        //this is for the purpose of working with single branch trees with potentially no left/right
//        if (node != null) {


//            if (node.isLeaf()) {
//
        TransitionWrap nodeWrap = transitionWrappers.get(node.getNr());
        DoubleMatrix leafPartials = zeros(nodeWrap.numStatuses + 1, 1);
        Integer observedStateKey = nodeWrap.statusMap.get(nodeWrap.leafState);
        leafPartials.put(observedStateKey, 1.0);
        setNodePartials(node.getNr(), leafPartials);


//            } else {
//                setPartials(node.getLeft());
//                setPartials(node.getRight());
//
//            }

    }
    //CHANGE block



    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];
    }



    public void setNodePartials(int nodeIndex, DoubleMatrix partial) {
        partials.put(nodeIndex + currentPartialsIndex[nodeIndex]*nodeIndex,partial);
    }


    public void setNodeMatrix(int nodeIndex, DoubleMatrix Mat) {
        ptMat.put(nodeIndex + currentPartialsIndex[nodeIndex] * nodeIndex,Mat);
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
        for(Node child: node.getChildren()) {
            int childNum = child.getNr();
            TransitionWrap childNodeWrap = transitionWrappers.get(childNum);
            DoubleMatrix Mat = ptMat.get(childNum + currentMatrixIndex[childNum] * childNum);

            //ptMats.put(child.getNr(),ptMat);
            //Log.info.println("PTMATRIX"+ptMat);

            //Get the probability for the data descended from the child node, assuming that the node
            //has a particular target tract repr.
            //These down probs are ordered according to the child node's numbering of the TTs states
            DoubleMatrix chOrderedDownProbs = Mat.mmul(partials.get(childNum + currentPartialsIndex[childNum] * childNum));
            DoubleMatrix downProbs = new DoubleMatrix();

            if (!node.isRoot()) {

                // Reorder summands according to node's numbering of tract_repr states
                downProbs = GeneralGestalt.reorderLikelihoods(chOrderedDownProbs, transitionWrappers.get(node.getNr()), childNodeWrap);

            }

            if (node.isRoot()) {
                //For the root node, we just want the probability where the root node is unmodified
                //No need to reorder
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




    @Override
    public boolean getUseScaling() {
        return true;}

    @Override
    public double getLogScalingFactor(int patternIndex_) {
        return 0;
    }

    @Override
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
    }

    public void setCurrentNodePartials(int nodeIndex, double[] partials) {
    }

//    @Override
//    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
//    }

    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
//        calculateIntegratePartials(partials[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
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
        //at this stage, partial sizes are not known!
//        this.partials[0][nodeIndex] = null;
//        this.partials[1][nodeIndex] = null;
    }

    @Override
    public void getNodePartials(int nodeIndex, double[] partialsOut) {
//        System.arraycopy(partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partialsOut, 0, partialsOut.length);
    }
    public DoubleMatrix getNodePartials(int nodeIndex) {
        return partials.get(nodeIndex + currentPartialsIndex[nodeIndex] *nodeIndex);
    }

    /**
     * Gets the partials for a particular node.
     *
     * @param nodeIndex   the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials) {
//        double[] partials1 = partials[currentPartialsIndex[nodeIndex]][nodeIndex];
//
//        System.arraycopy(partials1, 0, outPartials, 0, partialsSize);
    }


} // class BeerLikelihoodCore
