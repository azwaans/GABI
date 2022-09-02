package lineageTree.distributions;


import beast.evolution.alignment.TransitionWrap;
import beast.evolution.likelihood.LikelihoodCore;
import org.jblas.DoubleMatrix;

import java.util.Hashtable;

/**
 * standard likelihood core, uses no caching *
 */
public class GapmlLikelihoodCore extends LikelihoodCore {
    protected int nrOfStates;
    protected int nrOfNodes;



    protected static Hashtable<Integer, TransitionWrap> transitionWrappers;
    protected static Hashtable<Integer, DoubleMatrix> ptMats;
    protected static Hashtable<Integer, DoubleMatrix> partls;



    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;

    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;

    protected int[] currentWrapsIndex;
    protected int[] storedWrapsIndex;


    public GapmlLikelihoodCore() {
    } // c'tor







    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    @Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {

                sum += frequencies[i] * partials[v];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }

    @Override
    protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

    }


    /**
     * initializes partial likelihood arrays.
     *
     * @param nodeCount           the number of nodes in the tree
     * @param patternCount        the number of patterns
     * @param matrixCount         the number of matrices (i.e., number of categories)
     * @param integrateCategories whether sites are being integrated over all matrices
     */
	public static void init(int nodeCount) {

        this.nrOfNodes = nodeCount;
        partls = new Hashtable<>();
        ptMats = new Hashtable<>();

        currentMatrixIndex = new int[nodeCount];
        storedMatrixIndex = new int[nodeCount];

        currentPartialsIndex = new int[nodeCount];
        storedPartialsIndex = new int[nodeCount];

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

        partls = null;
        ptMats = null;
    }



    /**
     * Allocates partials for a node
     */
    @Override
	public void createNodePartials(int nodeIndex) {
        //at this stage, partial sizes are not known!
        this.partials[0][nodeIndex] = null;
        this.partials[1][nodeIndex] = null;
    }

    /**
     * Sets partials for a node
     */
    @Override
	public void setNodePartials(int nodeIndex, double[] partials) {

        if (this.partials[0][nodeIndex] == null) {
            createNodePartials(nodeIndex);
        }
//        if (partials.length < partialsSize) {
//            int k = 0;
//            for (int i = 0; i < nrOfMatrices; i++) {
//                System.arraycopy(partials, 0, this.partials[0][nodeIndex], k, partials.length);
//                k += partials.length;
//            }
//        } else {
            System.arraycopy(partials, 0, this.partials[0][nodeIndex], 0, partials.length);
        //}
    }

    @Override
    public void getNodePartials(int nodeIndex, double[] partialsOut) {
        System.arraycopy(partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partialsOut, 0, partialsOut.length);
    }

    /**
     * Allocates states for a node
     */
    public void createNodeStates(int nodeIndex) {

        this.states[nodeIndex] = new int[nrOfPatterns];
    }

    /**
     * Sets states for a node
     */
    @Override
	public void setNodeStates(int nodeIndex, int[] states) {

        if (this.states[nodeIndex] == null) {
            createNodeStates(nodeIndex);
        }
        System.arraycopy(states, 0, this.states[nodeIndex], 0, nrOfPatterns);
    }

    /**
     * Gets states for a node
     */
    @Override
	public void getNodeStates(int nodeIndex, int[] states) {
        System.arraycopy(this.states[nodeIndex], 0, states, 0, nrOfPatterns);
    }

    @Override
    public void setNodeMatrixForUpdate(int nodeIndex) {
        currentMatrixIndex[nodeIndex] = 1 - currentMatrixIndex[nodeIndex];

    }


//    public void setNodeWrapForUpdate(int nodeIndex) {
//        currentWrapIndex[nodeIndex] = 1 - currentWrapIndex[nodeIndex];
//
//    }


    /**
     * Sets probability matrix for a node
     */
	public static void setNodeMatrix(int nodeIndex, DoubleMatrix matrix) {
       ptMats.put( nodeIndex + nodeIndex*currentMatrixIndex[nodeIndex],matrix);
    }

    /**
     * Sets probability matrix for a node
     */
    @Override
    public void setNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
                matrixIndex * matrixSize, matrixSize);
    }

//    /**
//     * Sets Transition wrap for a node
//     */
//    public static void setNodeWrap(int nodeIndex, TransitionWrap wrap) {
//        transitionWrappers.put(nodeIndex, wrap);
//    }


//    /**
//     * Sets Transition wrap for a node
//     */
//    public static void getNodeWrap(int nodeIndex, TransitionWrap wrap) {
//        wrap = transitionWrappers.get(nodeIndex);
//    }

    public void setPaddedNodeMatrices(int nodeIndex, double[] matrix) {
        System.arraycopy(matrix, 0, matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
                0, nrOfMatrices * matrixSize);
    }


    /**
     * Gets probability matrix for a node
     * no usage!!! no need to change
     */
    @Override
	public void getNodeMatrix(int nodeIndex, int matrixIndex, double[] matrix) {
        System.arraycopy(matrices[currentMatrixIndex[nodeIndex]][nodeIndex],
                matrixIndex * matrixSize, matrix, 0, matrixSize);
    }

    @Override
    public void setNodePartialsForUpdate(int nodeIndex) {
        currentPartialsIndex[nodeIndex] = 1 - currentPartialsIndex[nodeIndex];
    }

    /**
     * Sets the currently updating node partials for node nodeIndex. This may
     * need to repeatedly copy the partials for the different category partitions
     */
    public void setCurrentNodePartials(int nodeIndex, double[] partials) {
        if (partials.length < partialsSize) {
            int k = 0;
            for (int i = 0; i < nrOfMatrices; i++) {
                System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], k, partials.length);
                k += partials.length;
            }
        } else {
            System.arraycopy(partials, 0, this.partials[currentPartialsIndex[nodeIndex]][nodeIndex], 0, partials.length);
        }
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
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                        states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            } else {
                calculateStatesPartialsPruning(states[nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            }
        } else {
            if (states[nodeIndex2] != null) {
                calculateStatesPartialsPruning(states[nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            } else {
                calculatePartialsPartialsPruning(partials[currentPartialsIndex[nodeIndex1]][nodeIndex1], matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                        partials[currentPartialsIndex[nodeIndex2]][nodeIndex2], matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3]);
            }
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }

//
//        int k =0;
//        for (int i = 0; i < patternCount; i++) {
//            double f = 0.0;
//
//            for (int j = 0; j < stateCount; j++) {
//                f += partials[currentPartialsIndices[nodeIndex3]][nodeIndex3][k];
//                k++;
//            }
//            if (f == 0.0) {
//                Logger.getLogger("error").severe("A partial likelihood (node index = " + nodeIndex3 + ", pattern = "+ i +") is zero for all states.");
//            }
//        }
    }

    @Override
    public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {

    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param matrixMap  a map of which matrix to use for each pattern (can be null if integrating over categories)
     */
    public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int[] matrixMap) {

    }


//    @Override
//	public void integratePartials(int nodeIndex, double[] proportions, double[] outPartials) {
//        calculateIntegratePartials(partials[currentPartialsIndex[nodeIndex]][nodeIndex], proportions, outPartials);
//    }




    /**
     * Gets the partials for a particular node.
     *
     * @param nodeIndex   the node
     * @param outPartials an array into which the partials will go
     */
    public void getPartials(int nodeIndex, double[] outPartials) {
        double[] partials1 = partials[currentPartialsIndex[nodeIndex]][nodeIndex];

        System.arraycopy(partials1, 0, outPartials, 0, partialsSize);
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

        int[] tmp3 = currentWrapsIndex;
        currentWrapsIndex = storedWrapsIndex;
        storedWrapsIndex = tmp3;
    }

    @Override
	public void unstore() {
        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
        System.arraycopy(storedWrapsIndex, 0, currentWrapsIndex, 0, nrOfNodes);

    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
        System.arraycopy(currentMatrixIndex, 0, storedMatrixIndex, 0, nrOfNodes);
        System.arraycopy(currentPartialsIndex, 0, storedPartialsIndex, 0, nrOfNodes);
        System.arraycopy(currentWrapsIndex, 0, storedWrapsIndex, 0, nrOfNodes);
    }



    
    @Override
    public boolean getUseScaling() {
        return useScaling;
    }

} // class BeerLikelihoodCore
