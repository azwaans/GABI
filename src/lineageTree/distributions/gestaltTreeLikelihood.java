package lineageTree.distributions;

import beast.evolution.likelihood.LikelihoodCore;
import org.jblas.DoubleMatrix;

import java.util.*;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.GeneralGestalt;

import static org.jblas.DoubleMatrix.*;
import static org.jblas.MatrixFunctions.*;

import static java.lang.Math.max;

import static com.google.common.math.DoubleMath.mean;


@Description("Generic tree likelihood for an alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model")
// Use this as base class to define any non-standard TreeLikelihood.
// Override Distribution.calculatLogP() to make this class functional.
//
// TODO: This could contain a generic traverse() method that takes dirty trees in account.
//
public class gestaltTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    protected GeneralGestalt substitutionModel;

    protected SiteModel.Base m_siteModel;

    protected BranchRateModel.Base branchRateModel;

    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;
    protected int[] partialsSizes;

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;

//    /**
//     * memory allocation for partials *
//     */
//    protected double[][][] Partials;

    protected Hashtable<Integer, DoubleMatrix> partials0;
//    protected Hashtable<Integer, DoubleMatrix> partials1;

//    protected Hashtable<Integer, Double> logScalingTerms ;

    protected Hashtable<Integer, DoubleMatrix> ptMat0;
    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    protected Hashtable<Integer, TransitionWrap> transitionWrappers;

    protected LikelihoodCore likelihoodCore;


    @Override
    public List<String> getArguments() {return null;}

    @Override
    public List<String> getConditions() {return null;}

    @Override
    public void sample(State state, Random random) {}

    public void initAndValidate() {
        Log.info.println("initializing tree likelihood");
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


        substitutionModel = (GeneralGestalt) m_siteModel.substModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];


        Alignment alignment = dataInput.get();

        //CHANGE BLOCK
//        likelihoodCore = new GapmlLikelihoodCore();
//        initCore();
        //CHANGE BLOCK

         partials0 = new Hashtable<>();
//         partials1 = new Hashtable<>();
//         logScalingTerms = new Hashtable<>();



    }

    //CHANGE BLOCK
//    protected void initCore() {
//        final int nodeCount = treeInput.get().getNodeCount();
//        likelihoodCore.initialize(
//                nodeCount,
//                dataInput.get().getPatternCount(),
//                m_siteModel.getCategoryCount(),
//                true, m_useAmbiguities.get()
//        );
//        //STEP 1
//        //ANTOINE PARTIALS is initialized as: partials = new double[2][nodeCount][];
//        final int extNodeCount = nodeCount / 2 + 1;
//        final int intNodeCount = nodeCount / 2;
//
//        /**
//         * set leaf partials in likelihood core *
//         */
//
//        setPartials(treeInput.get().getRoot());
//        hasDirt = Tree.IS_FILTHY;
//        for (int i = 0; i < intNodeCount; i++) {
//            likelihoodCore.createNodePartials(extNodeCount + i);
//            // ANTOINE: this is only to create an empty array of size PartialSize:
//            // the size is : partialsSize = patternCount * nrOfStates;
//            // in partials : [][][] -> [0/1][NodeIndex][double[partialsSize]]
//        }
//    }

    /**
     * set leaf partials in likelihood core *
     */
    //CHANGE BLOCK
    protected void setPartials(Node node) {
        //this is for the purpose of working with single branch trees with potentially no left/right
        if (node != null) {


            if (node.isLeaf()) {

                TransitionWrap nodeWrap = transitionWrappers.get(node.getNr());
                DoubleMatrix leafPartials = zeros(nodeWrap.numStatuses + 1, 1);
                Integer observedStateKey = nodeWrap.statusMap.get(nodeWrap.leafState);
                leafPartials.put(observedStateKey, 1.0);
                setNodePartials(node.getNr(), leafPartials);


            } else {
                setPartials(node.getLeft());
                setPartials(node.getRight());

            }
        }
    }
    //CHANGE block


    /**
     * Sets partials for a node
     */
    public void setNodePartials(int nodeIndex, DoubleMatrix partials) {

//        if (this.partials0.get(nodeIndex) == null) {
//            createNodePartials(nodeIndex);
//        }
        partials0.put(nodeIndex,partials);
    }

//    /**
//     * Sets partials for a node
//     */
//    public void setNodeMatrix(int nodeIndex, DoubleMatrix ptMat) {
//
////        if (this.partials0.get(nodeIndex) == null) {
////            createNodePartials(nodeIndex);
////        }
//        ptMat0.put(nodeIndex,ptMat);
//    }

    public double calculateLogP() {

        final TreeInterface tree = treeInput.get();

        try {
            //STEP2 :
            //CHANGE BLOCK
            //if (traverse(tree.getRoot()) != Tree.IS_CLEAN)
                calcLogP();
            //CHANGE BLOCK
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

        return logP;
    }

//    protected boolean wrapperRequiresRecalculation() {
//        if ( treeInput.get().getRoot().isDirty() == Tree.IS_FILTHY) {
//            return true;
//        }
//        else {
//            return false;
//        }

//    }


//    public DoubleMatrix getTransitionProbabilities(Node node, double parentheight, double endTime, double rate) {
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = child.getLength() * branchRate;
//        DoubleMatrix ptMat = expm((GeneralGestalt.createRateMatrix(childNodeWrap)).muli(branchTime));
//
//
//    }

//    /* Assumes there IS a branch rate model as opposed to traverse() */
//    protected int traverse(final Node node) {
//        //ANTOINE START traverse with node = root
//        // Has dirt is a general flag for calculation node
//        //when CLEAN=0, nothing needs to be recalculated for the node
//        //     * // when DIRTY=1 indicates a node partial needs to be recalculated
//        //     * // when FILTHY=2 indicates the indices for the node need to be recalculated
//        // BITWISE OR OPERATOR :
//        // 2|1 = 3
//        // 2|2 = 2
//        // 1|1 = 1
//        // 0|1 = 1
//        // 0|0 = 0
//
//        int update = (node.isDirty() | hasDirt);
//
//        final int nodeIndex = node.getNr();
//
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = node.getLength() * branchRate;
//
//        // First update the transition probability matrix(ices) for this branch
//        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
//        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
//            m_branchLengths[nodeIndex] = branchTime;
//            final Node parent = node.getParent();
//
//            //likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
//            DoubleMatrix ptMat = getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate);
//            setNodeMatrix(nodeIndex, ptMat);
//            //assign values of states for probability transition matrix for node with number nodeIndex
//
//            //ANTOINE read: update = update | Tree.IS_DIRTY
//            update |= Tree.IS_DIRTY;
//        }
//
//        // If the node is internal, update the partial likelihoods.
//        if (!node.isLeaf()) {
//
//            // Traverse down the two child nodes
//            final Node child1 = node.getLeft(); //Two children
//            final int update1 = traverse(child1);
//
//            final Node child2 = node.getRight();
//            final int update2 = traverse(child2);
//
//            // If either child node was updated then update this node too
//            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
//
//                final int childNum1 = child1.getNr();
//                final int childNum2 = child2.getNr();
//
//                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
//                update |= (update1 | update2);
//                if (update >= Tree.IS_FILTHY) {
//                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
//                }
//
//                likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
//
//
//                if (node.isRoot()) {
//                    // No parent this is the root of the beast.tree -
//                    // calculate the pattern likelihoods
//
//                    final double[] proportions = m_siteModel.getCategoryProportions(node);
//                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
//                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, patternLogLikelihoods);
//                }
//
//            }
//        }
//        return update;
//    } // traverseWithBRM



    void calcLogP() {
        //clean the hashtables:
        partials0 = new Hashtable<>();
//        partials1 = new Hashtable<>();
//        logScalingTerms = new Hashtable<>();

        final TreeInterface tree = treeInput.get();
        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree,dataInput.get(),substitutionModel.metaData.posSites,substitutionModel.metaData.nTargets);
        List<IndelSet.Singleton> list = substitutionModel.getAllSingletons(tree,statesDict);

        //update conditional probabilities
        substitutionModel.initSingletonProbs(list);

        //update target status rates
        substitutionModel.hazardAwayDict = substitutionModel.createHazardAwayDict();

        //empty targ_stat_transition_hazards_dict
        substitutionModel.targ_stat_transition_hazards_dict = new Hashtable<>();

        for(TargetStatus stat:substitutionModel.targStatTransitionsDict.keySet()) {
            Hashtable<TargetStatus,DoubleMatrix> empty = new Hashtable<>();
            substitutionModel.targ_stat_transition_hazards_dict.put(stat,empty);
        }

        //get transition wrappers
        long start1 = System.nanoTime();
        transitionWrappers = TransitionWrap.createTransitionWrappers(tree, dataInput.get(), substitutionModel.metaData, statesDict);
        long end1 = System.nanoTime();
        System.out.println("Elapsed Time in seconds: "+ (end1-start1)*0.000000001);
//        //initialize leaf partial likelihoods
//        setPartials(tree.getRoot());
//
//        //populate partial likelihoods with postorder traversal
//        traversal(tree.getRoot());
//
//        // take care of scaling terms
//        Collection<Double> logScalingTermsAll = logScalingTerms.values();
//        List<Double> terms = new ArrayList<>(logScalingTermsAll);
//        Double fullScaling = 0.0;
//        for(int i = 0; i < terms.size();++i) {
//            fullScaling = fullScaling + terms.get(i);
//        }
//
//        //get the likelihood from the root partial
//        int rootIndex = tree.getRoot().getNr();
//        DoubleMatrix logLikAlleles = logi(partials0.get(rootIndex)).addi(fullScaling);

        Double likelihood = createTopologyLogLikelihood(tree,transitionWrappers);
        Log.info.println("LOG LIKELIHOOD ALLELES"+likelihood);
        logP = likelihood;

    }


//    /**
//     * Allocates partials for a node
//     */
//    public void createNodePartials(int nodeIndex) {
//        //at this stage, partial sizes are not known!
//        partials0.put(nodeIndex,new DoubleMatrix(partialsSizes[nodeIndex])) ;
//        partials1.put(nodeIndex,new DoubleMatrix(partialsSizes[nodeIndex])) ;
//    }


//    public void traversal(Node node) {
//
//
//        if(node != null) {
//            if (!node.isLeaf()) {
//                // Traverse down the two child nodes
//                final Node child1 = node.getLeft(); //Two children
//                traversal(child1);
//
//                final Node child2 = node.getRight();
//                traversal(child2);
//
//
//                //DoubleMatrix nodeLogPartials = GeneralGestalt.initializeLowerLogProb(transitionWrappers.get(node.getNr()),node);
//                DoubleMatrix nodeLogPartials = DoubleMatrix.zeros(1);
//                DoubleMatrix hasPosProb = ones(1);
//
//
//                for (Node child : node.getChildren()) {
//
//
//                    TransitionWrap childNodeWrap = transitionWrappers.get(child.getNr());
//
//                    //Create the probability matrix exp(Qt)
//                    final double branchRate = branchRateModel.getRateForBranch(node);
//                    final double branchTime = child.getLength() * branchRate;
//                    //Log.info.println("Clock rate"+branchRate);
//                    DoubleMatrix rateM = substitutionModel.createRateMatrix(childNodeWrap);
//                    DoubleMatrix ptMat = expm(rateM.muli(branchTime));
//                    //ptMats.put(child.getNr(),ptMat);
//                    //Log.info.println("PTMATRIX"+ptMat);
//
//                    //Get the probability for the data descended from the child node, assuming that the node
//                    //has a particular target tract repr.
//                    //These down probs are ordered according to the child node's numbering of the TTs states
//                    DoubleMatrix chOrderedDownProbs = ptMat.mmul(partials0.get(child.getNr()));
//                    DoubleMatrix downProbs = new DoubleMatrix();
//
//                    if (!node.isRoot()) {
//
//                        // Reorder summands according to node's numbering of tract_repr states
//                        downProbs = GeneralGestalt.reorderLikelihoods(chOrderedDownProbs, transitionWrappers.get(node.getNr()), childNodeWrap);
//
//                    } else {
//                        //For the root node, we just want the probability where the root node is unmodified
//                        //No need to reorder
//                        Integer childUnmodifiedIndex = childNodeWrap.statusMap.get(new TargetStatus());
//                        double downProb = chOrderedDownProbs.get(childUnmodifiedIndex);
//                        downProb = max(downProb, 0.0);
//                        downProbs = new DoubleMatrix(1, 1, downProb);
//                    }
//                    //todo: implement node abundances
//                    //if (child.isLeaf()) {
//                    //Double leafAbundanceWeight = 1.0;
//                    hasPosProb = downProbs.ge(0);
//
//
//                    //protection against states with zero probability.
//                    nodeLogPartials = nodeLogPartials.addi(logi(downProbs.addi((hasPosProb.neg()).add(1).mul(1e-30))));
//
//
//                }
//                Double logScalingTerm = nodeLogPartials.max();
//
//
//                setNodePartials(node.getNr(), expi((nodeLogPartials.subi(logScalingTerm)).muli(hasPosProb)));
//
//                logScalingTerms.put(node.getNr(), logScalingTerm);
//            }
//
//
//        }
//
//    }


    public Double createTopologyLogLikelihood(TreeInterface tree, Hashtable<Integer, TransitionWrap> transitionWrappers) {
        Integer rootIndex = tree.getRoot().getNr();

        //Hashtable<Integer, DoubleMatrix> transMats = new Hashtable<>();
        //Hashtable<Integer, DoubleMatrix> Ddiags = new Hashtable<>();
        //Hashtable<Integer, DoubleMatrix> ptMats = new Hashtable<>();
        //Hashtable<Integer, DoubleMatrix> down_probs_dict= new Hashtable<>();
        //Hashtable<Integer, DoubleMatrix> partials = new Hashtable<>();

        Hashtable<Integer, Double> logScalingTerms = new Hashtable<>();

        for (Node node:tree.listNodesPostOrder(null,null)) {
            if(node.isLeaf()) {

                setPartials(node);

            }
            else {

                //DoubleMatrix nodeLogPartials = GeneralGestalt.initializeLowerLogProb(transitionWrappers.get(node.getNr()),node);
                DoubleMatrix nodeLogPartials = DoubleMatrix.zeros(1);
                DoubleMatrix hasPosProb = ones(1);


                for (Node child: node.getChildren()) {


                    TransitionWrap childNodeWrap = transitionWrappers.get(child.getNr());

                    //Create the probability matrix exp(Qt)
                    final double branchRate = branchRateModel.getRateForBranch(node);
                    final double branchTime = child.getLength() * branchRate;
                    DoubleMatrix ptMat = expm((substitutionModel.createRateMatrix(childNodeWrap)).muli(branchTime));
                    //ptMats.put(child.getNr(),ptMat);


                    //Get the probability for the data descended from the child node, assuming that the node
                    //has a particular target tract repr.
                    //These down probs are ordered according to the child node's numbering of the TTs states
                    DoubleMatrix chOrderedDownProbs = ptMat.mmul(partials0.get(child.getNr()));
                    DoubleMatrix downProbs = new DoubleMatrix();

                    if (! node.isRoot()) {

                        // Reorder summands according to node's numbering of tract_repr states
                        downProbs = GeneralGestalt.reorderLikelihoods(chOrderedDownProbs,transitionWrappers.get(node.getNr()),childNodeWrap);

                    }
                    else {
                        //For the root node, we just want the probability where the root node is unmodified
                        //No need to reorder
                        Integer childUnmodifiedIndex = childNodeWrap.statusMap.get(new TargetStatus());
                        double downProb = chOrderedDownProbs.get(childUnmodifiedIndex);
                        downProb = max(downProb,0.0);
                        downProbs = new DoubleMatrix(1,1,downProb);
                    }
                    //todo: implement node abundances
                    //if (child.isLeaf()) {
                    //Double leafAbundanceWeight = 1.0;
                    hasPosProb = downProbs.ge(0);


                    //protection against states with zero probability.
                    nodeLogPartials = nodeLogPartials.addi(logi(downProbs.addi((hasPosProb.neg()).add(1).mul( 1e-30 ))))  ;


                }
                Double logScalingTerm = nodeLogPartials.max();


                setNodePartials(node.getNr(),expi((nodeLogPartials.subi(logScalingTerm)).muli(hasPosProb)));

                logScalingTerms.put(node.getNr(), logScalingTerm);
            }



        }
        Collection<Double> logScalingTermsAll = logScalingTerms.values();


        List<Double> terms = new ArrayList<>(logScalingTermsAll);
        Double fullScaling = 0.0;
        for(int i = 0; i < terms.size();++i) {
            fullScaling = fullScaling + terms.get(i);
        }
        DoubleMatrix logLikAlleles = logi(partials0.get(rootIndex)).addi(fullScaling);

        /*partials = partials;
        down_probs_dict = down_probs_dict;
        pt_matrix = pt_matrix;
        trans_mats = trans_mats;
        trim_probs = trim_probs;*/


        return logLikAlleles.get(0);
    }

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



}
