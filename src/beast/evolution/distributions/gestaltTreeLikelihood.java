package beast.evolution.distributions;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Tree;
import org.jblas.DoubleMatrix;

import java.util.*;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.core.Log;
import beast.evolution.alignment.*;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.evolution.substitutionmodel.gestaltGeneral;

//import static beast.evolution.distributions.GapmlLikelihoodCore.ptMats;
//import static beast.evolution.distributions.GapmlLikelihoodCore.setNodeMatrix;
import static org.jblas.MatrixFunctions.*;

import static java.lang.Math.max;



@Description("Generic tree likelihood for an alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model")
// Use this as base class to define any non-standard TreeLikelihood.
// Override Distribution.calculatLogP() to make this class functional.
// TODO: This could contain a generic traverse() method that takes dirty trees in account.
//
public class gestaltTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    protected gestaltGeneral substitutionModel;

    protected SiteModel.Base m_siteModel;

    protected BranchRateModel.Base branchRateModel;

    protected double[] m_branchLengths;

    protected double[] storedBranchLengths;

    protected List<IndelSet.Singleton> singletonList;


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

        substitutionModel = (gestaltGeneral) m_siteModel.substModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        //CHANGE BLOCK
        likelihoodCore = new gestaltLikelihoodCore();
        initCore();



    }

    //CHANGE BLOCK
    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.init(nodeCount);

        /**
         * set leaf partials in likelihood core *
         */

        hasDirt = Tree.IS_FILTHY;

    }

    public double calculateLogP() {

        final TreeInterface tree = treeInput.get();

        //recording time
        long start1 = System.nanoTime();

        if(hasDirt == Tree.IS_FILTHY | likelihoodCore.transitionWrappers == null) {
            Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree,dataInput.get(),substitutionModel.metaData.posSites,substitutionModel.metaData.nTargets);
            singletonList = substitutionModel.getAllSingletons(tree,statesDict);
            Log.info.println("we are attempting to recalculate all of the transition wrappers, has dirt is:" + hasDirt);
            likelihoodCore.transitionWrappers = TransitionWrap.createTransitionWrappers(tree, dataInput.get(), substitutionModel.metaData, statesDict);
        }

        //if there is dirt, only
        if(hasDirt >= Tree.IS_DIRTY) {
            //update conditional probabilities
            substitutionModel.initSingletonProbs(singletonList);

            //update target status rates
            substitutionModel.hazardAwayDict = substitutionModel.createHazardAwayDict();

            //empty targ_stat_transition_hazards_dict
            substitutionModel.targ_stat_transition_hazards_dict = new Hashtable<>();

            for(TargetStatus stat:substitutionModel.targStatTransitionsDict.keySet()) {
                Hashtable<TargetStatus,DoubleMatrix> empty = new Hashtable<>();
                substitutionModel.targ_stat_transition_hazards_dict.put(stat,empty);
            }
        }

        long end1 = System.nanoTime();
        System.out.println("Elapsed Time in seconds: "+ (end1-start1)*0.000000001);

        try {

            //populate partial likelihoods with postorder traversal
            if (traverse(tree.getRoot()) != Tree.IS_CLEAN) {
                calcLogP();
            }

        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }

        return logP;
    }


   /* Assumes there IS a branch rate model as opposed to traverse() */
    protected int traverse(final Node node) {

        int update =  hasDirt;

        if(node != null) {
            update = (node.isDirty() | hasDirt);
            final int nodeIndex = node.getNr();
            Log.info.println("HAS DIRT in tree likelihood? " + hasDirt);
            Log.info.println("node is dirty? " + node.isDirty());

            final double branchRate = branchRateModel.getRateForBranch(node);
            final double branchTime = node.getLength() * branchRate;

            // First update the transition probability matrix(ices) for this branch
            //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
            if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
                m_branchLengths[nodeIndex] = branchTime;

                //likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
                TransitionWrap nodeWrap = likelihoodCore.transitionWrappers.get(node.getNr());
                DoubleMatrix ptMat = substitutionModel.getTransitionProbabilities(node.getParent(), nodeWrap, node.getLength(), branchRate);
                likelihoodCore.setNodeMatrix(nodeIndex, ptMat);


                update |= Tree.IS_DIRTY;
            }

            // Only update the partial likelihoods if the node is a leaf
            if(node.isLeaf()) {
                likelihoodCore.setLeafPartials(node);
            }
            else {

                // Traverse down the two child nodes
                final Node child1 = node.getLeft(); //Two children
                final int update1 = traverse(child1);

                final Node child2 = node.getRight();
                final int update2 = traverse(child2);

                // If either child node was updated then update this node too
                if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {


                    update |= (update1 | update2);
                    //likelihoodCore.setNodePartialsForUpdate(nodeIndex);
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
        Collection<Double> logScalingTermsAll = likelihoodCore.logScalingTerms.values();
        List<Double> terms = new ArrayList<>(logScalingTermsAll);
        Double fullScaling = 0.0;
        for(int i = 0; i < terms.size();++i) {
            fullScaling = fullScaling + terms.get(i);
        }

        //get the likelihood from the root partial
        int rootIndex = tree.getRoot().getNr();
        DoubleMatrix logLikAlleles = logi(likelihoodCore.getNodePartials(rootIndex)).addi(fullScaling);
        Log.info.println("LOG LIKELIHOOD ALLELES"+logLikAlleles.get(0));
        logP = logLikAlleles.get(0);

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



    /**
     * check state for changed variables and update temp results if necessary *
     */
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




}
