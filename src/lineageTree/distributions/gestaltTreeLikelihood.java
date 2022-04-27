package lineageTree.distributions;


import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;

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
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import lineageTree.substitutionmodel.GeneralGestalt;
import org.jblas.DoubleMatrix;
import static org.jblas.MatrixFunctions.*;

import static com.google.common.math.DoubleMath.mean;
import static org.jblas.MatrixFunctions.logi;


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

    protected Hashtable<Integer, TransitionWrap> transitionWrappers;


    @Override
    public List<String> getArguments() {return null;}

    @Override
    public List<String> getConditions() {return null;}

    @Override
    public void sample(State state, Random random) {}

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
        //substitutionModel = (GeneralGestalt) m_siteModel.substModelInput.get();
        //scarringStart = substitutionModel.getScarringHeight();
        //scarringStop = scarringStart - substitutionModel.getScarringDuration();

        substitutionModel = (GeneralGestalt) m_siteModel.substModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }

        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        Alignment alignment = dataInput.get();



    }

    protected boolean wrapperRequiresRecalculation() {
        if ( treeInput.get().getRoot().isDirty() == Tree.IS_FILTHY) {
            return true;
        }
        else {
            return false;
        }

    }

    @Override
    public double calculateLogP() {
//        ///test_branch_likelihood:
//        for (Node node : treeInput.get().listNodesPostOrder(null, null)) {
//            Log.info.println("node.isDirty() : "+ node.isDirty() + "\n");
//        }
       // int isROOTDIRTY = treeInput.get().getRoot().isDirty();
        // Log.info.println("ROOT.isDirty() : "+ isROOTDIRTY);

        //Log.info.println("wrapperrequiresRecalculation() : "+wrapperRequiresRecalculation() + "\n");
        //Hashtable<Integer, TransitionWrap> transitionWrappers = TransitionWrap.createTransitionWrappers(treeInput.get(),dataInput.get(),substitutionModel.metaData);
        if (wrapperRequiresRecalculation() == true) {
            transitionWrappers = TransitionWrap.createTransitionWrappers(treeInput.get(),dataInput.get(),substitutionModel.metaData);
        }
        //debugging for 222 seed
        ///  TransitionWrap eigth = transitionWrappers.get(8);
//        Log.info.println("node 8 node ID" + treeInput.get().getNode(8).getID());
//        Log.info.println("node 8 node number" + treeInput.get().getNode(8).getNr());
//        Log.info.println("node 8 is leaf? " + treeInput.get().getNode(8).isLeaf());
//        Log.info.println("node 8 status map size " + eigth.statusMap.size());

//        for (TargetStatus status : eigth.statusMap.keySet()) {
//            Log.info.println("in the eigth node status map" +Arrays.asList(status.getBinaryStatus(10)));
//        }
//        Log.info.println("TOTAL SIZE OF THE TRANSITION WRAPPER, whole tree"+ transitionWrappers.size());
//        //THAT?S CORRECT: 1 for the root and one for the leaf?
//        //check that the only possible state at the root is empty:
//        TransitionWrap rootWrap = transitionWrappers.get(treeInput.get().getRoot().getNr());
//        Log.info.println("root wrap size"+ rootWrap.numStatuses );
//        Log.info.println("root wrap ID"+ treeInput.get().getRoot().getNr() );
//        for(TargetStatus stat:rootWrap.transStatuses) {
//            Log.info.println("status in root wrap" + Arrays.toString(stat.getBinaryStatus(substitutionModel.metaData.nTargets)));
//        }
//        for(int i=0;i<2;i++) {
//            Log.info.println("non root wrap ID, is LEAF? T/F:"+ treeInput.get().getNode(i).isLeaf());
//            TransitionWrap leafWrap = transitionWrappers.get(i);
//            Log.info.println("non root wrap size" + leafWrap.numStatuses);
//            for (TargetStatus stat : leafWrap.transStatuses) {
//                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
//            }
//        }
        //THE WRAPPERS ARE CORRECT

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(treeInput.get(),dataInput.get(),substitutionModel.metaData.posSites,substitutionModel.metaData.nTargets);
        List<IndelSet.Singleton> list = substitutionModel.getAllSingletons(treeInput.get(),statesDict);
        //Log.info.println("GET ALL SINGLETONS SIZE" + list.size());
        substitutionModel.initSingletonProbs(list);
        Double likelihood = substitutionModel.createTopologyLogLikelihood(treeInput.get(),transitionWrappers);
        Log.info.println("likelihood TREE LIKELIHOOD" + likelihood + "\n");
        logP = (double) likelihood;

        //double pen_test = penalization();

        return logP;
    }

    protected double penalization() {

        for(int i =0;  i< treeInput.get().getNodeCount(); i++) {
            Node node = treeInput.get().getNode(i);
            final double branchRate = branchRateModel.getRateForBranch(node);
            final double branchTime = node.getLength() * branchRate;
            m_branchLengths[i] = branchTime;

            if (branchTime < 0.0) {
                    throw new RuntimeException("Negative branch length: " + branchTime);
            }

        }

        //Log.info.println("branch lengths"+ Arrays.toString(m_branchLengths));
        DoubleMatrix branchLengths = new DoubleMatrix(m_branchLengths);
        DoubleMatrix cutRates = new DoubleMatrix(Arrays.asList(substitutionModel.cutRates));
        DoubleMatrix logbranchLenghts = logi(branchLengths);
        DoubleMatrix logcutRates = logi(cutRates);
        Double branchLengPen = substitutionModel.branchLensPenalty*(powi(logbranchLenghts.subi(logbranchLenghts.mean()),2)).sum();
        Double cutRatesPen = substitutionModel.cutRatesPenalty*(powi(logcutRates.subi(logcutRates.mean()),2)).sum();

        return branchLengPen + cutRatesPen;

    }



}
