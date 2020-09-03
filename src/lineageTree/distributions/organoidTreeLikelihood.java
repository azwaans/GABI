package lineageTree.distributions;


import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;

@Description("Accepts or reject trees, given that they conform with assumptions on the scarring process")
public class organoidTreeLikelihood extends GenericTreeLikelihood {

    final public Input<Double> scarTimeInput = new Input<>("scarringTime", "Duration between the root and the onset of scarring", -1.0);
    //final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Input.Validate.REQUIRED);
    //final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Input.Validate.REQUIRED);


    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

        Alignment alignment = dataInput.get();
    }

    /**
     * check, that after scarring event every subtree has only leaf nodes of the same sequence.
     * @return -inf if tree is not allowed, 0 otw
     */
    public double calculateLogP(){

        final double scarTime = scarTimeInput.get();
        final Alignment alignment = dataInput.get();
        final TreeInterface tree = treeInput.get();
        double scarringHeight = 32 - scarTime; //tree.getRoot().getHeight()


        List<Node> nodesAfterScarring = getNodesAfterScarring(tree, scarringHeight);

        for  (int i=0; i < nodesAfterScarring.size(); i++){

            Node subTreeRoot = nodesAfterScarring.get(i);
            List<Node> leafNodes = subTreeRoot.getAllLeafNodes();
            String seqToCompare = alignment.getSequenceAsString(leafNodes.get(0).getID());

            for (int j=0; j < leafNodes.size(); j++){
                String leafNodeSeq = alignment.getSequenceAsString(leafNodes.get(j).getID());

                if( seqToCompare.compareTo(leafNodeSeq) != 0){
                    return (Double.NEGATIVE_INFINITY);
                }
            }
        }

        return(0);
    }



    /** Collects all nodes, that appear immediately after (~ toward the leaves) the scarring event using BFS
    **/
    public List<Node> getNodesAfterScarring(final TreeInterface tree, final double scarringHeight){

        List<Node> nodesAfterScarring = new ArrayList<>();

        // init queue and add children of root to it,
        // assuming that scarring happens after some rounds of cell cycles at least
        Queue<Node> queue = new ArrayDeque<>();
        queue.add(tree.getRoot().getChild(0));
        queue.add(tree.getRoot().getChild(1));

        while(!queue.isEmpty()){

            Node currentNode = queue.remove();
            // collect all internal nodes that initiate a subtree below the scarring time
            if (!(currentNode.isLeaf()) & currentNode.getHeight() <= scarringHeight & currentNode.getParent().getHeight() > scarringHeight){
                nodesAfterScarring.add(currentNode);
            }else if (!currentNode.isLeaf()){

                queue.add(currentNode.getChild(0));
                queue.add(currentNode.getChild(1));
            }
        }

        return  nodesAfterScarring;
    }


}
