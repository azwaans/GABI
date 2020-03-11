package lineageTree.distributions;


import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;


public class organoidTreeLikelihood extends GenericTreeLikelihood {

    final public Input<Double> scarTimeInput = new Input<>("scarringTime", "Duration between the root and the onset of scarring", -1.0);


    @Override
    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

    }

    /**
     * check, that after scarring event every subtree has only leaf nodes of the same sequence.
     * @return -inf if tree is not allowed, 0 otw
     */
    public double calculateLogP(){

        final TreeInterface tree = treeInput.get();
        final double scarTime = scarTimeInput.get();
        double scarringHeight = tree.getRoot().getHeight() - scarTime;
        Alignment alignment = dataInput.get();


        List<Node> nodesAfterScarring = getNodesAfterScarring(tree, scarringHeight);

        for  (int i=0; i < nodesAfterScarring.size(); i++){

            Node subTreeRoot = nodesAfterScarring.get(i);
            String subTreeRootSeq  = alignment.getSequenceAsString(subTreeRoot.getID());
            List<Node> leafNodes = subTreeRoot.getAllLeafNodes();

            for (int j=0; j < leafNodes.size(); j++){
                String leafNodeSeq = alignment.getSequenceAsString(leafNodes.get(j).getID());

                if( !( subTreeRootSeq == leafNodeSeq)){
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
            if (currentNode.getHeight() <= scarringHeight & currentNode.getParent().getHeight() > scarringHeight){
                nodesAfterScarring.add(currentNode);
            }else {
                queue.add(currentNode.getChild(0));
                queue.add(currentNode.getChild(1));
            }
        }

        return  nodesAfterScarring;
    }


}
