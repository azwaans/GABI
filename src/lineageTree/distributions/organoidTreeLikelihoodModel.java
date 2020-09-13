package lineageTree.distributions;


import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Queue;

@Description("Computes probability of the tree given a simple model")
public class organoidTreeLikelihoodModel extends GenericTreeLikelihood {

    final public Input<Double> scarringHeightInput = new Input<>("scarringHeight", "Duration between the onset of scarring and sampling of the cells", -1.0);

    double scarringHeight;

    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

        Alignment alignment = dataInput.get();
        scarringHeight = scarringHeightInput.get();
    }

    /**
     * check, that after scarring event every subtree has only leaf nodes of the same sequence.
     * @return -inf if tree is not allowed, 0 otw
     */
    public double calculateLogP(){

        // achtung ich suche nur nach subtreenodes, die unter der scarring height anfangen
        // damit "strafe ich keine nodes ab, die über der scarring height sind

        // füge funktion hinzu, die prüft, dass der rootnode aller identical sequences unter der scarringheight ist.
        // -> and <-


        final Alignment alignment = dataInput.get();
        TreeInterface tree = treeInput.get();


        // assign state to internal nodes
        //TODO check number of scars! Make state list available in metadastring to facilitate trouble shooting
        List<Integer> rootSeq = assign_internal_node_states(tree.getRoot(), alignment);


        //check that leaves below the subTreeRoot are identical (->)
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

        // check that all identical leave nodes coalesce before the scarringHeight (<-)
        List<Node> leaves = tree.getExternalNodes();

        outerloop:
        for (int i=0; i<leaves.size(); i++){
            for (int j=i+1; j<leaves.size(); j++){

                Node leaf_1 = leaves.get(i);
                Node leaf_2 = leaves.get(j);

                // if either leaf has no edits (is marked with cluster 0)
                if(leaf_1.getMetaData("cluster").equals(0)){
                    break;

                }else if(leaf_2.getMetaData("cluster").equals(0)){
                    continue;

                // else enforce condition
                }else{
                    String seq_1 = alignment.getSequenceAsString(leaf_1.getID());
                    String seq_2 = alignment.getSequenceAsString(leaf_2.getID());

                    //if the sequences are identical, their MRCA's height has to be <= scarringHeight
                    if (seq_1.compareTo(seq_2) == 0 ) {
                        // if the MRCA's height is above the scarring event, reject tree
                        if (!MRCA_before_scarring(scarringHeight, leaf_1, leaf_2)) {
                            return (Double.NEGATIVE_INFINITY);
                        }
                    }
                }
            }
        }

        return(0);
    }

    void calcLogP() {
        logP = 0.0;
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


    /**
     * Determines, whether the nodes 1 and 2 coalesce before (between the present and) the scarring event
     * @param scarringHeight
     * @return
     */
    public boolean MRCA_before_scarring(final double scarringHeight, Node leaf_1, Node leaf_2){

        Node parent = leaf_1.getParent();

        while (true){

            //System.out.println("Error from leaf1: "+leaf_1.getNr() + " leaf2 :  " + leaf_2.getNr());
            List<Node> childrenOfParent = parent.getAllChildNodesAndSelf();

            //if parent is the MRCA of leaf_1 and leaf_2
            if (childrenOfParent.contains(leaf_2)){
                // check whether MRCA is between leaves and scarringHeight
                MathContext m = new MathContext((3));
                BigDecimal parentHeight = new BigDecimal(parent.getHeight());
                BigDecimal roundedParentHeight = parentHeight.round(m);

                if( roundedParentHeight.doubleValue() <= scarringHeight){
                    return true;
                }else{
                    return false;
                }
            // if parent is not MRCA, choose its parent and test again
            }else{
                parent = parent.getParent();
            }
        }
    }

    /*
    assign scar sequence to internal nodes
    (1) the leaf nodes have the observed number of scars
    (2) the internal nodes up to the scarring event have the maximum number of scars of their children
    (3) the internal nodes directly after (above) the scarring event have M unedited scarring sites;
      the number of unedited scarring sites equals the sum over all scar types that would appear in this internal node by (2)
    (4) the internal nodes above 'the internal nodes directly above the scarring event' have the maximum number of
        unedited sites from either of their children
        */
    public List<Integer> assign_internal_node_states(Node node, Alignment alignment){
        DataType datT = new IntegerData();

        // given node is leaf, report leaf sequence
        if (node.isLeaf()){
            String leafSeq = alignment.getSequenceAsString(node.getID());
            List<Integer> intList = datT.stringToEncoding(leafSeq);

            // set scarring sequnce as metadata
            node.setMetaData("state", intList);
            node.metaDataString = (node.metaDataString + ",state="+
                    intList.toString().replace(",", "").replace("[", "").replace("]", "").trim());
            return(intList);

        // given internal node, check childrens sequences and ...
        }else{
            List<Integer> child1_seq = assign_internal_node_states(node.getChild(0), alignment);
            List<Integer> child2_seq = assign_internal_node_states(node.getChild(1), alignment);


            // ... create sequence for this node
            // use the same seq element from either child if they are equal
            // otw use the larger integer from either child
            List<Integer> node_seq =  new ArrayList<>();

            for (int i=0; i <child1_seq.size(); i++){
                Integer elem1 = child1_seq.get(i);
                Integer elem2 = child2_seq.get(i);

                if (elem1.equals(elem2)){
                        node_seq.add(elem1);

                }else if (elem1 > elem2){
                        node_seq.add(elem1);

                } else if (elem2 > elem1){
                        node_seq.add(elem2);
                } else{
                        System.err.println("Integer should be either equal, smaller or larger");
                }
            }

            // if node appears before scarring, then the number of scars have to be unedited scarring sites
            if (node.getHeight() > scarringHeight){
                int nScars = 0;
                for ( int i=0; i<node_seq.size(); i++){
                    int nScarsOfType = node_seq.get(i);
                    nScars += nScarsOfType;
                    // set all edited scars to 0, as they were aquired during the scarring event
                    node_seq.set(i, 0);
                }

                // set new sequence of unedited scars
                node_seq.set(0, nScars);
            }

            node.setMetaData("state", node_seq);
            node.metaDataString = (node.metaDataString + ",state="+
                    node_seq.toString().replace(",", "").replace("[", "").replace("]", "").trim());

            return node_seq;
        }
    }


}
