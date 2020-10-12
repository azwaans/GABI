package lineageTree.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;


@Description("Computes probability of the tree given a simple model")
public class organoidTreeLikelihoodModel extends GenericTreeLikelihood {

    final public Input<Double> scarringHeightInput = new Input<>("scarringHeight", "Duration between the onset of scarring and sampling of the cells", -1.0);
    final public Input<Double> scarringDurationInput = new Input<>("scarringDuration", "Duration of scarring");
    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SiteModel.Base m_siteModel;
    protected SubstitutionModel substitutionModel;

    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}


    double[] probabilities;
    double scarringHeight;
    double scarringDuration;

    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }

        //init site and substitution model related things
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();
        probabilities = new double[4];
        Arrays.fill(probabilities, 1.0);

        // init input arguments
        Alignment alignment = dataInput.get();
        scarringHeight = scarringHeightInput.get();
        scarringDuration = scarringDurationInput.get();
    }

    /**
     * check, that after scarring event every subtree has only leaf nodes of the same sequence.
     * @return -inf if tree is not allowed, 0 otw
     */
    public double calculateLogP(){

        // init
        final Alignment alignment = dataInput.get();
        TreeInterface tree = treeInput.get();


        // assign state to internal nodes
        List<Integer> rootSeq = assign_internal_node_states(tree.getRoot(), alignment);

        double logP = 0;
        List<Node> intNodes = tree.getInternalNodes();

        //loop over all branch lengths
        for (int i=0; i<intNodes.size(); i++){

            Node parent = intNodes.get(i);
            double ph = parent.getHeight();
            List<Integer> pSeq = (List<Integer>)parent.getMetaData("state");

            List<Node> children = parent.getChildren();

            for(int j=0; j<children.size(); j++){

                Node child = children.get(j);
                double ch = child.getHeight();
                List<Integer> cSeq = (List<Integer>)child.getMetaData("state");

                // test to which tree segment parent and child belong to get transition probabilities
                if (ph > scarringHeight){
                    if (ch > scarringHeight){
                        substitutionModel.getTransitionProbabilities(child, ph, ch, 1, probabilities);

                        logP += calcBranchLogPLoss(pSeq, cSeq, probabilities);

                    }else if(ch > (scarringHeight - scarringDuration)){
                        Node intermediate = new Node();
                        intermediate.setHeight(scarringHeight);

                        // add up segment 1 prob
                        substitutionModel.getTransitionProbabilities(intermediate, ph, scarringHeight, 1, probabilities);
                        // NOTE: Here I make an assumption on the assignment of internal node states as implemented in the
                        // function below. On the branch between a node in the scarring window and a node above the scarring window, only scarring occurs
                        // -> no loss -> therefore I set the seq here to the parent sequence
                        int nBarcodeP = pSeq.stream().mapToInt(Integer::intValue).sum();
                        int nBarcodeC = cSeq.stream().mapToInt(Integer::intValue).sum();

                        // no loss event
                        if( nBarcodeP == nBarcodeC){
                            logP += probabilities[0] * nBarcodeC; // not lost
                        // loss event needed
                        } else if(nBarcodeP > nBarcodeC){
                            int diff = nBarcodeP - nBarcodeC;
                            logP += probabilities[1]* diff;       // lost
                            logP += probabilities[0] * nBarcodeC; //not lost
                        }else {
                            throw new RuntimeException("Parent has " + nBarcodeP +" and child has " + nBarcodeC + " barcodes. Wrong internal node assignment!");
                        }

                        // segment 2 prob
                        substitutionModel.getTransitionProbabilities(child, scarringHeight, ch, 1, probabilities);
                        logP += calcBranchLogPScar(pSeq, cSeq, probabilities);

                    } else if (ch < (scarringHeight - scarringDuration)){
                        Node intermediate1 = new Node();
                        Node intermediate2 = new Node();
                        intermediate1.setHeight(scarringHeight);
                        intermediate2.setHeight(scarringHeight-scarringDuration);

                        // add up segment 1 prob
                        substitutionModel.getTransitionProbabilities(intermediate1, ph, scarringHeight, 1, probabilities);
                        logP += calcBranchLogPLoss(pSeq, pSeq, probabilities);

                        // segment 2 prob
                        substitutionModel.getTransitionProbabilities(intermediate2, scarringHeight, scarringHeight-scarringDuration, 1, probabilities);
                        logP += calcBranchLogPScar(pSeq, cSeq, probabilities);

                        //segment 3 prob
                        substitutionModel.getTransitionProbabilities(child, scarringHeight-scarringDuration, ch, 1, probabilities);
                        //logP += calcBranchLogPLoss();


                    }else{
                        throw new RuntimeException("Strange child height: "+ ch +"!");
                    }
                }else if(ph > (scarringHeight - scarringDuration)){
                    if(ch > (scarringHeight-scarringDuration)){

                        substitutionModel.getTransitionProbabilities(child, ph, ch, 1, probabilities);

                    } else if (ch < (scarringHeight - scarringDuration)) {
                        Node intermediate = new Node();
                        intermediate.setHeight(scarringHeight-scarringDuration);

                        //add segment 1
                        substitutionModel.getTransitionProbabilities(intermediate, ph, intermediate.getHeight(), 1, probabilities);
                        // and segment 2
                        substitutionModel.getTransitionProbabilities(child, intermediate.getHeight(), ch, 1, probabilities);
                    }else{
                        throw new RuntimeException("Strange child height: "+ ch +"!");
                    }
                }else if (ph < (scarringHeight - scarringDuration)){
                    substitutionModel.getTransitionProbabilities(child, ph, ch, 1, probabilities);
                } else{
                    throw new RuntimeException("Parent height: " + ph + " is not valid!");
                }



                // test how states relate to each other
                if (child.getMetaData("state").equals(parent.getMetaData("state")) ){
                    substitutionModel.getTransitionProbabilities(child, parent.getHeight(), child.getHeight(), 1, probabilities);

                }
            }

        }

        return(0);
    }

    //TODO potentially split this for before scarring window (then only diff of first array element needed) and after scarring
    double calcBranchLogPLoss(List<Integer> pSeq, List<Integer> cSeq, double[] probabilities) {

        double branchLogP = 0;
        // test for sequence similarity
        for (int scarType=0; scarType<pSeq.size(); scarType++){
            int pScar = pSeq.get(scarType);
            int cScar = cSeq.get(scarType);

            //identical barcodes
            if ( (pScar == cScar) & (pScar > 0)){
                 branchLogP += Math.log(probabilities[0]);

            // parent has more barcodes than child -> loss event
            }else if (pScar > cScar){
                int diff = pScar - cScar;
                branchLogP += Math.log(probabilities[1]) * diff;

            // parent has less barcodes than child -> impossible
            }else{
                return Double.NEGATIVE_INFINITY;
            }
        }

        return branchLogP;
    }

    double calcBranchLogPScar(List<Integer> pSeq, List<Integer> cSeq, double[] probabilities){

        double branchLogP =0;

        // in the first array element are the number of unedited barcodes
        int nNewScars = pSeq.get(0) - cSeq.get(0);


        branchLogP += Math.log(probabilities[1]) * nNewScars;
        branchLogP += Math.log(probabilities[0]) * cSeq.get(0);

        for (int scarType=1; scarType<pSeq.size(); scarType++) {
            int pScar = pSeq.get(scarType);
            int cScar = cSeq.get(scarType);

            //identical scars
            if ((pScar == cScar) & (pScar > 0)){
                // as this
                branchLogP += Math.log(probabilities[4]);

                // parent has more scars than child -> loss event which is not allowed during scarring
            }else if (pScar > cScar){
                return Double.NEGATIVE_INFINITY;
            }
        }

        return 0;

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

        // leaf node (1)
        if (node.isLeaf()){
            String leafSeq = alignment.getSequenceAsString(node.getID());
            List<Integer> intList = datT.stringToEncoding(leafSeq);

            // set scarring sequnce as metadata
            node.setMetaData("state", intList);
            node.metaDataString = (node.metaDataString + ",state="+
                    intList.toString().replace(",", "").replace("[", "").replace("]", "").trim());
            return(intList);

        // internal node
        }else{
            List<Integer> child1_seq = assign_internal_node_states(node.getChild(0), alignment);
            List<Integer> child2_seq = assign_internal_node_states(node.getChild(1), alignment);
            List<Integer> node_seq =  new ArrayList<>(child1_seq.size());

            // internal node above scarring event (3), (4)
            if (node.getHeight() > scarringHeight){
                int nScars_1 = 0;
                int nScars_2 = 0;

                for (int i=0; i<child1_seq.size(); i++){
                    nScars_1 += child1_seq.get(i);
                    nScars_2 += child2_seq.get(i);
                    node_seq.add(0);
                }

                node_seq.set(0, java.lang.Math.max(nScars_1, nScars_2)); // unedited scar type  == pos 0

            // internal node up to scarring event (2)
            }else {

                for (int i = 0; i < child1_seq.size(); i++) {
                    Integer elem1 = child1_seq.get(i);
                    Integer elem2 = child2_seq.get(i);

                    if (elem1.equals(elem2)) {
                        node_seq.add(elem1);

                    } else if (elem1 > elem2) {
                        node_seq.add(elem1);

                    } else if (elem2 > elem1) {
                        node_seq.add(elem2);
                    } else {
                        System.err.println("Integer should be either equal, smaller or larger");
                    }
                }
            }

            node.setMetaData("state", node_seq);
            node.metaDataString = (node.metaDataString + ",state="+
                    node_seq.toString().replace(",", "").replace("[", "").replace("]", "").trim());

            return node_seq;
        }
    }


}
