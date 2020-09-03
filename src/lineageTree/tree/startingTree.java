package lineageTree.tree;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;

import javax.management.RuntimeErrorException;
import java.util.List;
import java.util.Vector;
import java.util.Arrays;
import java.util.stream.IntStream;

@Description("This class provides the basic engine for generating a valid starting tree given information on the scarring experiment ")
public class startingTree extends Tree implements StateNodeInitialiser {
    final public Input<Alignment> taxaInput = new Input<>("taxa", "set of taxa to initialise tree specified by alignment");

    final public Input<Double> rootHeightInput = new Input<>("rootHeight", "Time from beginning of the experiment until sequencing");
    final public Input<Double> scarringHeightInput = new Input<>("scarringHeight", "Time from the onset of scarring until sequencing");
    final public Input<Integer> nClustersInput = new Input<>("nClusters", "Number of clusters, where each cluster consists of identical sequences");

    // set up useful parameters
    int[][] matchMatrix;
    int nTaxa;

    Node treeNode;

    @Override
    public void initAndValidate() {

        // init taxa
        Alignment taxa = taxaInput.get();
        List<String> taxanames = taxa.getTaxaNames();

        nTaxa = taxanames.size();

        double rootHeight = rootHeightInput.get();
        double scarringHeight = scarringHeightInput.get();

        if (scarringHeight >= rootHeight){
            throw new RuntimeException("ScarringHeight has to be strictly smaller than rootHeight");
        }

        int nClusters = nClustersInput.get();
        if (nClusters < 1){
            throw new RuntimeException("Cluster size has to be positive");
        }

        matchMatrix = set_match_matrix(taxa);


        Node treeNode = get_tree(rootHeight, scarringHeight, taxa, nClusters, matchMatrix);
        initArrays();
        assignFromWithoutID(new Tree(treeNode));

        initStateNodes();
        super.initAndValidate();
    }

    @Override
    public void initStateNodes() {
        //TODO which statenodes have to be initialised? -  the tree
        //TODO when does super.initAndValidate have to be called?

        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(m_initial.get());
    }

    public int[][] set_match_matrix(Alignment taxa){
    // calculate seq by seq matrix that shows identical (1) and non-identical (0) sequences

        matchMatrix = new int[nTaxa][nTaxa];
        List<String> taxanames = taxa.getTaxaNames();

        for (int i=0; i<nTaxa; i++){
            for (int j=i+1; j<nTaxa; j++){

                // get sequence
                String seq_i = taxa.getSequenceAsString(taxanames.get(i));
                String seq_j = taxa.getSequenceAsString(taxanames.get(j));

                matchMatrix[i][j] = (seq_i.equals(seq_j)) ? 1 : 0;
            }
        }
        return matchMatrix;
    }

    public Node get_cluster_tree(int iSeq, double scarringHeight, Alignment taxa){
        //build subtree for identical sequences with equally spaced divergence
        // times between the internal nodes up until the scarring event

        int nMatches = IntStream.of(matchMatrix[iSeq]).sum();
        double divTimes =  scarringHeight / nMatches;
        List<String> taxaNames = taxa.getTaxaNames();

        // generate left node
        Node nodeLeft = new Node();
        nodeLeft.setHeight(0);
        nodeLeft.setNr(iSeq);
        nodeLeft.setID(taxaNames.get(iSeq));

        for (int iMatch=1; iMatch<nMatches; iMatch++){
            //set up right node
            Node nodeRight = new Node();
            nodeRight.setHeight(0);
            nodeRight.setID(taxaNames.get(iSeq + iMatch));
            nodeRight.setNr(iSeq + iMatch);

            //set up parent node
            Node parent = new Node();
            // space internal nodes at equal intervals until scarringHeight
            parent.setHeight(divTimes * iMatch);
            parent.setNr(taxaNames.size() -1 + iSeq + iMatch);
            parent.addChild(nodeLeft);
            parent.addChild(nodeRight);

            nodeLeft = parent;
        }

        return nodeLeft;

    }

    public Node get_tree(double rootHeight, double scarringHeight, Alignment taxa, int nClusters, int[][] matchMatrix){

        //find starting sequence for cluster
        int iSeq = 0;
        int nMatches;
        int nTaxa = taxa.getTaxonCount(); // check that this is what I want

        // time between the internal nodes connecting the cluster trees
        double divTime = (rootHeight - scarringHeight) / nClusters;

        // get left subtree
        Node subtreeLeft = get_cluster_tree(iSeq, scarringHeight, taxa);
        // update iSeq to next cluster start
        nMatches = IntStream.of(matchMatrix[iSeq]).sum();
        iSeq += (nMatches + 1);


        for (int iCluster=1; iCluster < nClusters; iCluster++){

            // get right subtree
            Node subtreeRight = get_cluster_tree(iSeq, scarringHeight, taxa);
            nMatches = IntStream.of(matchMatrix[iSeq]).sum();
            iSeq += (nMatches + 1);

            // get parent
            Node parent = new Node();
            parent.setHeight(scarringHeight + iCluster * divTime);
            // Number of nodes in the final tree - number of nodes that have yet to be created
            parent.setNr(nTaxa + (nTaxa - 1)     - (nClusters-iCluster)); // check

            parent.addChild(subtreeLeft);
            parent.addChild(subtreeRight);

            // setup as left tree to end recursion
            subtreeLeft = parent;

        }

        if (iSeq > nTaxa){
            throw new RuntimeException("iSeq reached " + iSeq + " but it should never exceed the number of sequences: " + nTaxa + "!");
        }
        return subtreeLeft;
    }
}
