package lineageTree.distributions;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.substitutionmodel.GeneralScarringLoss;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class OTL2_Test {

    organoidTreeLikelihood2 likelihood1, likelihood2, likelihood3;
    Tree tree1, tree2, tree3;
    private organoidTreeLikelihood2 likelihood4;

    @Before
    public void setUp(){

        // init alignment
        Sequence a = new Sequence("0", "0,");
        Sequence b = new Sequence("1", "0,");
        Sequence b2 = new Sequence("1", "1,");
        Sequence b3 = new Sequence("1", "3,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "integer", "stateCount", 3);
        Alignment alignment2 = new Alignment();
        alignment2.initByName("sequence", a, "sequence", b2, "dataType", "integer", "stateCount", 3);
        Alignment alignment3 = new Alignment();
        alignment3.initByName("sequence", a, "sequence", b3, "dataType", "integer", "stateCount", 3);


        //init trees
        tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:28.5,1[&cluster=1]:28.5)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        tree2 = new TreeParser();
        tree2.initByName("IsLabelledNewick", true, "taxa", alignment2, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);
        tree3 = new TreeParser();
        tree3.initByName("IsLabelledNewick", true, "taxa", alignment3, "newick",
                "(0[&cluster=0]:25,1[&cluster=1]:25)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("0.2");
        RealParameter scarRates = new RealParameter("1.0 1.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        GeneralScarringLoss scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", 25.0,
                "scarringDuration", 2.0, "frequencies", frequencies);

        GeneralScarringLoss scarringModel3 = new GeneralScarringLoss();
        scarringModel3.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", 2.0,
                "scarringDuration", 2.0, "frequencies", frequencies);

        GeneralScarringLoss scarringModel4 = new GeneralScarringLoss();
        scarringModel4.initByName("scarringRates", new RealParameter("0.01 0.01"),
                "lossRate", new RealParameter("0.01"),
                "scarringHeight", 100.0,
                "scarringDuration", 100.0, "frequencies", frequencies);

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);
        SiteModel siteM4 = new SiteModel();
        siteM4.initByName("gammaCategoryCount", 0, "substModel", scarringModel4);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihood1 = new organoidTreeLikelihood2();
        likelihood1.initByName("data", alignment, "tree", tree1,
                "siteModel", siteM, "branchRateModel", clockModel);

        likelihood2 = new organoidTreeLikelihood2();
        likelihood2.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM, "branchRateModel", clockModel);

        likelihood3 = new organoidTreeLikelihood2();
        likelihood3.initByName("data", alignment3, "tree", tree3,
                "siteModel", siteM, "branchRateModel", clockModel);

        likelihood4 = new organoidTreeLikelihood2();
        likelihood4.initByName("data", alignment2, "tree", tree2,
                "siteModel", siteM4, "branchRateModel", clockModel);
    }

    @Test
    public void testLikelihood1() {
        // parent above scarring window, children below scarring window

        Node parent = tree1.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        // test correct likelihood calculation for branches that need helper nodes ( that is branches, that
        // cross points where the transition rate matrix changes)

        double[] partials = likelihood1.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 23.0, 25.0}, 1.0, 2);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.0001234098040866, 0, 0, 0}, 1e-15);

        partials = likelihood1.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                true, false);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{3.75566676593833e-9, 0, 0, 0}, 1e-15);

        double logP = likelihood1.calculateLogP();
        assertEquals(Math.log(3.75566676593833e-9), logP, 1e-13);
    }

    @Test
    public void testLikelihood2(){
        // parent within scarring window, children below scarring window
        Node parent = tree2.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        double[] partials = likelihood2.calculatePartialsBeforeParent(parent, child1, 0, 0,
                new double[]{0, 23.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.01005183574463, 0, 0, 0}, 1e-12);

        partials= likelihood2.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                false, true);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{4.081493696794465e-7, 0, 0, 0}, 1e-15);

        double logP = likelihood2.calculateLogP();
        assertEquals(Math.log(4.081493696794465e-7), logP, 1e-13);
    }

    @Test
    public void testLikelihood3(){
        // parent above scarring window, children within scarring window
        Node parent = tree3.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);

        double[] partials = likelihood3.calculatePartialsBeforeParent(parent, child2, 1, 1,
                new double[]{0, 2.0, Double.NEGATIVE_INFINITY}, 1.0, 1);
        assertArrayEquals("Assert correct likelihood at helper node before parent", partials,
                new double[]{0.329679953964361, 0.329679953964361, 0.329679953964361, 1.0}, 1e-12);

        partials= likelihood3.calculatePartialsForCrossBranches(partials, parent, child1, child2,
                false, true);
        assertArrayEquals("Assert correct likelihood at parent", partials,
                new double[]{1.225782753675767e-04, 0, 0, 0}, 1e-15);

        double logP = likelihood3.calculateLogP();
        assertEquals(Math.log(1.225782753675767e-04), logP, 1e-13);

    }

    @Test
    public void testLikelihood4(){
        //test likelihood calculation, when scarring window extends over the entire tree height
        // -> calculation does not require a changing rate matrix nor helper nodes

        // parent above scarring window, children within scarring window
        Node parent = tree3.getRoot();
        Node child1 = parent.getChild(0);
        Node child2 = parent.getChild(1);


        double logP = likelihood4.calculateLogP();
        assertEquals(Math.log(0.072374640511506), logP, 1e-14);

    }
}
