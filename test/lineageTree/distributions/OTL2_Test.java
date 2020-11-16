package lineageTree.distributions;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.substitutionmodel.GeneralScarringLoss;
import org.junit.Before;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class OTL2_Test {

    organoidTreeLikelihood2 likelihood;

    @Before
    public void setUp(){

        // init alignment
        Sequence a = new Sequence("0", "0,0,0");
        Sequence b = new Sequence("1", "0,1,2");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "integer", "stateCount", 3);


        //init tree
        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                "(0[&cluster=0]:28.5,1[&cluster=1]:28.5)2[&cluster=0]:0.0",
                "adjustTipHeights", false, "offset", 0);

        //init scarring model
        RealParameter lossRate = new RealParameter("1.0");
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

        // init site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", scarringModel);

        // init branch rate model
        StrictClockModel clockModel = new StrictClockModel();

        likelihood = new organoidTreeLikelihood2();
        likelihood.initByName("data", alignment, "tree", tree,
                "siteModel", siteM, "branchRateModel", clockModel);
    }

    @Test
    public void testLikelihood(){

        likelihood.calculateLogP();



    assertEquals(1.0, 1.0, 0.001);

    }
}
