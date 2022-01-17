package lineageTree.distributions;

import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.TargetStatus;
import beast.evolution.alignment.TransitionWrap;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import lineageTree.substitutionmodel.GeneralGestalt;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Hashtable;

import static org.junit.Assert.assertEquals;

public class GestaltLikelihoodTest {

    double runtime;

    /**
     * Basic test for the likelihood calculation with nothing sequenced
     * @throws Exception
     */
    @Test
    public void testLikelihoodSingleBranchNothing() {

        // Test for a single branch with no indels
        String newick = "(CHILD1:0.1):0.0;";

        Sequence a = new Sequence("CHILD1", "None,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-0.5602110979000002, logL, 1e-5);   // Reference GAPML python version


    }

    /**
     * Basic test for the likelihood calculation with a single indel
     * @throws Exception
     */
    @Test
    public void testLikelihoodSingleBranch() {

        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";

        Sequence a = new Sequence("CHILD1", "6_3_0_0_,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-37.798312496118861, logL, 1e-5);   // Reference GAPML python version


    }


    /**
     * Basic test for the likelihood calculation with a single indel
     * @throws Exception
     */
    @Test
    public void testLikelihoodSingleBranchIntertarget() {

        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";

        Sequence a = new Sequence("CHILD1", "6_28_0_2_ATC,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-11.804050, logL, 1e-5);   // Reference GAPML python version


    }


    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodDoubleBranchIntertarget() {

        // Test for a single branch
        String newick = "(CHILD2:5):10.0;";

        Sequence a = new Sequence("CHILD2", "6_28_0_2_ATC,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-11.804050, logL, 1e-5);   // Reference GAPML python version


    }


    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodTreeThreeLeaves() {

        // Test for a single branch
        String newick = "(((CHILD1:1,CHILD3:1)INTERNAL:1,CHILD2:2):2)0.0;";

        Sequence a = new Sequence("CHILD1", "None,");
        Sequence b = new Sequence("CHILD3", "6_10_0_0_,");
        Sequence c = new Sequence("CHILD2", "6_3_0_0_,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-49.686965041965479, logL, 1e-5);   // Reference GAPML python version


    }

    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodTree() {

        // Test for a single branch
        String newick = "(CHILD0:3,((CHILD4:1,COPY01:1)CHILD2:1,CHILD3:2)CHILD1:1);";
        Sequence a = new Sequence("CHILD0", "6_28_0_2_ATC,None,");
        Sequence b = new Sequence("CHILD3", "16_3_1_1_,6_3_0_0_,");
        Sequence c = new Sequence("CHILD4", "16_3_1_1_,6_10_0_0_,");
        Sequence d = new Sequence("COPY01", "16_3_1_1_,None,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"sequence", d,"dataType", "user defined");

        String barcodeSequence="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String cutSite="3";
        String crucialPos="3 3";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-46.02776497, logL, 1e-5);   // Reference GAPML python version


    }


    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodRealData() {


        // Test for a single branch
        String newick = "((3:4,4:4)1:1,1:5);";
        Sequence a = new Sequence("1", "36_3_0_0_aaGTATa,66_130_1_5_,251_16_8_8_,None,None,");
        Sequence b = new Sequence("3", "38_38_0_1_,82_47_2_3_,145_54_4_5_,251_11_8_8_,264_18_9_9_aa,");
        Sequence c = new Sequence("4", "38_1_0_0_,67_7_1_1_,115_144_3_8_,271_11_9_9_a,None,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-315.79117403, logL, 1e-5);   // Reference GAPML python version


    }

    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodRealData_1_branchEvents3132() {


        // Test for a single branch
        // events 3_1 and 3_2

        String newick = "(1:3);";
        Sequence a = new Sequence("1", "38_38_0_1_,82_47_2_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-69.70838033, logL, 1e-5);   // Reference GAPML python version


    }

    @Test
    public void testLikelihoodRealData_1_branchEventsTOGETHER() {


        // Test for a single branch
        // EVENTS 1_1 amd 1_2
        String newick = "(1:3);";
        Sequence a = new Sequence("1", "55_46_1_2_,118_5_3_3_,");
        /*Sequence b = new Sequence("3", "55_46_1_2,");
        Sequence c = new Sequence("4", ",");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-78.69923675, logL, 1e-5);   // Reference GAPML python version


    }

    @Test
    public void testLikelihoodRealData_1_branchEventsTOGETHER2() {


        // Test for a single branch
        // EVENTS 1_1 amd 1_2
        String newick = "(1:3);";
        Sequence a = new Sequence("1", "55_46_1_2_,108_98_3_6_,");
        /*Sequence b = new Sequence("3", "55_46_1_2,");
        Sequence c = new Sequence("4", ",");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-40.95627345, logL, 1e-5);   // Reference GAPML python version


    }


    @Test
    public void testLikelihoodRealData_1_branchEventsTOGETHER3() {


        // Test for a single branch
        // EVENTS 1_1 amd 1_2
        String newick = "(1:3);";
        Sequence a = new Sequence("1", "144_6_4_4_,167_81_5_7_");
        /*Sequence b = new Sequence("3", "55_46_1_2,");
        Sequence c = new Sequence("4", ",");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-81.08395031, logL, 1e-5);   // Reference GAPML python version


    }

    @Test
    public void testLikelihoodRealData_1_branchEventsALONE () {


        // Test for a single branch
        // EVENTS 1_1 amd 1_2
        String newick = "(1:3):0.0;";

        Sequence a = new Sequence("1", "36_3_0_0_aaGTATA,");
        /*Sequence b = new Sequence("3", "55_46_1_2,");
        Sequence c = new Sequence("4", ",");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-146.3612818, logL, 1e-5);   // Reference GAPML python version

        a = new Sequence("1", "66_130_1_5_,");
        alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");
        assertEquals(-62.48971707, logL, 1e-5);

        // Test for a single branch
        // EVENTS 1_1 amd 1_2

    }

    /**
     * Basic test for the likelihood calculation with a single indel, 2 branches
     * @throws Exception
     */
    @Test
    public void testLikelihoodRealData_1_branch2() {


        // Test for a single branch
        String newick = "(1:4):0.0;";
        Sequence a = new Sequence("1", "38_38_0_1_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-113.23692153, logL, 1e-5);   // Reference GAPML python version


    }

    @Test
    public void testLikelihoodRealData_1_branch3() {


        // Test for a single branch
        String newick = "(1:3):0.0;";
        Sequence a = new Sequence("1", "38_38_0_1,82_47_2_3_,145_54_4_5_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        String longTrimScaling="0.05 0.05";
        String doubleCutWeight="0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        GeneralGestalt gestaltModel = new GeneralGestalt();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScaling",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.substitutionModel = gestaltModel;
        double logL = likelihood.calculateLogP();

        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");

        assertEquals(-87.24370139, logL, 1e-5);   // Reference GAPML python version


    }





}
