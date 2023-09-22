package gestalt;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;

import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import gestalt.evolution.alignment.*;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import gestalt.evolution.likelihood.gestaltTreeLikelihood;

import gestalt.evolution.substitutionmodel.gestaltGeneral;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class GestaltPruningTest {

    double runtime;

    @Test
    public void test_transition_nothing_happened() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "None,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);


        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_one_thing_happened10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);
        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);
        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt = new ArrayList<>();
        attempt.add(new TargetDeactTract(0, 0));
        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_one_thing_happened0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt = new ArrayList<>();
        attempt.add(new TargetDeactTract(0, 0));
        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_two_things_happened10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,40_10_1_1_A,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);
        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));
        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));
        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt2)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_two_things_happened0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,40_10_1_1_A,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);
        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));
        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));
        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt2)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_intertarg_two_things0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "30_60_0_2_,110_10_3_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));

        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus()));

    }


    @Test
    public void test_transition_long_intertarget0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));

    }

    @Test
    public void test_transition_long_intertarget1Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(3, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(3, TransitionWraps.get(0).targetTractsTuples.size(), 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));

    }

    @Test
    public void test_transition_long_intertarget10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(3, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(3, TransitionWraps.get(0).targetTractsTuples.size(), 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget1Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(5, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(7, TransitionWraps.get(0).targetTractsTuples.size(), 1e-5);
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget2Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "2";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(5, TransitionWraps.get(0).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.get(0).equals(new TargetStatus()));


    }


    @Test
    public void test_transition_tree_no_extra() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,None,");
        Sequence b = new Sequence("CHILD2", "10_10_0_0_,40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));
        //checking the internal node wrapper
        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus(attempt1)));

        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);
        assertEquals(1, TransitionWraps.get(0).numStatuses, 1e-5);


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_tree_allow_extra() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,None,");
        Sequence b = new Sequence("CHILD2", "10_10_0_0_,40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));


        //checking the internal node wrapper
        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus(attempt1)));


        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(4, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);

        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt1)));


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_tree_allow_extra_intertarg() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_50_0_2_A,");
        Sequence b = new Sequence("CHILD2", "5_50_0_2_AA,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));


        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus(attempt2)));


        assertEquals(3, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(3, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);

        assertEquals(3, TransitionWraps.get(0).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt2)));


    }

    @Test
    public void test_transition_tree_no_extra_intertarg() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_50_0_2_A,");
        Sequence b = new Sequence("CHILD2", "40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);


        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));


        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus(attempt2)));


        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(1, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);

        assertEquals(2, TransitionWraps.get(0).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(0).transStatuses.contains(new TargetStatus(attempt2)));


    }

    @Test
    public void test_transition_start_middle() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "20_6_0_0_,40_5_1_1_,110_10_3_3_,");
        Sequence b = new Sequence("CHILD2", "40_5_1_1_,None,None,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));
        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);
        assertEquals(4, TransitionWraps.get(0).numStatuses, 1e-5);


    }

    @Test
    public void test_transition_start_cover() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "40_5_1_1_,");
        Sequence b = new Sequence("CHILD2", "20_100_0_3_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 4; i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0, 0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1, 1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0, 1));


        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus()));
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);
        assertEquals(1, TransitionWraps.get(0).numStatuses, 1e-5);


    }

    @Test
    public void testPruningREAL() {


        // Test for a single branch
        String newick = "(1:3);";
        Sequence a = new Sequence("1", "38_38_0_1_,82_47_2_3_,145_54_4_5_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "3000";
        String maxExtraSteps = "1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);


        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int i = 0; i < 2; i++) {
            Log.info.println(" wrap ID, is LEAF? T/F:" + tree1.getNode(i).isLeaf());
            Log.info.println(" wrap ID, is root? T/F:" + tree1.getNode(i).isRoot());
            TransitionWrap leafWrap = TransitionWraps.get(i);
            Log.info.println(" wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(0).numStatuses, 1e-5);   // Reference GAPML python version

    }

    @Test
    public void testPruning10() {


        // Test for a single branch
        String newick = "((7:1.0792180222784782,(2:0.49770756595224885,3:0.49770756595224885)17:0.5815104563262293)16:0.06657730305822662,(6:0.9556963887958418,(8:0.8906132932991206,(1:0.40186871143311353,(9:0.16422745979369024,(0:0.05947686862566007,(5:0.008303568372940242,4:0.008303568372940242)10:0.05117330025271982)11:0.10475059116803018)14:0.2376412516394233)15:0.4887445818660071)13:0.0650830954967212)12:0.19009893654086296)18:0.0;";
        Sequence a = new Sequence("0", "38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa,115_144_3_8_,271_11_9_9_a,None,None,None,");
        Sequence b = new Sequence("1", "37_2_0_0_ATaa,64_3_1_1_,114_14_3_3_,141_53_4_5_,197_4_6_6_a,244_34_8_8_,");
        Sequence c = new Sequence("2", "29_10_0_0_aaa,101_204_3_9_,None,None,None,None,");
        Sequence d = new Sequence("3", "36_3_0_0_aaGTATa,66_130_1_5_,251_16_8_8_,None,None,None,");
        Sequence e = new Sequence("4", "33_7_0_0_,63_3_1_1_TATGGAaaa,116_4_3_3_TTATCaaaaTTATGTTATTTGa,144_6_4_4_,167_81_5_7_,253_2_8_8_Caa,");
        Sequence f = new Sequence("5", "33_7_0_0_,58_9_1_1_,118_5_3_3_,128_44_4_4_,179_101_6_8_aa,None,");
        Sequence g = new Sequence("6", "34_12_0_0_,55_46_1_2_,108_98_3_6_,253_5_8_8_,None,None,");
        Sequence h = new Sequence("7", "34_12_0_0_,55_46_1_2_,108_98_3_6_,229_27_8_8_,280_2_9_9_a,None,");
        Sequence i = new Sequence("8", "38_38_0_1_,82_47_2_3_,174_27_5_6_,250_5_8_8_,None,None,");
        Sequence j = new Sequence("9", "34_5_0_0_,61_5_1_1_aaaaa,120_7_3_3_,132_33_4_4_,173_80_5_7_a,None,");


        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "sequence", f, "sequence", g, "sequence", h, "sequence", i, "sequence", j, "dataType", "user defined");

        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite = "6";
        String crucialPos = "6 6";
        String maxSumSteps = "20";
        String maxExtraSteps = "1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
        RealParameter longTrimScaling = new RealParameter("0.04 0.04");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("3.0 3.0");
        RealParameter trimLongParams = new RealParameter("4.0 4.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight = "0.3";

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);

        likelihood.populateStatesDict(tree1.getRoot());
        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1, gestaltModel.metaData, likelihood.statesDict, likelihood.currentStatesDictIndex);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for (int w = 0; w < 2; w++) {
            Log.info.println(" wrap ID, is LEAF? T/F:" + tree1.getNode(w).isLeaf());
            Log.info.println(" wrap ID, is root? T/F:" + tree1.getNode(w).isRoot());
            TransitionWrap leafWrap = TransitionWraps.get(w);
            Log.info.println(" wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }

        assertEquals(19, TransitionWraps.size(), 1e-5);   // Reference GAPML python version

        double logL = likelihood.calculateLogP();
        //
        // -1080.7329181626874
        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");
        assertEquals(logL, -1080.7329181626874, 1e-5);
    }

//    @Test
//    public void testPruning10_2() {
//
//
//        // Test for real data 10 leaves 2.7 to 3sec is baseline for the
//        String newick = "(5:2.17923259589162,((7:1.0792180222784782,(2:0.49770756595224885,3:0.49770756595224885)17:0.5815104563262293)16:0.06657730305822662,(6:0.9556963887958418,(8:0.8906132932991206,(1:0.40186871143311353,(9:0.16422745979369024,(0:0.05947686862566007,4:0.05947686862566007)11:0.10475059116803018)14:0.2376412516394233)15:0.4887445818660071)13:0.0650830954967212)12:0.19009893654086296)10:1.033437270554915)18:0.0;";
//        Sequence a = new Sequence("0", "38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa,115_144_3_8_,271_11_9_9_a,None,None,None,");
//        Sequence b = new Sequence("1", "37_2_0_0_ATaa,64_3_1_1_,114_14_3_3_,141_53_4_5_,197_4_6_6_a,244_34_8_8_,");
//        Sequence c = new Sequence("2", "29_10_0_0_aaa,101_204_3_9_,None,None,None,None,");
//        Sequence d = new Sequence("3", "36_3_0_0_aaGTATa,66_130_1_5_,251_16_8_8_,None,None,None,");
//        Sequence e = new Sequence("4", "33_7_0_0_,63_3_1_1_TATGGAaaa,116_4_3_3_TTATCaaaaTTATGTTATTTGa,144_6_4_4_,167_81_5_7_,253_2_8_8_Caa,");
//        Sequence f = new Sequence("5", "33_7_0_0_,58_9_1_1_,118_5_3_3_,128_44_4_4_,179_101_6_8_aa,None,");
//        Sequence g = new Sequence("6", "34_12_0_0_,55_46_1_2_,108_98_3_6_,253_5_8_8_,None,None,");
//        Sequence h = new Sequence("7", "34_12_0_0_,55_46_1_2_,108_98_3_6_,229_27_8_8_,280_2_9_9_a,None,");
//        Sequence i = new Sequence("8", "38_38_0_1_,82_47_2_3_,174_27_5_6_,250_5_8_8_,None,None,");
//        Sequence j = new Sequence("9", "34_5_0_0_,61_5_1_1_aaaaa,120_7_3_3_,132_33_4_4_,173_80_5_7_a,None,");
//
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "sequence", b, "sequence", c, "sequence", d, "sequence", e, "sequence", f, "sequence", g, "sequence", h, "sequence", i, "sequence", j, "dataType", "user defined");
//
//        String barcodeSequence = "CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
//        String cutSite = "6";
//        String crucialPos = "6 6";
//
//        //max sum steps is important. Some states will not be reached with too low of a maxsumsteps
//        String maxSumSteps = "40";
//        String maxExtraSteps = "1";
//        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
//        RealParameter longTrimScaling = new RealParameter("0.04 0.04");
//        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
//        RealParameter trimShortParams = new RealParameter("3.0 3.0");
//        RealParameter trimLongParams = new RealParameter("4.0 4.0");
//        String insertZeroProb = "0.5";
//        RealParameter insertParams = new RealParameter("2.0");
//        String doubleCutWeight = "0.3";
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        gestaltGeneral gestaltModel = new gestaltGeneral();
//        RealParameter freqs = new RealParameter("1.0 0 0");
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", freqs,
//                "estimate", false);
//        gestaltModel.initByName("barcodeSequence", barcodeSequence,
//                "cutSite", cutSite,
//                "crucialPos", crucialPos,
//                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps, "cutRates", cutRates, "longTrimScalingFactors", longTrimScaling, "doubleCutWeight", doubleCutWeight, "frequencies", frequencies, "insertZeroProb", insertZeroProb, "trimZeroProbs", trimZeroProbs, "trimShortParams", trimShortParams, "trimLongParams", trimLongParams, "insertParams", insertParams);
//
//        SiteModel siteM = new SiteModel();
//        siteM.initByName("gammaCategoryCount", 0, "substModel", gestaltModel);
//
//        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
//        likelihood.initByName("data", alignment, "tree", tree1, "siteModel", siteM);
//
//        //recording time
//        long start = System.nanoTime();
//        double logL = likelihood.calculateLogP();
//        long end = System.nanoTime();
//
//
//        for (int w = 0; w < 2; w++) {
//            Log.info.println(" wrap ID, is LEAF? T/F:" + tree1.getNode(w).isLeaf());
//            Log.info.println(" wrap ID, is root? T/F:" + tree1.getNode(w).isRoot());
//            TransitionWrap leafWrap = likelihood..get(w);
//            Log.info.println(" wrap size" + leafWrap.numStatuses);
//            for (TargetStatus stat : leafWrap.transStatuses) {
//                Log.info.println("status in wrap" + Arrays.toString(stat.getBinaryStatus(10)));
//            }
//        }
//
//        assertEquals(19, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
//
//
//        //
//        // -1080.7329181626874
//        System.out.println("gestaltLikelihood: " + logL + "\t- Test LikelihoodNothingHappens");
//
//    }


}