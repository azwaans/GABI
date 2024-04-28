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
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);



        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);// Reference GAPML python version
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_one_thing_happened10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt = new ArrayList<>();
        attempt.add(new TargetDeactTract(0,0));
        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_one_thing_happened0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt = new ArrayList<>();
        attempt.add(new TargetDeactTract(0,0));
        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_two_things_happened10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,40_10_1_1_A,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));
        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));
        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_two_things_happened0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,40_10_1_1_A,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));
        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));
        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus()));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }

    @Test
    public void test_transition_intertarg_two_things0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "30_60_0_2_,110_10_3_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));

        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus()));

    }


    @Test
    public void test_transition_long_intertarget0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));

    }

    @Test
    public void test_transition_long_intertarget1Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
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
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(3, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(3, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));

    }

    @Test
    public void test_transition_long_intertarget10Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_100_0_2_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(3, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(3, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget0Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget1Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
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
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(5, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(7, TransitionWraps.get(1).targetTractsTuples.size(), 1e-5);
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_super_long_intertarget2Steps() {


        // Test for a single branch
        String newick = "(CHILD1:10):0.0;";
        Sequence a = new Sequence("CHILD1", "10_130_0_3_,");
        /*Sequence b = new Sequence("3", "38_38_0_1_,");
        Sequence c = new Sequence("4", "38_1_0_0_,");*/

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="2";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,2));

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(5, TransitionWraps.get(1).numStatuses, 1e-5);// Reference GAPML python version
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus()));


    }


    @Test
    public void test_transition_tree_no_extra() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,None,");
        Sequence b = new Sequence("CHILD2", "10_10_0_0_,40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));
        //checking the internal node wrapper
        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(3).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.get(0).equals(new TargetStatus(attempt1)));

        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);



        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_tree_allow_extra() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_10_0_0_,None,");
        Sequence b = new Sequence("CHILD2", "10_10_0_0_,40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));



        //checking the internal node wrapper
        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(3).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus(attempt1)));



        assertEquals(4, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus(attempt1)));
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus(attempt3)));
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(4, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);

        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt1)));


        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));


    }

    @Test
    public void test_transition_tree_allow_extra_intertarg() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_50_0_2_A,");
        Sequence b = new Sequence("CHILD2", "5_50_0_2_AA,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="10";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));


        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));


        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(3).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus(attempt2)));


        assertEquals(3, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(3, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);

        assertEquals(3, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));





    }

    @Test
    public void test_transition_tree_no_extra_intertarg() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "10_50_0_2_A,");
        Sequence b = new Sequence("CHILD2", "40_10_1_1_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);


        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));


        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));


        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(3).targetTractsTuples.size(), 1e-5);
        assertTrue(TransitionWraps.get(3).transStatuses.get(0).equals(new TargetStatus(attempt2)));


        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(2).transStatuses.contains(new TargetStatus(attempt2)));
        assertEquals(1, TransitionWraps.get(2).targetTractsTuples.size(), 1e-5);

        assertEquals(2, TransitionWraps.get(1).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(1).transStatuses.contains(new TargetStatus(attempt2)));





    }

    @Test
    public void test_transition_start_middle() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "20_6_0_0_,40_5_1_1_,110_10_3_3_,");
        Sequence b = new Sequence("CHILD2", "40_5_1_1_,None,None,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));


        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));
        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertTrue(TransitionWraps.get(3).transStatuses.contains(new TargetStatus()));
        assertEquals(1, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);






    }

    @Test
    public void test_transition_start_cover() {


        // Test for a single branch
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "40_5_1_1_,");
        Sequence b = new Sequence("CHILD2", "20_100_0_3_,");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence",b,"dataType", "user defined");

        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="0";
        RealParameter cutRates = new RealParameter("1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452 1.03371947 1.02624685 1.09459452");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<4;i++) {
            Log.info.println("non root wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println("non root wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in non root wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }
        List<TargetDeactTract> attempt1 = new ArrayList<>();
        attempt1.add(new TargetDeactTract(0,0));

        List<TargetDeactTract> attempt2 = new ArrayList<>();
        attempt2.add(new TargetDeactTract(1,1));

        List<TargetDeactTract> attempt3 = new ArrayList<>();
        attempt3.add(new TargetDeactTract(0,1));


        assertTrue(TransitionWraps.get(4).transStatuses.get(0).equals(new TargetStatus()));
        assertTrue(TransitionWraps.get(3).transStatuses.contains(new TargetStatus()));
        assertEquals(2, TransitionWraps.get(3).numStatuses, 1e-5);
        assertEquals(2, TransitionWraps.get(2).numStatuses, 1e-5);
        assertEquals(1, TransitionWraps.get(1).numStatuses, 1e-5);


    }

    @Test
    public void testPruningREAL() {


        // Test for a single branch
        String newick = "(1:3);";
        Sequence a = new Sequence("1", "38_38_0_1_,82_47_2_3_,145_54_4_5_,");
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
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

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
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);

        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);

        Hashtable<Integer, AncStates> statesDict = TransitionWrap.createStatesDict(tree1,alignment,gestaltModel.metaData.posSites,gestaltModel.metaData.nTargets);

        Hashtable<Integer, TransitionWrap> TransitionWraps = TransitionWrap.createTransitionWraps(tree1,gestaltModel.metaData,statesDict);

        System.out.println("TransitionWraps size: " + TransitionWraps.size() + "\t- Test transition wrapper");

        for(int i=0;i<2;i++) {
            Log.info.println(" wrap ID, is LEAF? T/F:"+ tree1.getNode(i).isLeaf());
            Log.info.println(" wrap ID, is root? T/F:"+ tree1.getNode(i).isRoot());
            TransitionWrap leafWrap = TransitionWraps.get(i + 1);
            Log.info.println(" wrap size" + leafWrap.numStatuses);
            for (TargetStatus stat : leafWrap.transStatuses) {
                Log.info.println("status in wrap" + Arrays.toString(stat.getBinaryStatus(10)));
            }
        }

        assertEquals(2, TransitionWraps.size(), 1e-5);   // Reference GAPML python version
        assertEquals(4, TransitionWraps.get(1).numStatuses, 1e-5);   // Reference GAPML python version

    }


}