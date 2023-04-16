import beast.base.inference.parameter.RealParameter;
import gestalt.evolution.alignment.*;
import beast.base.evolution.substitutionmodel.Frequencies;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.jblas.DoubleMatrix;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.jblas.MatrixFunctions.powi;
import static org.junit.Assert.assertEquals;



public class gestaltGeneralTest {

    @Test
    public void test_hazard_1() {
        String barcodeSequence="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        String cutSite="6";
        String crucialPos="6 6";
        String maxSumSteps= "3000";
        String maxExtraSteps="1";
        RealParameter cutRates = new RealParameter("0.1 1.1 2.1 3.1 4.1 5.1 6.1 7.1 8.1 9.1");
              RealParameter longTrimScaling = new RealParameter("0.1 0.1");
        RealParameter trimZeroProbs = new RealParameter("0.5 0.5 0.5 0.5 0.5");
        RealParameter trimShortParams = new RealParameter("1.0 1.0");
        RealParameter trimLongParams = new RealParameter("1.0 1.0");
        String insertZeroProb = "0.5";
        RealParameter insertParams = new RealParameter("2.0");
        String doubleCutWeight="0.3";

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();
        //gestaltModel.trim_long_params = trim_long_params;
       // gestaltModel.trim_short_params =trim_short_params;
        //gestaltModel.longTrimFactors = new Double[]{0.05,0.05};
        IndelSet.TargetTract testing = new IndelSet.TargetTract(2,2,2,2);
        //1st, extract the indexc for target tract 2222 in the target tract dict:
        Integer testingIndex = gestaltModel.targetTractDict.get(testing);
        //2nd, extract value from the target tract hazards double matrix:
        Double testingValue = gestaltModel.targetTractHazards.get(testingIndex);
        //lastly, assert
        assertEquals(gestaltModel.cutRates.get(0).getValue(2), testingValue, 1e-5);
        testing = new IndelSet.TargetTract(0,2,0,2);
        //1st, extract the indexc for target tract 0001 in the target tract dict:
        testingIndex = gestaltModel.targetTractDict.get(testing);
        //2nd, extract value from the target tract hazards double matrix:
        testingValue = gestaltModel.targetTractHazards.get(testingIndex);
        //lastly, assert:
        assertEquals(gestaltModel.doubleCutWeight.getValue() * (gestaltModel.cutRates.get(0).getValue(0) + gestaltModel.cutRates.get(0).getValue(2)), testingValue, 1e-5);


    }

    @Test
    public void test_get_hazard_away_1() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest = new TargetDeactTract(0,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        assertEquals(0, hazawaynode, 1e-5);




    }


    @Test
    public void test_get_hazard_away_2() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest1 = new TargetDeactTract(0,1);
        TargetDeactTract tdttest2 = new TargetDeactTract(3,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest1);
        inputlist.add(tdttest2);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        Double myhazard = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(0)) * gestaltModel.cutRates.get(0).getValue(2);
        assertEquals(myhazard, hazawaynode, 1e-5);




    }

    @Test
    public void test_get_hazard_away_3() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest1 = new TargetDeactTract(0,1);
        TargetDeactTract tdttest2 = new TargetDeactTract(4,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest1);
        inputlist.add(tdttest2);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        Double trimlongshortboth = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(0));
        Double sum = gestaltModel.cutRates.get(0).getValue(2) + gestaltModel.cutRates.get(0).getValue(3);
        Double myhazard = (trimlongshortboth * sum) + (trimlongshortboth * sum * gestaltModel.doubleCutWeight.getValue()) ;
        assertEquals(myhazard, hazawaynode, 1e-5);




    }

    @Test
    public void test_get_hazard_away_4() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest1 = new TargetDeactTract(0,1);
        TargetDeactTract tdttest2 = new TargetDeactTract(3,5);
        TargetDeactTract tdttest3 = new TargetDeactTract(7,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest1);
        inputlist.add(tdttest2);
        inputlist.add(tdttest3);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        Double trimlongshortboth = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(1));
        Double sum = gestaltModel.cutRates.get(0).getValue(2) + gestaltModel.cutRates.get(0).getValue(6);
        Double myhazard = (trimlongshortboth * sum) * (1 + gestaltModel.doubleCutWeight.getValue()) ;
        assertEquals(myhazard, hazawaynode, 1e-5);




    }

    @Test
    public void test_get_hazard_away_5() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest1 = new TargetDeactTract(0,1);
        TargetDeactTract tdttest2 = new TargetDeactTract(4,4);
        TargetDeactTract tdttest3 = new TargetDeactTract(7,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest1);
        inputlist.add(tdttest2);
        inputlist.add(tdttest3);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        Double trimlongshortboth = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(0));
        Double sum1 = gestaltModel.cutRates.get(0).getValue(2) + gestaltModel.cutRates.get(0).getValue(3);
        Double sum2 = gestaltModel.cutRates.get(0).getValue(5) + gestaltModel.cutRates.get(0).getValue(6);
        Double myhazard = (trimlongshortboth * (sum1 + sum2)) + (gestaltModel.doubleCutWeight.getValue() * 3 * trimlongshortboth * (sum1 + sum2)) ;
        assertEquals(myhazard, hazawaynode, 1e-5);




    }

    @Test
    public void test_get_hazard_away_6() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();

        List<TargetStatus> testing = new ArrayList<>();
        TargetDeactTract tdttest1 = new TargetDeactTract(1,1);
        TargetDeactTract tdttest2 = new TargetDeactTract(3,9);
        List<TargetDeactTract> inputlist = new ArrayList<>();
        inputlist.add(tdttest1);
        inputlist.add(tdttest2);
        TargetStatus test = new TargetStatus(inputlist);
        testing.add(test);
        Double hazawaynode = gestaltModel.createHazardAwayTargetStatuses(testing).get(0);
        Double trimlongshortboth = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(0));
        Double myhazard = (trimlongshortboth * gestaltModel.cutRates.get(0).getValue(2)) + ((1 + gestaltModel.longTrimFactors.get(0).getValue(1))* gestaltModel.cutRates.get(0).getValue(0)) + (gestaltModel.cutRates.get(0).getValue(0) + gestaltModel.cutRates.get(0).getValue(2)) * gestaltModel.doubleCutWeight.getValue() * (1 + gestaltModel.longTrimFactors.get(0).getValue(1))  ;
        assertEquals(myhazard, hazawaynode, 1e-5);




    }

    @Test
    public void test_create_transition_matrix() {
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

         

        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps, "maxExtraSteps", maxExtraSteps,"cutRates",cutRates,"longTrimScalingFactors",longTrimScaling,"doubleCutWeight",doubleCutWeight,"frequencies",frequencies,"insertZeroProb",insertZeroProb, "trimZeroProbs",trimZeroProbs,"trimShortParams",trimShortParams,"trimLongParams",trimLongParams,"insertParams",insertParams);

        gestaltModel.initAndValidate();
        gestaltModel.hazardAwayDict = gestaltModel.createHazardAwayDict();

        //////////////////test_get_hazard_away 1
        TargetStatus target_stat_start = new TargetStatus();
        ArrayList<TargetDeactTract> temp1 = new ArrayList<>();
        temp1.add(new TargetDeactTract(0,1));
        TargetStatus target_stat1 = new TargetStatus(temp1);
        ArrayList<TargetDeactTract> temp2 = new ArrayList<>();
        temp2.add(new TargetDeactTract(0,3));
        TargetStatus target_stat2 = new TargetStatus(temp2);
        ArrayList<IndelSet> temp3 = new ArrayList<>();
        IndelSet.WildCard wc = new IndelSet.WildCard(0, 3);
        temp3.add(wc);
        AncStates anc_state = new AncStates(temp3);
        List<TargetStatus> transStatuses = new ArrayList<>();
        transStatuses.add(target_stat2);
        transStatuses.add(target_stat1);
        transStatuses.add(target_stat_start);


        List<List<IndelSet.TargetTract>> empty= new ArrayList();
        TransitionWrap transition_wrapper = new TransitionWrap(empty,anc_state,true);
        transition_wrapper.transStatuses=transStatuses;
        transition_wrapper.numStatuses = transStatuses.size();

        transition_wrapper.CreateStatusMap();

        DoubleMatrix qMat = gestaltModel.createRateMatrix(transition_wrapper);
        Double trimlongshortboth = (1 + gestaltModel.longTrimFactors.get(0).getValue(1))* (1 + gestaltModel.longTrimFactors.get(0).getValue(0));
        //FIRST CHECK: row sums = 0

        Double val1 = qMat.get(transition_wrapper.statusMap.get(target_stat_start),transition_wrapper.statusMap.get(target_stat1));
        Double hazard1 = gestaltModel.doubleCutWeight.getValue() * (gestaltModel.cutRates.get(0).getValue(0) + gestaltModel.cutRates.get(0).getValue(1)) + gestaltModel.cutRates.get(0).getValue(0) * gestaltModel.longTrimFactors.get(0).getValue(1) + gestaltModel.cutRates.get(0).getValue(1) * gestaltModel.longTrimFactors.get(0).getValue(0);
        assertEquals(hazard1, val1, 1e-5);


        Double val2 = qMat.get(transition_wrapper.statusMap.get(target_stat_start),transition_wrapper.statusMap.get(target_stat2));
        Double hazard2 = gestaltModel.doubleCutWeight.getValue() * (
                gestaltModel.cutRates.get(0).getValue(0) + gestaltModel.cutRates.get(0).getValue(3)
                 + (gestaltModel.cutRates.get(0).getValue(0) + gestaltModel.cutRates.get(0).getValue(2)) * gestaltModel.longTrimFactors.get(0).getValue(1)
                + (gestaltModel.cutRates.get(0).getValue(1) + gestaltModel.cutRates.get(0).getValue(3)) * gestaltModel.longTrimFactors.get(0).getValue(0)
                 + gestaltModel.longTrimFactors.get(0).getValue(0) * (gestaltModel.cutRates.get(0).getValue(1) + gestaltModel.cutRates.get(0).getValue(2)) * gestaltModel.longTrimFactors.get(0).getValue(1));
        assertEquals(hazard2, val2, 1e-5);

        Double val3 = qMat.get(transition_wrapper.statusMap.get(target_stat1),transition_wrapper.statusMap.get(target_stat2));
        Double hazard3 = (
                gestaltModel.doubleCutWeight.getValue() * (gestaltModel.cutRates.get(0).getValue(2) + gestaltModel.cutRates.get(0).getValue(3)) * (1 + gestaltModel.longTrimFactors.get(0).getValue(0))
                 + gestaltModel.cutRates.get(0).getValue(3) * gestaltModel.longTrimFactors.get(0).getValue(0)
                + gestaltModel.cutRates.get(0).getValue(2) * gestaltModel.longTrimFactors.get(0).getValue(1) * (1 + gestaltModel.longTrimFactors.get(0).getValue(0)));
        assertEquals(hazard3, val3, 1e-5);

        Double val4 = qMat.get(transition_wrapper.statusMap.get(target_stat2),transition_wrapper.statusMap.get(target_stat2));
        Double sumRates4_9 =  gestaltModel.cutRates.get(0).getValue(4) + gestaltModel.cutRates.get(0).getValue(5) + gestaltModel.cutRates.get(0).getValue(6) + gestaltModel.cutRates.get(0).getValue(7) + gestaltModel.cutRates.get(0).getValue(8);
        Double hazard4 = (
                gestaltModel.doubleCutWeight.getValue() * trimlongshortboth * sumRates4_9 * 4
                 + (1 + gestaltModel.longTrimFactors.get(0).getValue(0)) * (1 + gestaltModel.longTrimFactors.get(0).getValue(1)) * sumRates4_9)
                 + (1 + gestaltModel.longTrimFactors.get(0).getValue(0)) * gestaltModel.cutRates.get(0).getValue(9)
                 + gestaltModel.doubleCutWeight.getValue() * (1 + gestaltModel.longTrimFactors.get(0).getValue(0)) * (sumRates4_9 + 5 * gestaltModel.cutRates.get(0).getValue(9)) ;
        assertEquals(-hazard4, val4, 1e-5);

        assertEquals(0.0, qMat.get(transition_wrapper.statusMap.get(target_stat1),transition_wrapper.statusMap.get(target_stat_start)), 1e-5);
        assertEquals(0.0, qMat.get(transition_wrapper.statusMap.get(target_stat2),transition_wrapper.statusMap.get(target_stat_start)), 1e-5);

        DoubleMatrix zeros = new DoubleMatrix(qMat.rows);
        assertEquals(qMat.rowSums(),zeros);
        zeros = new DoubleMatrix(4);
        DoubleMatrix result = ((powi(qMat.getRow(3),2)).eqi(0.0));
        assertEquals(4.0,result.sum(),1e-5);

    }





}
