import beast.base.core.Log;
import beast.base.evolution.datatype.*;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import gestalt.evolution.alignment.BarcodeMeta;
import gestalt.evolution.alignment.GestaltEvent;
import gestalt.evolution.alignment.IndelSet;
import gestalt.evolution.alignment.TargetStatus;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;
import gestalt.evolution.simulation.SimulatedGestaltAlignment;
import gestalt.evolution.substitutionmodel.gestaltGeneral;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;

import static junit.framework.TestCase.assertEquals;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationTargetTractsonEmptyTargetStatus() {

        //initialise a tree
        Integer sequenceLength = 1;
        String outputFileName = "test/simAl.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("100");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate


        //initialise alignment
        SimulatedGestaltAlignment simAlignment = new SimulatedGestaltAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "outputFileName", outputFileName,
                "userDataType", integerData
        );
        simAlignment.initAndValidate();


        //call function:
        //1 Create an empty target status


        TargetStatus empty = new TargetStatus();

        double arbitraryBranchLength = 3;

        Pair<Double, IndelSet.TargetTract> outcome = simAlignment.raceTargetTracts(empty,arbitraryBranchLength);

        Log.info.println("Time at which the event occured: " + outcome.getFirst());

        Log.info.println("Hash of simulated TT : " + outcome.getSecond().hashCode());
        empty.addTarget(outcome.getSecond());
        Log.info.println("After adding one TT on an empty GESTALT barcode" + Arrays.toString(empty.getBinaryStatus(10)));
        String event = simAlignment.doRepair(outcome.getSecond());





        //AlignmentFromNexus expectedAlignment = new AlignmentFromNexus();
        //expectedAlignment.initByName("fileName",
        //        "test/sciphy/expectedAlignment.nexus",
        //        "userDataType", integerData);


        // assertEquals(expectedAlignment.getSequenceAsString("CHILD1"),
        //         simAlignment.getSequenceAsString("CHILD1"));

        // assertEquals(expectedAlignment.getSequenceAsString("CHILD2"),
        //         simAlignment.getSequenceAsString("CHILD2"));

    }

    // set up
    @Test
    public void testSimulationTargetStatuses() {

        //initialise a tree
        Integer sequenceLength = 1;
        String outputFileName = "test/simAlGESTALT.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:0.0):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        DataType Data = new UserDataType();

        // simulate


        //initialise alignment
        SimulatedGestaltAlignment simAlignment = new SimulatedGestaltAlignment();
        RealParameter originTime = new RealParameter("6");
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "outputFileName", outputFileName,
                "userDataType", Data, "origin", originTime
        );

        Log.info.println("absCutSites");
        Log.info.println(gestaltModel.metaData.absCutSites);
        //Take care of origin vs root!!


        //call function:
        //1 Create an empty target status
        //AlignmentFromNexus expectedAlignment = new AlignmentFromNexus();
        //expectedAlignment.initByName("fileName",
        //        "test/sciphy/expectedAlignment.nexus",
        //        "userDataType", integerData);


        // assertEquals(expectedAlignment.getSequenceAsString("CHILD1"),
        //         simAlignment.getSequenceAsString("CHILD1"));

        // assertEquals(expectedAlignment.getSequenceAsString("CHILD2"),
        //         simAlignment.getSequenceAsString("CHILD2"));

    }
    @Test
    public void testGestaltEventIntersection() {


        //simplest case, test for equality.
        GestaltEvent event1 = new GestaltEvent("38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");
        GestaltEvent event2 = new GestaltEvent("38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");

        GestaltEvent intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");
        Log.info.println("ev2: + " + "38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);


        //contiguous (this is theoretically not possible consecutively because it would require cutting at the same site)
        event1 = new GestaltEvent("38_1_0_0_TGGAGTCGAG");
        event2 = new GestaltEvent("40_3_0_0_AGCGCGCTCGTCGa");

        intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "38_1_0_0_TGGAGTCGAG");
        Log.info.println("ev2: + " + "40_3_0_0_AGCGCGCTCGTCGa");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);

        //overlapping
        event1 = new GestaltEvent("38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");
        event2 = new GestaltEvent("37_2_0_0_ATaa");

        intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "38_1_0_0_TGGAGTCGAGAGCGCGCTCGTCGa");
        Log.info.println("ev2: + " + "37_2_0_0_ATaa");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);




        // masking indels : we expect the outcome to be evt 2
        //151_6_5_5_GAGTTAA,125_85_4_7_TCTGAG,
        event1 = new GestaltEvent("151_6_5_5_GAGTTAA");
        event2 = new GestaltEvent("125_85_4_7_TCTGAG");

        intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "151_6_5_5_GAGTTAA");
        Log.info.println("ev2: + " + "125_85_4_7_TCTGAG");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);

        //totest
        // overlapping contiguous indels
        //97_162_3_9_,71_25_2_3_,
        event1 = new GestaltEvent("97_162_3_9_");
        event2 = new GestaltEvent("71_25_2_3_");

        intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "97_162_3_9_");
        Log.info.println("ev2: + " + "71_25_2_3_");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);

        //178_6_6_6_,149_31_5_6_TGCCAC
        event1 = new GestaltEvent("178_6_6_6_");
        event2 = new GestaltEvent("149_31_5_6_TGCCAC");

        intersection = GestaltEvent.intersect(event1,event2);
        Log.info.println("ev1: + " + "178_6_6_6_");
        Log.info.println("ev2: + " + "149_31_5_6_TGCCAC");
        Log.info.println("intersection of ev1 and ev2");
        Log.info.println("START"+intersection.startPos);
        Log.info.println("MINTARGET"+intersection.minTarg);
        Log.info.println("DELLEN"+intersection.delLen);
        Log.info.println("MAXTARG"+intersection.maxTarg);
        Log.info.println("INSSEQ"+intersection.insSeq);



    }

    String applyIndel(String barcodeSequence,String indel) {

        String[] strs = indel.split("_");
        int leftLen = Integer.parseInt(strs[0]);
        int leftCut = Integer.parseInt(strs[1]);
        String insSeq = strs[2];
        int rightCut = Integer.parseInt(strs[3]);
        int rightLen = Integer.parseInt(strs[4]);


        //first, split into left and right BARCODE strings (no regards to the cutsite
        int cutoffset = 6;
        String[] allele = barcodeSequence.split(" ");

        //left of the cut
        String[] alleleleft = Arrays.copyOfRange(allele, 0, leftCut * 2  +2);
        String leftStrings =  String.join(" ",alleleleft);

        //right of the cut
        String[] alleleright = Arrays.copyOfRange(allele, rightCut * 2 + 1, allele.length );
        String rightStrings = String.join(" ", alleleright);

        //from this, remove the offset for the left side to obtain a string until the cut
        String left_of_left_cut = leftStrings.substring(0,leftStrings.length() - cutoffset);

        //same to the right
        int length_cut_barcode = alleleright[0].length();
        String left_of_right_cut = rightStrings.substring(0, length_cut_barcode - cutoffset);
        String right_of_right_cut  = rightStrings.substring(length_cut_barcode - cutoffset,rightStrings.length());

        //if it is an intertarget deletion, there is also a central sequence
        String central = "";
        if (rightCut != leftCut) {
            String right_of_left_cut = leftStrings.substring(leftStrings.length() - cutoffset,leftStrings.length() );
            String[] allelecenter = Arrays.copyOfRange(allele, leftCut * 2 + 2, rightCut * 2 + 1);
            String centerStrings = String.join(" ", allelecenter);
            central = right_of_left_cut + " " + centerStrings + " " + left_of_right_cut;
        }

        //now create the deletions to the left:
        int deleted_left = 0;
        int index = left_of_left_cut.length();
        while(deleted_left != leftLen) {
            if(! left_of_left_cut.substring(index-1,index).equals(" ") ) {
                left_of_left_cut = left_of_left_cut.substring(0, index - 1) + "-" + left_of_left_cut.substring(index, left_of_left_cut.length());
                deleted_left += 1;
                index -= 1;
            }
            else {
                index -=1;
            }
        }
        Log.info.println(left_of_left_cut);

        int deleted_right = 0;
        index = 0;
        while(deleted_right != rightLen) {
            if(!right_of_right_cut.substring(index,index+1).equals(" ")) {
                right_of_right_cut = right_of_right_cut.substring(0, index) + "-" + right_of_right_cut.substring(index +1, right_of_right_cut.length());
                deleted_right += 1;
                index += 1;
            }
            else {
                index +=1;
            }
        }


        if(central != "") {
            central = central.replace("a", "");
            central = central.replace("c", "");
            central = central.replace("g", "");
            central = central.replace("t", "");
            central = central.replace("A", "-");
            central = central.replace("T", "-");
            central = central.replace("G", "-");
            central = central.replace("C", "-");
        }

        return left_of_left_cut + insSeq + central + right_of_right_cut;


    }


    @Test
    public void testAlleleMakingIntertargetNoInsertShort(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "3_6__7_3,7_1_TATGAC_5_7,5_8_GTCGCGCGCTTT_8_5,4_0_ACGAGATCT_0_4,4_9_CACAGCT_9_4,";

        String indel = "3_6__7_3";

        String finalAllele = applyIndel(barcodeSequence,indel);
        Log.info.println(finalAllele);
        assertEquals("CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCAC--------- ---- --------------------TGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA",finalAllele);
        ;
    }

    @Test
    public void testAlleleMakingIntertargetInsertShort(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "3_6__7_3,7_1_TATGAC_5_7,5_8_GTCGCGCGCTTT_8_5,4_0_ACGAGATCT_0_4,4_9_CACAGCT_9_4,";

        String indel = "7_1_tatgac_5_7";

        String finalAllele = applyIndel(barcodeSequence,indel);
        Log.info.println(finalAllele);
        String expected="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTC-------tatgac------ ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- -GTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";

        assertEquals(expected,finalAllele);
        ;
    }

    @Test
    public void testAlleleMakingFocalInsertShort(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "3_6__7_3,7_1_TATGAC_5_7,5_8_GTCGCGCGCTTT_8_5,4_0_ACGAGATCT_0_4,4_9_CACAGCT_9_4,";

        String indel = "5_8_gtcgcgcgcttt_8_5";

        String finalAllele = applyIndel(barcodeSequence,indel);
        Log.info.println(finalAllele);
        String expected="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTC-----gtcgcgcgcttt-----G AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAA";
        assertEquals(expected,finalAllele);

    }

    @Test
    public void testAlleleMakingFocalNoInsertLongLeft(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "3_6__7_3,7_1_TATGAC_5_7,5_8_GTCGCGCGCTTT_8_5,4_0_ACGAGATCT_0_4,4_9_CACAGCT_9_4,";

        String indel = "25_9__9_4";

        String finalAllele = applyIndel(barcodeSequence,indel);
        Log.info.println(finalAllele);
        String expected="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCT---- ---- ---------------------GG GAAAAAAAAAAAAAAA";
        assertEquals(expected,finalAllele);

    }

    @Test
    public void testAlleleMakingFocalNoInsert(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "3_6__7_3,7_1_TATGAC_5_7,5_8_GTCGCGCGCTTT_8_5,4_0_ACGAGATCT_0_4,4_9_CACAGCT_9_4,";

        String indel = "25_9__9_4";
        String finalAllele = applyIndel(barcodeSequence,indel);
        Log.info.println(finalAllele);
        String expected="CG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCT---- ---- ---------------------GG GAAAAAAAAAAAAAAA";
        assertEquals(expected,finalAllele);

    }

    @Test
    public void testAlleleMakingMultiple(){
        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);

        Log.info.println(gestaltModel.metaData.absCutSites);
        double[] cutsites = gestaltModel.metaData.absCutSites.toArray();
        int[] stringSites = new int[cutsites.length];
        stringSites[0] = (int) cutsites[0] + 1;
        for(int index=1;index<cutsites.length;index++) {
            stringSites[index] = (int) cutsites[index] + 1 + index*2;
        }
        Log.info.println(Arrays.toString(stringSites));

        //apply an indel to the allele:
        String AlleleIndels = "4_4__8_4,4_1__9_4,2_0_tcgagtaa_0_2";
        //iteratively apply indels:
        String finalAllele = barcodeSequence;
        for(String indel : AlleleIndels.split(",")) {
            finalAllele = applyIndel(finalAllele,indel);
        }


        Log.info.println(finalAllele);
        String expected="CG GATACGATACGCGCA--tcgagtaa--ATGG AGTC GACACGACTCGCG---------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ---------------------GG GAAAAAAAAAAAAAAA";
        assertEquals(expected,finalAllele);

    }
    @Test
    public void testAlleleObserve() {

        //initialise a substitution model
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
        gestaltModel.createTrimInsertDistributions(10);
        //initialise the site model with subst model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10.0");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", gestaltModel, "mutationRate", mutationRate);
        //            0           1                 2            3       4         5          6         7          8          9         10         11         12        13          14        15        16          17        18         19         20        21          22        23        24          25        26        27        28
        //////////////01 23456789012345678        901234 5678 90123456789012345678901 2345678901234567890123456789 0123 45678901234567890123456 7890 12345678901234567890123 4567 89012345678901234567890 1234 56789012345678901234567 8901 23456789012345678901234 5678 90123456789012345678901 2345 67890123456789012345678901234567890123456789
        String start="CG GATACGATACGCGCA--tcgagtaa--ATGG AGTC GACACGACTCGCG---------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ----------------------- ---- ---------------------GG GAAAAAAAAAAAAAAA";
        List<String> AlleleEvents = processAllele(start);
        String evt1 = AlleleEvents.get(0);
        String evt2 = AlleleEvents.get(1);
        assertEquals("17_21_4_tcgagtaa",evt1);
        assertEquals("42_266_224_",evt2);



        Log.info.println(AlleleEvents);
        Log.info.println(processEvents(AlleleEvents,gestaltModel.metaData));

    }

    List<String> processEvents(List<String> processedAllele, BarcodeMeta meta) {
        List<String> processedEvents = new ArrayList<>();
        for(String rawEvent : processedAllele) {
            String processedEvent = "";
            List<Integer> matchingTargets = new ArrayList<>();
            String[] splitevent = rawEvent.split("_");
            for(int targetindex= 0; targetindex < meta.absCutSites.length; ++targetindex) {
                double event0 = Double.parseDouble(splitevent[0]);
                double event1 = Double.parseDouble(splitevent[1]);
                if (event0 <=  meta.absCutSites.get(targetindex) && event1 >= meta.absCutSites.get(targetindex)) {
                    Log.info.println("Condition met");
                    matchingTargets.add(targetindex) ;
                }

            }
            //there is an insertion
            if(splitevent.length == 4) {
                processedEvent = splitevent[0] + "_" + splitevent[2] + "_" + matchingTargets.stream().mapToInt(v -> v).min().orElseThrow(NoSuchElementException::new) + "_" + matchingTargets.stream().mapToInt(v -> v).max().orElseThrow(NoSuchElementException::new) + "_" + splitevent[3];
                processedEvents.add(processedEvent);
            }
            //there is no insertion
            else {
                processedEvent = splitevent[0] + "_" + splitevent[2] + "_" + matchingTargets.stream().mapToInt(v -> v).min().orElseThrow(NoSuchElementException::new) + "_" + matchingTargets.stream().mapToInt(v -> v).max().orElseThrow(NoSuchElementException::new) ;
                processedEvents.add(processedEvent);
            }

        }

return processedEvents;

    }

    List<String> processAllele(String Allele) {

        List<String> indelEvents = new ArrayList<>();
        char[] allele = Allele.toCharArray();
        int unedited_index = 0;
        int edited_index = 0;
        while(edited_index != allele.length -1) {
            //spacer position
            if (allele[edited_index] == ' ') {
                edited_index += 1;
            }
            //unedited postion
            if (allele[edited_index] == 'A' || allele[edited_index] == 'T' || allele[edited_index] == 'G' || allele[edited_index] == 'C') {
                edited_index += 1;
                unedited_index +=1;
            }

            //edited position, create an indel
            if (allele[edited_index] == '-') {
                //the start position is inclusive
                int startPos = unedited_index;
                int delLen = 0;
                String insert = "";
                //iterate until end of indel
                while (allele[edited_index] != 'A' && allele[edited_index] != 'T' && allele[edited_index] != 'G' && allele[edited_index] != 'C') {
                    if (allele[edited_index] == '-') {
                        delLen += 1;
                        unedited_index += 1;
                        edited_index += 1;
                    }
                    if (allele[edited_index] == ' ') {
                        edited_index += 1;
                    }
                    if (allele[edited_index] == 'a' || allele[edited_index] == 'c' || allele[edited_index] == 'g' || allele[edited_index] == 't') {
                        insert = insert + allele[edited_index];
                        edited_index += 1;

                    }
                }
                //the end pos is exclusive
                int endPos = unedited_index ;
                String indelEvent = startPos + "_" + endPos + "_" + delLen + "_" + insert;
                indelEvents.add(indelEvent);
            }
        }

        //String eventInput = StartPos + "_" + DelLen + "_" + minDeac + "_" + maxDeac + "_" + insertSequence;
        return indelEvents;
    }

}
