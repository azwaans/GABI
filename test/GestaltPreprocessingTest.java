
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import gestalt.evolution.likelihood.gestaltTreeLikelihood;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.junit.Test;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

public class GestaltPreprocessingTest {

    @Test
    public void testAlignmentConversionGSMInput() throws IOException {
        // test the string conversion of sequences in the GSM format in GABI vs GAPML
        //metadata for V7 bcode with 20 bp padding and parameters
        String barcodeSequence="AAAAAAAAAAAAAAAAAAAACG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
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

        //substitution model
        gestaltGeneral gestaltModel = new gestaltGeneral();
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        gestaltModel.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps,
                "maxExtraSteps", maxExtraSteps,
                "cutRates",cutRates,
                "longTrimScalingFactors",longTrimScaling,
                "doubleCutWeight",doubleCutWeight,
                "frequencies",frequencies,
                "insertZeroProb",insertZeroProb,
                "trimZeroProbs",trimZeroProbs,
                "trimShortParams",trimShortParams,
                "trimLongParams",trimLongParams,
                "insertParams",insertParams);

        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", gestaltModel);
        //dummy tree and sequence for initialization
        String newick = "(1:1);";
        Sequence a = new Sequence("1", "None,None,");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"dataType", "user defined");
        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        //likelihood
        gestaltTreeLikelihood likelihood = new gestaltTreeLikelihood();
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);

        //read in processed alignment from GAPML
        FileReader input = new FileReader("test/processed_data_GSM2171788.txt");
        BufferedReader br = new BufferedReader(input);
        //ignore header
        br.readLine();

        //compare pre-processing output
        for(int index=1;index<361;index++) {
            String line = br.readLine();

            // split reads and processed sequence
            String[] strs = line.split(",");
            String unprocessed_read = strs[0];
            String unprocessed_read_comma = unprocessed_read.replace("_",",");
            String gapml_processed = strs[1];
            String gapml_processed_comma = gapml_processed.replace("=",",") + ",";
            String gabi_processed = likelihood.processSequence(unprocessed_read_comma);

            //empty sequences in gapml are encoded as "no_evts," in GAPML and "None," in GABI
            if(gapml_processed_comma.equals("no_evts,")) {
                assertEquals("None,",gabi_processed);
            }
            else {
                assertEquals(gapml_processed_comma, gabi_processed);
            }

        }
    }


    @Test
    public void testAlignmentConversionGSMInputLikelihood() throws IOException {
        // test that the likelihood calculation is unaffected by the GSM sequence conversion in GABI
        //metadata for V7 bcode with 20 bp padding and parameters
        String barcodeSequence="AAAAAAAAAAAAAAAAAAAACG GATACGATACGCGCACGCTATGG AGTC GACACGACTCGCGCATACGATGG AGTC GATAGTATGCGTATACGCTATGG AGTC GATATGCATAGCGCATGCTATGG AGTC GAGTCGAGACGCTGACGATATGG AGTC GCTACGATACACTCTGACTATGG AGTC GCGACTGTACGCACACGCGATGG AGTC GATACGTAGCACGCAGACTATGG AGTC GACACAGTACTCTCACTCTATGG AGTC GATATGAGACTCGCATGTGATGG GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
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

        //substitution model specifying the GSM input
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        gestaltGeneral gestaltModelGSMInput = new gestaltGeneral();
        gestaltModelGSMInput.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps,
                "maxExtraSteps", maxExtraSteps,
                "cutRates",cutRates,
                "longTrimScalingFactors", longTrimScaling,
                "doubleCutWeight",doubleCutWeight,
                "frequencies",frequencies,
                "insertZeroProb",insertZeroProb,
                "trimZeroProbs",trimZeroProbs,
                "trimShortParams",trimShortParams,
                "trimLongParams",trimLongParams,
                "insertParams",insertParams,
                "inputFormatGSM",true);

        SiteModel siteMGSMInput = new SiteModel();
        siteMGSMInput.initByName( "gammaCategoryCount", 0, "substModel", gestaltModelGSMInput);

        //substitution model specifying the Event input
        gestaltGeneral gestaltModelEventInput = new gestaltGeneral();
        gestaltModelEventInput.initByName("barcodeSequence", barcodeSequence,
                "cutSite", cutSite,
                "crucialPos", crucialPos,
                "maxSumSteps", maxSumSteps,
                "maxExtraSteps", maxExtraSteps,
                "cutRates",cutRates,
                "longTrimScalingFactors", longTrimScaling,
                "doubleCutWeight",doubleCutWeight,
                "frequencies",frequencies,
                "insertZeroProb",insertZeroProb,
                "trimZeroProbs",trimZeroProbs,
                "trimShortParams",trimShortParams,
                "trimLongParams",trimLongParams,
                "insertParams",insertParams,
                "inputFormatGSM",false);

        SiteModel siteMEventInput = new SiteModel();
        siteMEventInput.initByName( "gammaCategoryCount", 0, "substModel", gestaltModelEventInput);


        //read in processed alignment from GAPML
        FileReader input = new FileReader("test/processed_data_GSM2171788.txt");
        BufferedReader br = new BufferedReader(input);
        //ignore header
        br.readLine();

        //compare pre-processing output
        //using the first 40 sequences for testing
        for(int index=1;index<40;index++) {
            String line = br.readLine();

            // split reads and processed sequence
            String[] strs = line.split(",");
            String unprocessed_read = strs[0];
            String unprocessed_read_comma = unprocessed_read.replace("_",",");
            String gapml_processed = strs[1];
            String gapml_processed_comma = gapml_processed.replace("=",",") + ",";

            // CREATE GSM LIKELIHOOD //
            //create alignment
            Sequence GSM = new Sequence("TAXON1", unprocessed_read_comma);
            Alignment alignmentGSM = new Alignment();
            alignmentGSM.initByName("sequence", GSM,"dataType", "user defined");

            //create tree with alignment
            String newick = "(TAXON1:5);";
            Tree treeGSM = new TreeParser();
            treeGSM.initByName("IsLabelledNewick", true, "taxa", alignmentGSM, "newick",
                    newick,
                    "adjustTipHeights", false, "offset", 0);

            //calculate likelihood with the GSM
            gestaltTreeLikelihood likelihoodGSM = new gestaltTreeLikelihood();
            likelihoodGSM.initByName("data",alignmentGSM,"tree",treeGSM,"siteModel",siteMGSMInput);

            double logLGSM = likelihoodGSM.calculateLogP();

            // CREATE EVENT LIKELIHOOD //
            //create alignment
            Sequence eventFormat = new Sequence("TAXON1", gapml_processed_comma);
            Alignment alignmentEvent = new Alignment();
            alignmentEvent.initByName("sequence", eventFormat,"dataType", "user defined");

            //create tree with alignment
            Tree treeEvent = new TreeParser();
            treeEvent.initByName("IsLabelledNewick", true, "taxa", alignmentEvent, "newick",
                    newick,
                    "adjustTipHeights", false, "offset", 0);

            //calculate likelihood with the GSM
            gestaltTreeLikelihood likelihoodEvent = new gestaltTreeLikelihood();
            likelihoodEvent.initByName("data",alignmentEvent,"tree",treeEvent,"siteModel",siteMEventInput);

            double logLEvent = likelihoodEvent.calculateLogP();

            //assess equality of likelihood using the preprocessed sequences in Event format and GSM format
            assertEquals(logLEvent, logLGSM, 1e-8);
        }
    }




}
