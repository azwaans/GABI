package gestalt.evolution.alignment;

import beast.base.core.BEASTObject;
import org.apache.commons.math3.util.Pair;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.max;
import static java.lang.Math.min;

public class BarcodeMeta extends BEASTObject {

    /** Barcode metadata:
     *  User input of barcode design features and analysis metadata
     */

    /**
     * unedited barcode sequence, as a list of target sequences. Positions are indexed starting from 0.
     */
    protected List<String> uneditedBarcode;

    /**
     * number of targets
     */
    public int nTargets;

    /**
     * offset from 3' end of target for Cas9 cutting,
     * so a cut_site of 6 means that we start inserting
     * such that the inserted seq is 6 nucleotides from
     * the 3' end of the target
     */
    public int barcodeCutSite;


    /**
     * which positions to the left and right of the
     * cut site must not be disturbed for the target to
     * remain active
     */
    protected int[] crucialPosLen;


    public int maxSumSteps = 3000;
    public int maxExtraSteps = 1;

    //parameters computed with inputs:

    /**
     * Min length of a long trim for target i -- left
     */
    public List<Integer> leftLongTrimMin;

    /**
     * Min length of a long trim for target i -- right
     */
    public List<Integer> rightLongTrimMin;

    /**
     * Max length of any trim for target i -- left
     */
    public List<Integer> leftMaxTrim;

    /**
     * Max length of any trim for target i -- right
     */
    public List<Integer> rightMaxTrim;

    /**
     * absolute positions of cut locations
     */
    public DoubleMatrix absCutSites;


    /**
     * Range of positions for each target
     * regarding which positions must be unedited
     * for this target to still be active.
     */
    public List<Pair<Integer, Integer>> posSites;


    @Override
    public void initAndValidate() {

    }

    public BarcodeMeta(List<String> uneditedBcode, int cutSte, int[] crucialPsLn, int maxSumStps, int maxExtraStps) {
        uneditedBarcode = uneditedBcode;
        barcodeCutSite = cutSte;
        crucialPosLen = crucialPsLn;
        maxSumSteps = maxSumStps;
        maxExtraSteps = maxExtraStps;

        List<Integer> barcodeSubstringsLengths = new ArrayList<>();
        Integer origLength = 0;
        for (String str : uneditedBarcode) {
            barcodeSubstringsLengths.add(str.length());
            origLength = origLength + str.length();
        }
        nTargets = Math.floorDiv(uneditedBarcode.size() - 1, 2);

        absCutSites = new DoubleMatrix(nTargets);
        int[] cutSites = new int[nTargets];
        Arrays.fill(cutSites, barcodeCutSite);
        for (int i = 0; i < nTargets; i++) {

            List<Integer> temp = barcodeSubstringsLengths.subList(0, 2 * (i + 1));
            Integer sum = 0;
            for (Integer j : temp) {
                sum = sum + j;
            }
            absCutSites.put(i, sum - cutSites[i]);
        }
        posSites = new ArrayList<>();
        for (int i = 0; i < nTargets; i++) {
            int cutSite = (int) absCutSites.get(i);
            int right = min(origLength - 1, cutSite + crucialPosLen[1] - 1);
            int left = max(0, cutSite - crucialPosLen[0]);
            posSites.add(new Pair(left, right));
            assert left <= cutSite;
            assert right >= cutSite;
        }


        List<Integer> rightLongTrimMn = new ArrayList<>();
        for (int i = 0; i < nTargets - 1; ++i) {
            int site = (int) absCutSites.get(i);
            rightLongTrimMn.add(posSites.get(i + 1).getFirst() - site + 1);
        }
        rightLongTrimMin = rightLongTrimMn;

        List<Integer> leftMaxTrm = new ArrayList<>();
        leftMaxTrm.add((int) absCutSites.get(0));
        for (int i = 1; i < nTargets; ++i) {
            leftMaxTrm.add((int) absCutSites.get(i) - (int) absCutSites.get(i - 1) - 1);
        }

        leftMaxTrim = leftMaxTrm;

        List<Integer> leftLongTrimMn = new ArrayList<>();
        leftLongTrimMn.add(leftMaxTrm.get(0));
        for (int i = 1; i < nTargets; ++i) {
            int site = (int) absCutSites.get(i);
            leftLongTrimMn.add(site - posSites.get(i - 1).getSecond());
        }

        leftLongTrimMin = leftLongTrimMn;

        List<Integer> rightMaxTrm = new ArrayList<>();
        for (int i = 0; i < nTargets - 1; i++) {
            rightMaxTrm.add(((int) absCutSites.get(i + 1) - (int) absCutSites.get(i) - 1));
        }
        rightMaxTrm.add(origLength - (int) absCutSites.get(absCutSites.length - 1) - 1);
       /* rightMaxTrm.add((int) (origLength-absCutSites.get(absCutSites.length-1) -1));
        self.right_max_trim += [self.orig_length - self.abs_cut_sites[-1] - 1]*/

        rightLongTrimMn.add(rightMaxTrm.get(rightMaxTrm.size() - 1));

        rightMaxTrim = rightMaxTrm;


    }
}
