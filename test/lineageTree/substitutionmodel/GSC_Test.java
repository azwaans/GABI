package lineageTree.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.DoubleStream;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;

public class GSC_Test {

    GeneralScarringLoss scarringModel;

    @Before
    public void setup(){

        RealParameter lossRate = new RealParameter("1.0");
        RealParameter scarRates = new RealParameter("1.0 1.0");

        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);

        scarringModel = new GeneralScarringLoss();
        scarringModel.initByName("scarringRates", scarRates,
                "lossRate", lossRate,
                "scarringHeight", 25.0,
                "scarringDuration", 2.0, "frequencies", frequencies);
    }

    @Test
    public void testRateMatrix(){

        //System.out.println("relative rates :\n" +
        //        Arrays.toString(scarringModel.getRelativeRates()) + "\n");
        System.out.println("rate matrix :");
        double[][] rateM = scarringModel.getRateMatrix();
        for(int i = 0; i < rateM.length; i++)
            System.out.println(Arrays.toString(rateM[i]));


        double[][] correctMatrix = new double[rateM.length][rateM.length];
        correctMatrix[0][0] = 0;
        for (int i=1; i<rateM.length-2; i++){
            // insert scar rates into first row
            correctMatrix[0][i] = scarringModel.scarRates[i-1];

            // insert loss rates into last column
            correctMatrix[i-1][rateM.length-1] = scarringModel.lossRate;
        }
        correctMatrix[rateM.length-2][rateM.length-1] = scarringModel.lossRate;

        for (int i=0; i<rateM.length; i++){
            assertArrayEquals("Assert matrix entries", correctMatrix[i], rateM[i], 1e-15);
        }

    }

    @Test
    public void testTransitionProbabilities1(){
        /* Test transition probabilities for scarring and loss window
        */
        double startTime = 25;
        double endTime = 23;
        double rate = 1.0;

        System.out.println("freqs = \n" + Arrays.toString(scarringModel.getFrequencies()) + "\n");

        int len = scarringModel.getStateCount();
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        double[] correctMatrix = new double[]{
                0.002478752176666, 0.066428265529973, 0.066428265529973, 0.864664716763387,
                0, 0.135335283236613, 0, 0.864664716763387,
                0, 0, 0.135335283236613, 0.864664716763387,
                0, 0, 0, 1.0
        };

        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities",
                correctMatrix, prob, 1e-15);

    }

    @Test
    public void testTransitionProbabilities2(){
        /* Test transition probabilities for loss window
         */
        double startTime = 2;
        double endTime = 0;
        double rate = 1.0;

        System.out.println("freqs = \n" + Arrays.toString(scarringModel.getFrequencies()) + "\n");

        int len = scarringModel.getStateCount();
        double[] prob = new double[len*len];
        scarringModel.getTransitionProbabilities(new Node(), startTime, endTime, rate, prob);

        double[] correctMatrix = new double[]{
                0.135335283236613, 0, 0, 0.864664716763387,
                0, 0.135335283236613, 0, 0.864664716763387,
                0, 0, 0.135335283236613, 0.864664716763387,
                0, 0, 0, 1.0
        };

        System.out.println("\ntransition prob :\n" + Arrays.toString(prob));
        // P(t) row sum to 1
        for (int i=0; i < len; i++) {
            double[] row = new double[len];
            System.arraycopy(prob, i*len, row, 0, len);
            double sum = DoubleStream.of(row).sum();
            System.out.println("row " + i + " prob sum = " + sum);
            assertEquals("Assert row sum == 1",1, sum, 1e-15);
        }

        assertArrayEquals("Assert correct transition probabilities for loss",
                correctMatrix, prob, 0.01);

    }



}
