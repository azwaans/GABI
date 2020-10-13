package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GeneralScarringLoss extends SubstitutionModel.Base {

        final public Input<List<Double>> scarringRatesInput = new Input<>("scarringRates",
                "Rates, at which scars types are introduced into the scarring region",
                new ArrayList<>(), Input.Validate.REQUIRED);
        final public Input<Double> lossRateInput = new Input<>("lossRate",
                "Rate, at which scarring regions are lost");
        public Input<Double> scarringHeightInput = new Input<>("scarringHeight",
                "Duration between the onset of scarring and sampling of the cells");
        public Input<Double> scarringDurationInput = new Input<>("scarringDuration",
                "Duration of the scarring process");


        /**
         * Eigenvalue decomposition of rate matrix + its stored version *
         */
        private EigenDecomposition getEigenDecompositionScarring = null;
        private EigenDecomposition storedEigenDecompositionScarring = null;
        private EigenDecomposition getEigenDecompositionLoss = null;
        private EigenDecomposition storedEigenDecompositionLoss = null;
        /**
         * flag to indicate eigen decomposition is up to date *
         */
        private boolean updateEigenScar = true;
        private boolean updateEigenLoss = true;

        /**
         * flag to indicate matrix is up to date *
         */
        protected boolean updateMatrixScar = true;
        protected boolean updateMatrixLoss = true;


        double scarringHeight;
        double scarringDuration;
        double[] frequencies;
        //potentially delete
        double[][] rateMatrix;
        double[] scarRates;
        double lossRate;

    @Override
    public void initAndValidate() {

        List<Double> scarringRates = scarringRatesInput.get();

        // one state for each scar type + unedited + lost
        nrOfStates = scarringRates.size() + 2;
        rateMatrix = new double[nrOfStates][nrOfStates];

        scarRates = new double[nrOfStates-2];

        // assert positive rates; fill rate matrix
        double sumScarringRates = 0;
        for (int i =0; i <scarringRates.size(); i++){

            Double scarringRate = scarringRates.get(i);
            if (scarringRate <= 0) {
                throw new RuntimeException("All scarring rates must be positive!");
            }
            sumScarringRates += scarringRate;
            rateMatrix[0][i+1] = scarringRate;
            scarRates[i] = scarringRate;

        }

        lossRate = lossRateInput.get();
        if (lossRate <= 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }

        frequencies = new double[]{1, 0};

        scarringHeight = scarringHeightInput.get();
        scarringDuration = scarringDurationInput.get();
    }

    /*
    Calculate transition probability matrix for loss without scarring
    The transition probability matrix for loss with scarring differs only in the first row.
     */
    public void getLossProbabilities(double[] matrix, double expOfDeltaLoss){

        // fill diagonal and final column
        for (int i=0; i<nrOfStates; i++){
            for (int j=0; j<nrOfStates; j++){

                if ( i==j ){
                    matrix[i*nrOfStates + j] = expOfDeltaLoss;
                }else if(j == nrOfStates-1){
                    matrix[i*nrOfStates + j] = 1 - expOfDeltaLoss;
                }else{
                    matrix[i*nrOfStates + j] = 0;
                }
            }
        }
        // set final diagonal element to 1
        matrix[nrOfStates * nrOfStates - 1] = 1;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {

        double delta = startTime - endTime;
        double expOfDeltaLoss = Math.exp(-delta * lossRate);

        // calculate transition probabilities for loss process
        getLossProbabilities(matrix, expOfDeltaLoss);

        // for loss & scarring, add the scarring transition probabilities
        if (endTime >= (scarringHeight - scarringDuration)){
            double scarRateSum = Arrays.stream(scarRates).sum();

            // fill first row
            matrix[0] = Math.exp(-delta * (lossRate + scarRateSum));
            for (int i=0; i<nrOfStates-2; i++){
                matrix[i+1] = (scarRates[i] * expOfDeltaLoss - scarRates[i] * Math.exp(-delta * (lossRate + scarRateSum))) / scarRateSum;
            }
            matrix[nrOfStates-1] = 1 - expOfDeltaLoss;
        }
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {

        return null;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return false;
    }
}

