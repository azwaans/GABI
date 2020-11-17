package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

public class GeneralScarringLoss extends SubstitutionModel.Base {

        final public Input<List<RealParameter>> scarringRatesInput = new Input<>("scarringRates",
                "Rates, at which scars types are introduced into the scarring region",
                new ArrayList<>(), Input.Validate.REQUIRED);
        final public Input<RealParameter> lossRateInput = new Input<>("lossRate",
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
        //TODO potentially delete or add loss rates
        double[][] rateMatrix;
        Double[] scarRates;
        Double lossRate;

    @Override
    public void initAndValidate() {

        // one state for each scar type + unedited + lost
        nrOfStates = scarringRatesInput.get().get(0).getDimension() + 2;
        rateMatrix = new double[nrOfStates][nrOfStates];
        scarRates = scarringRatesInput.get().get(0).getValues();

        // assert positive rates
        double sumScarringRates = 0;

        // add scar rates to rate matrix
        for (int i=0; i<scarRates.length; i++){

            if (scarRates[i] <= 0) {
                throw new RuntimeException("All scarring rates must be positive!");
            }
            sumScarringRates += scarRates[i];
            rateMatrix[0][i+1] = scarRates[i];
        }

        lossRate = lossRateInput.get().getValue();
        if (lossRate <= 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }
        for (int i = 0; i<nrOfStates-1; i++){
            rateMatrix[i][nrOfStates-1] = lossRate;
        }

        // center root frequency on unedited state
        frequencies = new double[nrOfStates];
        frequencies[0] = 1;

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
        if ( (endTime >= (scarringHeight - scarringDuration)) & (endTime < scarringHeight)){

            Stream<Double> scarSum = Stream.of(scarRates);
            Double scarRateSum = scarSum.reduce(0.0, (subtotal, element) -> subtotal + element);
            //.sum();

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

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    public double getScarringHeight(){return scarringHeight;}

    public double getScarringDuration(){return scarringDuration;}

    public double[][] getRateMatrix(){return rateMatrix;}
}

