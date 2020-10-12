package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.DefaultEigenSystem;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Node;

public class MutationAndLoss extends GeneralSubstitutionModel{


    final public Input<RealParameter> scarringRateInput = new Input<>("scarringRate", "Rate, at which scars are introduced into the scarring region", Input.Validate.REQUIRED);
    final public Input<RealParameter> lossRateInput = new Input<>("lossRate", "Rate, at which scarred sites are lost");

    public MutationAndLoss() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
    }



    public void initAndValidate(){

        RealParameter scarringRate = scarringRateInput.get();
        RealParameter lossRate  = lossRateInput.get();
        frequencies = frequenciesInput.get();
        nrOfStates = frequencies.getFreqs().length;
        if (nrOfStates != 3) {
            throw new IllegalArgumentException("Frequencies has wrong size. Expected 3, but got " + nrOfStates);
        }

        rateMatrix = new double[(nrOfStates)][(nrOfStates)];

        eigenSystem = new DefaultEigenSystem(nrOfStates);
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        synchronized (this) {
            if (updateMatrix) {
                double [][] rateMatrix = getRateMatrix();

                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }

    public double [][] getRateMatrix(){

        if (nrOfStates == 2){
            double [] [] rateMatrix = new double[3][3];
            rateMatrix[0][0] = -scarringRateInput.get().getValue();
            rateMatrix[0][1] = scarringRateInput.get().getValue();
            rateMatrix[1][1] = 1;
            rateMatrix[1][0] = 0;
        }
        return rateMatrix;

    }


}
