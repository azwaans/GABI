package lineageTree.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;

/*
implements a model that allows scarring in a time window and a continuous loss process

Note: Here, I think of trees plotted with the root at the top and the leaves at the bottom.
 */
public class ScarringLoss extends SubstitutionModel.Base {
    final public Input<RealParameter> scarringRateInput = new Input<>("scarringRate", "Rate, at which scars are introduced into the scarring region", Input.Validate.REQUIRED);
    final public Input<RealParameter> lossRateInput = new Input<>("lossRate", "Rate, at which scarred sites are lost");
    public Input<Double> scarringHeightInput = new Input<>("scarringHeight", "Duration between the onset of scarring and sampling of the cells");
    public Input<Double> scarringDurationInput = new Input<>("scarringDuration", "Duration of the scarring process");


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



    @Override
    public void initAndValidate() {

        double s = scarringRateInput.get().getValue();
        if (s <= 0) {
            throw new RuntimeException("Scarring rate must be positive!");
        }

        // set up matrix for scar loss
        double l = lossRateInput.get().getValue();
        if (l <= 0) {
            throw new RuntimeException("Loss rate must be positive!");
        }


       /* if (frequenciesInput.get() != null) {
            throw new RuntimeException("Frequencies must not be specified in this model. No scars are assumed to exist in the beginning.");
        }*/

        frequencies = new double[]{1, 0};

        scarringHeight = scarringHeightInput.get();
        scarringDuration = scarringDurationInput.get();
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {

        // if node above scarring time, use scar loss matrix
        if (node.getHeight() >= scarringHeight | node.getHeight() < (scarringHeight - scarringDuration)) {
            synchronized (this) {
                double delta = (startTime - endTime);
                double l = lossRateInput.get().getValue();
                matrix[0] = Math.exp(-l * delta);
                matrix[1] = 1 - matrix[0];
                matrix[2] = 0;
                matrix[3] = 1;
            }

            // if node in scarring time, use scar aquisition matrix
        } else if (node.getHeight() >= (scarringHeight - scarringDuration)) {

            synchronized (this) {
                double delta = (startTime - endTime);
                double s = scarringRateInput.get().getValue();
                matrix[0] = Math.exp(-s * delta);
                matrix[1] = 1 - matrix[0];
                matrix[2] = 0;
                matrix[3] = 1;
            }

            // if node below scarring time window, use scar loss matrix
        } else {
            throw new RuntimeException("This node height is strange: " + node.getHeight());
        }
    }


    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {

        // if node above scarring time, no change happens
        if (node.getHeight() > scarringHeight | node.getHeight() < (scarringHeight - scarringDuration)) {

            if (updateEigenLoss) {
                // eigendecomposition for loss
                double l = lossRateInput.get().getValue();
                // eigendecomposition for loss
                double[] eval = new double[]{0.0, -l};
                double[] evec = new double[]{1.0, 1.0, 0, 1.0};
                double[] ivec = new double[]{0, 1.0, 1.0, -1.0};

                getEigenDecompositionLoss = new EigenDecomposition(evec, ivec, eval);
            }

            return getEigenDecompositionLoss;

        } else if (node.getHeight() >= (scarringHeight - scarringDuration)) {

            if (updateEigenScar) {

                // set up matrix for scar aquisition
                double s = scarringRateInput.get().getValue();

                // eigendecomposition for scarring
                double[] eval = new double[]{0.0, -s};
                double[] evec = new double[]{1.0, 1.0, 0, 1.0};
                double[] ivec = new double[]{0, 1.0, 1.0, -1.0};

                getEigenDecompositionScarring = new EigenDecomposition(evec, ivec, eval);

            }

            return getEigenDecompositionScarring;

        } else{
            throw new RuntimeException("Strange node height: " + node.getHeight());
        }
    }


    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof IntegerData;
    }


    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrixScar = true;
        updateEigenScar = true;
        updateEigenLoss = true;
        updateMatrixLoss = true;
        return true;
    }
    @Override
    protected void store() {
        if (getEigenDecompositionScarring != null) {
            storedEigenDecompositionScarring = getEigenDecompositionScarring.copy();
        }
        if (getEigenDecompositionLoss != null){
            storedEigenDecompositionLoss = getEigenDecompositionLoss.copy();
        }
        super.store();
    }

    @Override
    protected void restore() {
        updateMatrixScar = true;
        updateMatrixLoss = true;
        updateEigenScar = true;
        updateEigenLoss = true;

        if (storedEigenDecompositionScarring != null) {
            getEigenDecompositionScarring = storedEigenDecompositionScarring;
        }
        if (storedEigenDecompositionLoss != null){
            getEigenDecompositionLoss = storedEigenDecompositionLoss;
        }
        super.restore();
    }
}