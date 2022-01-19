package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.alignment.BarcodeMeta;
import beast.evolution.datatype.DataType.Base;

import java.util.Arrays;
import java.util.List;

/**
 todo: check codeMap examples.
 */
@Description("Datatype for string sequence representing barcodes states.")
public class GstaltData extends Base {


    public void initAndValidate(){
        Log.info.println("WE ARE IN INIT AND VALIDATE OF GESTALT DATA");
        //TODO split the string into list of strings
        stateCount  = -1;
        mapCodeToStateSet = null;
        codeLength = -1;
        codeMap = null;
    }

    @Override
    public String getTypeDescription() {
        return "gstaltData";
    }

    @Override
    public boolean isAmbiguousCode(int code) {
        return code < 0;
    }

    @Override
    public String getCharacter(int code) {
        if (code < 0) {
            Log.info.println("getCharacter is called and returns ? ");
            Log.info.println("CODE CODE " + code);

            return "?";
        }
        return code + "";
    }

    @Override
    public int[] getStatesForCode(int code) {

        return new int[]{code};
    }

    public boolean hasConstantCodeLength() {
        return false;
    }

}
