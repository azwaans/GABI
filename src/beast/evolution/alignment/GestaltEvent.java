package beast.evolution.alignment;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import org.antlr.v4.runtime.misc.Pair;

import java.util.Comparator;

@Description("Class representing gestalt event objects")
public class GestaltEvent extends BEASTObject implements Comparable<GestaltEvent> {

    final public Input<String> dataInput = new Input<>("value",
            "Event data, encoded as an underscore-separated list of integers and strings" +
                    "In either case, whitespace is ignored.", Input.Validate.REQUIRED);

    protected int startPos;
    protected int delLen;
    protected int minTarg;
    protected int maxTarg;
    protected String insSeq;
    protected int delEnd;

    public int compareTo(GestaltEvent ev2) {
        return this.startPos - ev2.startPos;
    }

    public GestaltEvent() {
    }
    /**
     * Constructor for testing.
     *
     * @param event
     */
    public GestaltEvent(String event) {
        dataInput.setValue(event, this);
        initAndValidate();
    }



    public void initAndValidate() {
        String data = dataInput.get();
        //string array for the event
        String[] strs = data.split("_");
        startPos =  Integer.parseInt(strs[0]);
        delLen = Integer.parseInt(strs[1]);
        minTarg = Integer.parseInt(strs[2]);
        maxTarg = Integer.parseInt(strs[3]);
        delEnd = startPos + delLen;

        if(strs.length == 5) {
            insSeq = strs[4];
        }
        else {insSeq ="";}





    } // initAndValidate


    /**
     * @return the data of this sequence as a string.
     */
    public final String getData() {
        return dataInput.get();
    }

    public org.apache.commons.math3.util.Pair<Integer,Integer> getMinMaxDeactTargets(java.util.List<org.apache.commons.math3.util.Pair<Integer,Integer>> posSites, int nTargets) {

        Integer minDeacTarg = minTarg;
        if((minTarg > 0) && (startPos < posSites.get(minTarg-1).getSecond())){
        minDeacTarg = minTarg -1;
        }
        Integer maxDeacTarg = maxTarg;
        //here is (n target -1)
        if((maxTarg < nTargets - 1) && (posSites.get(maxTarg+1).getFirst() < this.delEnd)){
            maxDeacTarg = maxTarg +1;
        }
        org.apache.commons.math3.util.Pair<Integer,Integer> MinMaxDeactTargets =new org.apache.commons.math3.util.Pair(minDeacTarg,maxDeacTarg);
        return MinMaxDeactTargets;
        }
        /*if self.min_target > 0 and self.start_pos <= bcode_meta.pos_sites[self.min_target - 1][1]:
        min_deact_target = self.min_target - 1
        else:
        min_deact_target = self.min_target

        if self.max_target < bcode_meta.n_targets - 1 and bcode_meta.pos_sites[self.max_target + 1][0] < self.del_end:
        max_deact_target = self.max_target + 1
        else:
        max_deact_target = self.max_target

        assert min_deact_target is not None
        assert max_deact_target is not None
        return min_deact_target, max_deact_target
*/



    }




