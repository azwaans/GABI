package gestalt.evolution.alignment;

import beast.base.core.BEASTObject;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class TargetDeactTract extends BEASTObject {

    /**
     * This class represents realisations of a deactivating event of a specific range of targets
     */

    protected int minTargDeac;
    protected int maxTargDeac;

    public TargetDeactTract(int minTrgtDeac, int maxTrgtDeac) {
        minTargDeac = minTrgtDeac;
        maxTargDeac = maxTrgtDeac;
        initAndValidate();
    }

    public void initAndValidate() {


    }

    /**
     * returns all possible statuses that could have lead to a DeacTract (all possible masked statuses)
     */

    public List<TargetStatus> getContainedStatuses() {
        //first create an array of 0 and 1
        int M = this.maxTargDeac - this.minTargDeac + 1;


        //creating all possible 0/1 cvombinations:
        int[][] allPossibleStat = new int[(int) (Math.pow(2, M))][M];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j <= (int) (Math.pow(2, i + 1) - 1); j++) {
                int startIndex = (int) ((int) (j) * (Math.pow(2, M - i - 1)));
                int endIndex = (int) ((j + 1) * (Math.pow(2, M - i - 1)) - 1);

                if (j % 2 == 0) {
                    //fill with 1
                    for (int h = startIndex; h <= endIndex; h++) {
                        allPossibleStat[h][i] = 1;
                    }
                }
                if (j % 2 == 1) {
                    //fill with 0
                    for (int h = startIndex; h <= endIndex; h++) {
                        allPossibleStat[h][i] = 0;
                    }
                }
            }
        }


        List<TargetStatus> targetStatuses = new ArrayList<>();
        for (int[] temp : allPossibleStat) {
            List<TargetDeactTract> deactTracts = new ArrayList<>();
            Integer currTrgtDeactStrt = null;
            int currentIndex = 0;
            for (int value : temp) {
                if (value == 1) {
                    if (currTrgtDeactStrt == null) {
                        currTrgtDeactStrt = currentIndex;
                    }
                } else {
                    if (currTrgtDeactStrt != null) {
                        //we just exited a target tract, need to save it! We also reset the current tract tracker back to nul
                        deactTracts.add(new TargetDeactTract(this.minTargDeac + currTrgtDeactStrt, this.minTargDeac + currentIndex - 1));
                        currTrgtDeactStrt = null;
                    }
                }
                currentIndex++;
            }
            /*We have reached the end of the Status encoding. If current deactivated tract tracker is still not null,
            it means the tract reaches the border of the tract of interest.
             */
            if (currTrgtDeactStrt != null) {
                deactTracts.add(new TargetDeactTract(this.minTargDeac + currTrgtDeactStrt, this.maxTargDeac));
            }
            targetStatuses.add(new TargetStatus(deactTracts));
        }

        //Sort the target statuses by number of deactivated targets (ASSUMING ASCENDING ORDER)
        Collections.sort(targetStatuses);

        return targetStatuses;
    }


}
