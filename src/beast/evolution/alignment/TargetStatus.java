package beast.evolution.alignment;

import beast.core.BEASTObject;
import beast.core.util.Log;
import org.antlr.v4.runtime.misc.Pair;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class TargetStatus extends BEASTObject implements Comparable < TargetStatus > {

    protected List<TargetDeactTract> targetDeacList;
    /*Target Status. The barcode target status consists of a collection of deactivation tracts. They result from applying
    a Target Tract indel to the barcode. The corresponding feature resulting from applying Target Tracts to a barcode
    are Target Deactivation Tracts. A Target Status can be converted into a binary format, as represented in J.Feng's paper.*/

    public TargetStatus() {
        targetDeacList = null;
        initAndValidate();
    }


    public TargetStatus(List<TargetDeactTract> list) {
        targetDeacList = list;
        initAndValidate();
    }

    public void initAndValidate() {
    }

    @Override
    public boolean equals(Object obj)
    {
        if (obj == this) {
            return true;
        }
        else {
            return this.hashCode() == obj.hashCode();
        }
    }
    @Override
    public int hashCode() {
        //here is ntarget
        Integer[] binary = this.getBinaryStatus(10);
        String concatenated = "";
        for(int i:binary) {
            concatenated = concatenated + i;
        }
        int code = Integer.parseInt(concatenated);
        return code;
    }


    public void addTarget(IndelSet.TargetTract targTrac) {
        //adding a target tract to an existing target status.
        //if the target status is empty, initialize it with the tract



        if (targetDeacList == null || targetDeacList.size() == 0 ) {

            List<TargetDeactTract> initList = new ArrayList<>();
            initList.add(new TargetDeactTract(targTrac.minTargDeac, targTrac.maxTargDeac));
            targetDeacList = initList;
        }
        //else, we find the max target already deactivated in the barcode
        else {

            int maxTargetCurrent = targetDeacList.get(targetDeacList.size() - 1).maxTargDeac + 1;

            int maxTargets = java.lang.Math.max(targTrac.maxTargDeac + 1, maxTargetCurrent);
        /*creating a binary representation of the barcode status, which is only as long as the final span of affected targets
        if the targets affected by this new target tract are already within the SPAN of previously affected targets, then the
        length of the binary representation does NOT change!*/
            Integer[] binaryStatus = this.getBinaryStatus(maxTargets);
            java.util.Arrays.fill(binaryStatus, targTrac.minTargDeac, targTrac.maxTargDeac + 1, 1);

            targetDeacList = convertBinaryToTarget(binaryStatus).targetDeacList;
        }
    }

    public List<Integer> getDeactTargets() {
        List<Integer> deactTargs = new ArrayList<>();
        if(targetDeacList != null) {
            for (TargetDeactTract tdt : targetDeacList) {
                for (int i = tdt.minTargDeac; i <= tdt.maxTargDeac; i++) {
                    deactTargs.add(i);
                }

            }
        }
        return deactTargs;
    }

    public static TargetStatus convertBinaryToTarget(Integer[] binaryStat) {
    /*converting a binary target status into a TargetStatus with TargetTractDeac attribute. We iterate through th
    binary representation. !The targets are indexed starting from 0! */
        int targetIndex = 0;
        Integer currentTrgtStart = null;
        List<TargetDeactTract> Deactracts = new ArrayList<>();
        for (int status : binaryStat) {
            if (status == 1) {
                //status is deactivated
                if (currentTrgtStart == null) {
                    //set current deactivated target to the index
                    currentTrgtStart = targetIndex;
                }

            } else if (currentTrgtStart != null) {
                /*status is back to being active, AND we left a deactivation tract (as recorded by current target start)
                we record the Deact Tract by appending the list of deac tracts. We use index -1 (we went past it)*/
                Deactracts.add(new TargetDeactTract(currentTrgtStart, targetIndex - 1));
                //we set the current target start back to none.
                currentTrgtStart = null;
            }
            targetIndex++;
        }
        /*Once we are done going through the binary target status representation, check if we have a target that extends until the
        end of the barcode*/
        if (currentTrgtStart != null) {
            Deactracts.add(new TargetDeactTract(currentTrgtStart, targetIndex - 1));
        }
        return new TargetStatus(Deactracts);
    }

    public Integer[] getBinaryStatus(int nTargets) {
        /*this creates a binary representation of the status of all targets <= n_targets. 0 if active, 1 if deactivated.
        should this be a boolean[]/a list?*/
        Integer[] binaryStatus = new Integer[nTargets];
        java.util.Arrays.fill(binaryStatus, 0);
        if (targetDeacList != null) {
            for (TargetDeactTract temp : targetDeacList) {

                java.util.Arrays.fill(binaryStatus, temp.minTargDeac, temp.maxTargDeac+1, 1);
            }
            return binaryStatus;
        } else {
            return binaryStatus;
        }
    }

    public void fillWithTargetTracts(List<IndelSet.TargetTract> targetTractList) {
        if (targetTractList != null) {
            for (IndelSet.TargetTract tt : targetTractList) {
                this.addTarget(tt);
            }
        }
    }

    public TargetStatus merge(TargetStatus otherTarget) {
        //merging 2 target statuses:
        if (this.targetDeacList.size() == 0) {
            return otherTarget;
        }
        if (otherTarget.targetDeacList.size() == 0) {
            return this;
        }
        /*Creating a new binary representation of the merged statuses. To do this, we first
        find out what will be the final length of the merged tract.
         */
        int maxThis = this.targetDeacList.get(targetDeacList.size() - 1).maxTargDeac + 1;
        int maxOther = otherTarget.targetDeacList.get(otherTarget.targetDeacList.size() - 1).maxTargDeac + 1;
        int maxTargets = java.lang.Math.max(maxThis, maxOther);
        /* create a new status based on the existing one, extending its size to max targets, the final status size (this just creates new empty
        targets in the binary representation if needed.
         */
        Integer[] binaryStatusExtended = this.getBinaryStatus(maxTargets);
        // fill the existing, extended status of *this with the other status
        for (TargetDeactTract deactTract : otherTarget.targetDeacList) {
            Arrays.fill(binaryStatusExtended, deactTract.minTargDeac, deactTract.maxTargDeac + 1, 1);
        }
        return convertBinaryToTarget(binaryStatusExtended);

    }

    public int compareTo(TargetStatus other) {
        return (this.targetDeacList.size()) - ((other.targetDeacList.size()));
    }

    public List<Integer> getInactiveTargets(int nTargets) {
        //returns edited target positions in the barcode.
        //TO DO: ADJUST TARGET NUMBER !!!!!!!!!!!!!!!!!!!!!!
        Integer[] binary = this.getBinaryStatus(nTargets);
        List<Integer> binaryStatus = Arrays.asList(binary);
        List<Integer> indices = new ArrayList<>();
        int index = binaryStatus.indexOf(1);
        if (index != -1) {
            while (index != -1) {
                binaryStatus.set(index, 0);
                indices.add(index);
                index = binaryStatus.indexOf(1);

            }
            return indices;
        } else {
            return indices;
        }
    }

    public List<Integer> getActiveTargets(int nTargets) {
        //returns unedited target positions in the barcode.
        //TO DO: ADJUST TARGET NUMBER !!!!!!!!!!!!!!!!!!!!!!
        Integer[] binary = this.getBinaryStatus(nTargets);
        List<Integer> binaryStatus = Arrays.asList(binary);
        List<Integer> indices = new ArrayList<>();
        int index = binaryStatus.indexOf(0);
        while (index != -1) {
            indices.add(index);
            binaryStatus.set(index, 1);
            index = binaryStatus.indexOf(0);
        }
        return indices;
    }


    public List<IndelSet.TargetTract> getPossibleTargetTracts(List<Integer> activeAny, int nTargets) {
        //given a specific target tract, returns all possible target tracts that could have generated it.
        //DO this by enumerating the start position and end positioon of the tt tract,
        //there are 2 edit types on each side of the deletion: short or long trim
        //the "basic" entry is a short trim

        if (activeAny.size() == 0) {
            activeAny = this.getActiveTargets(nTargets);
        }

        int nAnyTarg = activeAny.size();
        Hashtable<Integer, List<Pair<Integer, Integer>>> allStartOptions = new Hashtable<>();
        int i0Prime = 0;
        for (Integer t0prime : activeAny) {

            List<Pair<Integer, Integer>> basicEntry = new ArrayList<>();
            //adding a short left trim: cuts and deactivates the same target to the left
            basicEntry.add(new Pair(t0prime, t0prime));
            allStartOptions.put(i0Prime, basicEntry);
            //if the target tract does not start at the nfirst target, it can also affect the preceding target
            //the cut affects a target position to the left
            if (t0prime > 0) {
                List<Pair<Integer, Integer>> temp = allStartOptions.get(i0Prime);
                temp.add(new Pair(t0prime - 1, t0prime));
                allStartOptions.replace(i0Prime, temp);
            }
            ++i0Prime;
        }

        int i1Prime = 0;
        Hashtable<Integer, List<Pair<Integer, Integer>>> allEndOptions = new Hashtable<>();
        for (Integer t1prime : activeAny) {
            List<Pair<Integer, Integer>> basicEndEntry = new ArrayList<>();
            basicEndEntry.add(new Pair(t1prime, t1prime));

            allEndOptions.put(i1Prime, basicEndEntry);
            //change with metatdata in the future
            // here it is the last target index: (n tTarget -1)
            if (t1prime < nTargets -1) {
                List<Pair<Integer, Integer>> temp = allEndOptions.get(i1Prime);
                temp.add(new Pair(t1prime, t1prime + 1));
                allEndOptions.replace(i1Prime, temp);

            }
            ++i1Prime;
        }










        List<IndelSet.TargetTract> ttTracts = new ArrayList<>();
        for (Integer startKey : allStartOptions.keySet()) {

            List<Pair<Integer, Integer>> ttStarts = allStartOptions.get(startKey);
            //istart is correct value 0
            for (int k : IntStream.range(startKey, nAnyTarg).boxed().collect(Collectors.toList())) {
                List<Pair<Integer, Integer>> ttEnds = allEndOptions.get(k);
                for (Pair<Integer, Integer> ttStart : ttStarts) {
                    for (Pair<Integer, Integer> ttEnd : ttEnds) {
                        IndelSet.TargetTract ttevent = new IndelSet.TargetTract(ttStart.b, ttEnd.a, ttStart.a, ttEnd.b);

                        ttTracts.add(ttevent);
                    }
                }

            }
        }

        return ttTracts;

    }

    public Set<Integer> minus(TargetStatus originalTargetStatus) {
        Set<Integer> origTargs = new HashSet<>(originalTargetStatus.getDeactTargets());
        Set<Integer> selfTargs = new HashSet<>(this.getDeactTargets());
        if (selfTargs.containsAll(origTargs)) {
            selfTargs.removeAll(origTargs);
            return selfTargs;
        } else {
            return new HashSet<>();
        }


    }

    public static Hashtable<TargetStatus, Hashtable<TargetStatus, List<IndelSet.TargetTract>>> getAllTransitions(int nTargets) {
        Hashtable<TargetStatus, Hashtable<TargetStatus, List<IndelSet.TargetTract>>> targetStatusTransitionDict = new Hashtable<>();
        TargetDeactTract deacTargs = new TargetDeactTract(0, nTargets - 1);

        List<TargetStatus> targetStatuses = deacTargs.getContainedStatuses();
        for (TargetStatus targStat : targetStatuses) {
            Hashtable<TargetStatus, List<IndelSet.TargetTract>> targStatStartDict = new Hashtable<>();
            List<IndelSet.TargetTract> possibleTargetTracts = targStat.getPossibleTargetTracts(new ArrayList(),nTargets);
            for (IndelSet.TargetTract targetTract : possibleTargetTracts) {
                TargetStatus newTargStat = new TargetStatus(targStat.targetDeacList);
                newTargStat.addTarget(targetTract);
                if (targStatStartDict.containsKey(newTargStat)) {
                    List<IndelSet.TargetTract> toAppend = targStatStartDict.get(newTargStat);
                    toAppend.add(targetTract);
                    targStatStartDict.put(newTargStat, toAppend);

                } else {
                    List<IndelSet.TargetTract> entry = new ArrayList();
                    entry.add(targetTract);
                    targStatStartDict.put(newTargStat, entry);
                }
            }
            targetStatusTransitionDict.put(targStat,targStatStartDict);
        }


        return targetStatusTransitionDict;
    }

}