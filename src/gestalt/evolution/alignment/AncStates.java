package gestalt.evolution.alignment;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import org.apache.commons.math3.util.Pair;

import java.util.*;


public class AncStates extends BEASTObject {

    /**
     * AncStates is the set of likely ancestral states at a specific node in the tree, given leaf sequences. This set is characterized as any subset of alleles generated by
     * a list of indels, which are singletons or wild cards (AncState is then also defined as a union of disjoint WC/SWCS).
     */

    public List<IndelSet> indelSetList;

    /**
     * Constructor
     */
    public AncStates(List<IndelSet> list) {
        indelSetList = list;
        initAndValidate();
    }

    /**
     * No argument constructor
     */
    public AncStates() {
        indelSetList = new ArrayList<>();
        initAndValidate();
    }

    public void initAndValidate() {
    }


    /**
     * get all singletons from the indelset list
     */
    public List<IndelSet.Singleton> getSingletons() {
        List<IndelSet> singletonwcs = this.getSingletonWCs();
        List<IndelSet.Singleton> singletons = new ArrayList<>();
        if (singletonwcs != null) {
            for (IndelSet temp : singletonwcs) {
                singletons.add(temp.getSingleton());
            }
        }
        return singletons;
    }

    /**
     * get all singletonWCs from the indelset list
     */
    public List<IndelSet> getSingletonWCs() {
        List<IndelSet> singletonsWCs = new ArrayList<>();
        for (IndelSet temp : indelSetList) {
            if (temp.getClass() == IndelSet.SingletonWC.class) {
                singletonsWCs.add(temp);
            }
        }
        return singletonsWCs;
    }

    /**
     * Getting the intersection between two AncState sets. We first iterate through both indel set lists until we find
     * 2 overlapping indel sets. We return a new AncState set containing these intersections. (the intersection of 2
     * ancestral state sets is an ancestral set).
     */
    public static AncStates intersect(AncStates ancStates1, AncStates ancStates2) {

        //todo check that!
        if (ancStates1.indelSetList.size() == 0 || ancStates2.indelSetList.size() == 0) {
            return new AncStates();
        }
        //We can explore the indel sets:
        else {
            int indx1 = 0;
            int indx2 = 0;
            int leng1 = ancStates1.indelSetList.size();
            int leng2 = ancStates2.indelSetList.size();
            List<IndelSet> intersectionList = new ArrayList<>();
            while ((indx1 < leng1) && (indx2 < leng2)) {
                IndelSet indel1 = ancStates1.indelSetList.get(indx1);
                IndelSet indel2 = ancStates2.indelSetList.get(indx2);
                if (indel2.maxTarg < indel1.minTarg) {
                    //we move the index for indel list 1 attempting to get closer to the ones in indel2
                    indx2++;
                    continue;
                } else if (indel1.maxTarg < indel2.minTarg) {
                    //we move the index for indel list 1 attempting to get closer to the ones in indel2
                    indx1++;
                    continue;
                }

                //after this, we should be on overlapping indel sets. We find their intersection, and append the list:


                IndelSet interIndel = IndelSet.intersect(indel1, indel2);

                if (interIndel != null) {
                    intersectionList.add(interIndel);
                }

                if (indel1.maxTarg < indel2.maxTarg) {
                    indx1++;
                } else {
                    indx2++;
                }


            }
            //post processing of SGWC:

            intersectionList = processIntersection(intersectionList, ancStates1);
            intersectionList = processIntersection(intersectionList, ancStates2);

            return new AncStates(intersectionList);
        }
    }

    /**
     * returns true is the state defined by a given Target tract list is possible, against an ancestral state
     * in other words: is the input state a member of the AncStates list.
     */
    public boolean isPossible(List<IndelSet.TargetTract> state) {
        //1. Get the singleton list of the AncList
        List<IndelSet> sgs = this.getSingletonWCs();
        //2. Convert it into a target tract list
        List<IndelSet.TargetTract> sgsTTs = new ArrayList<>();
        for (IndelSet sg : sgs) {
            sgsTTs.add(sg.getTargetTract());
        }

        //go though the input TT list and construct a fake ancestral state based on these individual edits

        List<IndelSet> fakeAncState = new ArrayList<>();
        for (IndelSet.TargetTract tract : state) {
            //the target tract is already in the ANC State singleton tt list
            if (sgsTTs.contains(tract)) {
                fakeAncState.add(sgs.get(sgsTTs.indexOf(tract)));
            }
            //the target tract is not one of the singleton target tracts, create a wild card
            else {
                fakeAncState.add(new IndelSet.WildCard(tract.minTargDeac, tract.maxTargDeac));
            }
        }

        //intersect the fake ancestral state with the AncState


        List<IndelSet> intersectAncState = AncStates.intersect(this, new AncStates(fakeAncState)).indelSetList;

        //shortcut: if both ancestral states are of different size, already know that it is impossible
        if (intersectAncState.size() != fakeAncState.size()) {
            return false;

        }
        for (int i = 0; i < intersectAncState.size(); ++i) {
            if (!(intersectAncState.get(i).equals(fakeAncState.get(i)))) {
                return false;
            }
        }
        return true;

    }

    /**
     * Post-process intersection with a final check for singleton-wildcards
     * that deactivate the same target -- there is an ordering between these touching
     * singleton-wildcards that must be obeyed for irreversibility
     */
    private static List<IndelSet> processIntersection(List<IndelSet> intersectList, AncStates ancStates) {
        Hashtable<Integer, IndelSet> minTargetToIndelSet = new Hashtable<>();
        Hashtable<Integer, IndelSet> maxTargetToIndelSet = new Hashtable<>();
        //filling up dictionaries of indel sets where the min target == min target deac (clean cut)
        //or max target deact = max target deac
        for (IndelSet temp : ancStates.indelSetList) {
            if (temp.getminTargDeac() == temp.minTarg) {
                minTargetToIndelSet.put(temp.getminTargDeac(), temp);
            }
            if (temp.getmaxTargDeac() == temp.maxTarg) {
                maxTargetToIndelSet.put(temp.getmaxTargDeac(), temp);
            }
        }
        List<IndelSet> newIntersectionList = new ArrayList<>();
        int indx1 = 0;
        for (IndelSet temp : intersectList) {
            //check that the indel set is touching an indel set in the max dictionnary
            if ((temp.minTarg != temp.getminTargDeac()) && ((maxTargetToIndelSet.keySet()).contains(temp.getminTargDeac()))) {
                if (indx1 == 0 || (intersectList.get(indx1 - 1).getmaxTargDeac()) != temp.getminTargDeac()) {
                    indx1++;
                    continue;


                }
            }
            //check that the indel set is touching an indel set in the min dictionnary
            if ((temp.maxTarg != temp.getmaxTargDeac()) && ((minTargetToIndelSet.keySet()).contains(temp.getmaxTargDeac()))) {
                if (indx1 == intersectList.size() - 1 || (intersectList.get(indx1 + 1).getminTargDeac()) != temp.getmaxTargDeac()) {
                    indx1++;
                    continue;
                }
            }
            newIntersectionList.add(temp);
            indx1++;
        }
        return newIntersectionList;
    }

    public List<TargetStatus> getTargetSubStatuses(IndelSet indelSet) {
        /* we can only have WC and Singleton-WC.
        If the indelSet is a SGWC, the target sub statuses are the status generated by the Singleton + the substatuses generated by
        its inner wild card
         */

        if (indelSet instanceof IndelSet.SingletonWC) {
            List<TargetStatus> allStatuses = new ArrayList<>();
            // adding the singleton status
            List<TargetDeactTract> inputSGTractList = Arrays.asList(new TargetDeactTract(indelSet.getminTargDeac(), indelSet.getmaxTargDeac()));
            TargetStatus singletonStatus = new TargetStatus(inputSGTractList);
            allStatuses.add(singletonStatus);
            // adding the inner wild card statuses
            IndelSet.WildCard innerWC = indelSet.getInnerWC();
            if (innerWC != null) {
                TargetDeactTract deactTract = new TargetDeactTract(innerWC.minTarg, innerWC.maxTarg);
                allStatuses.addAll((deactTract.getContainedStatuses()));
            } else {
                //TO DO Check the necessity of explicitly adding an empty TargetStatus, as in GAPML ?????
                //subStatuses.add(null);
            }
            return allStatuses;
        } else {
            //the indel set is a wild card. we only have to get the substatuses generated by the wild card.
            TargetDeactTract deactTract = new TargetDeactTract(indelSet.getminTargDeac(), indelSet.getmaxTargDeac());
            return deactTract.getContainedStatuses();
        }
    }

    public List<TargetStatus> generatePossibleTargetStatuses() {
        if (indelSetList.size() == 0) {
            //TO DO Check null or empty target status?????
            return null;
        }

        //Start by creating a list of all TargetSubStatuses generated by all sets in the AncState.
        List<List<TargetStatus>> partitionTargetStatuses = new ArrayList<>();
        for (IndelSet set : indelSetList) {
            String className = set.getClass().getSimpleName();
            partitionTargetStatuses.add(getTargetSubStatuses(set));
        }

        // create the cartesian product of all sets
        SetOperations cartesian = new SetOperations();
        List<List<TargetStatus>> fullTargetStatuses = cartesian.cartesianProduct(partitionTargetStatuses);

        if (fullTargetStatuses == null) {
            return null;
        }
        //merge all resulting sets. For all lists in full target status, create the merged target status
        List<TargetStatus> mergedFullTargets = new ArrayList<>();
        for (List<TargetStatus> listToMerge : fullTargetStatuses) {

            TargetStatus merger = listToMerge.get(0);
            //TO DO remove this debuging print
            for (TargetStatus statusToMergein : listToMerge) {
                merger = merger.merge(statusToMergein);
            }
            mergedFullTargets.add(merger);
        }

        return mergedFullTargets;
    }

    static public AncStates createObservedAlleleSet(String sequenceInput, List<Pair<Integer, Integer>> posSites, Integer nTargets) {
        //create an ancestral state set for observed sequences at leaves
        //todo check if this split with ":" is necessary.
        String[] strings = sequenceInput.split(":");
        List<String> event_strings = new LinkedList<String>(Arrays.asList(strings[1].split(",")));
        event_strings.removeAll(Collections.singleton("None"));
        //create a SGWC list.
        List<GestaltEvent> eventList = new ArrayList<>();
        for (String i : event_strings) {

            GestaltEvent temp = new GestaltEvent(i);
            eventList.add(temp);

        }
        Collections.sort(eventList);
        List<IndelSet> singletonList = new ArrayList<>();


        for (GestaltEvent temp : eventList) {
            Pair<Integer, Integer> adjustedMinMax = temp.getMinMaxDeactTargets(posSites, nTargets);
            singletonList.add(new IndelSet.SingletonWC(temp.minTarg, temp.maxTarg, temp.delLen, temp.startPos, adjustedMinMax.getFirst(), adjustedMinMax.getSecond(), temp.insSeq));
        }
        //ORDER THE LIST HERE!!
        return new AncStates(singletonList);
    }

    public TargetStatus getMaxTargetStatus() {
        // returns the target status generated by applying all indel sets in the indel set list to the unedited barcode.

        if ((this.indelSetList).size() == 0) {
            return new TargetStatus();
        } else {
            int maxTarg = this.indelSetList.get(this.indelSetList.size() - 1).getmaxTargDeac();
            Integer[] binaryStatus = new Integer[maxTarg + 1];
            Arrays.fill(binaryStatus, 0);
            for (IndelSet temp : this.indelSetList) {
                Arrays.fill(binaryStatus, temp.getminTargDeac(), temp.getmaxTargDeac() + 1, 1);
            }
            TargetStatus maxStatus = TargetStatus.convertBinaryToTarget(binaryStatus);
            return maxStatus;
        }


    }

}