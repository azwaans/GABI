package beast.evolution.alignment;
import beast.evolution.tree.Node;
import beast.core.BEASTObject;
import beast.core.util.Log;
import org.antlr.v4.runtime.misc.Triple;
import org.apache.commons.math3.util.Pair;
import org.junit.Test;

import java.util.*;

public class TransitionWrap extends BEASTObject {


    public void initAndValidate() {

    } // initAndValidate

    public List<List<IndelSet.TargetTract>> ttTuple;
    public List<TargetStatus> transStatuses;
    public AncStates ancStates;
    public int numStatuses;
    //DO we really need this mapping?
    public Hashtable<TargetStatus, Integer> statusMap;
    public boolean isTreeLeaf;
    public TargetStatus leafState;

    public TransitionWrap(List<List<IndelSet.TargetTract>> ttupleIn, AncStates ancstateIn, boolean isLeafin) {
        ttTuple = ttupleIn;


        transStatuses = new ArrayList<>();
        Set<TargetStatus> noDuplicates = new HashSet<>();
        if (ttupleIn.size() != 0) {
            Set ohneDup = new LinkedHashSet();
            ohneDup.addAll(ttupleIn);
            ttupleIn.clear();
            ttupleIn.addAll(ohneDup);
            for (List<IndelSet.TargetTract> temp : ttupleIn) {
                TargetStatus tempStatus = new TargetStatus();
                tempStatus.fillWithTargetTracts(temp);
                noDuplicates.add(tempStatus);
            }

        }
        transStatuses.addAll(noDuplicates);
        statusMap = new Hashtable<>();
        //LATEST CHANGE
        for(int i=0;i<transStatuses.size();++i) {
            statusMap.put(transStatuses.get(i),i);
        }

        ancStates = ancstateIn;
        isTreeLeaf = isLeafin;
        numStatuses = transStatuses.size();
        leafState = null;
        if (isLeafin) {
            leafState = ancStates.getMaxTargetStatus();
        }
        initAndValidate();
    }



    public void CreateStatusMap() {
        statusMap = new Hashtable<>();
        for(int i=0;i<transStatuses.size();++i) {
            statusMap.put(transStatuses.get(i),i);
        }

    }

    public TransitionWrap() {
        ttTuple = new ArrayList<>();
        transStatuses = new ArrayList<>();
        ancStates = null;
        isTreeLeaf = false;
        numStatuses = 0;
        leafState = null;
        initAndValidate();
    }

    public static TransitionWrap getCloseTransitionWrapper(Node parNode,
                                                           Node childNode,
                                                           Hashtable<Integer, AncStates> ancStatesDict,
                                                           Hashtable<Integer, List<IndelSet.Singleton>> parsimonyStatesDict,
                                                           List<List<IndelSet.TargetTract>> parentTuples,
                                                           Integer maxSumStates,
                                                           int maxExtraSteps,
                                                           int nTargets
                                                           ) {

        TransitionWrap wrap = new TransitionWrap();
        AncStates nodeAncState = ancStatesDict.get(childNode.getNr());
        List<IndelSet.Singleton> parentSingletonState = parsimonyStatesDict.get(parNode.getNr());


        List<IndelSet.Singleton> childSingletonState = parsimonyStatesDict.get(childNode.getNr());

        //finding length of the set difference between child and parent to find to minimum number of steps:
        List<IndelSet.Singleton> differenceSteps = new ArrayList<>(childSingletonState);
        differenceSteps.removeAll(parentSingletonState);
        int minSteps = differenceSteps.size();

        //PREPARE THE DESIRED CHILD NODE STATUS IN THE PROPER FORMAT (FROM SINGLETON TO STATUS)

        TargetStatus requestedStatus = new TargetStatus();
        List<IndelSet.TargetTract> ttlist = new ArrayList<>();
        for (IndelSet.Singleton temp : childSingletonState) {
            ttlist.add(temp.getTargetTract());
        }
        requestedStatus.fillWithTargetTracts(ttlist);
        //CHILD NODE REQUESTED STATUS READY

        boolean statesTooMany = true;

        //TO DO check the extra steps calcultation rational?
        Integer maxExtraFinal = maxExtraSteps;


        if (maxSumStates != null) {
            if (!(java.lang.Math.pow(2, minSteps) <= maxSumStates)) {
                maxExtraFinal = java.lang.Math.max(maxExtraSteps - 1, 0);
            }
        }

        while (statesTooMany && (maxExtraFinal >= 0)) {

            List<List<IndelSet.TargetTract>> closeTTlist = getCloseStates((minSteps + maxExtraFinal),nTargets,parentTuples, nodeAncState, requestedStatus);



            if(closeTTlist.size() == 1 ) {
                List<IndelSet.TargetTract> empty = new ArrayList<>();
                if (closeTTlist.get(0) instanceof IndelSet.TargetTract) {
                    //only 1 element
                    IndelSet.TargetTract ttsingleton = (IndelSet.TargetTract) closeTTlist.get(0);
                    List<IndelSet.TargetTract> singleList = new ArrayList<IndelSet.TargetTract>();
                    singleList.add(ttsingleton);
                    List<List<IndelSet.TargetTract>> singleListList = new ArrayList<>();
                    singleListList.add(singleList);
                    wrap = new TransitionWrap(singleListList, nodeAncState, childNode.isLeaf());
                }
                else if (closeTTlist.get(0).size() == 0 ) {
                    //completely empty list
                    //null checking ?
                    List<IndelSet.TargetTract> singleList = new ArrayList<IndelSet.TargetTract>();
                    List<List<IndelSet.TargetTract>> singleListList = new ArrayList<>();
                    singleListList.add(singleList);
                    wrap = new TransitionWrap(singleListList, nodeAncState, childNode.isLeaf());
                }




            }
            else {

            wrap = new TransitionWrap(closeTTlist, nodeAncState, childNode.isLeaf());
            }

            statesTooMany = (maxExtraFinal != null && wrap.transStatuses.size() > maxSumStates);
            if (statesTooMany) {
                --maxExtraFinal;
            }

        }
        return wrap;

    }

    public static Hashtable<Integer, AncStates> createStatesDict(beast.evolution.tree.TreeInterface tree, Alignment alinmt, List<Pair<Integer,Integer>> posSites,int nTargets) {
        Hashtable<Integer, AncStates> statesDict = new Hashtable<>();
        for (Node node : tree.listNodesPostOrder(null, null)) {

            if (node.isLeaf()) {

                //for now, only 1 site without insertions
                //TO DO: find out why alignment is not containing the sequences properly
                String leafSeq = alinmt.sequences.get(node.getNr()).getData();

                //
                //leafSeq = leafSeq.substring(7, leafSeq.length() - 1);

                AncStates leafState = AncStates.createObservedAlleleSet(leafSeq, posSites,nTargets);


                statesDict.put(node.getNr(), leafState);

                //create anc state set based on the observed allele
            } else if (node.isRoot()) {
                //root node: create an empty AncState, as the root has no ancestral states (unedited barcode)
                statesDict.put(node.getNr(), new AncStates());

            } else {
                //it is an internal node. The node's ancestral state is the intersection of its children's ancestral states
                List<Node> childrn = node.getChildren();
                AncStates parntStat = statesDict.get((childrn.get(0)).getNr());
                for (int i = 1; i< childrn.size();i ++) {
                    AncStates otherChild = statesDict.get((childrn.get(i)).getNr());
                    parntStat = AncStates.intersect(parntStat, otherChild);
                }

                /// HERE IS ASSUMED BIFURCATING:
                statesDict.put(node.getNr(), parntStat);

            }

        }
        return statesDict;
    }

    public static Hashtable<Integer, TransitionWrap> createTransitionWrappers(beast.evolution.tree.TreeInterface tree, Alignment alinmt, BarcodeMeta metaData) {

        // annotate tree possible ancestral states HERE, by creating a dictionary with node indices.
        Hashtable<Integer, AncStates> statesDict = createStatesDict(tree,alinmt,metaData.posSites, metaData.nTargets);
        // annotate tree possible ancestral states HERE, by creating a dictionary with node indices.
        Hashtable<Integer, TransitionWrap> wrapDict = new Hashtable<>();
        //getting the list of all tree nodes in postorder: correct order! This is to annotate the tree nodes with
        //possible ancestral states before pruning down those as Target Statuses.

        for (Node node : tree.listNodesPostOrder(null, null)) {
            //creating an empty transition wrap, avoiding to use NULL as we need empty entries( not null)
            List<List<IndelSet.TargetTract>> initNull = new ArrayList<>();
            List<IndelSet.TargetTract> innerinitNull = new ArrayList<>();
            initNull.add(innerinitNull);
            TransitionWrap temp = new TransitionWrap(initNull, statesDict.get(node.getNr()), node.isLeaf());
            wrapDict.put(node.getNr(), temp);
        }
        //create a parsimony state dictionary
        Hashtable<Integer, List<IndelSet.Singleton>> parsimDict = new Hashtable<>();
        for (Integer key : statesDict.keySet()) {
            parsimDict.put(key, statesDict.get(key).getSingletons());
        }




        //Now can fill the dictionary:
        //preorder list of nodes :
        List<Node> preorderList = Arrays.asList(tree.listNodesPostOrder(null, null));
        for (int reverseIt = preorderList.size() - 1; reverseIt >= 0; reverseIt--) {
            Node parentNode = preorderList.get(reverseIt);
            //debugging code:

            List<List<IndelSet.TargetTract>> parentNodeTuples = wrapDict.get(parentNode.getNr()).ttTuple;
            List<List<IndelSet.TargetTract>> filteredTTTuples = parentNodeTuples;
            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap filterWrap = getCloseTransitionWrapper(parentNode, childNode, statesDict, parsimDict, parentNodeTuples, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);



                filteredTTTuples = IndelSet.TargetTract.intersect(filteredTTTuples,filterWrap.ttTuple);

                //ATTEMPT AT REMOVING DUPLICATES
                Set ohneDup = new LinkedHashSet();
                ohneDup.addAll(filteredTTTuples);
                filteredTTTuples.clear();
                filteredTTTuples.addAll(ohneDup);




            }


            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap finalWrap = getCloseTransitionWrapper(parentNode, childNode, statesDict, parsimDict, filteredTTTuples, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);
                //TODO check Attempting to reverse the tttupes
                List<TargetStatus> cleanTTUPLES = finalWrap.transStatuses;
                Collections.reverse(cleanTTUPLES);
                /*Set ohneDup = new LinkedHashSet();
                ohneDup.addAll(cleanTTUPLES);
                cleanTTUPLES.clear();
                cleanTTUPLES.addAll(ohneDup);*/
                finalWrap.transStatuses = cleanTTUPLES;

                wrapDict.put(childNode.getNr(), finalWrap);

            }
        }
        return wrapDict;
    }

    public static List<List<IndelSet.TargetTract>> getCloseStates(int maxSteps,
                                                                  int nTargets,
                                                                  List<List<IndelSet.TargetTract>> parentTtList,
                                                                  AncStates ancStates,
                                                                  TargetStatus requestStatus) {

        if (maxSteps == 0) {
            for (IndelSet temp : ancStates.indelSetList) {
                assert (temp.getClass() == IndelSet.SingletonWC.class);
            }
            List<IndelSet.TargetTract> noStepTuple = new ArrayList<>();
            for (IndelSet temp : ancStates.getSingletonWCs()) {
                noStepTuple.add(temp.getTargetTract());
            }

            return new ArrayList(noStepTuple);
        }

        List<Integer> requestedInactive = requestStatus.getInactiveTargets(nTargets);
        List<Triple<Integer, List<IndelSet.TargetTract>, Boolean>> stateQueue = new ArrayList<>();
        Hashtable<Pair<List<IndelSet.TargetTract>, Boolean>, Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>>> fullDict = new Hashtable<>();
        //fill the a dictionary of states for the parent states, and prepare the state quueue
        for (List<IndelSet.TargetTract> temp : parentTtList) {

            TargetStatus parentStat = new TargetStatus();
            parentStat.fillWithTargetTracts(temp);
            List<Integer> parentInactive = parentStat.getInactiveTargets(nTargets);
            //is the requested status a subset of the parent statuses?
            Boolean isPassed = parentInactive.containsAll(requestedInactive);

            Pair<List<IndelSet.TargetTract>, Boolean> parentEntry = new Pair(temp, isPassed);
            Triple<Integer, List<IndelSet.TargetTract>, Boolean> initialEntry = new Triple(0, temp, isPassed);
            stateQueue.add(initialEntry);
            //TO DO: CHECK initial entry put in there!
            //making an empty list to correctly initialize:
            List<IndelSet.TargetTract> empty = new ArrayList<>();
            fullDict.put(parentEntry, new Hashtable() {{
                put(0, new Pair(empty, false));
            }});

        }
        //Max targets
        List<Integer> maxTargets = ancStates.getMaxTargetStatus().getInactiveTargets(nTargets);

        //find all paths of at most max steps long from the parent states
        while (stateQueue.size() != 0) {

            //The priority queue is currently dequeued in the wrong order
            List<Integer> allDist = new ArrayList<>();
            for (int i = 0; i < stateQueue.size(); i++) {
                allDist.add(stateQueue.get(i).a);
            }
            int dequeueDist = Collections.min(allDist);
            int indexDequeue = 0;
            for(int i = 0;i<allDist.size();++i) {
                if (allDist.get(i)==dequeueDist) {
                    indexDequeue = i;
                }
            }

            Triple<Integer, List<IndelSet.TargetTract>, Boolean> currentTriple = stateQueue.remove(indexDequeue);

            if (currentTriple.a >= maxSteps) {
                continue;
            }

            Integer dist = currentTriple.a;
            List<IndelSet.TargetTract> state = currentTriple.b;
            boolean statePassedReq = currentTriple.c;


            TargetStatus targStat = new TargetStatus();
            targStat.fillWithTargetTracts(currentTriple.b);

            //getting the difference between AncState status and the obtained targStat:
            List<Integer> difference = new ArrayList<>(maxTargets);
            difference.removeAll(targStat.getInactiveTargets(nTargets));
            List<Integer> activeAnyTargs = difference;
            //HERE IS THE ISTAKE: WRONG INDICES INTRODUCED:
            List<IndelSet.TargetTract> possibleTTracts = targStat.getPossibleTargetTracts(activeAnyTargs,nTargets);

            //TO DO THE POSSIBLE TTRACT CONTAINS AN IMPOSSIBLE TTRACT SEE ABOVE ^
            for (IndelSet.TargetTract possibleTTract : possibleTTracts) {
                //EMPTY STATE, what to do?
                //attempting to merge the possible TTtracts in the state!
                //check the possible tract is it ok:
                //check that both merged states are the same
                List<IndelSet.TargetTract> newState = IndelSet.TargetTract.merge(possibleTTract, state);

                if (!ancStates.isPossible(newState)) {
                    //THE MERGED STATE IS IMPOSSIBLE
                    continue;
                }

                boolean newStatePassedReq = statePassedReq;
                if (!statePassedReq) {
                    //if the previous state was not passed, check if the new one is!
                    TargetStatus newStatus = new TargetStatus();
                    newStatus.fillWithTargetTracts(newState);
                    List<Integer> newInactiveTarg = newStatus.getInactiveTargets(nTargets);
                    // is the new status a superset of the ancestral state(max target) ?
                    newStatePassedReq = newInactiveTarg.containsAll(maxTargets);

                }
                //add the new state to queue
                Triple<Integer, List<IndelSet.TargetTract>, Boolean> childEntry = new Triple<>(dist + 1, newState, newStatePassedReq);

                stateQueue.add(childEntry);
                Pair<List<IndelSet.TargetTract>, Boolean> childKey = new Pair<>(newState, newStatePassedReq);
                //if no entry is present, just add a new one


                if (!fullDict.containsKey(childKey)) {
                    fullDict.put(childKey, new Hashtable() {{
                        put(dist + 1, new Pair(state, statePassedReq));
                    }});
                }
                //if there is already the key, append the nested hashtable:
                if (fullDict.containsKey(childKey)) {
                    // append the nested hashmap instead
                    Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> appendedHashtable = fullDict.get(childKey);
                    List<Pair<List<IndelSet.TargetTract>, Boolean>> replacementList = new ArrayList<>();
                    if (appendedHashtable.containsKey(dist+1)) {
                        // we append the list and change the hashtable
                        try {
                            replacementList = fullDict.get(childKey).get(dist + 1);
                        }
                        catch (Exception e) {
                            Pair<List<IndelSet.TargetTract>, Boolean> replacementPair = (Pair<List<IndelSet.TargetTract>, Boolean>) fullDict.get(childKey).get(dist + 1);
                            replacementList.add(replacementPair);

                        }
                        replacementList.add(new Pair(state, statePassedReq));
                        Set ohneDup = new LinkedHashSet();
                        ohneDup.addAll(replacementList);
                        replacementList.clear();
                        replacementList.addAll(ohneDup);
                        appendedHashtable.put(dist+1,replacementList);
                    }
                    else {
                        //we just add a dictionary entry in the hashtable
                        replacementList.add(new Pair(state, statePassedReq));
                        appendedHashtable.put(dist+1,replacementList);
                    }
                    fullDict.put(childKey,appendedHashtable);

                    /*if(replacementList == null ) {
                        replacementList = new ArrayList<Pair<List<IndelSet.TargetTract>, Boolean>>();
                        replacementList.add(new Pair(state, statePassedReq));
                        appendedHashtable.put(dist+1,replacementList);
                        fullDict.put(childKey,appendedHashtable);
                    }
                    else {
                        replacementList.add(new Pair(state, statePassedReq));
                        appendedHashtable.put(dist + 1, replacementList);
                        fullDict.put(childKey, appendedHashtable);
                        //TRYING TO SEE WHTHER IT IS AN APPENDING ISSUE
                    }*/


                }
            }
        }
        //backtracking the paths:
        List<Pair<Pair<List<IndelSet.TargetTract>, Boolean>, Integer>> closeBy = new ArrayList<>();


        for (Map.Entry<Pair<List<IndelSet.TargetTract>, Boolean>, Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>>> entry : fullDict.entrySet()) {

            Pair<List<IndelSet.TargetTract>, Boolean> metaStateKey = entry.getKey();
            if (!metaStateKey.getSecond()) {
                continue;
            }

            Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> subDict = entry.getValue();
            //DEBUGGING LOOPS
            //We explore the sub states obtained:

            for (Integer stepCount : subDict.keySet()) {

                Queue<Triple<List<IndelSet.TargetTract>, Boolean, Integer>> ancestorQueue = new LinkedList<>();
                Triple forqueue = new Triple(metaStateKey.getFirst(), metaStateKey.getSecond(), stepCount);
                ancestorQueue.add(forqueue);
                while (!ancestorQueue.isEmpty()) {

                    Triple<List<IndelSet.TargetTract>, Boolean, Integer> ancestorEntry = ancestorQueue.remove();

                    Pair dictExtractkey = new Pair(ancestorEntry.a, ancestorEntry.b);
                    Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> attempt = fullDict.get(dictExtractkey);

                    if (ancestorEntry.c >= 1) {
                        if ((!closeBy.contains(ancestorEntry)) && (fullDict.get(dictExtractkey).containsKey(ancestorEntry.c))) {
                            Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> subdictio = fullDict.get(dictExtractkey);
                            try {
                                for (Pair<List<IndelSet.TargetTract>, Boolean> temp : subdictio.get(ancestorEntry.c)) {
                                    ancestorQueue.add(new Triple(temp.getFirst(), temp.getSecond(), ancestorEntry.c - 1));
                                }
                            } catch (Exception e) {
                                Pair<List<IndelSet.TargetTract>, Boolean> temp = (Pair<List<IndelSet.TargetTract>, Boolean>) subdictio.get(ancestorEntry.c);
                                ancestorQueue.add(new Triple(temp.getFirst(), temp.getSecond(), ancestorEntry.c - 1));
                            } finally {
                            }

                        }

                    }
                    //TO DO: simpler way of appending closeby correct with single entry list!!!!

                    Pair<Pair<List<IndelSet.TargetTract>, Boolean>, Integer> closeByIn = new Pair(dictExtractkey, ancestorEntry.c);
                    closeBy.add(closeByIn);
                }
            }
        }
        assert closeBy.size() > 0;
        List<List<IndelSet.TargetTract>> finalOut = new ArrayList<>();
        for (int i = 0; i < closeBy.size(); ++i) {
            finalOut.add(closeBy.get(i).getFirst().getFirst());
        }

        return finalOut;
    }
}








