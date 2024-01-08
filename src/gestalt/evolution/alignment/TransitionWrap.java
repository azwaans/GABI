package gestalt.evolution.alignment;

import beast.base.core.BEASTObject;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class TransitionWrap extends BEASTObject {

    /**
     * A transition wrap stores information for a specific node to calculate transition probabilities.
     * It stores meta-states (allele states grouped) that can precede the observed data
     * The meta-states here are TargetStatuses (binary representations of the editing state)
     */

    public void initAndValidate() {

    }

    /**
     * The possible target tract tuples that can be introduced
     */
    public List<List<IndelSet.TargetTract>> targetTractsTuples;

    /**
     * The possible statuses that precede the observed data at the leaves, generated by target tract introduction
     */
    public List<TargetStatus> transStatuses;

    /**
     * the mapping from state to state index number for this node
     */
    public Hashtable<TargetStatus, Integer> statusMap;

    /**
     * ancestral state at this node
     */
    public AncStates ancStates;

    /**
     * If this node is a leaf, then also store the leaf state
     */
    public TargetStatus leafState;
    public int numStatuses;
    public boolean isTreeLeaf;


    //constructor
    public TransitionWrap(List<List<IndelSet.TargetTract>> ttupleIn,
                          AncStates ancstateIn,
                          boolean isLeafin) {
        //input list of target tract tuples
        targetTractsTuples = ttupleIn;

        //translate each of these into statuses
        //ensure that there are no duplicate entries
        transStatuses = new ArrayList<>();
        Set<TargetStatus> noDuplicates = new HashSet<>();
        if (ttupleIn.size() != 0) {
            Set noDup = new LinkedHashSet();
            noDup.addAll(ttupleIn);
            ttupleIn.clear();
            ttupleIn.addAll(noDup);
            for (List<IndelSet.TargetTract> temp : ttupleIn) {
                TargetStatus tempStatus = new TargetStatus();
                tempStatus.fillWithTargetTracts(temp);
                noDuplicates.add(tempStatus);
            }

        }
        transStatuses.addAll(noDuplicates);

        //generate a mapping of state to state index
        statusMap = new Hashtable<>();
        for (int i = 0; i < transStatuses.size(); ++i) {
            statusMap.put(transStatuses.get(i), i);
        }

        ancStates = ancstateIn;
        isTreeLeaf = isLeafin;
        numStatuses = transStatuses.size();
        leafState = null;

        //if the node is a leaf, generate the status
        if (isLeafin) {
            leafState = ancStates.getMaxTargetStatus();
        }
        initAndValidate();
    }


    public void CreateStatusMap() {
        statusMap = new Hashtable<>();
        for (int i = 0; i < transStatuses.size(); ++i) {
            statusMap.put(transStatuses.get(i), i);
        }

    }

    //no argument constructor todo check if needed?
    public TransitionWrap() {
        targetTractsTuples = new ArrayList<>();
        transStatuses = new ArrayList<>();
        ancStates = null;
        isTreeLeaf = false;
        numStatuses = 0;
        leafState = null;
        initAndValidate();
    }


    /**
     * This class helps prune the set of states that we need to calculate transition probabilities for.
     */

    public static TransitionWrap getCloseTransitionWrap(Node parNode,
                                                        Node childNode,
                                                        Hashtable<Integer, AncStates> ancStatesDict,
                                                        Hashtable<Integer, List<IndelSet.Singleton>> parsimonyStatesDict,
                                                        List<List<IndelSet.TargetTract>> parentTuples,
                                                        Integer maxSumStates,
                                                        int maxExtraSteps,
                                                        int nTargets) {

        TransitionWrap wrap = new TransitionWrap();
        AncStates nodeAncState = ancStatesDict.get(childNode.getNr());
        List<IndelSet.Singleton> parentSingletonState = parsimonyStatesDict.get(parNode.getNr());


        List<IndelSet.Singleton> childSingletonState = parsimonyStatesDict.get(childNode.getNr());

        //finding length of the set difference between child and parent to find to minimum number of steps:
        List<IndelSet.Singleton> differenceSteps = new ArrayList<>(childSingletonState);
        differenceSteps.removeAll(parentSingletonState);
        int minSteps = differenceSteps.size();

        //prepare the desired child node status in the proper format (from singleton to status)
        TargetStatus requestedStatus = new TargetStatus();
        List<IndelSet.TargetTract> ttlist = new ArrayList<>();
        for (IndelSet.Singleton temp : childSingletonState) {
            ttlist.add(temp.getTargetTract());
        }
        requestedStatus.fillWithTargetTracts(ttlist);

        boolean statesTooMany = true;
        Integer maxExtraFinal = maxExtraSteps;

        if (maxSumStates != null) {
            if (!(Math.pow(2, minSteps) <= maxSumStates)) {
                maxExtraFinal = Math.max(maxExtraSteps - 1, 0);
            }
        }

        while (statesTooMany && (maxExtraFinal >= 0)) {

            List<List<IndelSet.TargetTract>> closeTTlist = getCloseStates((minSteps + maxExtraFinal), nTargets, parentTuples, nodeAncState, requestedStatus);

            if (closeTTlist.size() == 1) {
                List<IndelSet.TargetTract> empty = new ArrayList<>();
                if (closeTTlist.get(0) instanceof IndelSet.TargetTract) {
                    //only 1 element
                    IndelSet.TargetTract ttsingleton = (IndelSet.TargetTract) closeTTlist.get(0);
                    List<IndelSet.TargetTract> singleList = new ArrayList<IndelSet.TargetTract>();
                    singleList.add(ttsingleton);
                    List<List<IndelSet.TargetTract>> singleListList = new ArrayList<>();
                    singleListList.add(singleList);
                    wrap = new TransitionWrap(singleListList, nodeAncState, childNode.isLeaf());
                } else if (closeTTlist.get(0).size() == 0) {
                    //completely empty list
                    List<IndelSet.TargetTract> singleList = new ArrayList<IndelSet.TargetTract>();
                    List<List<IndelSet.TargetTract>> singleListList = new ArrayList<>();
                    singleListList.add(singleList);
                    wrap = new TransitionWrap(singleListList, nodeAncState, childNode.isLeaf());
                }
            } else {

                wrap = new TransitionWrap(closeTTlist, nodeAncState, childNode.isLeaf());
            }

            statesTooMany = (maxExtraFinal != null && wrap.transStatuses.size() > maxSumStates);
            if (statesTooMany) {
                --maxExtraFinal;
            }

        }
        return wrap;

    }

    /**
     * in original implementation: annotate_ancestral_states
     * finds all possible Ancestral States at internal nodes
     */
    public static Hashtable<Integer, AncStates> createStatesDict(beast.base.evolution.tree.TreeInterface tree,
                                                                 Alignment alinmt,
                                                                 List<Pair<Integer, Integer>> posSites,
                                                                 int nTargets) {
        Hashtable<Integer, AncStates> statesDict = new Hashtable<>();
        for (Node node : tree.listNodesPostOrder(null, null)) {

            if (node.isLeaf()) {


                String leafSeq = alinmt.sequenceInput.get().get(node.getNr()).toString();
                AncStates leafState = AncStates.createObservedAlleleSet(leafSeq, posSites, nTargets);
                statesDict.put(node.getNr(), leafState);

                //create anc state set based on the observed allele
            } else if (node.isRoot()) {
                //root node: create an empty AncState, as the root has no ancestral states (unedited barcode)
                statesDict.put(node.getNr(), new AncStates());

            } else {
                //it is an internal node. The node's ancestral state is the intersection of its children's ancestral states
                List<Node> childrn = node.getChildren();
                AncStates parntStat = statesDict.get((childrn.get(0)).getNr());
                for (int i = 1; i < childrn.size(); i++) {
                    AncStates otherChild = statesDict.get((childrn.get(i)).getNr());
                    parntStat = AncStates.intersect(parntStat, otherChild);
                }

                statesDict.put(node.getNr(), parntStat);

            }

        }
        return statesDict;
    }

    /**
     * Creates a transition wrap for all nodes in the tree
     */

    public static Hashtable<Integer, TransitionWrap> createTransitionWraps(beast.base.evolution.tree.TreeInterface tree,
                                                                           BarcodeMeta metaData,
                                                                           Hashtable<Integer,
                                                                                   AncStates> statesDict) {

        Hashtable<Integer, TransitionWrap> wrapDict = new Hashtable<>();

        for (Node node : tree.listNodesPostOrder(null, null)) {
            //creating an empty transition wrap, avoiding to use NULL as we need empty entries( not null)
            List<List<IndelSet.TargetTract>> initNull = new ArrayList<>();
            List<IndelSet.TargetTract> innerinitNull = new ArrayList<>();
            initNull.add(innerinitNull);
            TransitionWrap temp = new TransitionWrap(initNull, statesDict.get(node.getNr()), node.isLeaf());
            wrapDict.put(node.getNr(), temp);
        }
        //dictionary of all singleton states (not node assigned)
        Hashtable<Integer, List<IndelSet.Singleton>> parsimDict = new Hashtable<>();
        for (Integer key : statesDict.keySet()) {
            parsimDict.put(key, statesDict.get(key).getSingletons());
        }

        //traverse the tree to fill up the dictionary of wraps
        List<Node> preorderList = Arrays.asList(tree.listNodesPostOrder(null, null));
        for (int reverseIt = preorderList.size() - 1; reverseIt >= 0; reverseIt--) {
            Node parentNode = preorderList.get(reverseIt);

            List<List<IndelSet.TargetTract>> parentNodeTuples = wrapDict.get(parentNode.getNr()).targetTractsTuples;
            List<List<IndelSet.TargetTract>> filteredtargetTractsTupless = parentNodeTuples;

            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap filterWrap = getCloseTransitionWrap(parentNode, childNode, statesDict, parsimDict, parentNodeTuples, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);
                filteredtargetTractsTupless = IndelSet.TargetTract.intersect(filteredtargetTractsTupless, filterWrap.targetTractsTuples);
                //removing duplicates
                Set noDup = new LinkedHashSet();
                noDup.addAll(filteredtargetTractsTupless);
                filteredtargetTractsTupless.clear();
                filteredtargetTractsTupless.addAll(noDup);

            }

            for (Node childNode : parentNode.getChildren()) {

                TransitionWrap finalWrap = getCloseTransitionWrap(parentNode, childNode, statesDict, parsimDict, filteredtargetTractsTupless, metaData.maxSumSteps, metaData.maxExtraSteps, metaData.nTargets);
                List<TargetStatus> cleanTTUPLES = finalWrap.transStatuses;
                Collections.reverse(cleanTTUPLES);
                finalWrap.transStatuses = cleanTTUPLES;
                wrapDict.put(childNode.getNr(), finalWrap);

            }
        }
        return wrapDict;
    }

    /**
     * Given possible ancestral states at a node, this returns target tract tuples that are within maxSteps of the parent statuses, while reaching a specific target status
     */
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

        //flag for old
        //List<Triple<Integer, List<IndelSet.TargetTract>, Boolean>> stateQueue = new ArrayList<>();
        //ATTEMPTNG TO USE THE priority queue implementation instead of the homemade list
        PriorityQueue<QueueObject> statequeue = new PriorityQueue<>();

        Hashtable<StateAlongPath, Hashtable<Integer, List<StateAlongPath>>> stateToParentDict = new Hashtable<>();
        //fill the a dictionary of states for the parent states, and prepare the state quueue


        //attempt at tthis loop with prioruty qwueue:

        for (List<IndelSet.TargetTract> temp : parentTtList) {
            TargetStatus parentStat = new TargetStatus();
            parentStat.fillWithTargetTracts(temp);
            List<Integer> parentInactive = parentStat.getInactiveTargets(nTargets);
            //is the requested status a subset of the parent statuses?
            Boolean isPassed = parentInactive.containsAll(requestedInactive);

            StateAlongPath parentEntry = new StateAlongPath(temp, isPassed);

            //flag for old
            //Triple<Integer, List<IndelSet.TargetTract>, Boolean> initialEntry = new Triple(0, temp, isPassed);
            QueueObject initialEntry = new QueueObject(0, temp, isPassed);
            statequeue.add(initialEntry);
            //TO DO: CHECK initial entry put in there!
            //making an empty list to correctly initialize:
            List<IndelSet.TargetTract> empty = new ArrayList<>();
            stateToParentDict.put(parentEntry, new Hashtable() {{
                put(0, new StateAlongPath(empty, false));
            }});

        }


//        for (List<IndelSet.TargetTract> temp : parentTtList) {
//
//            TargetStatus parentStat = new TargetStatus();
//            parentStat.fillWithTargetTracts(temp);
//            List<Integer> parentInactive = parentStat.getInactiveTargets(nTargets);
//            //is the requested status a subset of the parent statuses?
//            Boolean isPassed = parentInactive.containsAll(requestedInactive);
//
//            Pair<List<IndelSet.TargetTract>, Boolean> parentEntry = new Pair(temp, isPassed);
//            Triple<Integer, List<IndelSet.TargetTract>, Boolean> initialEntry = new Triple(0, temp, isPassed);
//            stateQueue.add(initialEntry);
//            //TO DO: CHECK initial entry put in there!
//            //making an empty list to correctly initialize:
//            List<IndelSet.TargetTract> empty = new ArrayList<>();
//            stateToParentDict.put(parentEntry, new Hashtable() {{
//                put(0, new Pair(empty, false));
//            }});
//
//        }
        //Max targets
        List<Integer> maxTargets = ancStates.getMaxTargetStatus().getInactiveTargets(nTargets);


        while (!statequeue.isEmpty()) {

//            //The priority queue is dequeued in the wrong order
//            List<Integer> allDist = new ArrayList<>();
//            for (int i = 0; i < stateQueue.size(); i++) {
//                allDist.add(stateQueue.get(i).a);
//            }
//            int dequeueDist = Collections.min(allDist);
//            int indexDequeue = 0;
//            for(int i = 0;i<allDist.size();++i) {
//                if (allDist.get(i)==dequeueDist) {
//                    indexDequeue = i;
//                }
//            }

            QueueObject current = statequeue.remove();

            if (current.distance >= maxSteps) {
                continue;
            }

            Integer dist = current.distance;
            List<IndelSet.TargetTract> state = current.state.targetTactList;
            boolean statePassedReq = current.state.statePassedRequested;


            TargetStatus targStat = new TargetStatus();
            targStat.fillWithTargetTracts(current.state.targetTactList);

            //getting the difference between AncState status and the obtained targStat:
            List<Integer> difference = new ArrayList<>(maxTargets);
            difference.removeAll(targStat.getInactiveTargets(nTargets));
            List<Integer> activeAnyTargs = difference;

            List<IndelSet.TargetTract> possibleTTracts = targStat.getPossibleTargetTracts(activeAnyTargs, nTargets);


            for (IndelSet.TargetTract possibleTTract : possibleTTracts) {
                //attempting to merge the possible TTtracts in the state!
                //check the possible tract is ok:
                //check that both merged states are the same
                List<IndelSet.TargetTract> newState = IndelSet.TargetTract.merge(possibleTTract, state);

                if (!ancStates.isPossible(newState)) {
                    //the merged state is impossible
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
                QueueObject childEntry = new QueueObject(dist + 1, newState, newStatePassedReq);

                statequeue.add(childEntry);
                StateAlongPath childKey = new StateAlongPath(newState, newStatePassedReq);
                //if no entry is present, just add a new one


                if (!stateToParentDict.containsKey(childKey)) {
                    stateToParentDict.put(childKey, new Hashtable() {{
                        put(dist + 1, new StateAlongPath(state, statePassedReq));
                    }});
                }
                //if there is already the key, append the nested hashtable:
                if (stateToParentDict.containsKey(childKey)) {
                    // append the nested hashmap instead
                    Hashtable<Integer, List<StateAlongPath>> appendedHashtable = stateToParentDict.get(childKey);
                    List<StateAlongPath> replacementList = new ArrayList<>();
                    if (appendedHashtable.containsKey(dist + 1)) {
                        // we append the list and change the hashtable
                        try {
                            replacementList = stateToParentDict.get(childKey).get(dist + 1);
                        } catch (Exception e) {
                            StateAlongPath replacementPair = (StateAlongPath) stateToParentDict.get(childKey).get(dist + 1);
                            replacementList.add(replacementPair);

                        }
                        replacementList.add(new StateAlongPath(state, statePassedReq));
                        Set noDup = new LinkedHashSet();
                        noDup.addAll(replacementList);
                        replacementList.clear();
                        replacementList.addAll(noDup);
                        appendedHashtable.put(dist + 1, replacementList);
                    } else {
                        //we just add a dictionary entry in the hashtable
                        replacementList.add(new StateAlongPath(state, statePassedReq));
                        appendedHashtable.put(dist + 1, replacementList);
                    }
                    stateToParentDict.put(childKey, appendedHashtable);


                }
            }
        }
//flag for old
        //find all paths of at most max steps long from the parent states
//        while (stateQueue.size() != 0) {
//
//            //The priority queue is dequeued in the wrong order
//            List<Integer> allDist = new ArrayList<>();
//            for (int i = 0; i < stateQueue.size(); i++) {
//                allDist.add(stateQueue.get(i).a);
//            }
//            int dequeueDist = Collections.min(allDist);
//            int indexDequeue = 0;
//            for(int i = 0;i<allDist.size();++i) {
//                if (allDist.get(i)==dequeueDist) {
//                    indexDequeue = i;
//                }
//            }
//
//            Triple<Integer, List<IndelSet.TargetTract>, Boolean> currentTriple = stateQueue.remove(indexDequeue);
//
//            if (currentTriple.a >= maxSteps) {
//                continue;
//            }
//
//            Integer dist = currentTriple.a;
//            List<IndelSet.TargetTract> state = currentTriple.b;
//            boolean statePassedReq = currentTriple.c;
//
//
//            TargetStatus targStat = new TargetStatus();
//            targStat.fillWithTargetTracts(currentTriple.b);
//
//            //getting the difference between AncState status and the obtained targStat:
//            List<Integer> difference = new ArrayList<>(maxTargets);
//            difference.removeAll(targStat.getInactiveTargets(nTargets));
//            List<Integer> activeAnyTargs = difference;
//
//            List<IndelSet.TargetTract> possibleTTracts = targStat.getPossibleTargetTracts(activeAnyTargs,nTargets);
//
//
//            for (IndelSet.TargetTract possibleTTract : possibleTTracts) {
//                //attempting to merge the possible TTtracts in the state!
//                //check the possible tract is ok:
//                //check that both merged states are the same
//                List<IndelSet.TargetTract> newState = IndelSet.TargetTract.merge(possibleTTract, state);
//
//                if (!ancStates.isPossible(newState)) {
//                    //the merged state is impossible
//                    continue;
//                }
//
//                boolean newStatePassedReq = statePassedReq;
//                if (!statePassedReq) {
//                    //if the previous state was not passed, check if the new one is!
//                    TargetStatus newStatus = new TargetStatus();
//                    newStatus.fillWithTargetTracts(newState);
//                    List<Integer> newInactiveTarg = newStatus.getInactiveTargets(nTargets);
//                    // is the new status a superset of the ancestral state(max target) ?
//                    newStatePassedReq = newInactiveTarg.containsAll(maxTargets);
//
//                }
//                //add the new state to queue
//                Triple<Integer, List<IndelSet.TargetTract>, Boolean> childEntry = new Triple<>(dist + 1, newState, newStatePassedReq);
//
//                stateQueue.add(childEntry);
//                Pair<List<IndelSet.TargetTract>, Boolean> childKey = new Pair<>(newState, newStatePassedReq);
//                //if no entry is present, just add a new one
//
//
//                if (!stateToParentDict.containsKey(childKey)) {
//                    stateToParentDict.put(childKey, new Hashtable() {{
//                        put(dist + 1, new Pair(state, statePassedReq));
//                    }});
//                }
//                //if there is already the key, append the nested hashtable:
//                if (stateToParentDict.containsKey(childKey)) {
//                    // append the nested hashmap instead
//                    Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> appendedHashtable = stateToParentDict.get(childKey);
//                    List<Pair<List<IndelSet.TargetTract>, Boolean>> replacementList = new ArrayList<>();
//                    if (appendedHashtable.containsKey(dist+1)) {
//                        // we append the list and change the hashtable
//                        try {
//                            replacementList = stateToParentDict.get(childKey).get(dist + 1);
//                        }
//                        catch (Exception e) {
//                            Pair<List<IndelSet.TargetTract>, Boolean> replacementPair = (Pair<List<IndelSet.TargetTract>, Boolean>) stateToParentDict.get(childKey).get(dist + 1);
//                            replacementList.add(replacementPair);
//
//                        }
//                        replacementList.add(new Pair(state, statePassedReq));
//                        Set noDup = new LinkedHashSet();
//                        noDup.addAll(replacementList);
//                        replacementList.clear();
//                        replacementList.addAll(noDup);
//                        appendedHashtable.put(dist+1,replacementList);
//                    }
//                    else {
//                        //we just add a dictionary entry in the hashtable
//                        replacementList.add(new Pair(state, statePassedReq));
//                        appendedHashtable.put(dist+1,replacementList);
//                    }
//                    stateToParentDict.put(childKey,appendedHashtable);
//
//
//                }
//            }
//        }
        //backtracking the paths:
//        List<Pair<Pair<List<IndelSet.TargetTract>, Boolean>, Integer>> closeBy = new ArrayList<>();
//
//
//        for (Map.Entry<Pair<List<IndelSet.TargetTract>, Boolean>, Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>>> entry : stateToParentDict.entrySet()) {
//
//            Pair<List<IndelSet.TargetTract>, Boolean> metaStateKey = entry.getKey();
//            if (!metaStateKey.getSecond()) {
//                continue;
//            }
//
//            Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> subDict = entry.getValue();
//
//            for (Integer stepCount : subDict.keySet()) {
//
//                Queue<Triple<List<IndelSet.TargetTract>, Boolean, Integer>> ancestorQueue = new LinkedList<>();
//                Triple forqueue = new Triple(metaStateKey.getFirst(), metaStateKey.getSecond(), stepCount);
//                ancestorQueue.add(forqueue);
//                while (!ancestorQueue.isEmpty()) {
//
//                    Triple<List<IndelSet.TargetTract>, Boolean, Integer> ancestorEntry = ancestorQueue.remove();
//
//                    Pair dictExtractkey = new Pair(ancestorEntry.a, ancestorEntry.b);
//                    Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> attempt = stateToParentDict.get(dictExtractkey);
//
//                    if (ancestorEntry.c >= 1) {
//                        if ((!closeBy.contains(ancestorEntry)) && (stateToParentDict.get(dictExtractkey).containsKey(ancestorEntry.c))) {
//                            Hashtable<Integer, List<Pair<List<IndelSet.TargetTract>, Boolean>>> subdictio = stateToParentDict.get(dictExtractkey);
//                            try {
//                                for (Pair<List<IndelSet.TargetTract>, Boolean> temp : subdictio.get(ancestorEntry.c)) {
//                                    ancestorQueue.add(new Triple(temp.getFirst(), temp.getSecond(), ancestorEntry.c - 1));
//                                }
//                            } catch (Exception e) {
//                                Pair<List<IndelSet.TargetTract>, Boolean> temp = (Pair<List<IndelSet.TargetTract>, Boolean>) subdictio.get(ancestorEntry.c);
//                                ancestorQueue.add(new Triple(temp.getFirst(), temp.getSecond(), ancestorEntry.c - 1));
//                            } finally {
//                            }
//
//                        }
//                    }
//
//                    Pair<Pair<List<IndelSet.TargetTract>, Boolean>, Integer> closeByIn = new Pair(dictExtractkey, ancestorEntry.c);
//                    closeBy.add(closeByIn);
//                }
//            }
//        }
//        assert closeBy.size() > 0;
//        List<List<IndelSet.TargetTract>> finalOut = new ArrayList<>();
//        for (int i = 0; i < closeBy.size(); ++i) {
//            finalOut.add(closeBy.get(i).getFirst().getFirst());
//        }
//
//        return finalOut;
//    }
        List<Pair<StateAlongPath, Integer>> closeBy = new ArrayList<>();


        for (Map.Entry<StateAlongPath, Hashtable<Integer, List<StateAlongPath>>> entry : stateToParentDict.entrySet()) {

            StateAlongPath metaStateKey = entry.getKey();
            if (!metaStateKey.statePassedRequested) {
                continue;
            }

            Hashtable<Integer, List<StateAlongPath>> subDict = entry.getValue();

            for (Integer stepCount : subDict.keySet()) {

                Queue<QueueObject> ancestorQueue = new LinkedList<>();
                QueueObject forqueue = new QueueObject(stepCount, metaStateKey.targetTactList, metaStateKey.statePassedRequested);
                ancestorQueue.add(forqueue);
                while (!ancestorQueue.isEmpty()) {

                    QueueObject ancestorEntry = ancestorQueue.remove();
                    StateAlongPath dictExtractkey = new StateAlongPath(ancestorEntry.state.targetTactList, ancestorEntry.state.statePassedRequested);
                    if (ancestorEntry.distance >= 1) {
                        if ((!closeBy.contains(ancestorEntry)) && ((stateToParentDict.get(dictExtractkey)).containsKey(ancestorEntry.distance))) {
                            Hashtable<Integer, List<StateAlongPath>> subdictio = stateToParentDict.get(dictExtractkey);
                            try {
                                for (StateAlongPath temp : subdictio.get(ancestorEntry.distance)) {
                                    ancestorQueue.add(new QueueObject(ancestorEntry.distance - 1, temp.targetTactList, temp.statePassedRequested));
                                }
                            } catch (Exception e) {
                                StateAlongPath temp = (StateAlongPath) subdictio.get(ancestorEntry.distance);
                                ancestorQueue.add(new QueueObject(ancestorEntry.distance - 1, temp.targetTactList, temp.statePassedRequested));
                            } finally {
                            }

                        }
                    }

                    Pair<StateAlongPath, Integer> closeByIn = new Pair(dictExtractkey, ancestorEntry.distance);
                    closeBy.add(closeByIn);
                }
            }
        }
        assert closeBy.size() > 0;
        List<List<IndelSet.TargetTract>> finalOut = new ArrayList<>();
        for (int i = 0; i < closeBy.size(); ++i) {
            finalOut.add(closeBy.get(i).getFirst().targetTactList);
        }

        return finalOut;
    }
}
