package gestalt.evolution.simulation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.pkgmgmt.BEASTClassLoader;
import beast.pkgmgmt.PackageManager;
import gestalt.evolution.alignment.GestaltEvent;
import gestalt.evolution.alignment.IndelSet;
import gestalt.evolution.alignment.TargetStatus;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.apache.commons.math3.util.Pair;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static gestalt.evolution.alignment.GestaltEvent.intersect;
import static java.lang.Math.max;
import static java.lang.Math.random;


@Description("A GESTALT alignment simulator adapted from Tim Vaughan's feast implementation")
public class SimulatedGestaltAlignment extends Alignment {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree down which to simulate sequence evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "Site model to use in simulation.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Start of the process, usually the experiment",
            Input.Validate.OPTIONAL);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "Name of file (if any) simulated alignment should be saved to.");

    private Tree tree;
    private SiteModel siteModel;
    private gestaltGeneral substModel;
    private DataType dataType;
    private double originHeight;
    private String[] alignment;

    public SimulatedGestaltAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        siteModel = siteModelInput.get();
        substModel = (gestaltGeneral) siteModel.getSubstitutionModel();
        sequences.clear();

        if (originInput.get() != null) {
            originHeight = originInput.get().getValue();
        }

        grabDataType();
        simulate();

        //find out length of longest sequence to adjust such that all are the same length
        int maxEdits = 0;
        for (int i = 0; i < tree.getLeafNodeCount(); ++i) {
            if (alignment[i] == "") {
                maxEdits = max(maxEdits, 0);
            } else {
                maxEdits = max(maxEdits, alignment[i].split(",").length);
            }
        }
        //append the shorter sequences with "None"
        for (int i = 0; i < tree.getLeafNodeCount(); ++i) {
            while (alignment[i].split(",").length < maxEdits) {
                if (alignment[i] == "") {
                    alignment[i] = "None";
                } else {
                    alignment[i] = alignment[i] + "," + "None";
                }
            }
        }
        //write alignment to file
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(outputFileNameInput.get(), StandardCharsets.UTF_8);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < tree.getLeafNodeCount(); ++i) {
            writer.println("<sequence id=\"" + tree.getTaxaNames()[i] + "\"" + " spec=\"Sequence\" taxon=\"" + tree.getTaxaNames()[i] + "\"  value=" + alignment[i] + ",\"/>");
        }
        writer.close();

    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        int nTaxa = tree.getLeafNodeCount();

        alignment = new String[nTaxa];

        Node root = tree.getRoot();

        TargetStatus rootStatus = new TargetStatus();
        List<IndelSet.TargetTract> rootTargetTracts = new ArrayList();
        String rootAllele = "";

        if (originHeight != 0) {
            // then parent sequence is sequence at origin and we evolve sequence first down to the root
            double deltaT = originHeight - root.getHeight();
            double clockRate = siteModel.getRateForCategory(0, root);
            while (deltaT > 0) {

                //Draw a cut and a time at which it happens
                Pair<Double, IndelSet.TargetTract> outcome = this.raceTargetTracts(rootStatus, clockRate * deltaT);

                //if the barcode is already saturated (outcome == null)
                //end simulation on that branch
                if (outcome == null) {
                    deltaT = 0;

                } else {

                    //time of that next cut
                    Double eventTime = outcome.getFirst();

                    //if the time drawn for the next barcode cut event exceeds the branch length
                    //end the simulation
                    if (eventTime >= deltaT) {
                        deltaT = 0;
                    }

                    //if (eventTime < deltaT) there is a viable cut on the branch.
                    else  {
                        //we update the remaining time, allowing for potentially more cuts
                        deltaT -= eventTime;
                        //we update the target status with the cut
                        //draw a repair event
                        String event = this.doRepair(outcome.getSecond());
                        if (rootAllele == "") {
                            rootAllele = event;

                        } else {
                            rootAllele = rootAllele + "," + event;
                        }
                        rootStatus.addTarget(outcome.getSecond());
                        rootTargetTracts.add(outcome.getSecond());

                    }


                }
            }
        }

        traverse(root, rootStatus, rootAllele, rootTargetTracts);

    }

    /**
     * Traverse a tree, simulating a sequence alignment down it.
     *
     * @param node         Node of the tree
     * @param parentStatus Sequence at the parent node in the tree
     */
    private void traverse(Node node,
                          TargetStatus parentStatus, String parentAllele, List<IndelSet.TargetTract> parentTracts) {

        for (Node child : node.getChildren()) {

            double deltaT = node.getHeight() - child.getHeight();
            double clockRate = siteModel.getRateForCategory(0, child);

            // Draw characters on child sequence
            TargetStatus childStatus = parentStatus;
            String childAllele = parentAllele;
            List<IndelSet.TargetTract> childTracts = parentTracts;

            while (deltaT > 0) {
                Pair<Double, IndelSet.TargetTract> outcome = this.raceTargetTracts(childStatus, clockRate * deltaT);
                if (outcome == null) {
                    //outcome is null when the target gets fully edited
                    break;
                }
                Double eventTime = outcome.getFirst();
                if (eventTime < deltaT) {
                    //there is a cut, we update the remaining time, allowing for potentially more cuts
                    deltaT -= eventTime;
                    //we update the target status with the cut
                    childStatus.addTarget(outcome.getSecond());
                    childTracts.add(outcome.getSecond());
                    //here, we want to check wether there is any overlap
                    String event = this.doRepair(outcome.getSecond());
                    if (childAllele == "") {
                        childAllele = event;

                    } else {
                        childAllele = childAllele + "," + event;
                    }
                }
                if (eventTime >= deltaT) {
                    deltaT = 0;
                    //there is no cut on that branch.
                }

            }
            if (child.isLeaf()) {
                //make a function that observes allele: this function will merge events that are contiguous together.
                //String correctedAllele = observeAllele(childAllele);
                alignment[child.getNr()] = childAllele;
            } else {
                traverse(child, childStatus, childAllele, childTracts);
            }
        }
    }

    //this function is here to postprocess a simulated allele into an "observed" allele.
    //Simulated alleles can contain contiguous indels/masked indels. masked events must be removed. In real data, contiguous/overlapping indels are observed as
    // single indels meaning if evt1 and evt2 are contiguous, we sequence a single indel: evt1+2 (intersection of both)
    String observeAllele(String simulatedAllele) {
        //extract events
        String[] stringInput = simulatedAllele.split(",");
        List<GestaltEvent> rawEvents = new ArrayList<>();
        for(String i : stringInput) {
            rawEvents.add(new GestaltEvent(i));
        }
        List<GestaltEvent> correctedEvents = new ArrayList<>();
        List<GestaltEvent> preEvents = rawEvents;

        //intersect list with itself until stable:
//        while(preEvents != correctedEvents ) {
//            preEvents = correctedEvents;
//            correctedEvents = intersectList(correctedEvents,correctedEvents);
//
//        }
        //postprocess
        List<GestaltEvent> processed = processSet(rawEvents);
        String finalAllele ="";
        for(GestaltEvent i : processed) {
            finalAllele = finalAllele + i.toString() + ",";

        }

        return finalAllele;




    }

    //process longitudinally the allele to clean up potential overlaps
    private List<GestaltEvent> processSet(List<GestaltEvent> correctedEvents) {
        List<GestaltEvent> processed = new ArrayList<>();
        for ( int i = 0; i < correctedEvents.size();++i) {
            GestaltEvent eventi = correctedEvents.get(i);
            GestaltEvent intersectioni_j = null;
            for (int j = i + 1; j < correctedEvents.size();++j) {
                intersectioni_j = intersect(eventi,correctedEvents.get(j));

            }
            if(intersectioni_j != null && (! processed.contains(intersectioni_j))) {
                processed.add(intersectioni_j);
            }
            else {
                processed.add(eventi);
            }
        }
        return processed;


    }

    //intersection of 2 lists of GESTALT events to create an Event set with all possible intersections
    public static List<GestaltEvent> intersectList(List<GestaltEvent> first, List<GestaltEvent> second) {
        List<GestaltEvent> intersectionList = new ArrayList<>();
        if(first.size()==0 && second.size()==0) {
            return intersectionList;
        }
        else if (first.size()==0 && second.size()!=0 || first.size()!=0 && second.size()==0) {
            return null;
        }
        for ( GestaltEvent event1: first) {
            for (GestaltEvent event2: second) {
                GestaltEvent inter = intersect(event1,event2);
                if(inter != null && (! intersectionList.contains(inter))) {
                    intersectionList.add(inter);
                }
            }
        }
        return intersectionList;
    }

    /**
     * HORRIBLE function to identify data type from given description.
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = PackageManager.find(DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                        | InstantiationException
                        | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }

    /**
     * Repairs allele per the target_tract
     * todo: handle overlapping deletions
     * todo: handle alignment
     *
     * @return a String of the indel, in the following format startPos_DelLen_MinTarg_MaxTarg_InsertSequence
     */
    public String doRepair(IndelSet.TargetTract targetTract) {
        //left cut site
        int target1 = targetTract.getminTarg();
        //right cut site
        int target2 = targetTract.getmaxTarg();
        boolean isIntertarg = (target1 != target2);

        boolean leftLong = targetTract.isLeftLong();
        boolean rightLong = targetTract.isRightLong();


// Keep running this until at  least a least a trim ot insertion occurs

        boolean didSomething = false;
        boolean doInsertion = false;
        boolean doDeletionLeft = false;
        boolean doDeletionRight = false;
        while (!didSomething) {
            doInsertion = random() > substModel.insertZeroProb.getValue();
            doDeletionLeft = random() > substModel.trimZeroProbsDict.get(0, isIntertarg ? 1 : 0);
            doDeletionRight = random() > substModel.trimZeroProbsDict.get(1, isIntertarg ? 1 : 0);

            if (isIntertarg == true || rightLong || leftLong) {
                didSomething = true;
            } else {
                didSomething = doInsertion || doDeletionLeft || doDeletionRight;
            }
        }
        //Determine insertion sequence
        int insertionLength = 0;
        if (doInsertion) {//drawv
            insertionLength = substModel.insertDist.nextInt();
        }
        //draw a string from ATCG
        String[] insertNucleotides = new String[insertionLength];
        //{0},  // A
        //{1},  // C
        //{2},  // G
        //{3},  // T
        double[] frequencies = {0.25, 0.25, 0.25, 0.25};
        for (int i = 0; i < insertNucleotides.length; i++) {
            int key = Randomizer.randomChoicePDF(frequencies);
            if (key == 0)
                insertNucleotides[i] = "a";
            if (key == 1)
                insertNucleotides[i] = "c";
            if (key == 2)
                insertNucleotides[i] = "g";
            if (key == 3)
                insertNucleotides[i] = "t";
        }
        String insertSequence = String.join("" +
                "", insertNucleotides);
        //Determine deletion lengths
        int leftDelLen = 0;
        if (doDeletionLeft || leftLong) {

            //the deletion distribution must be bounded (upper bound) such that it doesn't affect the previous target.
            // there are nTrimTypes which is why there are several delLongDists, They correspond to the trimShortParams_reshaped[] and trimLongParams_reshaped
            // they correspond to the means of  deletion length distribution, left [0] and right [1] of the cut
            //todo check that these length distributions/parameters are correctely bounded
            if (leftLong) {
                leftDelLen = substModel.delLongDist.get(0).nextInt() + substModel.metaData.leftLongTrimMin.get(target1);
                while (leftDelLen > substModel.metaData.leftMaxTrim.get(target1)) {
                    leftDelLen = substModel.delLongDist.get(0).nextInt() + substModel.metaData.leftLongTrimMin.get(target1);
                }
            } else {
                leftDelLen = substModel.delLongDist.get(0).nextInt() + 1;
                while (leftDelLen > substModel.metaData.leftLongTrimMin.get(target1) - 1) {
                    leftDelLen = substModel.delLongDist.get(0).nextInt() + 1;
                }
            }
        }
        int rightDelLen = 0;
        if (doDeletionRight || rightLong) {
            //the deletion distribution must be bounded (upper bound) such that it doesn't affect the previous target.
            // if it is a long deletion, there is also a lower bound needed to reach the previous target (lower bound)
            if (rightLong) {
                rightDelLen = substModel.delLongDist.get(1).nextInt() + substModel.metaData.rightLongTrimMin.get(target2);
                while (rightDelLen > substModel.metaData.rightMaxTrim.get(target2)) {
                    rightDelLen = substModel.delLongDist.get(1).nextInt() + substModel.metaData.rightLongTrimMin.get(target2);
                }
            } else {
                rightDelLen = substModel.delLongDist.get(1).nextInt() + 1;
                while (rightDelLen > substModel.metaData.rightLongTrimMin.get(target1) - 1) {
                    rightDelLen = substModel.delLongDist.get(1).nextInt() + 1;
                }


            }

        }
        //todo check for start pos whether there is aalso a shift by +/-
        // 1
        String StartPos = Integer.toString((int) (substModel.metaData.absCutSites.get(target1) - leftDelLen));
        //todo check the intertarget effect (if there is a +1 shift or not!, this might be wrong by +/-2 length)
        int interTargLen = (int) ((int) substModel.metaData.absCutSites.get(target2) - substModel.metaData.absCutSites.get(target1));

        String DelLen = Integer.toString(leftDelLen + rightDelLen + interTargLen);
        String minDeac = Integer.toString(targetTract.getminTargDeac());
        String maxDeac = Integer.toString(targetTract.getmaxTargDeac());
        String eventInput = StartPos + "_" + DelLen + "_" + minDeac + "_" + maxDeac + "_" + insertSequence;


        //note changing the event input to simulate alleles (vs the above)
        eventInput = leftDelLen + "_" + target1 + "_" + insertSequence + "_" + target2 + "_" + rightDelLen;
        return eventInput;
    }


    public Pair<Double, IndelSet.TargetTract> raceTargetTracts(TargetStatus startTS, double scaleHazard) {

        List<Integer> activeTargets = startTS.getActiveTargets(10);
        int numActive = activeTargets.size();
        //if the barcode is fully edited out, cannot edit anymore, return null
        if (numActive == 0) {
            return null;
        } else {
            //dummy activeAny parameter:
            List<Integer> activeAny = new ArrayList<>();
            List<IndelSet.TargetTract> targetTracts = startTS.getPossibleTargetTracts(activeAny, 10);
            double[] hazards = new double[targetTracts.size()];
            for (int i = 0; i <= targetTracts.size() - 1; i++) {
                hazards[i] = substModel.targetTractHazards.get(substModel.targetTractDict.get(targetTracts.get(i))) * scaleHazard;

            }
            double hazardSum = Arrays.stream(hazards).sum();
            double[] normalizedHazards = hazards;
            for (int i = 0; i <= hazards.length - 1; i++) {
                normalizedHazards[i] = normalizedHazards[i] / hazardSum;

            }

            Double minTime = Randomizer.nextExponential(hazardSum);
            int chosenIndex = Randomizer.randomChoicePDF(normalizedHazards);
            IndelSet.TargetTract chosenTt = targetTracts.get(chosenIndex);
            Log.info.println(chosenIndex);
            Log.info.println(minTime);
            Log.info.println(chosenTt.hashCode());
            return new Pair(minTime, chosenTt);
        }


    }
}
