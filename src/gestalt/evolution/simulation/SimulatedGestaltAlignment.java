package gestalt.evolution.simulation;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.pkgmgmt.BEASTClassLoader;
import beast.pkgmgmt.PackageManager;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import gestalt.evolution.alignment.IndelSet;
import gestalt.evolution.alignment.TargetStatus;
import gestalt.evolution.substitutionmodel.gestaltGeneral;
import org.apache.commons.math3.util.Pair;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

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

        //adjust sequence length to be sure all match
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
            if (alignment[i].equals("")) {
                alignment[i] = "None";}
            while (alignment[i].split(",").length < maxEdits) {
                 alignment[i] = alignment[i] + "," + "None";
            }
        }

        for (int leafIdx = 0; leafIdx < tree.getLeafNodeCount(); leafIdx++) {
            String seqString = alignment[leafIdx];
            String taxonName;
            if (tree.getNode(leafIdx).getID() != null)
                taxonName = tree.getNode(leafIdx).getID();
            else
                taxonName = "t" + leafIdx;

            sequenceInput.setValue(new Sequence(taxonName, seqString), this);
        }

        super.initAndValidate();

        // Write simulated alignment to disk if required
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(new TaxonSet(this)));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + outputFileNameInput.get() + ".");
            }
        }

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
        String rootAllele = substModel.uneditedBarcodeInput.get();

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
                        rootStatus.addTarget(outcome.getSecond());
                        rootTargetTracts.add(outcome.getSecond());

                        //draw a repair event
                        String indel = this.doRepair(outcome.getSecond());
                        //apply it to the allele
                        rootAllele = applyIndel(rootAllele,indel);



                    }


                }
            }
        }

        traverse(root, rootStatus, rootAllele, rootTargetTracts);


    }


    //apply an indel to the barcode sequence
    String applyIndel(String barcodeSequence,String indel) {

        String[] strs = indel.split("_");
        int leftLen = Integer.parseInt(strs[0]);
        int leftCut = Integer.parseInt(strs[1]);
        String insSeq = strs[2];
        int rightCut = Integer.parseInt(strs[3]);
        int rightLen = Integer.parseInt(strs[4]);


        //first, split into left and right BARCODE strings (no regards to the cutsite
        int cutoffset = substModel.metaData.barcodeCutSite;
        String[] allele = barcodeSequence.split(" ");

        //left of the cut
        String[] alleleleft = Arrays.copyOfRange(allele, 0, leftCut * 2  +2);
        String leftStrings =  String.join(" ",alleleleft);

        //right of the cut
        String[] alleleright = Arrays.copyOfRange(allele, rightCut * 2 + 1, allele.length );
        String rightStrings = String.join(" ", alleleright);

        //from this, remove the offset for the left side to obtain a string until the cut
        String left_of_left_cut = leftStrings.substring(0,leftStrings.length() - cutoffset);

        //same to the right
        int length_cut_barcode = alleleright[0].length();
        String left_of_right_cut = rightStrings.substring(0, length_cut_barcode - cutoffset);
        String right_of_right_cut  = rightStrings.substring(length_cut_barcode - cutoffset,rightStrings.length());

        //if it is an intertarget deletion, there is also a central sequence
        String central = "";
        if (rightCut != leftCut) {
            String right_of_left_cut = leftStrings.substring(leftStrings.length() - cutoffset,leftStrings.length() );
            String[] allelecenter = Arrays.copyOfRange(allele, leftCut * 2 + 2, rightCut * 2 + 1);
            String centerStrings = String.join(" ", allelecenter);
            central = right_of_left_cut + " " + centerStrings + " " + left_of_right_cut;
        }

        //now create the deletions to the left:
        int deleted_left = 0;
        int index = left_of_left_cut.length();
        while(deleted_left != leftLen) {
            if(! left_of_left_cut.substring(index-1,index).equals(" ") ) {
                left_of_left_cut = left_of_left_cut.substring(0, index - 1) + "-" + left_of_left_cut.substring(index, left_of_left_cut.length());
                deleted_left += 1;
                index -= 1;
            }
            else {
                index -=1;
            }
        }

        int deleted_right = 0;
        index = 0;
        while(deleted_right != rightLen) {
            if(!right_of_right_cut.substring(index,index+1).equals(" ")) {
                right_of_right_cut = right_of_right_cut.substring(0, index) + "-" + right_of_right_cut.substring(index +1, right_of_right_cut.length());
                deleted_right += 1;
                index += 1;
            }
            else {
                index +=1;
            }
        }


        if(central != "") {
            central = central.replace("a", "");
            central = central.replace("c", "");
            central = central.replace("g", "");
            central = central.replace("t", "");
            central = central.replace("A", "-");
            central = central.replace("T", "-");
            central = central.replace("G", "-");
            central = central.replace("C", "-");
        }

        return left_of_left_cut + insSeq + central + right_of_right_cut;


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
            TargetStatus childStatus = new TargetStatus(parentStatus);

            String childAllele = parentAllele;

            List<IndelSet.TargetTract> childTracts = parentTracts.stream().map(IndelSet.TargetTract::new).collect(Collectors.toList());

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

                    String indel = this.doRepair(outcome.getSecond());
                    childAllele = applyIndel(childAllele,indel);

                }
                if (eventTime >= deltaT) {
                    deltaT = 0;
                    //there is no cut on that branch.
                }

            }
            if (child.isLeaf()) {
                //make a function that observes allele: this function will merge events that are contiguous together.
                String correctedAllele = observeAllele(childAllele);
                alignment[child.getNr()] = correctedAllele;

            } else {
                traverse(child, childStatus, childAllele, childTracts);


            }


        }



    }

    //this function is here to postprocess a simulated allele into an "observed" allele.
    //Simulated alleles can contain contiguous indels/masked indels. masked events must be removed. In real data, contiguous/overlapping indels are observed as
    // single indels meaning if evt1 and evt2 are contiguous, we sequence a single indel: evt1+2 (intersection of both)
    public String observeAllele(String simulatedAllele) {

        List<String> alleleIndelFormat = processAllele(simulatedAllele);
        List<String> alleleEventFormat = processEvents(alleleIndelFormat);
        return String.join(",", alleleEventFormat);


    }

     List<String> processEvents(List<String> processedAllele) {
        List<String> processedEvents = new ArrayList<>();
        for(String rawEvent : processedAllele) {
            String processedEvent = "";
            List<Integer> matchingTargets = new ArrayList<>();
            String[] splitevent = rawEvent.split("_");

            for(int targetindex= 0; targetindex < substModel.metaData.absCutSites.length; ++targetindex) {
                double event0 = Double.parseDouble(splitevent[0]);
                double event1 = Double.parseDouble(splitevent[1]);
                if (event0 <=  substModel.metaData.absCutSites.get(targetindex) && event1 >= substModel.metaData.absCutSites.get(targetindex)) {
                    matchingTargets.add(targetindex) ;
                }

            }
            //there is an insertion
            if(splitevent.length == 4) {
                processedEvent = splitevent[0] + "_" + splitevent[2] + "_" + matchingTargets.stream().mapToInt(v -> v).min().orElseThrow(NoSuchElementException::new) + "_" + matchingTargets.stream().mapToInt(v -> v).max().orElseThrow(NoSuchElementException::new) + "_" + splitevent[3];
                processedEvents.add(processedEvent);
            }
            //there is no insertion
            else {
                processedEvent = splitevent[0] + "_" + splitevent[2] + "_" + matchingTargets.stream().mapToInt(v -> v).min().orElseThrow(NoSuchElementException::new) + "_" + matchingTargets.stream().mapToInt(v -> v).max().orElseThrow(NoSuchElementException::new) + "_" ;
                processedEvents.add(processedEvent);
            }

        }

        return processedEvents;

    }

     List<String> processAllele(String Allele) {

        List<String> indelEvents = new ArrayList<>();
        char[] allele = Allele.toCharArray();
        int unedited_index = 0;
        int edited_index = 0;
        while(edited_index != allele.length -1) {
            //spacer position
            if (allele[edited_index] == ' ') {
                edited_index += 1;
            }
            //unedited postion
            if (allele[edited_index] == 'A' || allele[edited_index] == 'T' || allele[edited_index] == 'G' || allele[edited_index] == 'C') {
                edited_index += 1;
                unedited_index +=1;
            }

            //edited position, create an indel
            if (allele[edited_index] == '-') {
                //the start position is inclusive
                int startPos = unedited_index;
                int delLen = 0;
                String insert = "";
                //iterate until end of indel
                while (allele[edited_index] != 'A' && allele[edited_index] != 'T' && allele[edited_index] != 'G' && allele[edited_index] != 'C') {
                    if (allele[edited_index] == '-') {
                        delLen += 1;
                        unedited_index += 1;
                        edited_index += 1;
                    }
                    if (allele[edited_index] == ' ') {
                        edited_index += 1;
                    }
                    if (allele[edited_index] == 'a' || allele[edited_index] == 'c' || allele[edited_index] == 'g' || allele[edited_index] == 't') {
                        insert = insert + allele[edited_index];
                        edited_index += 1;

                    }
                }
                //the end pos is exclusive
                int endPos = unedited_index ;
                String indelEvent = startPos + "_" + endPos + "_" + delLen + "_" + insert;
                indelEvents.add(indelEvent);
            }
        }

        //String eventInput = StartPos + "_" + DelLen + "_" + minDeac + "_" + maxDeac + "_" + insertSequence;
        return indelEvents;
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

        //note changing the event input to simulate alleles (vs the above)
        String indelInput = leftDelLen + "_" + target1 + "_" + insertSequence + "_" + target2 + "_" + rightDelLen;
        return indelInput;
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
            return new Pair(minTime, chosenTt);
        }


    }
}
