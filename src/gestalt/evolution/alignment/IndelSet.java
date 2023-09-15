package gestalt.evolution.alignment;

import beast.base.core.BEASTObject;

import java.math.BigInteger;
import java.util.*;

public abstract class IndelSet extends BEASTObject {

    /**
     * A superclass for things that represent sets of Indels
     */

    /**
     * min target cut by indel
     */
    protected int minTarg;
    /**
     * min target cut by indel
     */
    protected int maxTarg;

    /**
     * get the inner wildcard of any indel set
     */
    public abstract WildCard getInnerWC();

    @Override
    public boolean equals(Object other){
     if(this == other) {
         return true;
     }
     if(other == null || this.getClass() != other.getClass()) {
         return false;
     }

     IndelSet otherSet = (IndelSet) other;
    return  (this.minTarg == otherSet.minTarg) && (this.maxTarg == otherSet.maxTarg) && (this.getminTargDeac() == otherSet.getminTargDeac())
            && (this.getmaxTargDeac() == otherSet.getmaxTargDeac());
    }

    public int hashCode(){
        String concatenated = "" + minTarg + maxTarg + getminTargDeac() + getmaxTargDeac();
        int result = Integer.parseInt(concatenated);
        return result;

    }

    /**
     * get the singleton associated with any indel
     */
    public abstract Singleton getSingleton();

    /**
     * get the min target deactivated by said indel
     */
    public abstract int getminTargDeac();

    /**
     * get the min target cut by said indel
     */
    public int getminTarg() {
        return minTarg;
    }

    /**
     * get the max target cut by said indel
     */
    public int getmaxTarg(){
        return maxTarg;

    }

    /**
     * get the max target deactivated by said indel
     */
    public abstract int getmaxTargDeac();


    /**
     * get the target tract induced by said indel
     */
    public abstract TargetTract getTargetTract();

    /**
     * Obtain the intersection of 2 indel sets.
     */
    public static IndelSet intersect(IndelSet set1, IndelSet set2) {
        //We first check whether sets are identical

        if(set1.equals(set2)) {
            return set1;
        }
        /*if the indel sets are different, first, we check whether one is fully contained in the other (within the other
        set's inner wild-card. In that case the intersection is the indel set itself.
        */
        else {
            WildCard wild1 = set1.getInnerWC();
            WildCard wild2 = set2.getInnerWC();
            if ((wild1 != null) && (wild1.minTarg <= set2.getminTargDeac()) && (set2.getmaxTargDeac() <= wild1.maxTarg)) {
                return set2;
            }
            else if ((wild2 != null) && (wild2.minTarg <= set1.getminTargDeac()) && (set1.getmaxTargDeac() <= wild2.maxTarg)) {
                return set1;
            }
            else {
                //otherwise, we get the intersection of both WCs
                return IndelSet.WildCard.intersect(wild1,wild2);
            }

        }






    }


    public static class SingletonWC extends IndelSet {

        /**
         * A singleton wildcard is the union of a singleton and its inner wildcard. It represents all underlying possible indels
         * (inlc. masked) identifiable with an event in a sequence.
         */

        protected int startPos;
        protected int delLeng;
        protected int minTargDeac;
        protected int maxTargDeac;
        protected String insertSeq;

        @Override
        public boolean equals(Object other){
            if(this == other) {
                return true;
            }
            if(other == null || this.getClass() != other.getClass()) {
                return false;
            }

            SingletonWC otherSet = (SingletonWC) other;
            return  ((this.minTarg == otherSet.minTarg) && (this.maxTarg == otherSet.maxTarg) && (this.getminTargDeac() == otherSet.getminTargDeac())
                    && (this.getmaxTargDeac() == otherSet.getmaxTargDeac()) && (this.insertSeq == otherSet.insertSeq));
        }


        public SingletonWC(int minTrgt, int maxTrgt, int delLgth, int strtPos, int minTrgtDeac, int maxTrgtDeac, String insrtSeq) {
            startPos = strtPos;
            delLeng = delLgth;
            minTarg = minTrgt;
            maxTarg = maxTrgt;
            minTargDeac = minTrgtDeac;
            maxTargDeac = maxTrgtDeac;
            insertSeq = insrtSeq;
            initAndValidate();
        }


        public void initAndValidate() {
        } // initAndValidate

        public WildCard getInnerWC() {
            if (maxTarg - 1 >= minTarg + 1) {

                return new WildCard(minTarg + 1, maxTarg - 1);
            } else {
                return null;
            }
        }

        public TargetTract getTargetTract() {
            return new TargetTract(minTarg, maxTarg, minTargDeac, maxTargDeac);
        }


        public Singleton getSingleton() {
            return new Singleton(startPos, delLeng, minTarg, maxTarg, minTargDeac, maxTargDeac, insertSeq);

        }

        public int getminTargDeac() {
            return minTargDeac;
        }

        public int getmaxTargDeac() {
            return maxTargDeac;
        }




    }

    public static class TargetTract extends IndelSet {
        /**
         * The set of indel that cut and deactivate the same target(s) minTargDeac to maxTargDeac
         */
        protected int minTargDeac;
        protected int maxTargDeac;


        public TargetTract(int minTrgt, int maxTrgt, int minTrgtDeac, int maxTrgtDeac) {
            minTarg = minTrgt;
            maxTarg = maxTrgt;
            minTargDeac = minTrgtDeac;
            maxTargDeac = maxTrgtDeac;
            initAndValidate();

        }

        public void initAndValidate() {

        }

        public WildCard getInnerWC() {
            return null;
        }

        public TargetTract getTargetTract() {return this;}

        public Singleton getSingleton() {
            return null;
        }

        public int getminTargDeac() {
            return minTargDeac;
        }

        public int getmaxTargDeac() {
            return maxTargDeac;
        }

        public boolean isLeftLong() {return minTargDeac  != minTarg;}
        public boolean isRightLong() {return maxTarg  != maxTargDeac;}

        public static List<TargetTract> merge(TargetTract tractIn,List<TargetTract> listIn) {


            //returns flattened version of a list of tuples of target tract, taking care of masking cases
            if(listIn.size() == 0 ) {
                ArrayList<TargetTract> noMergeNeeded = new ArrayList<>();
                noMergeNeeded.add(tractIn);
                return noMergeNeeded;
            }
            List<TargetTract> reducedList = new ArrayList<>(listIn);
            reducedList.add(tractIn);
            reducedList.sort(Comparator.comparing(a -> a.minTarg));



            List<TargetTract> finalTTlist = new ArrayList<>();
            finalTTlist.add(0,(reducedList.get(0)));
            for(int i = 1; i < reducedList.size(); i++) {
                if (finalTTlist.get(finalTTlist.size()-1).maxTargDeac > reducedList.get(i).maxTargDeac) {
                    continue;
                }
                else {
                    assert (finalTTlist.get(finalTTlist.size()-1).maxTarg < reducedList.get(i).minTarg) ;
                    finalTTlist.add(reducedList.get(i));
                }
            }
            return finalTTlist;
        }

        public static List<List<TargetTract>> intersect(List<List<TargetTract>> firstListList,List<List<TargetTract>> secondListList) {

            Set<List<TargetTract>> set1 = new HashSet<>();
            set1.addAll(firstListList);
            Set<List<TargetTract>> set2 = new HashSet<>();
            set2.addAll(secondListList);
            set1.retainAll(set2);


            List<List<TargetTract>> set1List = new ArrayList<>();
            set1List.addAll(set1);
            return set1List;
        }

        public static List<TargetTract> intersectList(List<TargetTract> first, List<TargetTract> second) {
            List<TargetTract> intersectionList = new ArrayList<>();
            if(first.size()==0 && second.size()==0) {
                return intersectionList;
            }
            else if (first.size()==0 && second.size()!=0 || first.size()!=0 && second.size()==0) {
                return null;
            }
            for ( TargetTract target1: first) {
                for (TargetTract target2: second) {
                    IndelSet inter = intersect(target1,target2);
                    if(inter != null) {
                        intersectionList.add(inter.getTargetTract());

                    }
                }
            }
            return intersectionList;
        }




    }

    public static class Singleton extends IndelSet {
        /**
         * This represents an indel set containing a single indel event
         */

        protected int startPos;
        protected int delLen;
        protected int minTargDeac;
        protected int maxTargDeac;
        protected String insSeq;


        public Singleton(int strtPos, int delLn, int minTrgt, int maxTrgt, int minTrgtDeac, int maxTrgtDeac, String insrtSeq) {
            startPos = strtPos;
            delLen = delLn;
            minTarg = minTrgt;
            maxTarg = maxTrgt;
            minTargDeac = minTrgtDeac;
            maxTargDeac = maxTrgtDeac;
            insSeq = insrtSeq;
            initAndValidate();
        }

        public void initAndValidate() {

        }

        public WildCard getInnerWC() {
            return null;
        }

        public Singleton getSingleton() {
            return this;
        }

        public TargetTract getTargetTract() {
            return new TargetTract(minTarg, maxTarg, minTargDeac, maxTargDeac);
        }

        public int getminTargDeac() {
            return minTargDeac;
        }

        public int getmaxTargDeac() {
            return maxTargDeac;
        }

        public boolean isLeftLong() {return minTargDeac  != minTarg;}
        public boolean isRightLong() {return maxTarg  != maxTargDeac;}
        public boolean isIntertarget() {return minTarg != maxTarg;}
        public int getStartPos() {return startPos;}
        public int getEndPos() {
            return startPos+ delLen;}
        public int getInsertLength() {return insSeq.length();}

        public int hashCode(){
            String insertEncoded = "";
            for(int i=0;i<insSeq.length();i++) {
                if (insSeq.charAt(i) == 'A') {
                    insertEncoded += "1";
                }
                if (insSeq.charAt(i) == 'T') {
                    insertEncoded += "2";
                }
                if (insSeq.charAt(i) == 'G') {
                    insertEncoded += "3";
                }
                if (insSeq.charAt(i) == 'C') {
                    insertEncoded += "4";
                }
            }
            String concatenated = "" + minTarg + maxTarg + minTargDeac + maxTargDeac + delLen + insertEncoded ;
            BigInteger number = new BigInteger(concatenated);
            int result = number.hashCode();
            return result;

        }


    }


    public static class WildCard extends IndelSet {

        /**
         *  Set of all indel tracts that only deactivate targets within the range mintrgt and maxtargt, inclusive
         */


        public WildCard() {
        }

        public WildCard(int minTrgt, int maxTrgt) {
            minTarg = minTrgt;
            maxTarg = maxTrgt;
            initAndValidate();
        }


        public TargetTract getTargetTract() {
            return new TargetTract(minTarg,maxTarg,this.getminTargDeac(),this.getmaxTargDeac());
        }

        public void initAndValidate() {
        } // initAndValidate

        public WildCard getInnerWC() {
            return this;
        }

        public Singleton getSingleton() {
            return null;

        }

        public int getminTargDeac() {
            return minTarg;
        }

        public int getmaxTargDeac() {
            return maxTarg;
        }

        public static WildCard intersect(WildCard wildCard1, WildCard wildCard2) {
            if ((wildCard1 == null) || (wildCard2 == null)) {
                return null;
            } else {
                int minTrgt = Math.max(wildCard1.minTarg, wildCard2.minTarg);
                int maxTrgt = Math.min(wildCard1.maxTarg, wildCard2.maxTarg);
                if (minTrgt > maxTrgt) {
                    return null;
                }
                return new WildCard(minTrgt, maxTrgt);
            }
        }

    }

}