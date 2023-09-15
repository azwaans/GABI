package gestalt;

import gestalt.evolution.alignment.BarcodeMeta;
import gestalt.evolution.alignment.IndelSet;
import gestalt.evolution.alignment.TargetDeactTract;
import gestalt.evolution.alignment.TargetStatus;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class GestaltTargetStatusTest {



    @Test
    public void test_properties() {
        List<TargetDeactTract> deactTracts= Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1),new TargetDeactTract(3,3),new TargetDeactTract(5,7)});
        TargetStatus targStat = new TargetStatus(deactTracts);

        assertEquals((targStat.getDeactTargets()).toArray(),new Integer[]{0,1,3,5,6,7});


    }

    @Test
    public void test_merge() {

        List<TargetDeactTract> deactTracts1= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1),new TargetDeactTract(3,3),new TargetDeactTract(5,7)}));
        TargetStatus targStat1 = new TargetStatus(deactTracts1);

        List<TargetDeactTract> deactTracts2= Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(9,9)});
        TargetStatus targStat2 = new TargetStatus(deactTracts2);

        TargetStatus merge = targStat1.merge(targStat2);
        deactTracts1.addAll(deactTracts2);
        assertEquals(merge,new TargetStatus(deactTracts1));

    }

    @Test
    public void test_minus() {

        List<TargetDeactTract> deactTracts1= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)}));
        TargetStatus targStat1 = new TargetStatus(deactTracts1);
        Set<Integer> minus =  targStat1.minus(new TargetStatus());
        assertEquals(minus, new HashSet<>(Arrays.asList(new Integer[]{0, 1})));

        List<TargetDeactTract> deactTracts2= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1),new TargetDeactTract(3,3)}));
        TargetStatus targStat2 = new TargetStatus(deactTracts2);
        minus =  targStat2.minus(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)})));
        assertEquals(minus, new HashSet<>(Arrays.asList(new Integer[]{3})));

        List<TargetDeactTract> deactTracts3= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(3,6)}));
        TargetStatus targStat3 = new TargetStatus(deactTracts3);
        minus =  targStat3.minus(targStat2);
        assertEquals(minus, new HashSet<>());

    }

    @Test
    public void test_active_targets() {

        List<TargetDeactTract> deactTracts1= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1),new TargetDeactTract(4,4)}));
        TargetStatus targStat1 = new TargetStatus(deactTracts1);
        List<Integer> active = targStat1.getActiveTargets(10);
        assertEquals(active, Arrays.asList(new Integer[]{2,3,5,6,7,8,9}));

    }

    @Test
    public void test_possible_target_tracts() {
        ///WEIRD THAT GAPML TEST WORKS WITH A TOO SHORT BARCODE
        TargetStatus targetStatus = new TargetStatus();
        String uneditedBarcodeInput="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC ACTG TGACTAGC TT";
        String[] barcodeSplit = uneditedBarcodeInput.split(" ");
        int cutSite=3;
        int[] crucialPos= {3,3};
        int maxSumSteps= 3000;
        int maxExtraSteps=1;
        BarcodeMeta metaData = new BarcodeMeta(Arrays.asList(barcodeSplit),cutSite,crucialPos,maxSumSteps,maxExtraSteps);

        List<IndelSet.TargetTract> targetTracts = targetStatus.getPossibleTargetTracts(new ArrayList<Integer>(), metaData.nTargets);
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(0,2,0,2)));

        List<TargetDeactTract> deactTracts1= new ArrayList<TargetDeactTract>(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2),new TargetDeactTract(4,6),new TargetDeactTract(8,9)}));
        TargetStatus targStat1 = new TargetStatus(deactTracts1);
        targetTracts = targStat1.getPossibleTargetTracts(new ArrayList<Integer>(), metaData.nTargets);
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(3,3,2,4)));
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(3,3,3,3)));
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(3,7,3,8)));
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(3,7,2,7)));
        assertTrue(targetTracts.contains(new IndelSet.TargetTract(7,7,6,8)));



    }


    @Test
    public void test_get_contained_target_statuses() {

        List<TargetStatus> targetStatuses = (new TargetDeactTract(0,1)).getContainedStatuses();

        assertTrue(targetStatuses.contains(new TargetStatus()));
        assertTrue(targetStatuses.contains(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,0)}))));
        assertTrue(targetStatuses.contains(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(1,1)}))));
        assertTrue(targetStatuses.contains(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)}))));
        assertEquals( 4,targetStatuses.size());

    }

    @Test
    public void test_get_all_transitions() {


        String uneditedBarcodeInput="AA ATCGATCG ACTG ATCGATCG ACTG TGACTAGC TT";
        String[] barcodeSplit = uneditedBarcodeInput.split(" ");
        int cutSite=3;
        int[] crucialPos= {3,3};
        int maxSumSteps= 3000;
        int maxExtraSteps=1;
        BarcodeMeta metaData = new BarcodeMeta(Arrays.asList(barcodeSplit),cutSite,crucialPos,maxSumSteps,maxExtraSteps);

        Hashtable<TargetStatus, Hashtable<TargetStatus, List<IndelSet.TargetTract>>> allTransitions = TargetStatus.getAllTransitions(metaData.nTargets);
        assertEquals( 8,allTransitions.keySet().size());

        assertEquals(0,allTransitions.get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2)}))).size());

        Set<TargetStatus> set = new HashSet<TargetStatus>(Arrays.asList(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,0)})),
                new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)})),
                new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(1,1)})),
                new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(1,2)})),
                new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2)})),
                new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(2,2)}))));

        assertEquals(set,allTransitions.get(new TargetStatus()).keySet());

        Set<IndelSet.TargetTract> set2 = new HashSet<>(Arrays.asList(new IndelSet.TargetTract[]{
                new IndelSet.TargetTract(0, 1, 0, 1),
                new IndelSet.TargetTract(0, 0, 0, 1),
                new IndelSet.TargetTract(1, 1, 0, 1),
        }));



        assertTrue(set2.equals(new HashSet<>(allTransitions.get(new TargetStatus()).get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)}))))));

        Set<IndelSet.TargetTract> set3 = new HashSet<>(Arrays.asList(new IndelSet.TargetTract[]{
                new IndelSet.TargetTract(0, 2, 0, 2),
                new IndelSet.TargetTract(0, 1, 0, 2),
                new IndelSet.TargetTract(1, 1, 0, 2),
                new IndelSet.TargetTract(1, 2, 0, 2),
        }));



        assertTrue(set3.equals(new HashSet<>(allTransitions.get(new TargetStatus()).get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2)}))))));


        Set<IndelSet.TargetTract> set4 = new HashSet<>(Arrays.asList(new IndelSet.TargetTract[]{
                new IndelSet.TargetTract(2, 2, 2, 2),
                new IndelSet.TargetTract(2, 2, 1, 2),

        }));



        assertTrue(set4.equals(new HashSet<>(allTransitions.get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,1)}))).get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2)}))))));

        Set<IndelSet.TargetTract> set5 = new HashSet<>(Arrays.asList(new IndelSet.TargetTract[]{
                new IndelSet.TargetTract(0, 2, 0, 2),

        }));



        assertTrue(set5.equals(new HashSet<>(allTransitions.get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(1,1)}))).get(new TargetStatus(Arrays.asList(new TargetDeactTract[]{new TargetDeactTract(0,2)}))))));




    }




}
