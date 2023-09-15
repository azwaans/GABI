package gestalt;

import gestalt.evolution.alignment.AncStates;
import gestalt.evolution.alignment.IndelSet;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class GestaltIndelSetTest {


    @Test
    public void test_intersection() {

        IndelSet.WildCard w = new IndelSet.WildCard(2,4);
        IndelSet.SingletonWC s = new IndelSet.SingletonWC(1,1,10,1,1,1,"Asdf");
        assertNull(IndelSet.intersect(w,s));

         w = new IndelSet.WildCard(2,4);
         s = new IndelSet.SingletonWC(2,4,10,1,2,4,"Asdf");
        assertEquals(IndelSet.intersect(w,s),s);

        w = new IndelSet.WildCard(2,4);
        s = new IndelSet.SingletonWC(3,3,10,1,3,3,"Asdf");
        assertEquals(IndelSet.intersect(w,s),s);

        w = new IndelSet.WildCard(2,4);
        s = new IndelSet.SingletonWC(3,3,10,1,2,4,"Asdf");
        assertEquals(IndelSet.intersect(w,s),s);

        w = new IndelSet.WildCard(2,4);
        s = new IndelSet.SingletonWC(3,4,10,1,2,5,"Asdf");
        assertNull(IndelSet.intersect(w,s));

    }

    @Test
    public void test_intersect_ancstate() {

        List<IndelSet> l1 = new ArrayList<>();
        l1 = Arrays.asList(new IndelSet[]{new IndelSet.WildCard(1,1), new IndelSet.SingletonWC(2,2,30, 10,2,2, "asdf"),new IndelSet.SingletonWC(4, 4, 100,70,3,5,"")});
        List<IndelSet> l2 = new ArrayList<>();
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.WildCard(2,2), new IndelSet.WildCard(3,5)});
        AncStates ancState1 = new AncStates(l1);
        AncStates ancState2 = new AncStates(l2);
        AncStates intersect = AncStates.intersect(ancState1,ancState2);
        assertTrue(intersect.indelSetList.contains(l1.get(1)));
        assertTrue(intersect.indelSetList.contains(l1.get(2)));
        assertEquals(2,intersect.indelSetList.size());

        List<IndelSet> l3 = new ArrayList<>();
        l3 = Arrays.asList(new IndelSet[]{new IndelSet.WildCard(4,10)});
        AncStates ancState3 = new AncStates(l3);
        intersect = AncStates.intersect(ancState3,ancState2);
        assertEquals(intersect.indelSetList,Arrays.asList(new IndelSet[]{new IndelSet.WildCard(4,5)}));

        List<IndelSet> l4 = new ArrayList<>();
        l4 = Arrays.asList(new IndelSet[]{new IndelSet.WildCard(2,6)});
        AncStates ancState4 = new AncStates(l4);
        intersect = AncStates.intersect(ancState1,ancState4);
        assertTrue(intersect.indelSetList.contains(l1.get(1)));
        assertTrue(intersect.indelSetList.contains(l1.get(2)));
        assertEquals(2,intersect.indelSetList.size());

        List<IndelSet> l5 = new ArrayList<>();
        l5 = Arrays.asList(new IndelSet[]{new IndelSet.WildCard(8,10)});
        AncStates ancState5 = new AncStates(l5);
        intersect = AncStates.intersect(ancState5,ancState2);
        assertEquals(0,intersect.indelSetList.size());

        List<IndelSet> l6 = new ArrayList<>();
        l6 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(2,2,30, 10,2,2, "a")});
        AncStates ancState6 = new AncStates(l6);
        intersect = AncStates.intersect(ancState6,ancState1);
        assertEquals(0,intersect.indelSetList.size());

        List<IndelSet> l7 = l1.subList(1,1);
        AncStates ancState7 = new AncStates(l7);
        intersect = AncStates.intersect(ancState7,ancState1);
        assertEquals(ancState7.indelSetList,intersect.indelSetList);



    }

    @Test
    public void test_intersect_ancstate_offset_indel_sets() {

        List<IndelSet> l1 = new ArrayList<>();
        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,29, 17,0,1, "tggg"),new IndelSet.SingletonWC(2, 6, 116,73,2,6,"aggcga")});
        List<IndelSet> l2 = new ArrayList<>();
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(1,3,64, 42,1,3, "ac"),new IndelSet.SingletonWC(5, 5, 2,153,5,5,"a")});

        AncStates ancState1 = new AncStates(l1);
        AncStates ancState2 = new AncStates(l2);
        AncStates intersect = AncStates.intersect(ancState1,ancState2);
        assertTrue(intersect.indelSetList.contains(l2.get(1)));
        assertEquals(1,intersect.indelSetList.size());

        //////////////////

        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,49, 17,0,1, "tggg")});
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(1,3,64, 42,1,3, "ac")});

        ancState1 = new AncStates(l1);
        ancState2 = new AncStates(l2);
        intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(0,intersect.indelSetList.size());

        //////////////////


        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,3,60, 17,0,4, "")});
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(3,3,34, 100,2,3, "ac")});

        ancState1 = new AncStates(l1);
        ancState2 = new AncStates(l2);
        intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(0,intersect.indelSetList.size());


    }

    @Test
    public void test_intersect_ancstate_touching_sgwcs() {

        List<IndelSet> l1 = new ArrayList<>();
        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,29, 17,0,1, "tggg"),new IndelSet.SingletonWC(2, 2, 35,55,1,2,"a")});
        List<IndelSet> l2 = new ArrayList<>();
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,0,4, 20,0,0, "ac"),new IndelSet.SingletonWC(2, 2, 35,55,1,2,"a")});

        AncStates ancState1 = new AncStates(l1);
        AncStates ancState2 = new AncStates(l2);
        AncStates intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(0,intersect.indelSetList.size());

        //////////////////

        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,2,60, 17,0,2, "tggg"),new IndelSet.SingletonWC(3,3,35, 85,2,3, "a")});
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,5,100, 20,0,5, "ac")});

        ancState1 = new AncStates(l1);
        ancState2 = new AncStates(l2);
        intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(1,intersect.indelSetList.size());
        assertTrue(intersect.indelSetList.contains(new IndelSet.WildCard(1,1)));

        //////////////////

        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,45, 17,0,2, "tggg"),new IndelSet.SingletonWC(2,2,3, 70,2,2, "")});
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,45, 17,0,2, "tggg")});

        ancState1 = new AncStates(l1);
        ancState2 = new AncStates(l2);
        intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(0,intersect.indelSetList.size());

        ///////////////////

        l1 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,45, 17,0,2, "tggg"),new IndelSet.SingletonWC(2,2,3, 70,2,2, "")});
        l2 = Arrays.asList(new IndelSet[]{new IndelSet.SingletonWC(0,1,43, 19,0,2, ""),new IndelSet.SingletonWC(2,2,3, 70,2,2, "")});

        ancState1 = new AncStates(l1);
        ancState2 = new AncStates(l2);
        intersect = AncStates.intersect(ancState1,ancState2);
        assertEquals(1,intersect.indelSetList.size());
        assertTrue(intersect.indelSetList.contains(l2.get(1)));


    }

}
