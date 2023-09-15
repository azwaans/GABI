package gestalt.evolution.alignment;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Mahsa Sadi
 * @since 2020 - 03 - 25
 * <p>
 * License: Creative Commons
 * <p>
 * Copyright by Mahsa Sadi
 */


public class SetOperations {

    /**
     * Problem:
     * <p>
     * Calculate the Cartesian Products of the Given Sets.
     * <p>
     * <p>
     * <p>
     * Description:
     * <p>
     * Given a variable number of sets,
     * produce the Cartesian product of the given sets.
     * <p>
     * The sets are input as a list of lists.
     * <p>
     * <p>
     * <p>
     * <p>
     * Solution:
     * <p>
     * <p>
     * 1- Pick the first set from the list if the list is not empty
     * 2- For each element in the first list, add the element to the result.
     * 3- Remove the first set from the list.
     * 4- Repeat (Go to 1-).
     * 5- If the list is empty, add result to the list of answers.
     * 6- Remove the added element from the result in order to add the next element.
     * 7- Repeat (Go to 2-)
     */

    List<List<TargetStatus>> products;

    public SetOperations() {
        products = new ArrayList<List<TargetStatus>>();
    }


    public List<List<TargetStatus>> cartesianProduct(List<List<TargetStatus>> inputSets) {
        cartesianProduct(inputSets, new ArrayList<TargetStatus>());
        return products;
    }


    public void cartesianProduct(List<List<TargetStatus>> inputSets, List<TargetStatus> result) {
        if (inputSets.isEmpty()) {

            List<TargetStatus> copy = new ArrayList<TargetStatus>();

            for (TargetStatus e : result)
                copy.add(e);

            products.add(copy);

        } else {
            // 1- Pick the first set
            for (TargetStatus element : inputSets.get(0)) {
                // 2- Add the elements of the set to the result one by one
                result.add(element);

                // 3- Continue multiplication
                cartesianProduct(inputSets.subList(1, inputSets.size()), result);

                //4- Remove the element to be replaced by another one in the same set.
                result.remove(element);

            }
        }
    }


    public void printProducts() {
        for (int j = 0; j < products.size(); j++) {
            String s = "";

            for (int i = 0; i < products.get(j).size(); i++)
                s = s + products.get(j).get(i) + "  ";

            System.out.println(s);

        }

    }

}
