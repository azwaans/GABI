package beast.evolution.alignment;

import java.util.List;

public class QueueObject implements Comparable<QueueObject> {

    public int distance;
    public StateAlongPath state;

    public QueueObject(int dist, StateAlongPath State) {

        distance = dist;
        state = State;

    }

    public QueueObject(int dist, List<IndelSet.TargetTract> tttuple, boolean ispassed) {

        distance = dist;
        state = new StateAlongPath(tttuple,ispassed);

    }


    @Override
    public int compareTo(QueueObject o) {
        return this.distance - o.distance;
    }
}
