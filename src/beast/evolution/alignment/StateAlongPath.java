package beast.evolution.alignment;

import java.util.List;

public class StateAlongPath {

    public List<IndelSet.TargetTract> targetTactList;
    public boolean statePassedRequested;

    public StateAlongPath(List<IndelSet.TargetTract> ttTuple, boolean isPassed) {

        targetTactList = ttTuple;
        statePassedRequested = isPassed;

    }


    @Override
    public int hashCode() {
        return targetTactList.hashCode()*10 + ((this.statePassedRequested) ? 1 : 0);
    }


    public boolean equals(Object o)
    {
        if (o instanceof StateAlongPath) {
            return (this.targetTactList == ((StateAlongPath) o).targetTactList) && (this.statePassedRequested == ((StateAlongPath) o).statePassedRequested);
        }
        return false;
    }

}