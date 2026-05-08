package yimei.jss.simulation;

import yimei.jss.jobshop.Job;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.simulation.state.SystemState;

import java.util.ArrayList;
import java.util.List;

public class OrderingDecisionSituation extends DecisionSituation {

    private List<OperationOption> queue;

    private List<Job> jobQueue;
    private SystemState systemState;

    public OrderingDecisionSituation(List<OperationOption> operationOptions, SystemState systemState, String pair) {
        this.queue = operationOptions;
        this.systemState = systemState;
    }

    public List<OperationOption> getQueue() {
        return queue;
    }

/*    public OrderingDecisionSituation clone() {
        List<OperationOption> clonedQ = new ArrayList<>(queue);
        SystemState clonedState = systemState.clone();

        return new OrderingDecisionSituation(clonedQ, clonedState);
    }*/

    public OrderingDecisionSituation(List<Job> jobs, SystemState systemState) {
        this.jobQueue = jobs;
        this.systemState = systemState;
    }

    public List<Job> getJobQueue() {
        return jobQueue;
    }

    public SystemState getSystemState() {
        return systemState;
    }


    public OrderingDecisionSituation clone() {
        List<Job> clonedQ = new ArrayList<>(jobQueue);
        SystemState clonedState = systemState.clone();

        return new OrderingDecisionSituation(clonedQ, clonedState);
    }

}
