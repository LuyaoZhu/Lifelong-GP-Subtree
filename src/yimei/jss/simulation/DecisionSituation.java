package yimei.jss.simulation;

import yimei.jss.jobshop.Job;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.simulation.state.SystemState;

import java.util.List;

/**
 * Created by dyska on 22/09/17.
 */
public abstract class DecisionSituation {

    private List<OperationOption> queue = null;

    private List<Job> jobQueue = null;
    private SystemState systemState = null;

    public List<OperationOption> getQueue() {
        return queue;
    }

    public List<Job> getJobQueue() {
        return jobQueue;
    }

    public SystemState getSystemState() {
        return systemState;
    }

    //fzhang 2019.9.4
    //public abstract WorkCenter getWorkCenter();
}
