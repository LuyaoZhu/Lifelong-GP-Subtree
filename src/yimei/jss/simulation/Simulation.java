package yimei.jss.simulation;

import ec.gp.GPNode;
import yimei.jss.jobshop.*;
import yimei.jss.jobshop.Process;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.simulation.event.AbstractEvent;
import yimei.jss.simulation.event.ProcessStartEvent;
import yimei.jss.simulation.state.SystemState;

import java.util.*;

/**
 * The abstract simulation class for evaluating rules.
 * <p>
 * Created by yimei on 21/11/16.
 */
public abstract class Simulation {
    @Override
    public String toString() {
        return "Simulation{" +
                "sequencingRule=" + sequencingRule +
                ", routingRule=" + routingRule +
                ", systemState=" + systemState +
                ", eventQueue=" + eventQueue +
                ", numWorkCenters=" + numWorkCenters +
                ", numJobsRecorded=" + numJobsRecorded +
                ", warmupJobs=" + warmupJobs +
                ", numJobsArrived=" + numJobsArrived +
                ", throughput=" + throughput +
                '}';
    }

    protected AbstractRule sequencingRule;
    protected AbstractRule routingRule;
    protected AbstractRule orderingRule;
    protected SystemState systemState;
    protected PriorityQueue<AbstractEvent> eventQueue;

    protected int numWorkCenters;
    protected int numJobsRecorded;

    public int numBatchesRecorded;
    protected int warmupJobs;

    protected int warmupBatches;
    protected int numJobsArrived;

    protected int numBatchesArrived;
    public int throughput;
    //protected int[] jobStates;

    public AbstractRule defaultSequencingRule;
    public AbstractRule defaultRoutingRule;

    public static String optimisationObjective;

    //fzhang 3.6.2018  discard the individual(rule) can not complete the whole jobs well, take a long time (prefer to do part of each job)
    int beforeThroughput; //save the throughput value before updated (a job finished)
    int afterThroughput; //save the throughput value after updated (a job finished)
    int count = 0;


/*    public Simulation(AbstractRule sequencingRule,
                      AbstractRule routingRule,
                      int numWorkCenters,
                      int numJobsRecorded,
                      int warmupJobs) {
        this.sequencingRule = sequencingRule;
        this.routingRule = routingRule;
        this.numWorkCenters = numWorkCenters;
        this.numJobsRecorded = numJobsRecorded;
        this.warmupJobs = warmupJobs;

        systemState = new SystemState();
        eventQueue = new PriorityQueue<>();
//        int[] jobStates = new int[numJobsRecorded];
//        fill(jobStates, -1);
//        this.jobStates = jobStates;
    }*/

    public Simulation(AbstractRule sequencingRule,
                      AbstractRule routingRule,
                      int numWorkCenters,
                      int numBatchesRecorded,
                      int warmupBatches) {
        this.sequencingRule = sequencingRule;
        this.routingRule = routingRule;
        this.numWorkCenters = numWorkCenters;
        this.numBatchesRecorded = numBatchesRecorded;
        this.warmupBatches = warmupBatches;

        systemState = new SystemState();
        eventQueue = new PriorityQueue<>();
//        int[] jobStates = new int[numJobsRecorded];
//        fill(jobStates, -1);
//        this.jobStates = jobStates;
    }

    public AbstractRule getSequencingRule() {
        return sequencingRule;
    }

//    public int[] getJobStates() { return jobStates; }

    public AbstractRule getRoutingRule() {
        return routingRule;
    }

    public AbstractRule getDefaultSequencingRule() {
        return defaultSequencingRule;
    }

    public AbstractRule getDefaultRoutingRule() {
        return  defaultRoutingRule;
    }

    public AbstractRule getOrderingRule() {
        return orderingRule;
    }

    public SystemState getSystemState() {
        return systemState;
    }

    public PriorityQueue<AbstractEvent> getEventQueue() {
        return eventQueue;
    }

    public void setSequencingRule(AbstractRule sequencingRule) {
        this.sequencingRule = sequencingRule;
    }

    public Integer getNumJobsArrived() {
        return numJobsArrived;
    }

    public Integer getNumBatchesArrived() {
        return numBatchesArrived;
    }

    public Integer getWarmupBatches(){ return warmupBatches; }

    public Integer getBatchesRecorded(){ return numBatchesRecorded; }

//    public void setJobStates(int[] jobStates) { this.jobStates = jobStates; }

    public void setRoutingRule(AbstractRule routingRule) {
        this.routingRule = routingRule;
        //need to reset state as well, as the operationoptions associated
        //with workcenters are chosen using this routing rule, so current
        //values are outdated
        resetState();
    }

    public void setOrderingRule(AbstractRule orderingRule) {
        this.orderingRule = orderingRule;
        //need to reset state as well, as the operationoptions associated
        //with workcenters are chosen using this routing rule, so current
        //values are outdated
//        resetState();
    }

    public void setBothRuleWithoutReset(AbstractRule sequencingRule, AbstractRule routingRule) {
        this.sequencingRule = sequencingRule;
        this.routingRule = routingRule;
    }

    public double getClockTime() {
        return systemState.getClockTime();
    }

    public void addEvent(AbstractEvent event) {
        eventQueue.add(event);
    }

    public boolean canAddToQueue(Process process) {
        Iterator<AbstractEvent> e = eventQueue.iterator();
        if (e.hasNext()) {
            AbstractEvent a = e.next();
            if (a instanceof ProcessStartEvent) {
                if (((ProcessStartEvent) a).getProcess().getWorkCenter().getId() ==
                        process.getWorkCenter().getId()) {
                    return false;
                }
            }
        }
        return true;
    }

    //int countBadrun =0;
    public void run() {

        while (!eventQueue.isEmpty() && throughput < numBatchesRecorded) { //numJobsRecorded == 5000

            AbstractEvent nextEvent = eventQueue.poll(); // the head of this queue, or null if this queue is empty

//            if(eventQueue.isEmpty())
//                System.out.println(1);

//            System.out.println("EventQueue's size: " + eventQueue.size());
            //fzhang 3.6.2018  fix the stuck problem
            beforeThroughput = throughput; //save the throughput value before updated (a job finished)

            systemState.setClockTime(nextEvent.getTime());
            nextEvent.trigger(this); //nextEvent includes many different types of events

            afterThroughput = throughput; //save the throughput value after updated (a job finished)

            if (throughput > 0 & afterThroughput - beforeThroughput == 0) { //if the value was not updated
                count++;
            }

            //otherwise the fitness of test rules may
//
            //System.out.println("count "+count);
            if (count > 150000) {
                count = 0;
                systemState.setClockTime(Double.MAX_VALUE);
                eventQueue.clear();
            }

            //===================ignore busy machine here==============================
            //when nextEvent was done, check the numOpsInQueue
            for (WorkCenter w : systemState.getWorkCenters()) {
                if (w.numOpsInQueue() > 150) {
                    systemState.setClockTime(Double.MAX_VALUE);
                    eventQueue.clear();
                    //countBadrun++;
                }
            }
//            }
        }
        //modified by fzhang 18.04.2018
        /*if(countBadrun>0) {
        	 System.out.println("The number of badrun grasped in simulation: "+ countBadrun);
         }*/

        if (!systemState.getJobsInSystem().isEmpty() && !(this instanceof DynamicSimulation)) {
            System.out.println("Event queue is empty but simulation is not complete.");
            System.out.println("Makespan is garbage - cannot continue.");
            System.exit(0);
        }

    }


//    private boolean eventIsDuplicate(AbstractEvent event) {
//        if (event instanceof ProcessFinishEvent) {
//            Process p = ((ProcessFinishEvent) event).getProcess();
//            //want to check whether this operation has already been performed
//            int jobId = p.getOperationOption().getJob().getId();
//            if (jobId >= 0) {
//                int jobState = jobStates[jobId];
//                int opNum = p.getOperationOption().getOperation().getId();
//                if ((jobState+1) != opNum) {
//                    //upcoming event should only be the next job in the sequence,
//                    //not a job we've already done, or one ahead of the next one
//                    return true;
//                }
//            }
//        }
//        return false;
//    }

    public void rerun() {
        //original
        //fzhang 2018.11.5 this is used for generate different instances in a generation.
        //if the replications is 1, does not have influence
        resetState();

        //reset(): reset seed value, will get the same instance
        //reset();
        run();
    }

    public void completeJob(Job job) {
//        if (numBatchesArrived > warmupBatches && job.batch.getId() >= 0
//                && job.batch.getId() < numBatchesRecorded + warmupBatches) {
        if (job.batch.getId() >=  warmupBatches              // Collect data from the fixed ID (51-300)
                && job.batch.getId() < numBatchesRecorded + warmupBatches) {

            if (job.batch.completedJobs == job.batch.numOfJobs) {

                job.getBatch().setCompletionTime(job.getCompletionTime());

                throughput++;  //before only have this line

                count = 0;

                systemState.addCompletedBatch(job.getBatch());
            }

            systemState.addCompletedJob(job);
//            int a = systemState.getJobsCompleted().size();
//            System.out.println("The number of completed jobs: "+systemState.getJobsCompleted().size());
        }
        systemState.removeJobFromSystem(job);
        if (job.batch.completedJobs == job.batch.numOfJobs) {
            systemState.removeBatchFromSystem(job.batch);
        }
    }

    public double makespan() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            double tmp = job.getCompletionTime();
            if (value < tmp)
                value = tmp;
        }

        return value;
    }

    public double meanFlowtime() {
        double value = 0.0;

        for (Job job : systemState.getJobsCompleted()) {
            value += job.flowTime();
        }

        return value / systemState.getJobsCompleted().size();

    }

    public double batchMeanFlowtime() {
        double value = 0.0;

        for (Batch batch : systemState.getBatchesCompleted()) {
            value += batch.flowTime();
        }

        return value / numBatchesRecorded;

    }

    public double maxFlowtime() {
        double value = 0.0;

        for (Job job : systemState.getJobsCompleted()) {
            double tmp = job.flowTime();
            if (value < tmp)
                value = tmp;
        }

        return value;

    }

    public double batchMaxFlowtime() {
        double value = 0.0;

        for (Batch batch : systemState.getBatchesCompleted()) {
            double tmp = batch.flowTime();
            if (value < tmp)
                value = tmp;
        }

        return value;

    }

    public double meanWeightedFlowtime() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            value += job.weightedFlowTime();
        }

        return value / systemState.getJobsCompleted().size();
    }

    public double batchMeanWeightedFlowtime() {
        double value = 0.0;
        for (Batch batch : systemState.getBatchesCompleted()) {
            value += batch.weightedFlowTime();
        }

        return value / numBatchesRecorded;
    }

    public double maxWeightedFlowtime() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            double tmp = job.weightedFlowTime();
            if (value < tmp)
                value = tmp;
        }

        return value;
    }

    public double meanTardiness() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
//            System.out.println(batch.tardiness());
            value += job.tardiness();
        }

        return value / systemState.getJobsCompleted().size();

    }

    public double batchMeanTardiness() {
        double value = 0.0;
        for (Batch batch : systemState.getBatchesCompleted()) {
//            System.out.println(batch.tardiness());
            value += batch.tardiness();
        }

        return value / numBatchesRecorded;

    }

    public double maxTardiness() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            double tmp = job.tardiness();

            if (value < tmp)
                value = tmp;
        }

        return value;
    }

    public double meanWeightedTardiness() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            value += job.weightedTardiness();
        }

        return value / systemState.getJobsCompleted().size();
    }

    public double batchMeanWeightedTardiness() {
        double value = 0.0;
        for (Batch batch : systemState.getBatchesCompleted()) {
            value += batch.weightedTardiness();
        }

        return value / numBatchesRecorded;
    }

    public double batchMaxWeightedTardiness() {
        double value = 0.0;
        for (Batch batch : systemState.getBatchesCompleted()) {
            double tmp = batch.weightedTardiness();

            if (value < tmp)
                value = tmp;

        }

        return value;
    }

    public double maxWeightedTardiness() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            double tmp = job.weightedTardiness();

            if (value < tmp)
                value = tmp;
        }
        return value;
    }

    public double jobProfit() {
        double value = 0.0;            //This is the profit = revenue - penalty
        double revenue = 0.0;
        double penalty = 0.0;

        //calculate the revenue, completed jobs
        for (Job job : systemState.getJobsCompleted()) {
            double profit = job.getRevenue() * (1 -  job.getWeight()/2 * (job.tardiness() / (job.getTotalProcTime()*2)));  //        value = revenue - penalty;
            value += profit;
        }
        return value;

    }

    public double batchProfit() {
        double value = 0.0;            //This is the profit = revenue - penalty
        double revenue = 0.0;
        double penalty = 0.0;

        //calculate the revenue, completed jobs
        for (Batch batch : systemState.getBatchesCompleted()) {
            double profit = batch.revenue * (1 -  batch.getWeight() * (batch.tardiness() / (4*batch.estimatedProcTime)));
//            double profit = batch.revenue * (1 - batch.getWeight() *  (batch.tardiness() / (8 * batch.estimatedProcTime)));
            value += profit;
        }
        return value;

    }

    public double propTardyJobs() {
        double value = 0.0;
        for (Job job : systemState.getJobsCompleted()) {
            if (job.getCompletionTime() > job.getDueDate())
                value++;
        }

        return value / numJobsRecorded;
    }

    //2018.12.20 define rule size as an objective
    public int rulesize() {
        int value = 0;
        GPRule seqRule = null;
        GPRule routRule = null;

        seqRule = (GPRule) this.getSequencingRule();
        routRule = (GPRule) this.getRoutingRule();
        int seqRuleSize = seqRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);
        int routRuleSize = routRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);

        value = seqRuleSize + routRuleSize;
    	/*System.out.println("==========================");
    	System.out.println("RuleSize "+value);*/
        return value;
    }

    //2019.2.26 define routing rule size as an objective
    public int rulesizeR() {
        int value = 0;
        GPRule routRule = null;
        routRule = (GPRule) this.getRoutingRule();
        int routRuleSize = routRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);
//    	System.out.println("routRuleSize "+routRuleSize);
        value = routRuleSize;
        return value;
    }

    public int rulesizeS() {
        int value = 0;
        GPRule seqRule = null;
        seqRule = (GPRule) this.getSequencingRule();
        int seqRuleSize = seqRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);
//    	System.out.println("seqRuleSize "+seqRuleSize);
        value = seqRuleSize;
        return value;
    }

    public double objectiveValue(Objective objective) {
        switch (objective) {
            case MAKESPAN:
                return makespan();
            case MEAN_FLOWTIME:
                return meanFlowtime();
            case MAX_FLOWTIME:
                return maxFlowtime();
            case MEAN_WEIGHTED_FLOWTIME:
                return meanWeightedFlowtime();
            case MAX_WEIGHTED_FLOWTIME:
                return maxWeightedFlowtime();
            case MEAN_TARDINESS:
                return meanTardiness();
            case MAX_TARDINESS:
                return maxTardiness();
            case MEAN_WEIGHTED_TARDINESS:
                return meanWeightedTardiness();
            case MAX_WEIGHTED_TARDINESS:
                return maxWeightedTardiness();
            case PROP_TARDY_JOBS:
                return propTardyJobs();
            case RULESIZE:
                return rulesize();
            case RULESIZER:
                return rulesizeR();
            case RULESIZES:
                return rulesizeS();
            case JOB_PROFIT:
                return jobProfit();

            case BATCH_MEAN_FLOWTIME:
                return batchMeanFlowtime();
            case BATCH_MAX_FLOWTIME:
                return batchMaxFlowtime();
            case BATCH_MEAN_TARDINESS:
                return batchMeanTardiness();
            case BATCH_MEAN_WEIGHTED_TARDINESS:
                return batchMeanWeightedTardiness();
            case BATCH_MAX_WEIGHTED_TARDINESS:
                return batchMaxWeightedTardiness();
            case BATCH_PROFIT:
                return batchProfit();
        }

        return -1.0;
    }

    //modified by mengxu for one instance multi case evaluation.
    public double[] objectiveValueMultiCase(Objective objective, int numCase) {
        switch (objective) {

            case MEAN_FLOWTIME:
                return multiCaseMeanFlowtime(numCase);
            case MAX_FLOWTIME:
                return multiCaseMaxFlowtime(numCase);
            /*case MEAN_WEIGHTED_FLOWTIME:
                return multiCaseMeanWeightedFlowtime(numCase);
            case MAX_WEIGHTED_FLOWTIME:
                return multiCaseMaxWeightedFlowtime(numCase);
            case MEAN_TARDINESS:
                return multiCaseMeanTardiness(numCase);
            case MAX_TARDINESS:
                return multiCaseMaxTardiness(numCase);*/
            case MEAN_WEIGHTED_TARDINESS:
                return multiCaseMeanWeightedTardiness(numCase);
            case MAX_WEIGHTED_TARDINESS:
                return multiCaseMaxWeightedTardiness(numCase);
//            case RULESIZE:
//                return rulesize();
//            case RULESIZER:
//                return rulesizeR();
//            case RULESIZES:
//                return rulesizeS();
            case BATCH_MEAN_FLOWTIME:
                return multiCaseBatchMeanFlowtime(numCase);
            case BATCH_MAX_FLOWTIME:
                return multiCaseBatchMaxFlowtime(numCase);
            case BATCH_MEAN_WEIGHTED_TARDINESS:
                return multiCaseBatchMeanWeightedTardiness(numCase);
            case BATCH_MAX_WEIGHTED_TARDINESS:
                return multiCaseBatchMaxWeightedTardiness(numCase);
        }

        double[] numCaseFitness = new double[numCase];
        Arrays.fill(numCaseFitness,-1);

        return numCaseFitness;
    }

    public double[] multiCaseMeanFlowtime(int numCase) {
        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double sumFlowTime = 0.0;
            int jobsInThisCase = 0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    sumFlowTime = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    sumFlowTime = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    for (Job job : batch.getJobs()) {
                        sumFlowTime += job.flowTime();
                        jobsInThisCase++;
                    }
                }
            }

            numCaseFitness[caseIdx] = sumFlowTime/jobsInThisCase;
        }

        return numCaseFitness;
    }

    public double[] multiCaseBatchMeanFlowtime(int numCase) {
        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double sumFlowTime = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    sumFlowTime = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    sumFlowTime = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    sumFlowTime += batch.flowTime();
                }
            }

            numCaseFitness[caseIdx] = sumFlowTime/totalBatches;
        }

        return numCaseFitness;
    }

    public double[] multiCaseMaxFlowtime(int numCase) {
        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double maxFlowTime = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    maxFlowTime = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    maxFlowTime = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    for (Job job : batch.getJobs()) {
                        maxFlowTime = Math.max(maxFlowTime, job.flowTime());
                    }
                }
            }

            numCaseFitness[caseIdx] = maxFlowTime;
        }

        return numCaseFitness;
    }

    public double[] multiCaseBatchMaxFlowtime(int numCase) {
        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double maxFlowTime = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    maxFlowTime = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    maxFlowTime = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    maxFlowTime = Math.max(maxFlowTime, batch.flowTime());
                }
            }

            numCaseFitness[caseIdx] = maxFlowTime;
        }

        return numCaseFitness;
    }

    public double[] multiCaseMeanWeightedFlowtime(int numCase) {
        int totalJobs = systemState.getJobsCompleted().size();
        double[] numCaseFitness = new double[numCase];

        // 出错保护：如果 throughput 不正常，全部返回 MAX
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.MAX_VALUE);
            return numCaseFitness;
        }

        // 按照 Job ID 顺序排序
        Collections.sort(systemState.getJobsCompleted(), Comparator.comparingInt(Job::getId));

        // 均匀分配 job 数
        int baseJobsPerCase = totalJobs / numCase;
        int remainder = totalJobs % numCase;

        int jobIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int jobsInThisCase = baseJobsPerCase + (caseIdx < remainder ? 1 : 0);
            double sumWeightedFlowTime = 0.0;

            for (int j = 0; j < jobsInThisCase; j++) {
                if (jobIndex >= totalJobs) {
                    sumWeightedFlowTime = Double.POSITIVE_INFINITY;
                    break;
                }

                Job job = systemState.getJobsCompleted().get(jobIndex++);
                if (job == null) {
                    sumWeightedFlowTime = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    sumWeightedFlowTime += job.weightedFlowTime();
                }
            }

            numCaseFitness[caseIdx] = (sumWeightedFlowTime == Double.POSITIVE_INFINITY)
                    ? Double.MAX_VALUE
                    : sumWeightedFlowTime / jobsInThisCase;
        }

        return numCaseFitness;
    }


    public double[] multiCaseMaxWeightedFlowtime(int numCase) {
        int totalJobs = systemState.getJobsCompleted().size();
        double[] numCaseFitness = new double[numCase];

        // 异常情况：throughput不一致，返回无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 按照 Job ID 排序
        Collections.sort(systemState.getJobsCompleted(), Comparator.comparingInt(Job::getId));

        // 均匀分配 job 数
        int baseJobsPerCase = totalJobs / numCase;
        int remainder = totalJobs % numCase;

        int jobIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int jobsInThisCase = baseJobsPerCase + (caseIdx < remainder ? 1 : 0);
            double maxWeightedFlow = Double.NEGATIVE_INFINITY;

            for (int j = 0; j < jobsInThisCase; j++) {
                if (jobIndex >= totalJobs) {
                    maxWeightedFlow = Double.POSITIVE_INFINITY;
                    break;
                }

                Job job = systemState.getJobsCompleted().get(jobIndex++);
                if (job == null) {
                    maxWeightedFlow = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    maxWeightedFlow = Math.max(maxWeightedFlow, job.weightedFlowTime());
                }
            }

            numCaseFitness[caseIdx] = (maxWeightedFlow == Double.POSITIVE_INFINITY)
                    ? Double.POSITIVE_INFINITY
                    : maxWeightedFlow;
        }

        return numCaseFitness;
    }


    public double[] multiCaseMeanTardiness(int numCase) {
        int totalJobs = systemState.getJobsCompleted().size();
        double[] numCaseFitness = new double[numCase];

        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 按 Job ID 升序排序（确保顺序一致）
        Collections.sort(systemState.getJobsCompleted(), Comparator.comparingInt(Job::getId));

        // 计算每个 case 的 job 数量
        int baseJobsPerCase = totalJobs / numCase;
        int remainder = totalJobs % numCase;

        int jobIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int jobsInThisCase = baseJobsPerCase + (caseIdx < remainder ? 1 : 0);
            double totalTardiness = 0.0;

            for (int j = 0; j < jobsInThisCase; j++) {
                if (jobIndex >= totalJobs) {
                    totalTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Job job = systemState.getJobsCompleted().get(jobIndex++);
                if (job == null) {
                    totalTardiness = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    totalTardiness += job.tardiness();
                }
            }

            if (totalTardiness == Double.POSITIVE_INFINITY || jobsInThisCase == 0) {
                numCaseFitness[caseIdx] = Double.POSITIVE_INFINITY;
            } else {
                numCaseFitness[caseIdx] = totalTardiness / jobsInThisCase;
            }
        }

        return numCaseFitness;
    }


    public double[] multiCaseMaxTardiness(int numCase) {
        int totalJobs = systemState.getJobsCompleted().size();
        double[] numCaseFitness = new double[numCase];

        if (throughput != numJobsRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 计算每个 case 应该包含的 job 数（尽可能均匀）
        int baseJobsPerCase = totalJobs / numCase;
        int remainder = totalJobs % numCase;

        int jobIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int jobsInThisCase = baseJobsPerCase + (caseIdx < remainder ? 1 : 0);
            double maxTardiness = 0.0;

            for (int j = 0; j < jobsInThisCase; j++) {
                if (jobIndex >= totalJobs) {
                    maxTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Job job = systemState.getJobsCompleted().get(jobIndex++);
                if (job == null) {
                    maxTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                double t = job.tardiness();
                if (t > maxTardiness) {
                    maxTardiness = t;
                }
            }

            numCaseFitness[caseIdx] = maxTardiness;
        }

        return numCaseFitness;
    }


    public double[] multiCaseMeanWeightedTardiness(int numCase) {

        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double sumWeightedTardiness = 0.0;
            int jobsInThisCase = 0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    sumWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    sumWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    for (Job job : batch.getJobs()) {
                        sumWeightedTardiness += job.weightedTardiness();
                        jobsInThisCase++;
                    }
                }
            }

            numCaseFitness[caseIdx] = sumWeightedTardiness/jobsInThisCase;
        }

        return numCaseFitness;
    }

    public double[] multiCaseBatchMeanWeightedTardiness(int numCase) {

        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double sumWeightedTardiness = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    sumWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    sumWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    sumWeightedTardiness += batch.weightedTardiness();
                }
            }

            numCaseFitness[caseIdx] = sumWeightedTardiness/totalBatches;
        }

        return numCaseFitness;
    }

    public double[] multiCaseMaxWeightedTardiness(int numCase) {


        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double maxWeightedTardiness = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    maxWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    maxWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                } else {
                    for (Job job : batch.getJobs()) {
                        maxWeightedTardiness = Math.max(maxWeightedTardiness, job.weightedTardiness());
                    }
                }
            }

            numCaseFitness[caseIdx] = maxWeightedTardiness;
        }

        return numCaseFitness;
    }

    public double[] multiCaseBatchMaxWeightedTardiness(int numCase) {


        int totalBatches = systemState.getBatchesCompleted().size();

        double[] numCaseFitness = new double[numCase];

        // 如果加工数不等于记录数，直接返回全为无穷大
        if (throughput != numBatchesRecorded) {
            Arrays.fill(numCaseFitness, Double.POSITIVE_INFINITY);
            return numCaseFitness;
        }

        // 排序以确保 job 顺序一致
        Collections.sort(systemState.getBatchesCompleted(), Comparator.comparingInt(Batch::getId));

        // 计算每个 case 应该包含的 job 数
        int baseBatchesPerCase = totalBatches / numCase;
        int remainder = totalBatches % numCase; // 有多少 case 需要多分一个 job

        int batchIndex = 0;
        for (int caseIdx = 0; caseIdx < numCase; caseIdx++) {
            int batchesInThisCase = baseBatchesPerCase;

            double maxWeightedTardiness = 0.0;

            for (int j = 0; j < batchesInThisCase; j++) {
                if (batchIndex >= systemState.getBatchesCompleted().size()) {
                    maxWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                }

                Batch batch = systemState.getBatchesCompleted().get(batchIndex++);
                if (batch == null) {
                    maxWeightedTardiness = Double.POSITIVE_INFINITY;
                    break;
                } else {
                        maxWeightedTardiness = Math.max(maxWeightedTardiness, batch.weightedTardiness());
                }
            }

            numCaseFitness[caseIdx] = maxWeightedTardiness;
        }

        return numCaseFitness;
    }

    public double workCenterUtilLevel(int idx) {
        return systemState.getWorkCenter(idx).getBusyTime() / getClockTime();
    }

    public String workCenterUtilLevelsToString() {
        String string = "[";
        for (int i = 0; i < systemState.getWorkCenters().size(); i++) {
            string += String.format("%.3f ", workCenterUtilLevel(i));
        }
        string += "]";

        return string;
    }

    public abstract void setup();

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Simulation that = (Simulation) o;

        if (numWorkCenters != that.numWorkCenters) return false;
        if (numJobsRecorded != that.numJobsRecorded) return false;
        if (warmupJobs != that.warmupJobs) return false;
        if (numJobsArrived != that.numJobsArrived) return false;
        if (throughput != that.throughput) return false;
        if (sequencingRule != null ? !sequencingRule.equals(that.sequencingRule) : that.sequencingRule != null)
            return false;
        if (routingRule != null ? !routingRule.equals(that.routingRule) : that.routingRule != null) return false;
        if (systemState != null ? !systemState.equals(that.systemState) : that.systemState != null) return false;
        return eventQueue != null ? eventQueue.equals(that.eventQueue) : that.eventQueue == null;
    }

    @Override
    public int hashCode() {
        int result = sequencingRule != null ? sequencingRule.hashCode() : 0;
        result = 31 * result + (routingRule != null ? routingRule.hashCode() : 0);
        result = 31 * result + (systemState != null ? systemState.hashCode() : 0);
        result = 31 * result + (eventQueue != null ? eventQueue.hashCode() : 0);
        result = 31 * result + numWorkCenters;
        result = 31 * result + numJobsRecorded;
        result = 31 * result + warmupJobs;
        result = 31 * result + numJobsArrived;
        result = 31 * result + throughput;
        return result;
    }

    public abstract void resetState();

    public abstract void reset();

    public abstract void rotateSeed();

    public abstract void generateJob();

    public abstract Simulation surrogate(int numWorkCenters, int numJobsRecorded,
                                         int warmupJobs);

    public abstract Simulation surrogateBusy(int numWorkCenters, int numJobsRecorded,
                                             int warmupJobs);

    public abstract void generateBatch();
}