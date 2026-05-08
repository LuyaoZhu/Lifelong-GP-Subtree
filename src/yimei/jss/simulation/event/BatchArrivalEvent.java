package yimei.jss.simulation.event;

import yimei.jss.gp.CalcPriorityProblem;
import yimei.jss.gp.data.DoubleData;
import yimei.jss.jobshop.*;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.simulation.OrderingDecisionSituation;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;
import yimei.jss.simulation.Simulation;
import yimei.jss.simulation.state.SystemState;

import java.util.*;

/**
 * Created by yimei on 22/09/16.
 */
public class BatchArrivalEvent extends AbstractEvent {

    private Job job;

    private Batch batch;

    public BatchArrivalEvent(double time, Batch batch) {
        super(time);
        this.batch = batch;
    }

    public Batch getBatch() {return batch;}

    public BatchArrivalEvent(Batch batch) {
        this(batch.getArrivalTime(), batch);
    }

    @Override
    public void trigger(Simulation simulation) {

//        GPRuleEvolutionState.sumOccurrance++;

/*        if(simulation.getWarmupBatches() < simulation.getNumBatchesArrived()){
        int idleMachine = 0;

        for (WorkCenter w : simulation.getSystemState().getWorkCenters()) {
            if (w.getReadyTime() <= simulation.getSystemState().getClockTime()) {
                idleMachine++;
//                GPRuleEvolutionState.sumIdleMachines++;
            }
        }

        System.out.println("The number of idle machines is " + idleMachine);}*/

        if(simulation.getNumBatchesArrived() > simulation.getWarmupBatches() + 1){

//        RandomOrdering(simulation);

//          SPTJobs(simulation);

//        OrderingJobs(simulation);

//        OrderingJobsUpdate(simulation);

//        startOrdering(simulation);

        orderingPairsUpdated(simulation);

//        MixedOrdering(simulation);

//        MixedOrderingV1(simulation);
        } else {
            RandomOrdering(simulation);
        }

        //------------------------------------------end: prioritise pairs --------------------------------------------------------


        //generate next batch, One batch by one batch.
        simulation.generateBatch();

    }

    private void SPTJobs(Simulation simulation) {

        //------------------------------------begin:prioritise jobs----------------------------------------------------------
        //calculate the priority values of each jobs in this batch by ordering rule
        OrderingDecisionSituation orderingDecisionSituation =
                new OrderingDecisionSituation(batch.getJobs(), simulation.getSystemState());


        orderingDecisionSituation.getJobQueue().sort(Comparator.comparingDouble(job -> job.getTotalProcTime()));

        for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
            Job job = orderingDecisionSituation.getJobQueue().get(j);
            new JobArrivalEvent(job).trigger(simulation);
        }

    }

    private void RandomOrdering(Simulation simulation) {

        //------------------------------------begin:prioritise jobs----------------------------------------------------------
        //calculate the priority values of each jobs in this batch by ordering rule
       OrderingDecisionSituation orderingDecisionSituation =
                new OrderingDecisionSituation(batch.getJobs(), simulation.getSystemState());

        //use ordering rule to determine the priority
//        nextJob(orderingDecisionSituation, simulation.getOrderingRule());
//
        //manual rules SPT/LPT
/*        for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
            Job job = orderingDecisionSituation.getJobQueue().get(j);
//            job.setPriority(job.getId()); //Manual rule: FIFS
            job.setPriority(job.getTotalProcTime()); //Manual rule: SPT
        }
        orderingDecisionSituation.getJobQueue().sort(Comparator.comparingDouble(job -> job.getPriority()));*/

        //Randomly

        Random random = new Random(Math.round(simulation.getClockTime()));
        Collections.shuffle(orderingDecisionSituation.getJobQueue(), random);

        for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
            Job job = orderingDecisionSituation.getJobQueue().get(j);
            new JobArrivalEvent(job).trigger(simulation);
        }

    }

    private void OrderingJobs(Simulation simulation) {

        //------------------------------------begin:prioritise jobs----------------------------------------------------------
        //calculate the priority values of each jobs in this batch by ordering rule
        OrderingDecisionSituation orderingDecisionSituation =
                new OrderingDecisionSituation(batch.getJobs(), simulation.getSystemState());

        //use ordering rule to determine the priority
        nextJob(orderingDecisionSituation, simulation.getOrderingRule());

        orderingDecisionSituation.getJobQueue().sort(Comparator.comparingDouble(job -> job.getPriority()));

        for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
            Job job = orderingDecisionSituation.getJobQueue().get(j);
            new JobArrivalEvent(job).trigger(simulation);
        }

    }

    private void OrderingJobsUpdate(Simulation simulation) {

        //------------------------------------begin:prioritise jobs----------------------------------------------------------
        //calculate the priority values of each jobs in this batch by ordering rule
        OrderingDecisionSituation orderingDecisionSituationOrigin =
                new OrderingDecisionSituation(batch.getJobs(), simulation.getSystemState());

        OrderingDecisionSituation orderingDecisionSituation = orderingDecisionSituationOrigin.clone();

        while (orderingDecisionSituation.getJobQueue().size() > 0) {

            //use ordering rule to determine the priority
            Job selectedJob = nextJob(orderingDecisionSituation, simulation.getOrderingRule());

            new JobArrivalEvent(selectedJob).trigger(simulation);

            orderingDecisionSituation.getJobQueue().remove(selectedJob);

        }
    }

    private void orderingPairsUpdated(Simulation simulation) {

        List<OperationOption> candidatePair = new ArrayList<>();
        for (Job job : batch.getJobs()) {
            for (OperationOption opOption : job.getOperation(0).getOperationOptions()) {
                candidatePair.add(opOption);
            }

            Operation operation = job.getOperation(0);

            for (OperationOption op : operation.getOperationOptions())
                op.setReadyTime(job.getReleaseTime());

        }

        OrderingDecisionSituation orderingDecisionSituation =
                new OrderingDecisionSituation(candidatePair,
                        simulation.getSystemState(), new String("pair"));

        //-------------------------------update the machine state-----------------------------------------------

        List<OperationOption> unUpdatedPairs = new ArrayList<>();

        while (orderingDecisionSituation.getQueue().size() > 0 || unUpdatedPairs.size() > 0) {

            //1.assign one pair, and remove all pairs with the same operation, update related pairs
            if(orderingDecisionSituation.getQueue().size() > 0) {

/*                //update the terminal
                for (WorkCenter workCenter:simulation.getSystemState().getWorkCenters()) {
                    workCenter.numberOfOperationsCanProcess = 0;
                    for (OperationOption operationOption: orderingDecisionSituation.getQueue()){
                        if(workCenter.getId() == operationOption.getWorkCenter().getId())
                            workCenter.numberOfOperationsCanProcess++;
                    }
                    for (OperationOption operationOption: unUpdatedPairs){
                        if(workCenter.getId() == operationOption.getWorkCenter().getId())
                            workCenter.numberOfOperationsCanProcess++;
                    }
                }*/

//                simulation.getOrderingRule().nextPairOption(orderingDecisionSituation); //use ordering rule to assign priority
//        simulation.getSequencingRule().nextPairOption(orderingDecisionSituation);
                simulation.getRoutingRule().nextPairOption(orderingDecisionSituation);
            }

            orderingDecisionSituation.getQueue().addAll(unUpdatedPairs);

            unUpdatedPairs.clear();

            orderingDecisionSituation.getQueue().sort(Comparator.comparingDouble(operationOption -> operationOption.getPriority()));

            OperationOption opOptionSelected = orderingDecisionSituation.getQueue().get(0);

            Job job = opOptionSelected.getJob();
            new OperationVisitEvent(job.getReleaseTime(), opOptionSelected).trigger(simulation);

            //remove all pairs with the same operation
            for (int o = 0; o < orderingDecisionSituation.getQueue().size(); o++) {
                OperationOption opOption = orderingDecisionSituation.getQueue().get(o);
                if(job.getId() == opOption.getJob().getId()){
                    orderingDecisionSituation.getQueue().remove(o);
                    o--;
                } else if (!(opOptionSelected.getWorkCenter().getId() == opOption.getWorkCenter().getId())){
                    unUpdatedPairs.add(opOption);
                    orderingDecisionSituation.getQueue().remove(o);
                    o--;
                }
            }

        }

    }



    private void MixedOrdering(Simulation simulation) {

        List[] sequencingInMachines = new List[simulation.getSystemState().getWorkCenters().size()];

        for (int i = 0; i < sequencingInMachines.length; i++) {
            List<OperationOption> candidateInOneMachine = new ArrayList<>();
            sequencingInMachines[i] = candidateInOneMachine;
        }


        List<OperationOption> candidatePair = new ArrayList<>();


        for (Job job : batch.getJobs()) {
            for (OperationOption opOption : job.getOperation(0).getOperationOptions()) {
                candidatePair.add(opOption);
            }

            Operation operation = job.getOperation(0);

            for (OperationOption op : operation.getOperationOptions())
                op.setReadyTime(job.getReleaseTime());

        }

        for (OperationOption op : candidatePair) {

            sequencingInMachines[op.getWorkCenter().getId()].add(op);

        }


        for (int m = 0; m < sequencingInMachines.length; m++) {

            Arrays.sort(sequencingInMachines, Comparator.comparingInt(List::size));

            List sequencingDecision = sequencingInMachines[m];
            //if there are some operations can be processed by this idle machine
            if (sequencingDecision.size() >= 1 && ((OperationOption) sequencingDecision.get(0)).getWorkCenter().getReadyTime() <= simulation.getSystemState().getClockTime()) {

                SequencingDecisionSituation sequencingDecisionSituation = new SequencingDecisionSituation(sequencingDecision, ((OperationOption) sequencingDecision.get(0)).getWorkCenter(), simulation.getSystemState());
                OperationOption opOption =
                        simulation.getSequencingRule().priorOperation(sequencingDecisionSituation);
                Job job = opOption.getJob();
//                job.arrivalTime += 0.001 * count;
//                count++;
//                job.releaseTime = job.arrivalTime;
                OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
//                simulation.addEvent(operationVisitEvent);
                operationVisitEvent.trigger(simulation);
                AbstractEvent nextEvent = simulation.getEventQueue().poll();
//                simulation.getSystemState().setClockTime(nextEvent.getTime());
                nextEvent.trigger(simulation);

                //remove options in the remaining queues
                for (int n = m + 1; n < sequencingInMachines.length; n++) {
                    for (int o = 0; o < sequencingInMachines[n].size(); o++) {
                        OperationOption op = (OperationOption) sequencingInMachines[n].get(o);
                        if (op.getOperation().equals(opOption.getOperation())) {
                            sequencingInMachines[n].remove(op);
                            o--;
                        }
                    }
                }
                //remove options related to this job
                candidatePair.removeAll(opOption.getOperation().getOperationOptions());
            }
            sequencingInMachines[m].clear();
        }

        //and then do routing
        List<OperationOption> routingDecision = new ArrayList<>();
        for (int c = 0; c < candidatePair.size(); c++) {
            OperationOption option = candidatePair.get(c);
            routingDecision.add(option);
            candidatePair.remove(option);
            c--;
            for (int n = c + 1; c < candidatePair.size(); n++) {
                if (candidatePair.isEmpty())
                    break;
                OperationOption optionNext = candidatePair.get(n);
                if (optionNext.getOperation().equals(option.getOperation())) {
                    routingDecision.add(optionNext);
                    candidatePair.remove(optionNext);
                    n--;
                } else
                    break;
            }

            RoutingDecisionSituation routingDecisionSituation = new RoutingDecisionSituation(routingDecision, simulation.getSystemState());
            OperationOption opOption =
                    simulation.getRoutingRule().nextOperationOption(routingDecisionSituation);
            Job job = opOption.getJob();
//                job.arrivalTime += 0.001 * count;
//                count++;
//                job.releaseTime = job.arrivalTime;
            OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
//                simulation.addEvent(operationVisitEvent);
            operationVisitEvent.trigger(simulation);

            routingDecision.clear();
        }


    }

    //First assign to idle machines , and then assign to busy machines based on ready time
    /*private void MixedOrderingV1(Simulation simulation) {


        List[] sequencingInMachines = new List[simulation.getSystemState().getWorkCenters().size()];

        for (int i = 0; i < sequencingInMachines.length; i++) {
            List<OperationOption> candidateInOneMachine = new ArrayList<>();
            sequencingInMachines[i] = candidateInOneMachine;

        }


        for (Job job : batch.getJobs()) {
            for (OperationOption op : job.getOperation(0).getOperationOptions()) {
                sequencingInMachines[op.getWorkCenter().getId()].add(op);
            }

            Operation operation = job.getOperation(0);

            for (OperationOption op : operation.getOperationOptions())
                op.setReadyTime(job.getReleaseTime());
        }


        for (int m = 0; m < sequencingInMachines.length; m++) {

            Arrays.sort(sequencingInMachines, Comparator.comparingInt(List::size));

            List sequencingDecision = sequencingInMachines[m];
            //if there are some operations can be processed by this idle machine
            if (sequencingDecision.size() >= 1 && ((OperationOption) sequencingDecision.get(0)).getWorkCenter().getReadyTime() <= simulation.getSystemState().getClockTime()) {

                SequencingDecisionSituation sequencingDecisionSituation = new SequencingDecisionSituation(sequencingDecision, ((OperationOption) sequencingDecision.get(0)).getWorkCenter(), simulation.getSystemState());
                OperationOption opOption =
                        simulation.getSequencingRule().priorOperation(sequencingDecisionSituation);
                Job job = opOption.getJob();
                OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
                operationVisitEvent.trigger(simulation);
                AbstractEvent nextEvent = simulation.getEventQueue().poll();
                nextEvent.trigger(simulation);

                //remove options in the remaining queues
                for (int n = 0; n < sequencingInMachines.length; n++) {
                    for (int o = 0; o < sequencingInMachines[n].size(); o++) {
                        OperationOption op = (OperationOption) sequencingInMachines[n].get(o);
                        if (op.getOperation().equals(opOption.getOperation())) {
                            sequencingInMachines[n].remove(op);
                            o--;
                        }
                    }
                }
            }
        }


        while (Arrays.stream(sequencingInMachines).allMatch(List::isEmpty)) {

            WorkCenter earliestMachine = simulation.getSystemState().getWorkCenters().stream().min(Comparator.comparingDouble(WorkCenter::getReadyTime)).get();

            SequencingDecisionSituation sequencingDecisionSituation = new SequencingDecisionSituation(sequencingInMachines[],
            ((OperationOption) sequencingDecision.get(0)).getWorkCenter(), simulation.getSystemState());
            OperationOption opOption =
                    simulation.getSequencingRule().priorOperation(sequencingDecisionSituation);
            Job job = opOption.getJob();
            OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
            operationVisitEvent.trigger(simulation);
            AbstractEvent nextEvent = simulation.getEventQueue().poll();
            nextEvent.trigger(simulation);

            //remove options in the remaining queues
            for (int n = 0; n < sequencingInMachines.length; n++) {
                for (int o = 0; o < sequencingInMachines[n].size(); o++) {
                    OperationOption op = (OperationOption) sequencingInMachines[n].get(o);
                    if (op.getOperation().equals(opOption.getOperation())) {
                        sequencingInMachines[n].remove(op);
                        MachineQueue.get(op.getWorkCenter());
                        MachineQueue.get(op.getWorkCenter()).remove(op);
                        o--;
                    }
                }
            }
        }
    }*/





public void startOrdering(Simulation simulation) {


/*        ArrayList<Job> unProcessedJobs = new ArrayList<>();

        for (Job job : systemState.getJobsInSystem()) {
            if (job.getArrivalTime() >= systemState.getClockTime()) {
                unProcessedJobs.add(job);
                if (eventQueue.contains(new JobArrivalEvent(job)))
                    eventQueue.remove(new JobArrivalEvent(job)); // remove unprocessedJobs
            }
        }*/

/*                for (int j = 0; j < unProcessedJobs.size(); j++) {
                    Job job = unProcessedJobs.get(j);
                    job.arrivalTime += 0.001 * j;
                    job.releaseTime = job.arrivalTime;
                    eventQueue.add(new JobArrivalEvent(job));
                }*/


    //add (operation,machine) to the candidate queue
    List<OperationOption> candidatePair = new ArrayList<>();
    for (Job job : batch.getJobs()) {
        for (OperationOption opOption : job.getOperation(0).getOperationOptions()) {
            candidatePair.add(opOption);
        }

        Operation operation = job.getOperation(0);

        for (OperationOption op : operation.getOperationOptions())
            op.setReadyTime(job.getReleaseTime());

    }

    OrderingDecisionSituation orderingDecisionSituation =
            new OrderingDecisionSituation(candidatePair,
                    simulation.getSystemState(), new String("pair"));

    Set<Integer> filterPair = new HashSet<>();
//        simulation.getOrderingRule().nextPairOption(orderingDecisionSituation);
//        simulation.getSequencingRule().nextPairOption(orderingDecisionSituation);
            simulation.getRoutingRule().nextPairOption(orderingDecisionSituation);

    orderingDecisionSituation.getQueue().sort(Comparator.comparingDouble(OperationOption -> OperationOption.getPriority()));

    //we need to remove all option belong to one operation
    for (int o = 0; o < orderingDecisionSituation.getQueue().size(); o++) {
        OperationOption opOption = orderingDecisionSituation.getQueue().get(o);
        if (filterPair.contains(opOption.getJob().getId()) ) {
            orderingDecisionSituation.getQueue().remove(o);
            o--;
        } else {

          // -------------------directly trigger---------------------------------------
            filterPair.add(opOption.getJob().getId());

            Job job = opOption.getJob();

            new OperationVisitEvent(job.getReleaseTime(), opOption).trigger(simulation);

//            int before = simulation.getEventQueue().size();
//            OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
//            operationVisitEvent.trigger(simulation);
//            int after = simulation.getEventQueue().size();
//            if(after - before > 0) {
//                AbstractEvent nextEvent1 = simulation.getEventQueue().poll();
//                nextEvent1.trigger(simulation);
//            }

//            job.arrivalTime += o*0.0000000001;
//            job.releaseTime = job.arrivalTime;
//            OperationVisitEvent operationVisitEvent = new OperationVisitEvent(job.getReleaseTime(), opOption);
//            simulation.addEvent(operationVisitEvent);

        }

    }

}

public Job nextJob(OrderingDecisionSituation orderingDecisionSituation, AbstractRule orderingRule) {

    List<Job> queue = orderingDecisionSituation.getJobQueue();
    SystemState systemState = orderingDecisionSituation.getSystemState();
    //================original=================
    //==================start==================
    Job bestJob = queue.get(0);
    bestJob.setPriority(priority(bestJob, systemState, orderingRule));
    // loop all the options, save the best one as "selected" one
    for (int i = 1; i < queue.size(); i++) {
        Job job = queue.get(i);
        job.setPriority(priority(job, systemState, orderingRule));

        if (job.priorTo(bestJob)) {
            bestJob = job;
        }
    }
    return bestJob;// this links which job will be chosen.
}

public double priority(Job job,
                       SystemState systemState, AbstractRule orderingRule) {
    CalcPriorityProblem calcPrioProb =
            new CalcPriorityProblem(job, systemState);

    DoubleData tmp = new DoubleData();

    ((GPRule) orderingRule).gpTree.child.eval(null, 0, tmp, null, null, calcPrioProb);

    return tmp.value;
}

@Override
public void addSequencingDecisionSituation(Simulation simulation,
                                           List<SequencingDecisionSituation> situations,
                                           int minQueueLength) {
    trigger(simulation);

}

@Override
public void addRoutingDecisionSituation(Simulation simulation,
                                        List<RoutingDecisionSituation> situations,
                                        int minQueueLength) {
    trigger(simulation);
}

public void addOrderingDecisionSituation(Simulation simulation,
                                            List<OrderingDecisionSituation> situations,
                                            int minQueueLength) {


    OrderingDecisionSituation orderingDecisionSituation =
            new OrderingDecisionSituation(batch.getJobs(), simulation.getSystemState());

    if(orderingDecisionSituation.getJobQueue().size() >= minQueueLength) {
        situations.add(orderingDecisionSituation.clone());
    }

    //use ordering rule to determine the priority

        for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
            Job job = orderingDecisionSituation.getJobQueue().get(j);
            job.setPriority(job.getId()); //Manual rule: FIFS
//            job.setPriority(job.getOperation(0).getMedianProcessTime()); //Manual rule: SPT
        }

    //learned ordering rule
    orderingDecisionSituation.getJobQueue().sort(Comparator.comparingDouble(job -> job.getPriority()));

    for (int j = 0; j < orderingDecisionSituation.getJobQueue().size(); j++) {
        Job job = orderingDecisionSituation.getJobQueue().get(j);

        job.arrivalTime += j*0.0001;
        job.releaseTime = job.arrivalTime;
        simulation.addEvent(new JobArrivalEvent(job));

    }

    simulation.generateBatch();

}

@Override
public String toString() {
    return String.format("%.1f: batch %d arrives.\n", time, batch.getId());
}

@Override
public int compareTo(AbstractEvent other) {
    if (time < other.time)
        return -1;

    if (time > other.time)
        return 1;

    if (other instanceof BatchArrivalEvent) {
        BatchArrivalEvent otherJAE = (BatchArrivalEvent) other;

        if (job.getId() < otherJAE.job.getId())
            return -1;

        if (job.getId() > otherJAE.job.getId())
            return 1;
    }

    return -1;
}
}
