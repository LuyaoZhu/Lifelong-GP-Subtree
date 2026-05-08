package yimei.jss.simulation;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import yimei.jss.gp.CalcPriorityProblem;
import yimei.jss.gp.data.DoubleData;
import yimei.jss.jobshop.*;
import yimei.jss.jobshop.Process;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.SPT;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.simulation.event.*;
import yimei.jss.simulation.state.SystemState;
import yimei.util.random.*;

import java.io.*;
import java.util.*;

import static yimei.jss.algorithm.downSampling.GPRuleEvolutionStateDownSampling.*;

/**
 * The dynamic simulation -- discrete event simulation
 * <p>
 * Created by yimei on 22/09/16.
 */
public class DynamicSimulation extends Simulation {

    public final static int SEED_ROTATION = 10000;

    private static final Logger log = LoggerFactory.getLogger(DynamicSimulation.class);

    private long seed;
    private RandomDataGenerator randomDataGenerator;

    private final int minNumOperations;
    private final int maxNumOperations;

    public final int minBatchSize;
    public final int maxBatchSize;
    private final double utilLevel;
    private final double dueDateFactor;
    private final boolean revisit;

    private AbstractIntegerSampler numOperationsSampler;
    //modified by fzhang, 17.04.2018  in order to set options from 2 to 10
    //private AbstractIntegerSampler numOptionsSampler;

    private AbstractRealSampler procTimeSampler;
    private AbstractRealSampler interArrivalTimeSampler;
    private AbstractRealSampler jobWeightSampler;

    private AbstractRealSampler dueDateFactorSampler;

    private AbstractIntegerSampler numJobSampler;

    private AbstractRealSampler jobRevenueSampler;

    public long getSeed(){return seed;}

    public enum ShopType {
        UNIFORM,
        NORMAL,
        EXPONENTIAL,
        GAMMA,
        UNIFORM1,
        NORMAL1,
        EXPONENTIAL1,
        GAMMA1

    }


    /*private DynamicSimulation(long seed,
                              AbstractRule sequencingRule,
                              AbstractRule routingRule,
                              int numWorkCenters,
                              int numJobsRecorded,
                              int warmupJobs,
                              int minNumOperations,
                              int maxNumOperations,
                              double utilLevel,
                              double dueDateFactor,
                              boolean revisit,
                              AbstractIntegerSampler numOperationsSampler,
                              //modified by fzhang, 17.04.2018
                              //AbstractIntegerSampler numOptionsSampler,

                              AbstractRealSampler procTimeSampler,
                              AbstractRealSampler interArrivalTimeSampler,
                              AbstractRealSampler jobWeightSampler,
                              AbstractIntegerSampler numJobSampler) {
        super(sequencingRule, routingRule, numWorkCenters, numJobsRecorded, warmupJobs);

        this.seed = seed;
        this.randomDataGenerator = new RandomDataGenerator();
        this.randomDataGenerator.reSeed(seed);

        this.minNumOperations = minNumOperations;
        this.maxNumOperations = maxNumOperations;
        this.utilLevel = utilLevel;
        this.dueDateFactor = dueDateFactor;
        this.revisit = revisit;

        this.numOperationsSampler = numOperationsSampler;
        //modified by fzhang 17.04.2018
        //this.numOptionsSampler = numOptionsSampler;

        this.procTimeSampler = procTimeSampler;
        this.interArrivalTimeSampler = interArrivalTimeSampler;
        this.jobWeightSampler = jobWeightSampler;
        this.numJobSampler = numJobSampler;

        setInterArrivalTimeSamplerMean();

        // Create the work centers, with empty queue and ready to go initially.
        for (int i = 0; i < numWorkCenters; i++) {
            systemState.addWorkCenter(new WorkCenter(i));
        }

        setup();
    }*/

/*    public DynamicSimulation(long seed,
                             AbstractRule sequencingRule,
                             AbstractRule routingRule,
                             int numWorkCenters,
                             int numJobsRecorded,
                             int warmupJobs,
                             int minNumOperations,
                             int maxNumOperations,
                             double utilLevel,
                             double dueDateFactor,
                             boolean revisit) {
        this(seed, sequencingRule, routingRule, numWorkCenters, numJobsRecorded, warmupJobs,
                minNumOperations, maxNumOperations, utilLevel, dueDateFactor, revisit,
                //here, specifiy the range of UniformIntegerSample to (1,10)
                new UniformIntegerSampler(minNumOperations, maxNumOperations), //these two values will be changed during the evolutionary process, because different models are called.
                //the surrogate model will set them to 1 and 5, but the original model will set them to 1 and 10
                //when calculate the phenotype, in this code, full simulation is used, they will be set to 10 and 10.
                //modified by fzhang 17.04.2018
                //new UniformIntegerSampler(1, numWorkCenters), //in this way, whether add this parameter or not is the same
                //new UniformIntegerSampler(1, 5), //one operation only can be processed at 5 machines

                new UniformSampler(1, 99),
                new ExponentialSampler(),
                new TwoSixTwoSampler(),
                new UniformIntegerSampler(10, 50));
    }*/

    private DynamicSimulation(long seed,
                              AbstractRule sequencingRule,
                              AbstractRule routingRule,
                              int numWorkCenters,
                              int numBatchesRecorded,
                              int warmupBatches,
                              int minNumOperations,
                              int maxNumOperations,
                              int minBatchSize,
                              int maxBatchSize,
                              double utilLevel,
//                              double dueDateFactor,
                              double dueDateFactor,
                              boolean revisit,
                              AbstractIntegerSampler numOperationsSampler,
                              //modified by fzhang, 17.04.2018
                              //AbstractIntegerSampler numOptionsSampler,

                              AbstractRealSampler procTimeSampler,
                              AbstractRealSampler interArrivalTimeSampler,
                              AbstractRealSampler jobWeightSampler,
                              AbstractIntegerSampler numJobSampler,
                              AbstractRealSampler jobRevenueSampler) {
        super(sequencingRule, routingRule, numWorkCenters, numBatchesRecorded, warmupBatches);

        this.seed = seed;
        this.randomDataGenerator = new RandomDataGenerator();
        this.randomDataGenerator.reSeed(seed);

        this.minNumOperations = minNumOperations;
        this.maxNumOperations = maxNumOperations;
        this.minBatchSize = minBatchSize;
        this.maxBatchSize = maxBatchSize;
        this.utilLevel = utilLevel;
        this.dueDateFactor = dueDateFactor;
        this.revisit = revisit;

        this.numOperationsSampler = numOperationsSampler;
        //modified by fzhang 17.04.2018
        //this.numOptionsSampler = numOptionsSampler;
        this.procTimeSampler = procTimeSampler;
        this.interArrivalTimeSampler = interArrivalTimeSampler;
        this.jobWeightSampler = jobWeightSampler;
        this.numJobSampler = numJobSampler;
        this.jobRevenueSampler = jobRevenueSampler;

        setInterArrivalTimeSamplerMean();

        defaultSequencingRule = new SPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        if (!(optimisationObjective == null) && optimisationObjective.contains("weight")) {
            defaultSequencingRule = new WSPT(RuleType.SEQUENCING);
        }
        defaultRoutingRule = new WIQ(RuleType.ROUTING);


        // Create the work centers, with empty queue and ready to go initially.
        for (int i = 0; i < numWorkCenters; i++) {
            systemState.addWorkCenter(new WorkCenter(i));
        }

        setup();
    }

    public DynamicSimulation(long seed,
                             AbstractRule sequencingRule,
                             AbstractRule routingRule,
                             int numWorkCenters,
                             int numBatchesRecorded,
                             int warmupBatches,
                             int minNumOperations,
                             int maxNumOperations,
                             int minBatchSize,
                             int maxBatchSize,
                             double utilLevel,
                             double dueDateFactor,
                             boolean revisit) {
        this(seed, sequencingRule, routingRule, numWorkCenters, numBatchesRecorded, warmupBatches,
                minNumOperations, maxNumOperations, minBatchSize, maxBatchSize, utilLevel, dueDateFactor, revisit,
                //here, specifiy the range of UniformIntegerSample to (1,10)
                new UniformIntegerSampler(minNumOperations, maxNumOperations), //traditional process form
//                new BinomialIntegerSampler(minNumOperations, maxNumOperations,0.7), //small product manufactory
//                new NormalIntegerSampler( (int)((minNumOperations+maxNumOperations)*0.5), 1, minNumOperations, maxNumOperations), //customised form

                //the surrogate model will set them to 1 and 5, but the original model will set them to 1 and 10
                //when calculate the phenotype, in this code, full simulation is used, they will be set to 10 and 10.
                //modified by fzhang 17.04.2018
                //new UniformIntegerSampler(1, numWorkCenters), //in this way, whether add this parameter or not is the same
                //new UniformIntegerSampler(1, 5), //one operation only can be processed at 5 machines

                new UniformSampler(1, 99),           //traditional process form
//                new GammaSampler(5,6,1,99),
//                new ExponentialSampler(9,1,99),//small product manufactory
//                new NormalSampler(50,10,1,99),    //customised form
//                new LogNormalSampler(4,0.7,1,99),

                new ExponentialSampler(),
                new TwoSixTwoSampler(),
//                new PoissonSampler(),
                new UniformIntegerSampler(minBatchSize, maxBatchSize),
//                new NormalIntegerSampler( (int)((minBatchSize+maxBatchSize)*0.5), 1, minBatchSize, maxBatchSize), //customised form
//                new BinomialIntegerSampler(minBatchSize, maxBatchSize,0.7), //small product manufactory
                new LogNormalSampler(0, 1));
    }

    public DynamicSimulation(long seed,
                             AbstractRule sequencingRule,
                             AbstractRule routingRule,
                             int numWorkCenters,
                             int numBatchesRecorded,
                             int warmupBatches,
                             int minNumOperations,
                             int maxNumOperations,
                             int minBatchSize,
                             int maxBatchSize,
                             double utilLevel,
                             double dueDateFactor,
                             boolean revisit,
                             ShopType type) {

        this(seed, sequencingRule, routingRule,
                numWorkCenters, numBatchesRecorded, warmupBatches,
                minNumOperations, maxNumOperations,
                minBatchSize, maxBatchSize,
                utilLevel, dueDateFactor, revisit,

                // 所有逻辑都写在参数表达式里，不会违反 “this 必须第一行”
                new UniformIntegerSampler(minNumOperations, maxNumOperations),
                createProcessingTimeSampler(type),
                new ExponentialSampler(),       // 你原来的到达/间隔分布
                new TwoSixTwoSampler(),         // 你原来的 due date 分布
                new UniformIntegerSampler(minBatchSize, maxBatchSize),
                new LogNormalSampler(0, 1));    // 你原来的权重分布

    }

    private static AbstractRealSampler createProcessingTimeSampler(ShopType type) {
        switch (type) {
            case UNIFORM:
                return new UniformSampler(1, 99);

            case NORMAL:
                return new NormalSampler(80,10,1,99);

            case EXPONENTIAL:
                return new ExponentialSampler(20,1,99);

            case GAMMA:
                return new GammaSampler(2,25,1,99);
                //return new GammaSampler(5,17,1,99);

            case UNIFORM1:
                return new UniformSampler(50, 300);

            case NORMAL1:
                return new NormalSampler(100,30,1,300);

            case EXPONENTIAL1:
                return new ExponentialSampler(20,1,300);

            case GAMMA1:
                return new GammaSampler(5,40,1,300);

            default:
                return new UniformSampler(1, 99);
        }
    }



    public int getNumWorkCenters() {
        return numWorkCenters;
    }

    public int getNumJobsRecorded() {
        return numJobsRecorded;
    }

    public int getWarmupJobs() {
        return warmupJobs;
    }

    public int getMinNumOperations() {
        return minNumOperations;
    }

    public int getMaxNumOperations() {
        return maxNumOperations;
    }

    public double getUtilLevel() {
        return utilLevel;
    }

/*    public double getDueDateFactor() {
        return dueDateFactor;
    }*/

    public boolean isRevisit() {
        return revisit;
    }

    public RandomDataGenerator getRandomDataGenerator() {
        return randomDataGenerator;
    }

    public AbstractIntegerSampler getNumOperationsSampler() {
        return numOperationsSampler;
    }

    public AbstractRealSampler getProcTimeSampler() {
        return procTimeSampler;
    }

    public AbstractRealSampler getInterArrivalTimeSampler() {
        return interArrivalTimeSampler;
    }

    public AbstractRealSampler getJobWeightSampler() {
        return jobWeightSampler;
    }

    public AbstractIntegerSampler getNumJobSampler() {
        return numJobSampler;
    }

    public AbstractRealSampler getJobRevenueSampler() {
        return jobRevenueSampler;
    }

    @Override
    public void setup() {
        numJobsArrived = 0;
        numBatchesArrived = 0;
        throughput = 0;
//        generateJob();
//        setInterArrivalTimeSamplerMean(1);
        generateBatch();
    }

    @Override
    public void resetState() {
        systemState.reset();
        eventQueue.clear();
        count = 0;
        setup();
    }

    @Override
    public void reset() {
        reset(seed);
    }

    public void reset(long seed) {
        reseed(seed);
        resetState();
    }

    public void reseed(long seed) {
        this.seed = seed;
        randomDataGenerator.reSeed(seed);
    }

    @Override
    public void rotateSeed() {//this is use for changing seed value in next generation
        //this only relates to generation
        seed += SEED_ROTATION;
        reset();
        //System.out.println(seed);//when seed=0, after Gen0, the value is 10000, after Gen1, the value is 20000....
    }

    @Override
    public void generateJob() {
        //runExperiments();
        //modified by fzhang 15.5.2018  to avoid negative time  finallly decide to keep double type: to avoid same arrival time
        double arrivalTime = getClockTime()
                + interArrivalTimeSampler.next(randomDataGenerator);

        double weight = jobWeightSampler.next(randomDataGenerator);
        double revenue = jobRevenueSampler.next(randomDataGenerator);
        Job job = new Job(numJobsArrived, new ArrayList<>(),
                arrivalTime, arrivalTime, 0, weight, revenue);
        int numOperations = numOperationsSampler.next(randomDataGenerator);

        for (int i = 0; i < numOperations; i++) {
            Operation o = new Operation(job, i);
            //modified by fzhang 17.04.2018
            //int numOptions = numOptionsSampler.next(randomDataGenerator);
            int numOptions = numOperationsSampler.next(randomDataGenerator);
            //System.out.println("numOptions: "+numOptions);

            int[] route = randomDataGenerator.nextPermutation(numWorkCenters, numOptions);
            //nextPermutation(n,k)
            //Generates an integer array of length k whose entries are selected randomly, without repetition, from the integers 0, ..., n - 1 (inclusive).

            //modified by fzhang  14.5.2018  in order to avoid negative or positive time(equal = 0)  finallly decide to keep double type
            double procTime = procTimeSampler.next(randomDataGenerator); //use same proc time for all options for now
            //================start==========
/*            for (int j = 0; j < numOptions; j++) {//9
                double procTime = procTimeSampler.next(randomDataGenerator);
                procTime = (procTime - 1) / 98 * 40 + 30;
                o.addOperationOption(new OperationOption(o, j, procTime, systemState.getWorkCenter(route[j])));
            }*/
            //==========end===========

            //modified by fzhang  29.5.2018  set different processing time for different machines
           /* for (int j = 0; j < numOptions; j++) {
            	double procTime = procTimeSampler.next(randomDataGenerator);
                o.addOperationOption(new OperationOption(o,j,procTime,systemState.getWorkCenter(route[j])));
            }
*/

            //fzhang 2019.6.22 set different processtime to each machine
            //============================start========================================================
            /*double ptmean =  procTimeSampler.next(randomDataGenerator);// set processtime of each option
            AbstractRealSampler ptnsampler=new NormalSampler(ptmean, ptmean/10);

            for (int j = 0; j < numOptions; j++) {
                double procTime= ptnsampler.next(randomDataGenerator);
                o.addOperationOption(new OperationOption(o,j,procTime,systemState.getWorkCenter(route[j])));
            }*/
            //==============================end=================================================

            job.addOperation(o);
        }

        job.linkOperations();
        //just set totalProcTime to average value, as we don't know which option will be chosen
        //this is just used to define dueDate value
        double totalProcTime = numOperations * procTimeSampler.getMean();
        double dueDate = job.getReleaseTime() + dueDateFactor * totalProcTime;

        job.setDueDate(dueDate);
//        if (job.getId() > 501) {
//            int a  = 1;
//        }

        systemState.addJobToSystem(job);
        numJobsArrived++;

        eventQueue.add(new JobArrivalEvent(job));
    }

    // luyao on 20/12/24
    public void generateBatch() {

        //runExperiments();
        //modified by fzhang 15.5.2018  to avoid negative time  finally decide to keep double type: to avoid same arrival time
        double arrivalTime = 0;

        arrivalTime = getClockTime()
                + interArrivalTimeSampler.next(randomDataGenerator);

        double batchWeight = jobWeightSampler.next(randomDataGenerator);

        int currentBatchNumJob = (int) numJobSampler.next(randomDataGenerator);


//        System.out.print(currentBatchNumJob + " , ");

        //adaptively adjust the inner arrival time
//        setInterArrivalTimeSamplerMean(currentBatchNumJob);
//        setInterArrivalTimeSamplerMean(30);

        Batch batch = new Batch(numBatchesArrived, currentBatchNumJob, 0, arrivalTime, arrivalTime, 0, batchWeight);


        for (int jobIndex = 0; jobIndex < currentBatchNumJob; jobIndex++) {

            double weight;

            if (optimisationObjective == null || !optimisationObjective.startsWith("batch"))
                weight = jobWeightSampler.next(randomDataGenerator);
            else
                weight = batchWeight;

            double revenue = jobRevenueSampler.next(randomDataGenerator);

//            arrivalTime += 0.1;
            Job job = new Job(numJobsArrived, new ArrayList<>(),
                    arrivalTime, arrivalTime, 0, weight, revenue);
            int numOperations = numOperationsSampler.next(randomDataGenerator);

            job.operationRemaining = numOperations;

            for (int i = 0; i < numOperations; i++) {
                Operation o = new Operation(job, i);
                //modified by fzhang 17.04.2018
                //int numOptions = numOptionsSampler.next(randomDataGenerator);
                int numOptions = numOperationsSampler.next(randomDataGenerator);
                //System.out.println("numOptions: "+numOptions);


                int[] route = randomDataGenerator.nextPermutation(numWorkCenters, numOptions);
                //nextPermutation(n,k)
                //Generates an integer array of length k whose entries are selected randomly, without repetition, from the integers 0, ..., n - 1 (inclusive).

                //modified by fzhang  14.5.2018  in order to avoid negative or positive time(equal = 0)  finallly decide to keep double type
                //modified by luyao 1.10.2024 the processing time of one operation on candidate machines should differ.
                double procTimeBase = procTimeSampler.next(randomDataGenerator); //use same proc time for all options for now
                double noise = 10;
                double lowerBound = Math.max(procTimeSampler.getLower(), procTimeBase - noise);
                double upperBound = Math.min(procTimeSampler.getUpper(), procTimeBase + noise);
                //================start==========
                for (int j = 0; j < numOptions; j++) {//9

                    // change around the base-procTime i.e.,base-procTime is 50, then the processing time on candidate machines may from 40 to 60.
                    double procTime = randomDataGenerator.nextUniform(lowerBound,upperBound);
//                    double procTime = procTimeBase;
//                    double perturbationFactor = randomDataGenerator.nextUniform(0.8,1.2);
//                    double procTime = Math.max(1.0, Math.min(99.0, procTimeBase*perturbationFactor));

                    o.addOperationOption(new OperationOption(o, j, procTime, systemState.getWorkCenter(route[j])));
                }
                //==========end===========

                //modified by fzhang  29.5.2018  set different processing time for different machines
           /* for (int j = 0; j < numOptions; j++) {
            	double procTime = procTimeSampler.next(randomDataGenerator);
                o.addOperationOption(new OperationOption(o,j,procTime,systemState.getWorkCenter(route[j])));
            }
*/

                //fzhang 2019.6.22 set different processtime to each machine
                //============================start========================================================
            /*double ptmean =  procTimeSampler.next(randomDataGenerator);// set processtime of each option
            AbstractRealSampler ptnsampler=new NormalSampler(ptmean, ptmean/10);

            for (int j = 0; j < numOptions; j++) {
                double procTime= ptnsampler.next(randomDataGenerator);
                o.addOperationOption(new OperationOption(o,j,procTime,systemState.getWorkCenter(route[j])));
            }*/
                //==============================end=================================================

                job.addOperation(o);
            }

            job.linkOperations();
            //just set totalProcTime to average value, as we don't know which option will be chosen
            //this is just used to define dueDate value
//            double totalProcTime = numOperations * procTimeSampler.getMean();

            //set job's due date
            double dueDate = job.getReleaseTime() + dueDateFactor * job.getTotalProcTime();
            job.setDueDate(dueDate);
            job.workLoadRemaining = job.getTotalProcTime();

            //set batch's due date
            batch.totalProcTime += job.getTotalProcTime();
            batch.revenue += job.getRevenue();


            systemState.addJobToSystem(job);
            numJobsArrived++;
            batch.addJob(job);
            job.batch = batch;

//            eventQueue.add(new JobArrivalEvent(job));

        }
        systemState.addBatchToSystem(batch);
        numBatchesArrived++;

        double maxProcessTimeJob = batch.getJobs().stream()
                .mapToDouble(Job::getTotalProcTime)
                .max()
                .orElseThrow(() -> new RuntimeException("Empty job list"));
        double parallelProcessTime = batch.totalProcTime / numWorkCenters;
        batch.estimatedProcTime = Math.max(maxProcessTimeJob, parallelProcessTime);

        double dueDate = batch.getReleaseTime() + dueDateFactor * batch.estimatedProcTime;

        batch.dueDate = dueDate;

        //this is for the objectives related to batch
        if (!(optimisationObjective == null) && optimisationObjective.startsWith("batch")) {
            for (Job job : batch.getJobs()) {
                job.setDueDate(dueDate);
            }
        }

        eventQueue.add(new BatchArrivalEvent(batch));

    }

    public void runRecordSimulation() {

        while (!eventQueue.isEmpty() && throughput < numBatchesRecorded) { //numJobsRecorded == 5000

            //if the next event is batchArrival (the beginning of each case), we will add it to the recorded simulation
            if(eventQueue.iterator().next() instanceof BatchArrivalEvent &&
               startBatchID.contains(((BatchArrivalEvent)eventQueue.iterator().next()).getBatch().getId())) {
                recordedSimulations.add(this.deepCloneRuntimeNoReflect());
            }

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
            if (count > 100000) {
                count = 0;
                systemState.setClockTime(Double.MAX_VALUE);
                eventQueue.clear();
            }

            //===================ignore busy machine here==============================
            //when nextEvent was done, check the numOpsInQueue
            for (WorkCenter w : systemState.getWorkCenters()) {
                if (w.numOpsInQueue() > 100) {
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

/*        for (Batch batch:systemState.getBatchesCompleted()){
            if(batch.getId() == ){
                eventQueue.clear();
            }
        }*/

    }

    public void assignJobsPriority(OrderingDecisionSituation orderingDecisionSituation) {

        List<Job> queue = orderingDecisionSituation.getJobQueue();
        SystemState systemState = orderingDecisionSituation.getSystemState();
        //================original=================
        //==================start==================
/*        Job bestJobOption = queue.get(0);
        bestJobOption
                .setPriority(priority(bestJobOption, systemState));
        // loop all the options, save the best one as "selected" one
        for (int i = 1; i < queue.size(); i++) {
            Job jobOption = queue.get(i);
            jobOption.setPriority(priority(jobOption, systemState));
            if (jobOption.priorTo(bestJobOption)) {
                bestJobOption = jobOption;
            }
        }*/

        for (int j = 0; j < queue.size(); j++) {
            queue.get(j).setPriority(priority(queue.get(j), systemState));
        }

    }

    public double priority(Job job,
                           SystemState systemState) {
        CalcPriorityProblem calcPrioProb =
                new CalcPriorityProblem(job, systemState);

        DoubleData tmp = new DoubleData();

        ((GPRule) this.getOrderingRule()).gpTree.child.eval(null, 0, tmp, null, null, calcPrioProb);

        return tmp.value;
    }

//    private void runExperiments() {
//        double interArrivalSum = 0.0;
//        double numOperationsSum = 0.0;
//        double numOptionsSum = 0.0;
//        double procTimeSum = 0.0;
//        int numRuns = 5000000;
//
//        for (int i = 0; i < numRuns; ++i) {
//            interArrivalSum += interArrivalTimeSampler.next(randomDataGenerator);
//            numOperationsSum += numOperationsSampler.next(randomDataGenerator);
//            numOptionsSum += numOperationsSampler.next(randomDataGenerator);
//            procTimeSum += procTimeSampler.next(randomDataGenerator);
//        }
//        System.out.println("Average interarrival time: "+interArrivalSum/numRuns);
//        System.out.println("Average num operations: "+numOperationsSum/numRuns);
//        System.out.println("Average num options: "+numOptionsSum/numRuns);
//        System.out.println("Average procedure time: "+procTimeSum/numRuns);
//        System.out.println();
//    }

    //control the inter time of job arrival
    public double interArrivalTimeMean(int numWorkCenters,
                                       int minNumOps,
                                       int maxNumOps,
                                       double utilLevel) {
        double meanNumOps = 0.5 * (minNumOps + maxNumOps); //(1+10)/2=5.5 average operations for a job is 5.5
        double meanProcTime = procTimeSampler.getMean(); //(1+99)/2=50   average processing time for an operation is 50
        double meanNumJobs = numJobSampler.getMean();

        //for machines with same capacity, this return value is the same.
        //for machines with different capacities, this return value is different because utilLevel is dynamic
        return (meanNumOps * meanProcTime * meanNumJobs) / (utilLevel * numWorkCenters); // the time to processing a job on each workcenter
    }

    public void setInterArrivalTimeSamplerMean() {
        double mean = interArrivalTimeMean(numWorkCenters, minNumOperations, maxNumOperations, utilLevel);
        interArrivalTimeSampler.setMean(mean);
    }

    public List<SequencingDecisionSituation> sequencingDecisionSituations(int minQueueLength) {
        List<SequencingDecisionSituation> sequencingDecisionSituations = new ArrayList<>();
        //numJobsRecorded

        while (!eventQueue.isEmpty() && throughput < numBatchesRecorded) {
            AbstractEvent nextEvent = eventQueue.poll();
//            System.out.println("throughput "+throughput);
            systemState.setClockTime(nextEvent.getTime());
            nextEvent.addSequencingDecisionSituation(this, sequencingDecisionSituations, minQueueLength);
        }

        resetState();

        return sequencingDecisionSituations;
    }

    public List<OrderingDecisionSituation> orderingDecisionSituations(int minQueueLength) {
        List<OrderingDecisionSituation> orderingDecisionSituations = new ArrayList<>();
        //numJobsRecorded

        while (!eventQueue.isEmpty() && throughput < numBatchesRecorded) {
            AbstractEvent nextEvent = eventQueue.poll();
//            System.out.println("throughput "+throughput);
            systemState.setClockTime(nextEvent.getTime());
            nextEvent.addOrderingDecisionSituation(this, orderingDecisionSituations, minQueueLength);
        }

        resetState();

        return orderingDecisionSituations;
    }

    public List<RoutingDecisionSituation> routingDecisionSituations(int minQueueLength) {
        List<RoutingDecisionSituation> routingDecisionSituations = new ArrayList<>();

        while (!eventQueue.isEmpty() && throughput < numBatchesRecorded) {
            AbstractEvent nextEvent = eventQueue.poll();

            systemState.setClockTime(nextEvent.getTime());
            nextEvent.addRoutingDecisionSituation(this, routingDecisionSituations, minQueueLength);
        }

        resetState();

        return routingDecisionSituations;
    }

    @Override
    public Simulation surrogate(int numWorkCenters, int numJobsRecorded,
                                int warmupJobs) {
        int surrogateMaxNumOperations = maxNumOperations;

        AbstractIntegerSampler surrogateNumOperationsSampler = numOperationsSampler.clone();
        AbstractIntegerSampler surrogateNumOptionsSampler = numOperationsSampler.clone();
        AbstractRealSampler surrogateInterArrivalTimeSampler = interArrivalTimeSampler.clone();

        if (surrogateMaxNumOperations > numWorkCenters) {
            surrogateMaxNumOperations = numWorkCenters;
            surrogateNumOperationsSampler.setUpper(surrogateMaxNumOperations);

            surrogateInterArrivalTimeSampler.setMean(interArrivalTimeMean(numWorkCenters,
                    minNumOperations, surrogateMaxNumOperations, utilLevel));
        }

        Simulation surrogate = new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters,
                numJobsRecorded, warmupJobs, minNumOperations, surrogateMaxNumOperations, minBatchSize, maxBatchSize,
                utilLevel, dueDateFactor, revisit, surrogateNumOperationsSampler,
                procTimeSampler, surrogateInterArrivalTimeSampler, jobWeightSampler, numJobSampler, jobRevenueSampler);

        //modified by fzhang 17.04.2018
       /* Simulation surrogate = new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters,
                numJobsRecorded, warmupJobs, minNumOperations, surrogateMaxNumOperations,
                utilLevel, dueDateFactor, revisit, surrogateNumOperationsSampler,
                numOptionsSampler, procTimeSampler, surrogateInterArrivalTimeSampler, jobWeightSampler);*/

        return surrogate;
    }

    @Override
    public Simulation surrogateBusy(int numWorkCenters, int numJobsRecorded,
                                    int warmupJobs) {
        double utilLevel = 1;
        int surrogateMaxNumOperations = maxNumOperations;

        AbstractIntegerSampler surrogateNumOperationsSampler = numOperationsSampler.clone();
        AbstractRealSampler surrogateInterArrivalTimeSampler = interArrivalTimeSampler.clone();

        if (surrogateMaxNumOperations > numWorkCenters) {
            surrogateMaxNumOperations = numWorkCenters;
            surrogateNumOperationsSampler.setUpper(surrogateMaxNumOperations);

            surrogateInterArrivalTimeSampler.setMean(interArrivalTimeMean(numWorkCenters,
                    minNumOperations, surrogateMaxNumOperations, utilLevel));
        }

        Simulation surrogate = new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters,
                numJobsRecorded, warmupJobs, minNumOperations, surrogateMaxNumOperations, minBatchSize, maxBatchSize, utilLevel,
                dueDateFactor, revisit, surrogateNumOperationsSampler, procTimeSampler,
                surrogateInterArrivalTimeSampler, jobWeightSampler, numJobSampler, jobRevenueSampler);

        //modified by fzhang 17.04.2018
     /*   Simulation surrogate = new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters,
                numJobsRecorded, warmupJobs, minNumOperations, surrogateMaxNumOperations, utilLevel,
                dueDateFactor, revisit, surrogateNumOperationsSampler, numOptionsSampler, procTimeSampler,
                surrogateInterArrivalTimeSampler, jobWeightSampler);*/

        return surrogate;
    }

    public static DynamicSimulation standardFull(
            long seed,
            AbstractRule sequencingRule,
            AbstractRule routingRule,
            int minBatchSize,
            int maxBatchSize,
            int numWorkCenters,
            int numJobsRecorded,
            int warmupJobs,
            double utilLevel,
            double dueDateFactor) {
        return new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters, numJobsRecorded,
                warmupJobs, numWorkCenters, numWorkCenters, minBatchSize, maxBatchSize, utilLevel, dueDateFactor, false);
    }


    public static DynamicSimulation standardMissing(
            long seed,
            AbstractRule sequencingRule,
            AbstractRule routingRule,
            int minBatchSize,
            int maxBatchSize,
            int numWorkCenters,
            int numJobsRecorded,
            int warmupJobs,
            double utilLevel,
            double dueDateFactor) {
        return new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters, numJobsRecorded,
                warmupJobs, 1, numWorkCenters, minBatchSize, maxBatchSize, utilLevel, dueDateFactor, false);
    }

    public static DynamicSimulation standardMissing(
            long seed,
            AbstractRule sequencingRule,
            AbstractRule routingRule,
            int minBatchSize,
            int maxBatchSize,
            int numWorkCenters,
            int numJobsRecorded,
            int warmupJobs,
            double utilLevel,
            double dueDateFactor,
            ShopType type) {
        return new DynamicSimulation(seed, sequencingRule, routingRule, numWorkCenters, numJobsRecorded,
                warmupJobs, 1, numWorkCenters, minBatchSize, maxBatchSize, utilLevel, dueDateFactor, false, type);
    }

    public DynamicSimulation deepCloneRuntimeNoReflect() { //this will change something, such as warmupBatches, numBatchRecorded
        // 1) 创建壳子
        DynamicSimulation copy = new DynamicSimulation(
                this.seed,
                this.sequencingRule, // 若 GP 规则有状态，建议另写 deepCopy
                this.routingRule,
                this.numWorkCenters,
                batchInOneCase,
                this.numBatchesArrived-1,
                this.minNumOperations,
                this.maxNumOperations,
                this.minBatchSize,
                this.maxBatchSize,
                this.utilLevel,
                this.dueDateFactor,
                this.revisit,
                this.getNumOperationsSampler().clone(),
                this.getProcTimeSampler().clone(),
                this.getInterArrivalTimeSampler().clone(),
                this.getJobWeightSampler().clone(),
                this.getNumJobSampler().clone(),
                this.getJobRevenueSampler().clone()
        );

        // 2) 拷贝全局状态
        copy.getSystemState().setClockTime(this.getSystemState().getClockTime());
        copy.throughput = 0;
        copy.beforeThroughput = 0;
        copy.afterThroughput = 0;
        copy.count = 0;
        copy.numBatchesArrived = this.numBatchesArrived;
        copy.numJobsArrived = this.numJobsArrived;

        // 3) 映射表
        Map<WorkCenter, WorkCenter> wcMap = new IdentityHashMap<>();
        Map<Job, Job> jobMap = new IdentityHashMap<>();
        Map<Operation, Operation> opMap = new IdentityHashMap<>(); // 这是 operationMap
        Map<OperationOption, OperationOption> optMap = new IdentityHashMap<>(); // 这是 optionMap
        Map<Batch, Batch> batchMap = new IdentityHashMap<>();

        // 4) 深拷 WorkCenter
        copy.getSystemState().getWorkCenters().clear();
        for (WorkCenter srcWc : this.getSystemState().getWorkCenters()) {
            List<Double> mrtCopy = new ArrayList<>(srcWc.getMachineReadyTimes());
            LinkedList<OperationOption> newQueue = new LinkedList<>();
            WorkCenter newWc = new WorkCenter(
                    srcWc.getId(),
                    srcWc.getNumMachines(),
                    newQueue,
                    mrtCopy,
                    srcWc.getWorkInQueue(),
                    srcWc.getBusyTime()
            );
            wcMap.put(srcWc, newWc);
            copy.getSystemState().addWorkCenter(newWc);
        }

        // 5) 系统内 Job 深拷
        copy.getSystemState().getJobsInSystem().clear();
        copy.getSystemState().getBatchesInSystem().clear();
        for (Job srcJob : this.getSystemState().getJobsInSystem()) {
            Job newJob = cloneOneJobGraph(srcJob, wcMap, jobMap, opMap, optMap);
            if (srcJob.batch != null) {
                Batch newBatch = cloneOneBatch(srcJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                newJob.batch = newBatch;
            }
            copy.getSystemState().addJobToSystem(newJob);
            if(!copy.getSystemState().getBatchesInSystem().contains(newJob.batch)) {
                copy.getSystemState().addBatchToSystem(newJob.batch);
            }
        }


// 6) 预扫描事件队列，提前克隆未来对象（Job / Batch / OpOption / Process）
        for (AbstractEvent evt : this.eventQueue) {

            if (evt instanceof BatchArrivalEvent) {
                // 先确保 Batch 及批内 Job 都进入映射表
                Batch oldBatch = ((BatchArrivalEvent) evt).getBatch();
                if (!batchMap.containsKey(oldBatch)) {
                    cloneOneBatch(oldBatch, wcMap, jobMap, opMap, optMap, batchMap);
                }

            } else if (evt instanceof JobArrivalEvent) {
                Job oldJob = ((JobArrivalEvent) evt).getJob();
                if (!jobMap.containsKey(oldJob)) {
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof OperationVisitEvent) {
                OperationOption oldOp = ((OperationVisitEvent) evt).getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    // 通过克隆 job 图补齐 op/opt 映射
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof ProcessStartEvent) {
                Process oldP = ((ProcessStartEvent) evt).getProcess();
                OperationOption oldOp = oldP.getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof ProcessFinishEvent) {
                Process oldP = ((ProcessFinishEvent) evt).getProcess();
                OperationOption oldOp = oldP.getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }
            }
        }

        // 7) 恢复各 WC 队列
        for (WorkCenter srcWc : this.getSystemState().getWorkCenters()) {
            WorkCenter newWc = wcMap.get(srcWc);
            for (OperationOption srcOpt : srcWc.getQueue()) {
                OperationOption newOpt = optMap.get(srcOpt);
                newWc.getQueue().add(newOpt);
            }
        }

        // 8) 克隆事件队列
        copy.eventQueue.clear();
        for (AbstractEvent evt : this.eventQueue) {
            AbstractEvent clonedEvt = cloneEvent(evt, jobMap, opMap,optMap, wcMap, batchMap);
            if (clonedEvt != null) {
                copy.eventQueue.add(clonedEvt);
            }
        }

        // 9) 随机源
        copy.randomDataGenerator = deepCopy(randomDataGenerator);

        return copy;
    }

    public DynamicSimulation deepCloneAllSame() { //this will not change anything
        // 1) 创建壳子
        DynamicSimulation copy = new DynamicSimulation(
                this.seed,
                this.sequencingRule, // 若 GP 规则有状态，建议另写 deepCopy
                this.routingRule,
                this.numWorkCenters,
                numBatchesRecorded,
                warmupBatches,
                this.minNumOperations,
                this.maxNumOperations,
                this.minBatchSize,
                this.maxBatchSize,
                this.utilLevel,
                this.dueDateFactor,
                this.revisit,
                this.getNumOperationsSampler().clone(),
                this.getProcTimeSampler().clone(),
                this.getInterArrivalTimeSampler().clone(),
                this.getJobWeightSampler().clone(),
                this.getNumJobSampler().clone(),
                this.getJobRevenueSampler().clone()
        );

        // 2) 拷贝全局状态
        copy.getSystemState().setClockTime(this.getSystemState().getClockTime());
        copy.throughput = 0;
        copy.beforeThroughput = 0;
        copy.afterThroughput = 0;
        copy.count = 0;
        copy.numBatchesArrived = this.numBatchesArrived;
        copy.numJobsArrived = this.numJobsArrived;

        // 3) 映射表
        Map<WorkCenter, WorkCenter> wcMap = new IdentityHashMap<>();
        Map<Job, Job> jobMap = new IdentityHashMap<>();
        Map<Operation, Operation> opMap = new IdentityHashMap<>(); // 这是 operationMap
        Map<OperationOption, OperationOption> optMap = new IdentityHashMap<>(); // 这是 optionMap
        Map<Batch, Batch> batchMap = new IdentityHashMap<>();

        // 4) 深拷 WorkCenter
        copy.getSystemState().getWorkCenters().clear();
        for (WorkCenter srcWc : this.getSystemState().getWorkCenters()) {
            List<Double> mrtCopy = new ArrayList<>(srcWc.getMachineReadyTimes());
            LinkedList<OperationOption> newQueue = new LinkedList<>();
            WorkCenter newWc = new WorkCenter(
                    srcWc.getId(),
                    srcWc.getNumMachines(),
                    newQueue,
                    mrtCopy,
                    srcWc.getWorkInQueue(),
                    srcWc.getBusyTime()
            );
            wcMap.put(srcWc, newWc);
            copy.getSystemState().addWorkCenter(newWc);
        }

        // 5) 系统内 Job 深拷
        copy.getSystemState().getJobsInSystem().clear();
        copy.getSystemState().getBatchesInSystem().clear();
        for (Job srcJob : this.getSystemState().getJobsInSystem()) {
            Job newJob = cloneOneJobGraph(srcJob, wcMap, jobMap, opMap, optMap);
            if (srcJob.batch != null) {
                Batch newBatch = cloneOneBatch(srcJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                newJob.batch = newBatch;
            }
            copy.getSystemState().addJobToSystem(newJob);
            if(!copy.getSystemState().getBatchesInSystem().contains(newJob.batch)) {
                copy.getSystemState().addBatchToSystem(newJob.batch);
            }
        }

// 6) 预扫描事件队列，提前克隆未来对象（Job / Batch / OpOption / Process）
        for (AbstractEvent evt : this.eventQueue) {

            if (evt instanceof BatchArrivalEvent) {
                // 先确保 Batch 及批内 Job 都进入映射表
                Batch oldBatch = ((BatchArrivalEvent) evt).getBatch();
                if (!batchMap.containsKey(oldBatch)) {
                    cloneOneBatch(oldBatch, wcMap, jobMap, opMap, optMap, batchMap);
                }

            } else if (evt instanceof JobArrivalEvent) {
                Job oldJob = ((JobArrivalEvent) evt).getJob();
                if (!jobMap.containsKey(oldJob)) {
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof OperationVisitEvent) {
                OperationOption oldOp = ((OperationVisitEvent) evt).getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    // 通过克隆 job 图补齐 op/opt 映射
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof ProcessStartEvent) {
                Process oldP = ((ProcessStartEvent) evt).getProcess();
                OperationOption oldOp = oldP.getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }

            } else if (evt instanceof ProcessFinishEvent) {
                Process oldP = ((ProcessFinishEvent) evt).getProcess();
                OperationOption oldOp = oldP.getOperationOption();
                if (!optMap.containsKey(oldOp)) {
                    Job oldJob = oldOp.getJob();
                    Job newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
                    if (oldJob.batch != null && newJob.batch == null) {
                        Batch newBatch = batchMap.containsKey(oldJob.batch)
                                ? batchMap.get(oldJob.batch)
                                : cloneOneBatch(oldJob.batch, wcMap, jobMap, opMap, optMap, batchMap);
                        newJob.batch = newBatch;
                    }
                }
            }
        }

        // 7) 恢复各 WC 队列
        for (WorkCenter srcWc : this.getSystemState().getWorkCenters()) {
            WorkCenter newWc = wcMap.get(srcWc);
            for (OperationOption srcOpt : srcWc.getQueue()) {
                OperationOption newOpt = optMap.get(srcOpt);
                newWc.getQueue().add(newOpt);
            }
        }

        // 8) 克隆事件队列
        copy.eventQueue.clear();
        for (AbstractEvent evt : this.eventQueue) {
            AbstractEvent clonedEvt = cloneEvent(evt, jobMap, opMap,optMap, wcMap, batchMap);
            if (clonedEvt != null) {
                copy.eventQueue.add(clonedEvt);
            }
        }

        // 9) 随机源
        copy.randomDataGenerator = deepCopy(randomDataGenerator);

        return copy;
    }

    private AbstractEvent cloneEvent(AbstractEvent evt,
                                     Map<Job, Job> jobMap,
                                     Map<Operation, Operation> opMap,
                                     Map<OperationOption, OperationOption> optMap,
                                     Map<WorkCenter, WorkCenter> wcMap,
                                     Map<Batch, Batch> batchMap) {

        if (evt instanceof JobArrivalEvent) {
            Job oldJob = ((JobArrivalEvent) evt).getJob();
            Job newJob = jobMap.get(oldJob);
            if (oldJob.batch != null && newJob.batch == null) {
                newJob.batch = batchMap.get(oldJob.batch);
            }
            return new JobArrivalEvent(evt.getTime(), newJob);

        } else if (evt instanceof BatchArrivalEvent) {
            Batch oldBatch = ((BatchArrivalEvent) evt).getBatch();
            Batch newBatch = batchMap.get(oldBatch);
            if (newBatch == null) {
                newBatch = cloneOneBatch(oldBatch, wcMap, jobMap, opMap, optMap, batchMap);
            }
            return new BatchArrivalEvent(evt.getTime(), newBatch);

        } else if (evt instanceof OperationVisitEvent) {
            OperationOption oldOp = ((OperationVisitEvent) evt).getOperationOption();
            OperationOption newOp = optMap.get(oldOp);
            if (newOp.getJob().batch == null && oldOp.getJob().batch != null) {
                newOp.getJob().batch = batchMap.get(oldOp.getJob().batch);
            }
            return new OperationVisitEvent(evt.getTime(), newOp);

        } else if (evt instanceof ProcessStartEvent) {
            Process oldP = ((ProcessStartEvent) evt).getProcess();
            WorkCenter newWC = wcMap.get(oldP.getWorkCenter());
            OperationOption newOp = optMap.get(oldP.getOperationOption());
            if (newOp.getJob().batch == null && oldP.getOperationOption().getJob().batch != null) {
                newOp.getJob().batch = batchMap.get(oldP.getOperationOption().getJob().batch);
            }
            Process newP = new Process(newWC, oldP.getMachineId(), newOp, oldP.getStartTime());
            return new ProcessStartEvent(evt.getTime(), newP);

        } else if (evt instanceof ProcessFinishEvent) {
            Process oldP = ((ProcessFinishEvent) evt).getProcess();
            WorkCenter newWC = wcMap.get(oldP.getWorkCenter());
            OperationOption newOp = optMap.get(oldP.getOperationOption());
            if (newOp.getJob().batch == null && oldP.getOperationOption().getJob().batch != null) {
                newOp.getJob().batch = batchMap.get(oldP.getOperationOption().getJob().batch);
            }
            Process newP = new Process(newWC, oldP.getMachineId(), newOp, oldP.getStartTime());
            return new ProcessFinishEvent(evt.getTime(), newP);
        }

        throw new IllegalArgumentException("Unsupported event type: " + evt.getClass());
    }



    private Batch cloneOneBatch(Batch oldBatch,
                                Map<WorkCenter, WorkCenter> wcMap,
                                Map<Job, Job> jobMap,
                                Map<Operation, Operation> opMap,
                                Map<OperationOption, OperationOption> optMap,
                                Map<Batch, Batch> batchMap) {
        Batch existed = batchMap.get(oldBatch);
        if (existed != null) return existed;

        Batch newBatch = new Batch(
                oldBatch.getId(),
                oldBatch.getNumOfJobs(),
                oldBatch.getTotalProcTime(),
                oldBatch.getArrivalTime(),
                oldBatch.getReleaseTime(),
                oldBatch.getDueDate(),
                oldBatch.getWeight()
        );

        for (Job oldJob : oldBatch.getJobs()) {
            Job newJob = jobMap.get(oldJob);
            if (newJob == null) {
                newJob = cloneOneJobGraph(oldJob, wcMap, jobMap, opMap, optMap);
            }
            newBatch.addJob(newJob);
            newJob.batch = newBatch; // 保证 batch 引用
        }

        batchMap.put(oldBatch, newBatch);
        return newBatch;
    }

    private Job cloneOneJobGraph(Job srcJob,
                                 Map<WorkCenter, WorkCenter> wcMap,
                                 Map<Job, Job> jobMap,
                                 Map<Operation, Operation> opMap,
                                 Map<OperationOption, OperationOption> optMap) {
        Job existed = jobMap.get(srcJob);
        if (existed != null) return existed;

        Job newJob = new Job(
                srcJob.getId(),
                new ArrayList<>(),
                srcJob.getArrivalTime(),
                srcJob.getReleaseTime(),
                srcJob.getDueDate(),
                srcJob.getWeight(),
                srcJob.getRevenue()
        );
        newJob.setCompletionTime(srcJob.getCompletionTime());
        newJob.setPriority(srcJob.getPriority());
        newJob.workLoadRemaining = srcJob.workLoadRemaining;
        newJob.operationRemaining = srcJob.operationRemaining;

        jobMap.put(srcJob, newJob);

        // Operation
        for (Operation srcOp : srcJob.getOperations()) {
            Operation newOp = new Operation(newJob, srcOp.getId());
            opMap.put(srcOp, newOp);
            newJob.addOperation(newOp);
        }

        // Option
        for (Operation srcOp : srcJob.getOperations()) {
            Operation newOp = opMap.get(srcOp);
            for (OperationOption srcOpt : srcOp.getOperationOptions()) {
                WorkCenter newWc = wcMap.get(srcOpt.getWorkCenter());
                OperationOption newOpt = new OperationOption(
                        newOp,
                        srcOpt.getOptionId(),
                        srcOpt.getProcTime(),
                        newWc
                );
                newOp.addOperationOption(newOpt);
                optMap.put(srcOpt, newOpt);
            }
        }

        newJob.linkOperations();
        return newJob;
    }

    @SuppressWarnings("unchecked")
    public static <T extends Serializable> T deepCopy(T obj) {
        try (ByteArrayOutputStream bos = new ByteArrayOutputStream();
             ObjectOutputStream out = new ObjectOutputStream(bos)) {
            out.writeObject(obj);
            out.flush();
            try (ObjectInputStream in =
                         new ObjectInputStream(new ByteArrayInputStream(bos.toByteArray()))) {
                return (T) in.readObject();
            }
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException("Deep copy failed", e);
        }
    }

}
