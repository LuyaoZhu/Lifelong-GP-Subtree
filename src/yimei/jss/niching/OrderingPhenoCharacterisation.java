package yimei.jss.niching;

import yimei.jss.jobshop.FlexibleStaticInstance;
import yimei.jss.jobshop.Job;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.SPT;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.OrderingDecisionSituation;
import yimei.jss.simulation.StaticSimulation;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class OrderingPhenoCharacterisation extends PhenoCharacterisation {
    public List<OrderingDecisionSituation> decisionSituations;
    private int[] referenceIndexes;

    public OrderingPhenoCharacterisation(AbstractRule orderingReferenceRule,
                                         List<OrderingDecisionSituation> decisionSituations) {
        super(orderingReferenceRule);
        this.decisionSituations = decisionSituations;
        this.referenceIndexes = new int[decisionSituations.size()];

        calcReferenceIndexes();
    }

    public int[] characterise(AbstractRule rule) {
        int[] charList = new int[decisionSituations.size()];

        for (int i = 0; i < decisionSituations.size(); i++) {
            OrderingDecisionSituation situation = decisionSituations.get(i);
            List<Job> queue = situation.getJobQueue();

            int refIdx = referenceIndexes[i]; //[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  20 values---20 situations

            // Calculate the priority for all the operations.
            for (Job job : queue) {
                job.setPriority(rule.priorityJob(
                        job,situation.getSystemState(),rule));
            }

            // get the rank of the processing chosen by the reference rule.
            //the operation chosen by reference rule, rank to where by examined rules
            int rank = 1;
            for (int j = 0; j < queue.size(); j++) {
                if (queue.get(j).priorTo(queue.get(refIdx))) {
                    rank ++;
                }
            }

            charList[i] = rank;
        }

        return charList;
    }

    public int[][] characteriseRanks(AbstractRule rule) {
        int[][] charList = new int[decisionSituations.size()][];

        for (int i = 0; i < decisionSituations.size(); i++) {
            OrderingDecisionSituation situation = decisionSituations.get(i);

        }

        return charList;
    }

    protected void calcReferenceIndexes() {
        for (int i = 0; i < decisionSituations.size(); i++) {
            OrderingDecisionSituation situation = decisionSituations.get(i);

            int startID = situation.getJobQueue().get(0).getId();

            for (int j = 0; j < situation.getJobQueue().size(); j++) {
                Job job = situation.getJobQueue().get(j);
//            job.setPriority(job.getId()); //Manual rule: FIFS
                job.setPriority(job.getOperation(0).getMedianProcessTime()); //Manual rule: SPT
            }

            //learned ordering rule
            situation.getJobQueue().sort(Comparator.comparingDouble(job -> job.getPriority()));

            int index = situation.getJobQueue().get(0).getId() - startID + 1;
            referenceIndexes[i] = index;
        }
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation() {
        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        //the larger the weight, the smaller the WSPT value
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);
        AbstractRule defaultOrderingRule = new SPT(RuleType.ORDERING);

        //fzhang 2019.6.22 original
        //int minQueueLength = 8;
        //int numDecisionSituations = 20;

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
        int minQueueLength = 10;
        int numDecisionSituations = 100;//used for measuring the behavior of different rules
        long shuffleSeed = 8295342;

        /*DynamicSimulation simulation = DynamicSimulation.standardFull(0, defaultSequencingRule,
                defaultRoutingRule, 10, 500, 0,
                0.95, 4.0); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state
*/
        DynamicSimulation simulation = DynamicSimulation.standardFull(10, defaultSequencingRule,
                defaultRoutingRule, 10,25,10, 150, 0,
                0.95, 1.2); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state

        simulation.setOrderingRule(defaultOrderingRule);

        List<OrderingDecisionSituation> situations = simulation.orderingDecisionSituations(minQueueLength); //situations have 20 elements
        Collections.shuffle(situations, new Random(shuffleSeed)); //Randomly permute the specified list using the specified source of randomness.
        //randomly change the sorting of list of situations

        situations = situations.subList(0, numDecisionSituations); //Returns a view of the portion of this list between the specified fromIndex,
        //inclusive, and toIndex, exclusive. (If fromIndex and toIndex are equal, the returned list is empty.)
        return new OrderingPhenoCharacterisation(defaultOrderingRule, situations);
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation(String filePath) {
        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING);
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);
        AbstractRule defaultOrderingRule = new WIQ(RuleType.ORDERING);
        FlexibleStaticInstance flexibleStaticInstance = FlexibleStaticInstance.readFromAbsPath(filePath);
        StaticSimulation simulation = new StaticSimulation(defaultSequencingRule, defaultRoutingRule,
                flexibleStaticInstance);

        //int minQueueLength = 8; // some flexible static instances will have short queues
        //int numDecisionSituations = 20;

        int minQueueLength = 7; // some flexible static instances will have short queues
        int numDecisionSituations = 100;
        long shuffleSeed = 8295342;

        //the number of sequencing decisions available of a given queue length will vary
        //greatly for different statics instances, so we'll start with 8 and decrease
        //the min queue length until we can get at least 20
        List<OrderingDecisionSituation> situations = simulation.orderingDecisionSituations(minQueueLength);
        while (situations.size() < numDecisionSituations && minQueueLength > 2) {
            minQueueLength--;
            situations = simulation.orderingDecisionSituations(minQueueLength);
        }

        if (minQueueLength == 2 && situations.size() < 100) {
            //no point going to queue length of 1, as this will only have 1 outcome
            System.out.println("Only "+situations.size() +" instances available for sequencing pheno characterisation.");
            numDecisionSituations = situations.size();
        }

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
        return new OrderingPhenoCharacterisation(defaultOrderingRule, situations);
    }

    public List<OrderingDecisionSituation> getDecisionSituations() {
        return decisionSituations;
    }

    public int[] getReferenceIndexes() {
        return referenceIndexes;
    }
}