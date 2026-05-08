package yimei.jss.niching;

import ec.gp.GPIndividual;
import smile.clustering.KMeans;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.FlexibleStaticInstance;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.LPT;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.SequencingDecisionSituation;
import yimei.jss.simulation.StaticSimulation;

import java.util.*;
import java.util.stream.IntStream;

public class SequencingPhenoCharacterisation extends PhenoCharacterisation {
    public List<SequencingDecisionSituation> decisionSituations;
    private int[][] referenceIndexes;

    public SequencingPhenoCharacterisation(AbstractRule sequencingReferenceRule,
                                           List<SequencingDecisionSituation> decisionSituations) {
        super(sequencingReferenceRule);
        this.decisionSituations = decisionSituations;
        this.referenceIndexes = new int[decisionSituations.size()][];

        calcReferenceIndexes();
    }

    public int[] characterise(AbstractRule rule) {
        int[] charList = new int[decisionSituations.size()];

        for (int i = 0; i < decisionSituations.size(); i++) {
            SequencingDecisionSituation situation = decisionSituations.get(i);
            List<OperationOption> queue = situation.getQueue();


            // Calculate the priority for all the operations.
            for (OperationOption op : queue) {
                op.setPriority(rule.priority(
                        op, situation.getWorkCenter(), situation.getSystemState()));
            }

            //get the index of best operation
            int idxBestOp = 0;
            for (int j = 0; j < queue.size(); j++) {
                if (queue.get(j).priorTo(queue.get(idxBestOp))) {
                    idxBestOp = j;
                }
            }
            charList[i] = referenceIndexes[i][idxBestOp];
        }

        return charList;
    }

    public int[][] characteriseRanks(AbstractRule rule) {
        int[][] charList = new int[decisionSituations.size()][];

        for (int i = 0; i < decisionSituations.size(); i++) {
            SequencingDecisionSituation situation = decisionSituations.get(i);
            charList[i] = rule.priorValueOperation(situation);
        }

        return charList;
    }

    protected void calcReferenceIndexes() {
        for (int i = 0; i < decisionSituations.size(); i++) {
            SequencingDecisionSituation situation = decisionSituations.get(i);
            int[] ranks = referenceRule.priorValueOperation(situation);

            //int index = situation.getQueue().indexOf(op);
            referenceIndexes[i] = ranks;
        }
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation() {

        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();

        if (DynamicSimulation.optimisationObjective.contains("batch")) {
            defaultSequencingRule = new LPT(RuleType.SEQUENCING);
        }

        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);


        //fzhang 2019.6.22 original
        //int minQueueLength = 8;
        //int numDecisionSituations = 20;

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
        int minQueueLength = 6;
        int numDecisionSituations = 20;//used for measuring the behavior of different rules
        long shuffleSeed = 8295342;

        /*DynamicSimulation simulation = DynamicSimulation.standardFull(0, defaultSequencingRule,
                defaultRoutingRule, 10, 500, 0,
                0.95, 4.0); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state
*/
        DynamicSimulation simulation = DynamicSimulation.standardFull(10, defaultSequencingRule,
                defaultRoutingRule, 10,20,5, 100, 10,
                0.85, 1.2); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state

        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        List<SequencingDecisionSituation> situations = simulation.sequencingDecisionSituations(minQueueLength); //situations have 20 elements

        //after get situations,we need to select them based on ranks[], if same, delete.
/*        int[][] referenceIndexes = new int[situations.size()][];
        for (int i = 0; i < situations.size(); i++) {
            SequencingDecisionSituation situation = situations.get(i);
            int[] ranks = defaultSequencingRule.priorValueOperation(situation);
            referenceIndexes[i] = ranks;
        }

        double[][] data = new double[referenceIndexes.length][referenceIndexes[0].length];
        for (int i = 0; i < referenceIndexes.length; i++) {
            for (int j = 0; j < referenceIndexes[0].length; j++) {
                data[i][j] = referenceIndexes[i][j];
            }
        }

        KMeans kMeans = KMeans.lloyd(data, numDecisionSituations);

        int[] centerIndices = getClusterCenterIndices(kMeans,data);

        situations = selectByIndices(situations, centerIndices);*/

        Collections.shuffle(situations, new Random(shuffleSeed)); //Randomly permute the specified list using the specified source of randomness.
        //randomly change the sorting of list of situations

        situations = situations.subList(0, numDecisionSituations); //Returns a view of the portion of this list between the specified fromIndex,
        //inclusive, and toIndex, exclusive. (If fromIndex and toIndex are equal, the returned list is empty.)
        return new SequencingPhenoCharacterisation(defaultSequencingRule, situations);
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation(String filePath) {
        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING);
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);
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
        List<SequencingDecisionSituation> situations = simulation.sequencingDecisionSituations(minQueueLength);
        while (situations.size() < numDecisionSituations && minQueueLength > 2) {
            minQueueLength--;
            situations = simulation.sequencingDecisionSituations(minQueueLength);
        }

        if (minQueueLength == 2 && situations.size() < 100) {
            //no point going to queue length of 1, as this will only have 1 outcome
            System.out.println("Only " + situations.size() + " instances available for sequencing pheno characterisation.");
            numDecisionSituations = situations.size();
        }

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
        return new SequencingPhenoCharacterisation(defaultSequencingRule, situations);
    }

    public List<SequencingDecisionSituation> getDecisionSituations() {
        return decisionSituations;
    }

    public int[][] getReferenceIndexes() {
        return referenceIndexes;
    }


    // 获取每个簇中距离质心最近的原始样本索引
    public static int[] getClusterCenterIndices(KMeans kmeans, double[][] data) {
        int k = kmeans.k;
        double[][] centroids = kmeans.centroids;
        int[] labels = kmeans.y;
        int[] centerIndices = new int[k];

        for (int clusterId = 0; clusterId < k; clusterId++) {
            double minDist = Double.MAX_VALUE;
            int bestIndex = -1;

            for (int i = 0; i < data.length; i++) {
                if (labels[i] == clusterId) {
                    double dist = euclideanDistance(centroids[clusterId], data[i]);
                    if (dist < minDist) {
                        minDist = dist;
                        bestIndex = i;
                    }
                }
            }

            centerIndices[clusterId] = bestIndex;
        }

        return centerIndices;
    }

    // 欧几里得距离计算
    private static double euclideanDistance(double[] a, double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    public static <T> List<T> selectByIndices(List<T> list, int[] indices) {
        List<T> result = new ArrayList<>();
        for (int index : indices) {
            result.add(list.get(index));
        }
        return result;
    }

    public static PhenoCharacterisation currentGenPhenoCharacterisation(GPRuleEvolutionState state, RuleOptimizationProblem problem) {

        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);

//        AbstractRule defaultSequencingRule = new GPRule(RuleType.SEQUENCING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[0]);
//        AbstractRule defaultRoutingRule = new GPRule(RuleType.ROUTING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[1]);


        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel)problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
        int minQueueLength = simulation.getNumWorkCenters();
        int numDecisionSituations = 100;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        List<SequencingDecisionSituation> situations = simulation.sequencingDecisionSituations(minQueueLength); //situations have 20 elements

        Collections.shuffle(situations, new Random(shuffleSeed)); //Randomly permute the specified list using the specified source of randomness.
        //randomly change the sorting of list of situations

//        situations = selectSituations(situations,state);

        situations = situations.subList(0, numDecisionSituations); //Returns a view of the portion of this list between the specified fromIndex,
        //inclusive, and toIndex, exclusive. (If fromIndex and toIndex are equal, the returned list is empty.)

        simulation.reset();

        return new SequencingPhenoCharacterisation(defaultSequencingRule, situations);
    }

    private static List<SequencingDecisionSituation> selectSituations(List<SequencingDecisionSituation> situations, GPRuleEvolutionState state) {

        //1. randomly select 200 cases to save resources
        int selectedCases = Math.min(500, situations.size());
        situations = situations.subList(0, selectedCases);

        //2. select 10 individuals from 1%, 10%, 20%, 30%, ..., 90%
        ArrayList<GPIndividual> selectedIndividuals = new ArrayList<>();

        // 百分位列表（注意 1% 的索引为 0.01 * total，向上取整）
//        int total = state.population.subpops[0].individuals.length;
//        int[] percentiles = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
//
//        for (int p : percentiles) {
//            int index = (int) Math.round(p / 100.0 * (total - 1));
//            index = Math.max(0, Math.min(index, total - 1)); // 保证不越界
//            GPIndividual ind = (GPIndividual) state.population.subpops[0].individuals[index];
//            selectedIndividuals.add(ind);
//        }

        for (int index=0;index<state.population.subpops[0].individuals.length;index++) {
            GPIndividual ind = (GPIndividual) state.population.subpops[0].individuals[index];
            selectedIndividuals.add(ind);
        }

        //3. get the decision behaviour of these 10 individuals
        double[] entropies = new double[selectedCases];

        for (int i = 0; i < situations.size(); i++) {
            SequencingDecisionSituation situation = situations.get(i);
            int [][] behaviourVector = new int[selectedIndividuals.size()][];
            for (int j = 0; j < selectedIndividuals.size(); j++) {
                GPIndividual selectedIndividual = selectedIndividuals.get(j);
                GPRule sequencingRule = new GPRule(RuleType.SEQUENCING,selectedIndividual.trees[0]);
                behaviourVector[j] = sequencingRule.priorValueOperation(situation);
            }
            entropies[i] = PopulationUtils.entropy(behaviourVector);
        }

//        KneePointFinder.KneeResult result = KneePointFinder.findKneePoint(entropies,0);
//
//        List<SequencingDecisionSituation> finalSituations = new ArrayList<>();
//        for (int i = 0; i < entropies.length; i++) {
//            if (entropies[i] < result.kneeValue) {
//               finalSituations.add(situations.get(i));
//            }
//        }

        //4. select final 20 cases based on the entropy of decision situations
        situations = selectTopKSituations(entropies,situations,100);

        return situations;
    }

    public static int[] getTopKIndices(double[] values, int k) {
        Integer[] indices = IntStream.range(0, values.length).boxed().toArray(Integer[]::new);
        Arrays.sort(indices, (i, j) -> Double.compare(values[j], values[i]));
        return Arrays.stream(indices).limit(k).mapToInt(i -> i).toArray();
    }

    public static List<SequencingDecisionSituation> selectTopKSituations(
            double[] entropies, List<SequencingDecisionSituation> situations, int k) {

        int[] topIndices = getTopKIndices(entropies, k);
        List<SequencingDecisionSituation> selected = new ArrayList<>();

        for (int idx : topIndices) {
            selected.add(situations.get(idx));
        }

        return selected;
    }

    public static PhenoCharacterisation currentBestRulePhenoCharacterisation(GPRuleEvolutionState state, RuleOptimizationProblem problem) {

//        AbstractRule defaultSequencingRule = new SPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
//        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);

        AbstractRule defaultSequencingRule = new GPRule(RuleType.SEQUENCING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[0]);
        AbstractRule defaultRoutingRule = new GPRule(RuleType.ROUTING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[1]);

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
        int minQueueLength = 2;
        int numDecisionSituations = 20;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        DynamicSimulation simulationCurrent = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel)problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);


        DynamicSimulation simulation = DynamicSimulation.standardMissing(0, defaultSequencingRule,
                defaultRoutingRule, simulationCurrent.minBatchSize ,simulationCurrent.maxBatchSize,simulationCurrent.getNumWorkCenters(),simulationCurrent.getBatchesRecorded() , simulationCurrent.getWarmupBatches(),
                simulationCurrent.getUtilLevel(), 1.2); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state


        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        List<SequencingDecisionSituation> situations = simulation.sequencingDecisionSituations(minQueueLength); //situations have 20 elements

        Collections.shuffle(situations, new Random(shuffleSeed)); //Randomly permute the specified list using the specified source of randomness.
        //randomly change the sorting of list of situations

        HashMap hashMap = new HashMap();

        for (SequencingDecisionSituation situation: situations){
            hashMap.put(((SequencingDecisionSituation)situation).getQueue().size(),situation);
        }

        situations = selectSituations(situations,state);

//        situations = situations.subList(0, numDecisionSituations); //Returns a view of the portion of this list between the specified fromIndex,
        //inclusive, and toIndex, exclusive. (If fromIndex and toIndex are equal, the returned list is empty.)

        simulation.reset();

        return new SequencingPhenoCharacterisation(defaultSequencingRule, situations);
    }

    public static PhenoCharacterisation currentTaskPhenoCharacterisation(DynamicSimulation simulation) {

        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);

//        AbstractRule defaultSequencingRule = new GPRule(RuleType.SEQUENCING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[0]);
//        AbstractRule defaultRoutingRule = new GPRule(RuleType.ROUTING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[1]);

        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
//        int minQueueLength = simulation.getNumWorkCenters();
        int minQueueLength = 2;
        int numDecisionSituations = 100;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        List<SequencingDecisionSituation> situations = simulation.sequencingDecisionSituations(minQueueLength); //situations have 20 elements

        Collections.shuffle(situations, new Random(shuffleSeed)); //Randomly permute the specified list using the specified source of randomness.
        //randomly change the sorting of list of situations

//        situations = selectSituations(situations,state);

        situations = situations.subList(0, numDecisionSituations); //Returns a view of the portion of this list between the specified fromIndex,
        //inclusive, and toIndex, exclusive. (If fromIndex and toIndex are equal, the returned list is empty.)

        simulation.reset();

        return new SequencingPhenoCharacterisation(defaultSequencingRule, situations);
    }

}