package yimei.jss.niching;

import ec.gp.GPIndividual;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.FlexibleStaticInstance;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.SPT;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.StaticSimulation;

import java.util.*;
import java.util.stream.IntStream;

public class RoutingPhenoCharacterisation extends PhenoCharacterisation {
    //private List<RoutingDecisionSituation> decisionSituations;
    public List<RoutingDecisionSituation> decisionSituations;
    private int[][] referenceIndexes;

    public RoutingPhenoCharacterisation(AbstractRule routingReferenceRule,
                                        List<RoutingDecisionSituation> decisionSituations) {
        super(routingReferenceRule);
        this.decisionSituations = decisionSituations;
        this.referenceIndexes = new int[decisionSituations.size()][];

        calcReferenceIndexes();
    }

    public int[] characterise(AbstractRule rule) {
        int[] charList = new int[decisionSituations.size()];

        for (int i = 0; i < decisionSituations.size(); i++) {
            //this is for routing rule
            RoutingDecisionSituation situation = decisionSituations.get(i);
            List<OperationOption> queue = situation.getQueue();

            //int refIdx = referenceIndexes[i];

            // Calculate the priority for all the operations.
            for (OperationOption op : queue) {
                op.setPriority(rule.priority(
                        op, op.getWorkCenter(), situation.getSystemState()));
            }
            // get the rank of the processing chosen by the reference rule.
            int idxBestOption = 0;
            for (int j = 0; j < queue.size(); j++) {
                if (queue.get(j).priorTo(queue.get(idxBestOption))) {
                    idxBestOption = j;
                }
            }
            charList[i] = referenceIndexes[i][idxBestOption];
        }
        return charList;
    }

    public int[][] characteriseRanks(AbstractRule rule) {
        int[][] charList = new int[decisionSituations.size()][];

        for (int i = 0; i < decisionSituations.size(); i++) {
            //this is for routing rule
            RoutingDecisionSituation situation = decisionSituations.get(i);

            charList[i] = rule.priorValueOperationMahcine(situation);
        }
        return charList;
    }

    protected void calcReferenceIndexes() {
        for (int i = 0; i < decisionSituations.size(); i++) {
            RoutingDecisionSituation situation = decisionSituations.get(i);
            int[] ranks = referenceRule.priorValueOperationMahcine(situation);
            referenceIndexes[i] = ranks;
        }
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation() {
        AbstractRule defaultSequencingRule = new SPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);

        //fzhang 2019.6.22 original code
        //int minQueueLength = 8; //original setting

        //fzhang 2019.6.22
        int minQueueLength = 6;
        int numDecisionSituations = 20;//used for measuring the behavior of different rules
        long shuffleSeed = 8295342;

        /*DynamicSimulation simulation = DynamicSimulation.standardFull(0, defaultSequencingRule,
                defaultRoutingRule, 10, 500, 0,
                0.95, 4.0); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state
*/
        DynamicSimulation simulation = DynamicSimulation.standardFull(10, defaultSequencingRule,
                defaultRoutingRule, 10,20,6, 100, 10,
                0.85, 1.2); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state


        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        List<RoutingDecisionSituation> situations = simulation.routingDecisionSituations(minQueueLength);

        //after get situations,we need to select them based on ranks[], if same, delete.
/*        int[][] referenceIndexes = new int[situations.size()][];
        for (int i = 0; i < situations.size(); i++) {
            RoutingDecisionSituation situation = situations.get(i);
            int[] ranks = defaultRoutingRule.priorValueOperationMahcine(situation);
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

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
        return new RoutingPhenoCharacterisation(defaultRoutingRule, situations);
    }

    public static PhenoCharacterisation defaultPhenoCharacterisation(String filePath) {
        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING);
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);
        FlexibleStaticInstance flexibleStaticInstance = FlexibleStaticInstance.readFromAbsPath(filePath);
        StaticSimulation simulation = new StaticSimulation(defaultSequencingRule, defaultRoutingRule,
                flexibleStaticInstance);

        int minQueueLength = 8;
        int numDecisionSituations = 100;
        long shuffleSeed = 8295342;

        List<RoutingDecisionSituation> situations = simulation.routingDecisionSituations(minQueueLength);
        while (situations.size() < numDecisionSituations && minQueueLength > 2) {
            minQueueLength--;
            situations = simulation.routingDecisionSituations(minQueueLength);
        }

        if (minQueueLength == 2 && situations.size() < numDecisionSituations) {
            //no point going to queue length of 1, as this will only have 1 outcome
            System.out.println("Only "+situations.size() +" instances available for routing pheno characterisation.");
            numDecisionSituations = situations.size();
        }

        Collections.shuffle(situations, new Random(shuffleSeed));


        situations = situations.subList(0, numDecisionSituations);
        return new RoutingPhenoCharacterisation(defaultRoutingRule, situations);
    }

    public List<RoutingDecisionSituation> getDecisionSituations() {
        return decisionSituations;
    }

    public int[][] getReferenceIndexes() {
        return referenceIndexes;
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
        int numDecisionSituations = 20;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        List<RoutingDecisionSituation> situations = simulation.routingDecisionSituations(minQueueLength);

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
//        situations = selectSituations(situations, state);

        simulation.reset();
//
        return new RoutingPhenoCharacterisation(defaultRoutingRule, situations);

    }

    private static List<RoutingDecisionSituation> selectSituations(List<RoutingDecisionSituation> situations, GPRuleEvolutionState state) {

        //1. randomly select 200 cases to save resources
        int selectedCases = Math.min(500, situations.size());
        situations = situations.subList(0, selectedCases);

        //2. select 10 individuals from 1%, 10%, 20%, 30%, ..., 90%, 100%
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
            RoutingDecisionSituation situation = situations.get(i);
            int [][] behaviourVector = new int[selectedIndividuals.size()][];
            for (int j = 0; j < selectedIndividuals.size(); j++) {
                GPIndividual selectedIndividual = selectedIndividuals.get(j);
                GPRule routingRule = new GPRule(RuleType.ROUTING,selectedIndividual.trees[1]);
                behaviourVector[j] = routingRule.priorValueOperationMahcine(situation);
            }
            entropies[i] = PopulationUtils.entropy(behaviourVector);
        }

        //4. select final 20 cases based on the entropy of decision situations
        situations = selectTopKSituations(entropies,situations,100);

        return situations;
    }

    public static int[] getTopKIndices(double[] values, int k) {
        Integer[] indices = IntStream.range(0, values.length).boxed().toArray(Integer[]::new);
        Arrays.sort(indices, (i, j) -> Double.compare(values[j], values[i]));
        return Arrays.stream(indices).limit(k).mapToInt(i -> i).toArray();
    }

    public static List<RoutingDecisionSituation> selectTopKSituations(
            double[] entropies, List<RoutingDecisionSituation> situations, int k) {

        int[] topIndices = getTopKIndices(entropies, k);
        List<RoutingDecisionSituation> selected = new ArrayList<>();

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
        int minQueueLength = 5;
        int numDecisionSituations = 20;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        DynamicSimulation simulationCurrent = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel)problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);

        DynamicSimulation simulation = DynamicSimulation.standardMissing(0, defaultSequencingRule,
                defaultRoutingRule, simulationCurrent.minBatchSize ,simulationCurrent.maxBatchSize,simulationCurrent.getNumWorkCenters(),simulationCurrent.getBatchesRecorded() , simulationCurrent.getWarmupBatches(),
                simulationCurrent.getUtilLevel(), 1.2); //use this simulation, no warmup jobs? --- here, just to measure the behavior of rule, so need to get a steady state

        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        List<RoutingDecisionSituation> situations = simulation.routingDecisionSituations(minQueueLength);

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
//        situations = selectSituations(situations, state);

        simulation.reset();

        return new RoutingPhenoCharacterisation(defaultRoutingRule, situations);

    }

    public static PhenoCharacterisation currentTaskPhenoCharacterisation(DynamicSimulation simulation) {

        AbstractRule defaultSequencingRule = new WSPT(RuleType.SEQUENCING); //op.getProcTime() / op.getJob().getWeight();
        AbstractRule defaultRoutingRule = new WIQ(RuleType.ROUTING);

//        AbstractRule defaultSequencingRule = new GPRule(RuleType.SEQUENCING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[0]);
//        AbstractRule defaultRoutingRule = new GPRule(RuleType.ROUTING,((GPIndividual)state.population.subpops[0].individuals[0]).trees[1]);


        simulation.setSequencingRule(defaultSequencingRule);
        simulation.setRoutingRule(defaultRoutingRule);

        //fzhang 2019.6.22 change to 7, otherwise, can not get this kinds of simulations---because the simulation can not enough queue size as 7
        int minQueueLength = 2;
        int numDecisionSituations = 100;//used for measuring the behavior of different rules
        long shuffleSeed = 19980818;

        List<RoutingDecisionSituation> situations = simulation.routingDecisionSituations(minQueueLength);

        Collections.shuffle(situations, new Random(shuffleSeed));

        situations = situations.subList(0, numDecisionSituations);
//        situations = selectSituations(situations, state);

        simulation.reset();
//
        return new RoutingPhenoCharacterisation(defaultRoutingRule, situations);

    }

}