package yimei.jss.algorithm.lifelongGP;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.SortComparatorL;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import smile.clustering.HierarchicalClustering;
import smile.clustering.linkage.CompleteLinkage;
import yimei.jss.algorithm.surrogateCaseLS.SilhouetteScore;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.AbstractRuleHelper;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static yimei.jss.algorithm.surrogateAccuracy.surrogateClearingMultitreeEvaluatorV1.SamePCNum;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * in each task, select individuals with good quality and diversity as the initial population of next task
 * in addition, the final population will be evaluated in 10 unseen instances to build surrogate
 *
 * @author luyao
 */

public class GPRuleEvolutionStateLifelongGPV1 extends GPRuleEvolutionStateLifelongGPV3 {

    ArrayList<Integer> uniquePCNum = new ArrayList<>();
    ArrayList<Integer> samePCNumInCases = new ArrayList<>();

    Map<List<Integer>, List<Integer>> pcToIndexMap = new HashMap<>();

    public int pcDistance;

    public static final String P_REPLICATIONS = "num-Rep";
    public int numRep;

    public static final String P_ELITES = "num-elites";
    public static int numElites;

    public static int simulationsPerTask;

    public double coefficientCurrentTask;

    public static int generationPerTask;

    ArrayList<GPIndividual> QDIndividuals = new ArrayList<>();

    HashMap<List<Integer>, GPIndividual> PCIndividualMap = new HashMap<>();

    ArrayList<Individual> savedIndividuals = new ArrayList<>();

    ArrayList<Integer> archiveSampleNumber = new ArrayList<>();


    @Override
    public void setup(EvolutionState state, Parameter base) {
        Parameter p;
        //fzhang 2018.11.8 I need to do this to be able to load seed values in the AbtractRule class.
        AbstractRuleHelper.state = this;

        // Get the job seed.
        p = new Parameter("seed").push("" + 0);
        jobSeed = parameters.getLongWithDefault(p, null, 0);

        setupTerminals();

        super.setup(this, base);

        phenoCharacterisation = new PhenoCharacterisation[2];

        pcDistance = parameters.getIntWithDefault(new Parameter("pcDistance"), null, 1);

        simulationsPerTask = 10;

        coefficientCurrentTask = parameters.getDoubleWithDefault(new Parameter("coefficientCurrentTask"), null, 0.5);

        generationPerTask = parameters.getIntWithDefault(new Parameter("generationPerTask"), null, 50);

        phenoCharacterisation[0] =
                SequencingPhenoCharacterisation.defaultPhenoCharacterisation();
        phenoCharacterisation[1] =
                RoutingPhenoCharacterisation.defaultPhenoCharacterisation();


    }

    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        //in each generation, calculate phenoCharacterisation
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(generation / generationPerTask);

        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);

        for (int s = 0; s < indsCharListsMultiTree.length; s++) {
            ((GPIndividual) population.subpops[0].individuals[s]).PC = indsCharListsMultiTree[s];
        }

        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
//        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);


        // EVALUATION
        statistics.preEvaluationStatistics(this);

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

        double refFit = calculateReferenceRuleFitness(simulation);

        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            if (individual.fitness.fitness() < Double.MAX_VALUE) {
                double[] objective = new double[1];
                objective[0] = individual.fitness.fitness() / refFit;
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
            }
        }

        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind].clone();
            if (individual.fitness.fitness() < 1) {
                List<Integer> key = Arrays.stream(individual.PC)
                        .boxed()
                        .collect(Collectors.toList());
                PCIndividualMap.put(key, individual);
            }
        }

        System.out.println(PCIndividualMap.size());

        archiveSampleNumber.add(PCIndividualMap.size());

        //then the fitness of one individual should be normalised raw fitness + normalised estimated fitness
        if (generation >= generationPerTask) {

            double[][] estimatedFitness = new double[surrogateFitness.size()][population.subpops[0].individuals.length];

            for (int t = 0; t < surrogateSamples.size(); t++) {
                estimatedFitness[t] = evaluatePopulation(this, surrogateSamples.get(t), surrogateFitness.get(t), 10);
            }

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {

                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];

                double combinedFitness = 0;

                for (int t = 0; t < estimatedFitness.length; t++) {
                    combinedFitness += estimatedFitness[t][ind];
                }

                double[] objective = new double[1];
                objective[0] = coefficientCurrentTask * individual.fitness.fitness() + (1 - coefficientCurrentTask) * combinedFitness / estimatedFitness.length;

                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
            }

        }

        if (generation == generationPerTask - 1 || generation == generationPerTask * 2 - 1 || generation == generationPerTask * 3 - 1) {

            // -------------------------------
            // first select QD individuals (value) AND PCs (key) with guaranteed alignment
            // -------------------------------
            List<Map.Entry<List<Integer>, GPIndividual>> qdEntries =
                    new ArrayList<>(PCIndividualMap.entrySet());

            QDIndividuals = new ArrayList<>(qdEntries.size());
            double[][] PCs = new double[qdEntries.size()][];

            for (int i = 0; i < qdEntries.size(); i++) {
                Map.Entry<List<Integer>, GPIndividual> e = qdEntries.get(i);

                // value -> individual
                QDIndividuals.add(e.getValue());

                // key -> PC (int[])
                List<Integer> key = e.getKey();
                double[] pc = new double[key.size()];
                for (int j = 0; j < key.size(); j++) {
                    pc[j] = key.get(j);
                }
                PCs[i] = pc;
            }

            int[][] PC;

            if(PCs.length >= 500) {

                //then do k-means to select 500 individuals
                SimpleKMedoids.Result r = SimpleKMedoids.fit(PCs, population.subpops[0].individuals.length, 20250101L, 15);

                PC = new int[population.subpops[0].individuals.length][PCs[0].length];

                for (int a = 0; a < r.medoids.length; a++) {

                    int m = r.medoids[a];
                    savedIndividuals.add(QDIndividuals.get(m));   // 每个簇中心个体

                    for (int i = 0; i < PCs[m].length; i++) {
                        PC[a][i] = (int) PCs[m][i];
                    }
                }
            } else {
                PC = new int[PCs.length][PCs[0].length];

                for (int a = 0; a < QDIndividuals.size(); a++) {
                    savedIndividuals.add(QDIndividuals.get(a));   // 每个簇中心个体
                }

                for (int i = 0; i < PCs.length; i++) {
                    for (int j = 0; j < PCs[i].length; j++) {
                        PC[i][j] = (int) PCs[i][j];
                    }
                }
            }


            //then evaluate them

            double[] fitnessOneSurrogate = new double[savedIndividuals.size()];

            for (int s = 0; s < simulationsPerTask; s++) {

                simulation.reset(608 + s);
                double refFitness = calculateReferenceRuleFitness(simulation);

                // （可选但推荐）避免 refFit1[1] 为 0 导致除零
                if (refFitness == 0.0) {
                    throw new IllegalStateException("refFit1[1] is 0, cannot normalise ObjValue / refFit1[1].");
                }

                for (int ind = 0; ind < savedIndividuals.size(); ind++) {

                    GPIndividual indi = (GPIndividual) savedIndividuals.get(ind);
                    GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
                    GPRule routingRule = new GPRule(RuleType.ROUTING, indi.trees[1]);
                    simulation.setSequencingRule(sequencingRule);
                    simulation.setRoutingRule(routingRule);

                    simulation.run();
                    String objectiveName = parameters.getStringWithDefault(
                            new Parameter("eval.problem.eval-model.objectives.0"), null, "");
                    Objective objective = Objective.get(objectiveName);

                    double ObjValue = simulation.objectiveValue(objective);

                    for (WorkCenter w : simulation.getSystemState().getWorkCenters()) {
                        if (w.numOpsInQueue() > 100) {
                            if (objective.getName().endsWith("profit"))
                                ObjValue = -Double.MAX_VALUE;
                            else
                                ObjValue = Double.MAX_VALUE;
                            break;
                        }
                    }

                    simulation.reset();

                    // normalise
                    fitnessOneSurrogate[ind] += ObjValue / refFitness;
                }
            }

            for (int ind = 0; ind < fitnessOneSurrogate.length; ind++) {
                fitnessOneSurrogate[ind] /= simulationsPerTask;
            }

            surrogateFitness.add(fitnessOneSurrogate);
            surrogateSamples.add(PC);   // ✅ 这里 PCs[ind] 与 QDIndividuals.get(ind) 严格一一对应

        }

        statistics.postEvaluationStatistics(this);

        //After evaluate all individuals, we record feature information about 5 elites
        // SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }


        // SHOULD WE QUIT?
        if (generation == numGenerations - 1) {

            writeDiversityToFile(entropyDiversity);

//            writeAccuracyToFile(MSEGen, SpearmanCorrelationGen, SamePCNum, AveragePCDistance);

            writeArchiveSampleNumToFile(archiveSampleNumber);

//            writeTerimalOccuranceToFile(seqFeatureOccurrences, 0, "seq");
//            writeTerimalOccuranceToFile(rouFeatureOccurrences, 1, "rou");
//            if (ordFeatureOccurrences.size() >= 1)
//                writeTerimalOccuranceToFile(ordFeatureOccurrences, 2, "ord");

//            writeSpearmanToFile(spearmanEstimatedFitnessCurrentGen, "estimatedCurrent");

//            writeSpearmanToFile(spearmanRealFitnessCurrentGen, "realCurrent");

//            writeSpearmanToFile(spearmanEstimatedRealFitness, "estimatedReal");

//            writeFirstSelectedCasesNumToFile(firstSelectedCasesNum);

            generation++; // in this way, the last generation value will be printed properly.  fzhang 28.3.2018
            return R_FAILURE;
        }

        // PRE-BREEDING EXCHANGING
        statistics.prePreBreedingExchangeStatistics(this);

        population = exchanger.preBreedingExchangePopulation(this);  /** Simply returns state.population. */
        statistics.postPreBreedingExchangeStatistics(this);

        String exchangerWantsToShutdown = exchanger.runComplete(this);  /** Always returns null */
        if (exchangerWantsToShutdown != null) {
            output.message(exchangerWantsToShutdown);
            /*
             * Don't really know what to return here.  The only place I could
             * find where runComplete ever returns non-null is
             * IslandExchange.  However, that can return non-null whether or
             * not the ideal individual was found (for example, if there was
             * a communication error with the server).
             *
             * Since the original version of this code didn't care, and the
             * result was initialized to R_SUCCESS before the while loop, I'm
             * just going to return R_SUCCESS here.
             */

            return R_SUCCESS;
        }

        // BREEDING
        statistics.preBreedingStatistics(this);

        if (generation == generationPerTask - 1 || generation == generationPerTask * 2 - 1 || generation == generationPerTask * 3 - 1) {


            ArrayList<Individual> savedTopIndividuals = new ArrayList<>();
            PopulationUtils.sort(population);

            for (int i=0; i<population.subpops[0].individuals.length * seedingRatio; i++ ) {
                savedTopIndividuals.add(population.subpops[0].individuals[i]);
            }

            population.clear();

            population = initializer.initialPopulation(this, 0);

            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                for (int replace = 0; replace<savedTopIndividuals.size(); replace++) {
                    population.subpops[sub].individuals[replace] = savedTopIndividuals.get(replace);
                }
            }

            savedIndividuals.clear();
            PCIndividualMap.clear();

        } else {
            population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N
        }//        population = preselection();


// POST-BREEDING EXCHANGING
        statistics.postBreedingStatistics(this);   //position 1  here, a new pop has been generated.

// POST-BREEDING EXCHANGING
        statistics.prePostBreedingExchangeStatistics(this);

        population = exchanger.postBreedingExchangePopulation(this);   /** Simply returns state.population. */
        statistics.postPostBreedingExchangeStatistics(this);  //position 2


// Generate new instances if needed
        if (problem.getEvaluationModel().isRotatable()) {
            problem.rotateEvaluationModel();
        }


// INCREMENT GENERATION AND CHECKPOINT
        generation++;
        if (checkpoint && generation % checkpointModulo == 0) {
            output.message("Checkpointing");
            statistics.preCheckpointStatistics(this);
            Checkpoint.setCheckpoint(this);
            statistics.postCheckpointStatistics(this);
        }

        return R_NOTDONE;
    }

    private double calculateReferenceRuleFitness(DynamicSimulation simulation) {
        // 两个参考 rule 组合：SPT+WIQ 和 WSPT+WIQ

        AbstractRule referenceSeqRule2 = new WSPT(RuleType.SEQUENCING);
        AbstractRule referenceRouRule2 = new WIQ(RuleType.ROUTING);

        double fitness;

        // 统一从参数中取一次 objective
        String objectiveName = parameters.getStringWithDefault(
                new Parameter("eval.problem.eval-model.objectives.0"), null, "");
        Objective objective = Objective.get(objectiveName);

        // 计算两个参考组合的 fitness
        fitness = evaluateRulePair(simulation, referenceSeqRule2, referenceRouRule2, objective);

        return fitness;
    }

    /**
     * 使用给定的排序/路由规则运行一次仿真，返回该组合在指定 objective 下的值。
     * 如果仿真没有正常结束（clockTime == Double.MAX_VALUE），抛出异常。
     */
    private double evaluateRulePair(DynamicSimulation simulation,
                                    AbstractRule sequencingRule,
                                    AbstractRule routingRule,
                                    Objective objective) {

        simulation.setSequencingRule(sequencingRule);
        simulation.setRoutingRule(routingRule);
        simulation.run();

        try {
            if (simulation.getSystemState().getClockTime() == Double.MAX_VALUE) {
                // 说明这个 rule 组合导致系统跑挂了，这里直接报错
                throw new IllegalStateException(
                        "Reference rule pair " + sequencingRule.getClass().getSimpleName() +
                                " + " + routingRule.getClass().getSimpleName() +
                                " failed to finish the simulation (clockTime == Double.MAX_VALUE).");
            }

            return simulation.objectiveValue(objective);
        } finally {
            // 无论成功还是失败，都保证重置 simulation 状态
            simulation.reset();
        }
    }


    private void setupTerminals() {
        Parameter p;

        //Need to know how many populations we're expecting here, as will need
        //one terminal set per population
        int numSubPops = parameters.getInt(new Parameter("pop.subpops"), null);
        int numTrees = parameters.getInt(new Parameter("pop.subpop.0.species.ind.numtrees"), null);


        int num = Math.max(numSubPops, numTrees);

        if (num == 1) {

            p = new Parameter(P_TERMINALS_FROM);

            terminalsFrom = new String[]{parameters.getStringWithDefault(p,
                    null, "relative")};

            p = new Parameter(P_INCLUDE_ERC);
            //includeErc seems like does not have influence.
            includeErc = new boolean[]{parameters.getBoolean(p, null, false)};
            initTerminalSet();
        } else if (num == 2) {
            terminalsFrom = new String[num];
            includeErc = new boolean[num];
            int subPopNum = 0;

            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop1TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop1TerminalSet == null) {
                //might have provided other value by mistake, we should check for this
                p = new Parameter(P_TERMINALS_FROM);
                subPop1TerminalSet = parameters.getStringWithDefault(p,
                        null, "relative");
                output.warning("No terminal set for subpopulation 1 specified - using " + subPop1TerminalSet + ".");

            }
//            terminalsFrom[subPopNum] = subPop1TerminalSet;
            terminalsFrom[subPopNum] = subPop1TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop2TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop2TerminalSet == null) {
                //use whatever we settled on for first population
                subPop2TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 2 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop2TerminalSet;
            //TODO: Add support for erc - will be false by default

            initTerminalSet(); //right
        } else {
            terminalsFrom = new String[num];
            includeErc = new boolean[num];
            int subPopNum = 0;

            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop1TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop1TerminalSet == null) {
                //might have provided other value by mistake, we should check for this
                p = new Parameter(P_TERMINALS_FROM);
                subPop1TerminalSet = parameters.getStringWithDefault(p,
                        null, "relative");
                output.warning("No terminal set for subpopulation 1 specified - using " + subPop1TerminalSet + ".");

            }
            terminalsFrom[subPopNum] = subPop1TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop2TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop2TerminalSet == null) {
                //use whatever we settled on for first population
                subPop2TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 2 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop2TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop3TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop3TerminalSet == null) {
                //use whatever we settled on for first population
                subPop3TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 3 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop3TerminalSet;

            //TODO: Add support for erc - will be false by default

            initTerminalSet(); //right
        }
    }


    //2021.4.16 fzhang save the diversity value to csv
    public void writeDiversityToFile(ArrayList<double[]> entropyDiversity) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".diversity.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!entropyDiversity.isEmpty()) {
                for (int j = 0; j < entropyDiversity.get(0).length; j++) {
                    header.append(",diversitySubpop").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < entropyDiversity.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                for (double value : entropyDiversity.get(i)) {
                    line.append(", ").append(value);
                }
                writer.write(line.toString());
                writer.newLine();
            }

            entropyDiversity.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static class EliteComparator implements SortComparatorL {
        Individual[] inds;

        public EliteComparator(Individual[] inds) {
            super();
            this.inds = inds;
        }

        public boolean lt(long a, long b) {
            return inds[(int) b].fitness.betterThan(inds[(int) a].fitness);
        }

        public boolean gt(long a, long b) {
            return inds[(int) a].fitness.betterThan(inds[(int) b].fitness);
        }
    }

    //calculate the occurrance of one terminal
    public int countSubstring(String rule, String terminal) {
        int count = 0;
        int index = 0;

        while ((index = rule.indexOf(terminal, index)) != -1) {
            count++;
            index += terminal.length();
        }

        if (terminal.equals("W")) {
            return count
                    - countSubstring(rule, "WIQ")
                    - countSubstring(rule, "WKR")
                    - countSubstring(rule, "MWT")
                    - countSubstring(rule, "OWT")
                    - countSubstring(rule, "AWIS")
                    - countSubstring(rule, "MWIS")
                    - countSubstring(rule, "BWKR")
                    - countSubstring(rule, "AMWTF")
                    - countSubstring(rule, "AWIQF");
        } else if (terminal.equals("PT")) {
            return count
                    - countSubstring(rule, "NPT")
                    - countSubstring(rule, "TPT")
                    - countSubstring(rule, "APTF");
        } else if (terminal.equals("R")) {
            return count
                    - countSubstring(rule, "WKR")
                    - countSubstring(rule, "NOR")
                    - countSubstring(rule, "RDD")
                    - countSubstring(rule, "BNOR")
                    - countSubstring(rule, "BWKR");
        } else if (terminal.equals("WKR")) {
            return count
                    - countSubstring(rule, "BWKR");  // 防止 BWKR 统计两次 WKR
        } else if (terminal.equals("NOR")) {
            return count
                    - countSubstring(rule, "BNOR");  // 防止 BNOR 统计两次 NOR
        } else if (terminal.equals("MWT")) {
            return count
                    - countSubstring(rule, "AMWTF");  // 防止 BWKR 统计两次 WKR
        } else if (terminal.equals("WIQ")) {
            return count
                    - countSubstring(rule, "AWIQF");  // 防止 BNOR 统计两次 NOR
        } else if (terminal.equals("NIQ")) {
            return count
                    - countSubstring(rule, "ANIQF");  // 防止 BNOR 统计两次 NOR
        } else {
            return count;
        }
    }

    //2021.4.16 fzhang save the diversity value to csv
    public void writeTerimalOccuranceToFile(ArrayList<int[]> featureOccurance, int t, String s) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + s + "FeatureOccurance.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < terminals[t].length; i++) {
                sb.append(terminals[t][i].toStringForHumans());
                if (i < terminals[t].length - 1) {
                    sb.append(","); // 添加逗号分隔符
                }
            }

            writer.write("Gen" + "," + sb.toString());
            writer.newLine();
            for (int i = 0; i < featureOccurance.size(); i++) { //every two into one generation
                //writer.newLine();
                String[] stringArray = Arrays.stream(featureOccurance.get(i))
                        .mapToObj(String::valueOf)
                        .toArray(String[]::new);

                String content = String.join(",", stringArray);

                writer.write(i + "," + content + "\n");
            }

            featureOccurance.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // 将 int[] 转换为 List<Integer>，确保可以用作 HashMap 的键
    private static List<Integer> arrayToList(int[] arr) {
        List<Integer> list = new ArrayList<>();
        for (int num : arr) {
            list.add(num);
        }
        return list;
    }

    public static int findKneePoint(List<Double> values) {
        if (values.size() < 3) return -1; // 至少需要 3 个点

        // 找到最小值和最大值的索引
        int minIndex = 0, maxIndex = 0;
        for (int i = 0; i < values.size(); i++) {
            if (values.get(i) < values.get(minIndex)) minIndex = i;
            if (values.get(i) > values.get(maxIndex)) maxIndex = i;
        }

        // 直线参数 (P1: minIndex, P2: maxIndex)
        double x1 = minIndex, y1 = values.get(minIndex);
        double x2 = maxIndex, y2 = values.get(maxIndex);

        // 计算点到直线的最大垂直距离
        int inflectionIndex = -1;
        double maxDistance = -1;

        for (int i = 0; i < values.size(); i++) {
            if (i == minIndex || i == maxIndex) continue; // 跳过最小和最大点

            double x0 = i, y0 = values.get(i);
            double distance = pointToLineDistance(x1, y1, x2, y2, x0, y0);

            if (distance > maxDistance) {
                maxDistance = distance;
                inflectionIndex = i;
            }
        }

        return inflectionIndex;
    }

    public static double pointToLineDistance(double x1, double y1, double x2, double y2, double x0, double y0) {
        double numerator = Math.abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1);
        double denominator = Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
        return numerator / denominator;
    }

    public void writeFirstSelectedCasesNumToFile(ArrayList<Integer> casesNum) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".firstSelectedCasesNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,casesNum");
            writer.newLine();
            for (int i = 0; i < casesNum.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + casesNum.get(i) + "\n");
            }
            casesNum.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] computeRanks(double[] doubleArray) {
        int n = doubleArray.length;
        Integer[] indices = new Integer[n];

// 初始化索引数组
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }

// 按 doubleArray 的值排序索引数组
        Arrays.sort(indices, Comparator.comparingDouble(i -> doubleArray[i]));

// 生成排名数组
        double[] rankArray = new double[n];
        for (int rank = 0; rank < n; rank++) {
            rankArray[indices[rank]] = rank + 1; // 排名从 1 开始
        }

        return rankArray;
    }

    // 计算 Spearman 相关性矩阵
    public static double[][] computeSpearmanMatrix(double[][] A) {
        int numVectors = A.length;
        double[][] spearmanMatrix = new double[numVectors][numVectors];

        SpearmansCorrelation spearman = new SpearmansCorrelation();

        for (int i = 0; i < numVectors; i++) {
            for (int j = 0; j < numVectors; j++) {
                if (i == j) {
                    spearmanMatrix[i][j] = 1.0; // 自相关性 = 1
                } else {
                    spearmanMatrix[i][j] = spearman.correlation(A[i], A[j]);
                }
            }
        }
        return spearmanMatrix;
    }

    // 计算最佳 K（使用轮廓系数）
    public static int findBestK(double[][] distanceMatrix, int minK, int maxK) {
        int bestK = minK;
        double bestScore = -1;

        for (int k = minK; k <= maxK; k++) {
            int[] labels = hierarchicalClustering(distanceMatrix, k);
            double score = SilhouetteScore.computeSilhouette(distanceMatrix, labels);
            System.out.println("K = " + k + ", 轮廓系数 = " + score);

            if (score > bestScore) {
                bestScore = score;
                bestK = k;
            }
        }
        return bestK;
    }

    // 进行层次聚类
    public static int[] hierarchicalClustering(double[][] distanceMatrix, int numClusters) {
        // 使用 Complete Linkage 进行层次聚类
        HierarchicalClustering hc = new HierarchicalClustering(new CompleteLinkage(distanceMatrix));

        // 分割成 numClusters 个簇
        return hc.partition(numClusters);
    }

    // 计算距离矩阵 (1 - 相关性)
    public static double[][] computeDistanceMatrix(double[][] similarityMatrix) {
        int n = similarityMatrix.length;
        double[][] distanceMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distanceMatrix[i][j] = 1 - similarityMatrix[i][j]; // 相关性转距离
            }
        }
        return distanceMatrix;
    }

    // 计算均值
    private static double mean(List<Double> values) {
        return values.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
    }

    // 计算标准差
    private static double stddev(List<Double> values, double mean) {
        double variance = values.stream()
                .mapToDouble(v -> Math.pow(v - mean, 2))
                .average()
                .orElse(0.0);
        return Math.sqrt(variance);
    }

    // Z-score 归一化
    private static List<Double> zScoreNormalize(List<Double> values) {
        double mean = mean(values);
        double std = stddev(values, mean);
        double epsilon = 1e-9;  // 避免标准差为 0
        return values.stream()
                .map(v -> (v - mean) / (std + epsilon))
                .collect(Collectors.toList());
    }

    private static List<Double> minMaxNormalize(List<Double> values) {
        double min = values.stream().min(Double::compareTo).orElse(0.0);
        double max = values.stream().max(Double::compareTo).orElse(1.0);
        double epsilon = 1e-9;
        return values.stream()
                .map(v -> (v - min) / (max - min + epsilon))
                .collect(Collectors.toList());
    }

    // 对 fitness 进行 Log 变换
    private static List<Double> logTransform(List<Double> values) {
        double epsilon = 1e-9;
        return values.stream()
                .map(v -> Math.log(v + epsilon))  // 避免 log(0)
                .collect(Collectors.toList());
    }

    public void writeSpearmanToFile(ArrayList<double[]> spearman, String s) {
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + s + "Spearman.csv");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // 写入标题行（header）
            StringBuilder header = new StringBuilder("Gen");
            if (!spearman.isEmpty()) {
                for (int j = 0; j < spearman.get(0).length; j++) {
                    header.append(",").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // 写入每一行数据
            for (int i = 0; i < spearman.size(); i++) {
                StringBuilder line = new StringBuilder();
                line.append(i + 10); // 第 i 代（generation）

                double[] row = spearman.get(i);
                for (double value : row) {
                    line.append(",").append(String.format("%.4f", value)); // 保留4位小数
                }

                writer.write(line.toString());
                writer.newLine();
            }

            spearman.clear(); // 清空列表
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeAccuracyToFile(ArrayList<Double> mseGen, ArrayList<Double> spearmanCorrelationGen, ArrayList<Integer> samePCNum, ArrayList<Double> avePCDistance) {

        File weightFile = new File(out_dir + "/job." + jobSeed + "." + "SurrogateAccuracy.csv");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // 写入标题行（header）
            StringBuilder header = new StringBuilder("Gen");
            header.append(",").append("MSE").append(",").append("SpearmanCorrelation").append(",").append("SamePCNum").append(",").append("AvePCDistance");
            writer.write(header.toString());
            writer.newLine();

            // 写入每一行数据
            for (int i = 0; i < mseGen.size(); i++) {
                StringBuilder line = new StringBuilder();
                line.append(i); // 第 i 代（generation）

                line.append(",").append(String.format("%.4f", mseGen.get(i))).append(",").append(String.format("%.4f", spearmanCorrelationGen.get(i)))
                        .append(",").append(String.format("%d", SamePCNum.get(i))).append(",").append(String.format("%.4f", avePCDistance.get(i)));


                writer.write(line.toString());
                writer.newLine();
            }
            mseGen.clear();
            spearmanCorrelationGen.clear(); // 清空列表
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static ArrayList<Integer> getTopIndices(double[] array, double threshold, int maxSize) {
        // 使用一个临时列表保存满足条件的值和索引
        ArrayList<int[]> filteredList = new ArrayList<>();

        for (int i = 0; i < array.length; i++) {
            if (array[i] > threshold) {
                filteredList.add(new int[]{i, Double.valueOf(array[i]).hashCode()});  // 保存索引和对应值
            }
        }

        // 根据值排序，降序（值越大排越前）
        filteredList.sort((a, b) -> Double.compare(array[b[0]], array[a[0]]));

        // 最多只保留前 maxSize 个索引
        ArrayList<Integer> result = new ArrayList<>();
        for (int i = 0; i < Math.min(maxSize, filteredList.size()); i++) {
            result.add(filteredList.get(i)[0]);
        }

        return result;
    }

    // 计算每一行的中位数，返回 double[] 结果
    public static double[] rowWiseMedians(double[][] data) {
        if (data == null || data.length == 0) {
            throw new IllegalArgumentException("Input array is null or empty");
        }

        double[] medians = new double[data.length];

        for (int i = 0; i < data.length; i++) {
            double[] row = data[i];

            if (row == null || row.length == 0) {
                throw new IllegalArgumentException("Row " + i + " is null or empty");
            }

            double[] sortedRow = Arrays.copyOf(row, row.length);
            Arrays.sort(sortedRow);
            int len = sortedRow.length;

            if (len % 2 == 0) {
                medians[i] = (sortedRow[len / 2 - 1] + sortedRow[len / 2]) / 2.0;
            } else {
                medians[i] = sortedRow[len / 2];
            }
        }

        return medians;
    }

    public static double[] getTop100Min(double[] arr) {
        // 创建一个副本，避免修改原数组
        double[] copy = Arrays.copyOf(arr, arr.length);

        // 升序排序
        Arrays.sort(copy);

        // 如果长度小于100，返回全部；否则返回前100个
        return Arrays.copyOf(copy, Math.min(100, copy.length));
    }

    public Population preselection() {

        //ensure the elites can be saved in next generation
        ArrayList<Individual[]> elites = new ArrayList<Individual[]>(this.population.subpops.length);
        PopulationUtils.sort(this.population);

        numElites = this.parameters.getIntWithDefault(new Parameter(P_ELITES), null, 5);

        for (int pop = 0; pop < this.population.subpops.length; pop++) {
            List<Individual> tempElites = new ArrayList<>();
            for (int e = 0; e < numElites; e++) {
                tempElites.add(this.population.subpops[pop].individuals[e]);
            }
            if (elites.size() == 0 || elites.size() == 1) {
                elites.add(pop, tempElites.toArray(new Individual[tempElites.size()]));
            } else {
                Individual[] combineElites = ArrayUtils.addAll(elites.get(pop), tempElites.toArray(new Individual[tempElites.size()]));
                elites.set(pop, combineElites);
            }
        }

        Population newPop = (Population) this.population.emptyClone();//save the population with k*populationsize individuals
        Population tempNewPop; //save the population with populationsize individuals for combining them together to newPop
        numRep = this.parameters.getIntWithDefault(new Parameter(P_REPLICATIONS), null, 1);
        for (int i = 0; i < numRep; i++) {
            tempNewPop = breeder.breedPopulation(this);
            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                //combinedInds = new Individual[subpopsLength];
                //System.arraycopy(tempNewPop.subpops[sub].individuals, 0, combinedInds, 0, tempNewPop.subpops[sub].individuals.length);
                if (i == 0) {
                    newPop.subpops[sub].individuals = tempNewPop.subpops[sub].individuals;
                } else {
                    newPop.subpops[sub].individuals = ArrayUtils.addAll(newPop.subpops[sub].individuals, tempNewPop.subpops[sub].individuals);
                }
            }
        }

        population = newPop;

        //evaluate the population based on surrogate model
        evaluator.evaluatePopulation(this);
        //2021.7.29 when sorting the inds in the intermediate pop, mapping with the distance
/*        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            distanceFitness = new ArrayList<>();
            for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++){
                distanceFitness.add(new ImmutablePair<>(pcDistance.get(ind),this.population.subpops[0].individuals[ind].fitness.fitness()));
            }
            distanceFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better
            pcDistance = new ArrayList<>();
        }*/
        //================================end===========================

        PopulationUtils.sort(population); //sort the inds in the intermediate pop by estimated fitness


        //2021.7.26 this is for investigating the effectiveness of the surrogate. To see whether if there are good inds are ignored based on surrogate
        //=====================begin==================
/*        ArrayList<Pair<Integer, Double>> indexEstimatedFitness;
        indexEstimatedFitness = new ArrayList<>();
        for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
            indexEstimatedFitness.add(new ImmutablePair<>(ind, this.population.subpops[0].individuals[ind].fitness.fitness()));
        }
        indexEstimatedFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better
        int[] indexSelectedIndsEstimated = new int[population.subpops[0].individuals.length / numRep];
        for (int i = 0; i < population.subpops[0].individuals.length / numRep; i++) {
            indexSelectedIndsEstimated[i] = indexEstimatedFitness.get(i).getKey();
        }
        //2021.7.29 only consider the top 30% inds to investigate whether the selected inds by real evaluation are selected by surrogate or not.
        int[] indexSelectedIndsEstimatedTop30Percent = new int[(int)((population.subpops[0].individuals.length / numRep) * 0.3)];
        for (int i = 0; i < indexSelectedIndsEstimatedTop30Percent.length; i++) {
            indexSelectedIndsEstimatedTop30Percent[i] = indexEstimatedFitness.get(i).getKey();
        }

        //real evaluation
//        realEvaluationIntermediate = true;
//        evaluator.evaluatePopulation(this);
//        realEvaluationIntermediate = false;
//        ArrayList<Pair<Integer, Double>> indexRealFitness;
//        indexRealFitness = new ArrayList<>();
//        for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
//            indexRealFitness.add(new ImmutablePair<>(ind, this.population.subpops[0].individuals[ind].fitness.fitness()));//this is real fitness
//        }

        //2021.7.29 get the distance VS fitnessGap
        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            ArrayList<Double> distanceFinal  = new ArrayList<>();
            ArrayList<Double> fitnessGap = new ArrayList<>();
            ArrayList<Double> realFitness = new ArrayList<>();
            ArrayList<Double> surrogateFitness = new ArrayList<>();
            ArrayList<Integer> overOrUnderEstimated = new ArrayList<>();
            for (int i = 0; i < distanceFitness.size() / numRep; i++) {
                if (distanceFitness.get(i).getKey() != Double.MAX_VALUE) { //means this inds are evaluated, 1. cleared inds 2. generate by reproductions
                    distanceFinal.add(distanceFitness.get(i).getKey());

//                    System.out.println(i + "," + distanceFitness.get(i).getKey() + "," + Math.abs(this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue()));
                    fitnessGap.add(Math.abs(this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue()));
                    realFitness.add(this.population.subpops[0].individuals[i].fitness.fitness());
                    surrogateFitness.add(distanceFitness.get(i).getValue());
                    if (this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue() > 0) {
                        overOrUnderEstimated.add(-1); //estimated a smaller value --- bad becomes good
                    } else if (this.population.subpops[0].individuals[i].fitness.fitness() - distanceFitness.get(i).getValue() < 0) {
                        overOrUnderEstimated.add(1); //estimated a larger value --- treat good as bad (serious)
                    } else {
                        overOrUnderEstimated.add(0); //exactly the same
                    }
                }
            }
            writeDistanceFitnessGap(distanceFinal, realFitness, surrogateFitness, fitnessGap, overOrUnderEstimated);
        }

        //2021.7.27 the inds in the intermediate population is sorted based on the estimated fitness. We also get the real fitness of inds---ordered
        //by the estimated fitness
        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            writeEstimatedRealFitnessIndsIntermediatePop(indexEstimatedFitness, indexRealFitness);
        }

        indexRealFitness.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better

        int[] indexSelectedIndsReal = new int[population.subpops[0].individuals.length / numRep];
        for (int i = 0; i < population.subpops[0].individuals.length / numRep; i++) {
            indexSelectedIndsReal[i] = indexRealFitness.get(i).getKey();
        }

        //2021.7.29 only consider the top 30% inds to investigate whether the selected inds by real evaluation are selected by surrogate or not.
        int[] indexSelectedIndsReal30Percent = new int[(int)((population.subpops[0].individuals.length / numRep) * 0.3)];
        for (int i = 0; i < indexSelectedIndsReal30Percent.length; i++) {
            indexSelectedIndsReal30Percent[i] = indexRealFitness.get(i).getKey();
        }

        int tempSelectedRealUnselectedEstimated = 0;
        for (int selectedReal = 0; selectedReal < indexSelectedIndsReal.length; selectedReal++) {
            if (ArrayUtils.contains(indexSelectedIndsReal, indexSelectedIndsEstimated[selectedReal]) == false) {//if the selected inds
                //by real evaluation is not selected by estimated fitness (surrogate)
                tempSelectedRealUnselectedEstimated++;
            }
        }
        selectedRealUnselectedEstimated.add(tempSelectedRealUnselectedEstimated);

        //2021.7.29
        int tempSelectedRealUnselectedEstimatedtop30Percent = 0;
        for (int selectedReal = 0; selectedReal < indexSelectedIndsReal30Percent.length; selectedReal++) {
            if (ArrayUtils.contains(indexSelectedIndsReal30Percent, indexSelectedIndsEstimatedTop30Percent[selectedReal]) == false) {//if the selected inds
                //by real evaluation is not selected by estimated fitness (surrogate)
                tempSelectedRealUnselectedEstimatedtop30Percent++;
            }
        }
        selectedRealUnselectedEstimatedTop30Percent.add(tempSelectedRealUnselectedEstimatedtop30Percent);
//        System.out.println(selectedRealUnselectedEstimated);

        if (this.generation == 0 || this.generation == 10 || this.generation == 20 || this.generation == 30
                || this.generation == 40 || this.generation == 49) {
            writeRealEstimatedFitness(indexEstimatedFitness, indexRealFitness); //ordered by index
        }*/
        //=====================================end=====================================

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            population.subpops[sub].resize(population.subpops[sub].individuals.length / numRep);
        }

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            int e = 0;
            for (int replace = population.subpops[sub].individuals.length - 1; replace >= population.subpops[sub].individuals.length - numElites; replace--) {
                population.subpops[sub].individuals[replace] = elites.get(sub)[e];
                e++;
            }
        }

        return population;
    }

    public double[] evaluatePopulation(final EvolutionState state, int[][][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int[][][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulationRanks(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        double[] estimatedFitness = new double[population.subpops[0].individuals.length];

        for (int sub = 0; sub < state.population.subpops.length; sub++) {

            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;

                    int[][] pcIntermediate = indsCharListsIntermediatePop[i];
                    //calculate the fitness based on surrogate model
                    for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                        int[][] pcModel = indsCharListsMultiTree[pc];
                        double d = PhenoCharacterisation.spearmanDistance(pcIntermediate, pcModel);
                        if (d == 0) {
                            dMin = d;
                            index = pc;
                            break;
                        }

                        if (d < dMin) {
                            dMin = d;
                            index = pc;
                        }

                    }

                    estimatedFitness[i] = fitnessesForModel[index];
//                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel[index]);
                    //individual.evaluated = true;
                    //individual.fitness.trials = new ArrayList();
                    //individual.fitness.trials.add(individual.fitness.fitness());


                } else {
                    estimatedFitness[i] = Double.MAX_VALUE;
                }
            }

        }

        return estimatedFitness;

    }

    public double[] evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel, double threshold) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        double[] estimatedFitness = new double[population.subpops[0].individuals.length];

        for (int sub = 0; sub < state.population.subpops.length; sub++) {

            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;

                    int[] pcIntermediate = indsCharListsIntermediatePop[i];
                    //calculate the fitness based on surrogate model
                    for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                        int[] pcModel = indsCharListsMultiTree[pc];
                        double d = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                        if (d == 0) {
                            dMin = d;
                            index = pc;
                            break;
                        }

                        if (d < dMin) {
                            dMin = d;
                            index = pc;
                        }

                    }

                    if (dMin <= threshold) {
                        estimatedFitness[i] = fitnessesForModel[index];
                    } else {
                        estimatedFitness[i] = Double.MAX_VALUE;
                    }


//                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel[index]);
                    //individual.evaluated = true;
                    //individual.fitness.trials = new ArrayList();
                    //individual.fitness.trials.add(individual.fitness.fitness());


                } else {
                    estimatedFitness[i] = Double.MAX_VALUE;
                }
            }

        }

        return estimatedFitness;

    }

    private void writeArchiveSampleNumToFile(ArrayList<Integer> archiveSampleNumber) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".archiveSampleNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!archiveSampleNumber.isEmpty()) {
                header.append(",archiveSampleNum");

            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < archiveSampleNumber.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                line.append(", ").append(archiveSampleNumber.get(i));
                writer.write(line.toString());
                writer.newLine();
            }

            archiveSampleNumber.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}


