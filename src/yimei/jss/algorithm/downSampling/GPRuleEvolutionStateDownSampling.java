package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.ParameterDatabase;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.Simulation;
import yimei.jss.simulation.event.BatchArrivalEvent;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;

import static yimei.jss.algorithm.downSampling.SimpleStatisticsSaveRulesizeGen.writeAveGenRulesizeToFile;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * in each generation. we use the simulation in this generation, so inds will have different PC
 * from generation 10, we record the simulation
 *
 * @author yimei
 */

public class GPRuleEvolutionStateDownSampling extends GPRuleEvolutionState {

    ArrayList<Integer> uniquePCNum = new ArrayList<>();
    ArrayList<Integer> samePCNumInCases = new ArrayList<>();

    Map<List<Integer>, List<Integer>> pcToIndexMap = new HashMap<>();

    public int pcDistance;

    public static final String P_REPLICATIONS = "num-Rep";
    public int numRep;

    public static final String P_ELITES = "num-elites";
    public static int numElites;

    public double totalTime;
    public double endTime;
    public double duration;
    public int caseInOneInstance;
    public static List<Simulation> recordedSimulations = new ArrayList<>();


    public static List<Simulation> finalSelectedSimulations = new ArrayList<>();

    public static ArrayList<Integer> startBatchID = new ArrayList<>();
    public static int batchInOneCase;

    ArrayList<double[]> fitnessInDifferentCases = new ArrayList<>();
    ArrayList<Individual> topIndividuals = new ArrayList<>();
    public ArrayList<Integer> finalSelectedCasesIndex = new ArrayList<>();
    double correlationThreshold;
    public int switchGen;
    public double[][] spearmanMatrix;

    // ===== 放在类字段位置（新增的状态变量）=====
    // 滑窗相关 / Sliding window for last 10 generations
    public Deque<DynamicSimulation> windowSimulations = new ArrayDeque<>();
    public ArrayList<double[]> windowFitness = new ArrayList<>(); // one double[] per case in window
    public int lastEliteUpdateGen = 0;   // 上次精英更新的代（10 的倍数）/ last elite refresh gen
    public int lastCaseEvaluatedIndex = 0; // recordedSimulations 已评估到的 index / last evaluated idx

    // 配置项：每多少代更新一次精英/滑窗宽度（与题述一致为 10）
    private static final int ELITE_PERIOD = 10;
    private static final int WINDOW_SPAN_GEN = 10;
    ArrayList<Integer> selectedCaseNum = new ArrayList<>();

    ArrayList<int[]> selectedCasesIndex = new ArrayList<>();
    ArrayList<double[]> selectedCasesStage = new ArrayList<>();

    int wholeBatchNum;

    int triggerEnd = 0; //record the time that bigger than end time

    @Override
    public void setup(EvolutionState state, Parameter base) {
        super.setup(this, base);

        phenoCharacterisation = new PhenoCharacterisation[2];

        endTime = parameters.getDoubleWithDefault(new Parameter("endTime"), null, 7200);
        caseInOneInstance = parameters.getIntWithDefault(new Parameter("caseInOneInstance"), null, 25);

        RuleOptimizationProblem problem = (RuleOptimizationProblem) state.evaluator.p_problem;
        phenoCharacterisation[0] = SequencingPhenoCharacterisation.currentGenPhenoCharacterisation((GPRuleEvolutionState) state, problem);
        phenoCharacterisation[1] = RoutingPhenoCharacterisation.currentGenPhenoCharacterisation((GPRuleEvolutionState) state, problem);

        correlationThreshold = parameters.getDoubleWithDefault(new Parameter("correlationThreshold"), null, 0.6);
        switchGen = parameters.getIntWithDefault(new Parameter("switchGen"), null, 10);

        int warmupBatch = ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0).getWarmupBatches();
        batchInOneCase = ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0).getBatchesRecorded() / caseInOneInstance;
        for (int i = 0; i < caseInOneInstance; i++) {
            startBatchID.add(warmupBatch + i * batchInOneCase);
        }

        DynamicSimulation simulation = ((DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0));

        wholeBatchNum = 6000 / ((simulation.maxBatchSize + simulation.minBatchSize) / 2);

    }

    public void run(int condition) {
        totalTime = 0;

        if (condition == C_STARTED_FRESH) {
            startFresh();
        } else {
            startFromCheckpoint();
        }

        int result = R_NOTDONE; //2, means not finished, continue to do
        while (result == R_NOTDONE) {
            //fzhang 21.7.2018  after startFresh(), we will in this loop
            //here, to reFresh breeder only
            //startFreshResetOperatorProb();

            long start = System.currentTimeMillis();//yimei.util.Timer.getCpuTime();

            result = evolve();//System.out.println(result);

            long finish = System.currentTimeMillis();// yimei.util.Timer.getCpuTime();
            double duration = (finish - start) / 1000;//000000;
            genTimes.add(duration);
            totalTime += duration;

            output.message("Generation " + (generation - 1) + " elapsed " + duration + " seconds.");// time used for each generation
        }

        output.message("The whole program elapsed " + totalTime + " seconds."); // time used for total program

//		File timeFile = new File("job." + jobSeed + ".time.csv"); //jobSeed = 0
        File timeFile = new File(out_dir + "/job." + jobSeed + ".time.csv"); //jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(timeFile));
            writer.write("Gen,Time");
            writer.newLine();
            for (int gen = 0; gen < genTimes.size(); gen++) {
                writer.write(gen + "," + genTimes.get(gen));
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        finish(result);
    }

    @Override
    public int evolve() {

        long start = System.currentTimeMillis();

        if (generation > 0)
            output.message("Generation " + generation);

        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        //in each generation, calculate phenoCharacterisation
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;

        //-------1. use one reference rule to run the simulation and record the sub-simulation cases-----------
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
        simulation.setRoutingRule(simulation.defaultRoutingRule);
        simulation.setSequencingRule(simulation.defaultSequencingRule);
/*        if(generation == 0) {
            simulation.setRoutingRule(simulation.defaultRoutingRule);
            simulation.setSequencingRule(simulation.defaultSequencingRule);
        } else {
            simulation.setSequencingRule(new GPRule(RuleType.SEQUENCING,((GPIndividual)topIndividuals.get(0)).trees[0]));
            simulation.setRoutingRule(new GPRule(RuleType.ROUTING,((GPIndividual)topIndividuals.get(0)).trees[1]));
        }*/
        simulation.runRecordSimulation();
        simulation.reset();


        // 2) 每 10 代更新精英并重置“滑窗的起点”
        // ====== 10 的倍数代：刷新精英 & 重评整个滑窗 ======
/*        if (generation > 0 && generation % ELITE_PERIOD == 0) {

            for (int idx = lastCaseEvaluatedIndex; idx < recordedSimulations.size(); idx++) {
                DynamicSimulation recSim = (DynamicSimulation) recordedSimulations.get(idx);
                windowSimulations.addLast(recSim);
            }

            // 1) 刷新精英
            PopulationUtils.sort(population.subpops[0].individuals);
            int numIndividuals = population.subpops[0].individuals.length;

            List<Integer> sorted = new ArrayList<>();
            for (int i = 0; i < numIndividuals; i++) {
                sorted.add(i);
            }
            int cutoff = (int) Math.ceil(numIndividuals * 0.3);
            List<Integer> top30 = sorted.subList(0, cutoff);

            // 从前 30% 随机抽 10 个（不放回）
            Collections.shuffle(top30, new Random(seed));
            List<Integer> chosen = top30.subList(0, Math.min(10, top30.size()));

            for (int i : chosen) {
                topIndividuals.add(population.subpops[0].individuals[i]);
            }

//            for (int i = 0; i < 10; i++) topIndividuals.add(population.subpops[0].individuals[i]);

            // 3) 用“新精英”重评整个滑窗
            windowFitness.clear();
            if (generation > 10) {
                for (int k = 0; k < 200 && !windowSimulations.isEmpty(); k++) {
                    windowSimulations.pollFirst();  // 删除最前面的
                }
            }
            for (DynamicSimulation recSim : windowSimulations) {
                windowFitness.add(evaluateOneCaseWithElites(recSim, topIndividuals, parameters));
                // 如果你的 evaluate 版本需要 objective/problem，可用你已有的重载
                // windowFitness.add(evaluateOneCaseWithElites(recSim, topIndividuals, problem, parameters));
            }

            // 4) 标记
            lastEliteUpdateGen = generation;

            // ====== 非 10 的倍数代：只评估“新增 case”，并追加到滑窗 ======
        } else {
            boolean elitesReady = (topIndividuals != null && !topIndividuals.isEmpty());

            // 找到新增 cases 的范围：从 lastCaseEvaluatedIndex 到 recordedSimulations.size()-1
            for (int idx = lastCaseEvaluatedIndex; idx < recordedSimulations.size(); idx++) {
                DynamicSimulation recSim = (DynamicSimulation) recordedSimulations.get(idx);
                windowSimulations.addLast(recSim);

                // 只有在“已经有精英”的情况下才评估新增 case；
                // 否则（gen=1–9）先累计 case，等到 gen=10 统一重评。
                if (elitesReady) {
                    windowFitness.add(evaluateOneCaseWithElites(recSim, topIndividuals, parameters));
                }
            }
        }*/
        windowSimulations.clear();
        windowFitness.clear();
        for (int idx = lastCaseEvaluatedIndex; idx < recordedSimulations.size(); idx++) {
            DynamicSimulation recSim = (DynamicSimulation) recordedSimulations.get(idx);
            windowSimulations.addLast(recSim);
        }

        // 更新“已评估到哪里”的全局指针
        lastCaseEvaluatedIndex = recordedSimulations.size();

        // 3) 达到 switchGen 之后，基于“滑窗”的 fitness 做 spearman + 代表性筛选
        if (generation >= switchGen) {

/*            for (int idx = 0; idx < windowSimulations.size(); idx++) {
                DynamicSimulation recSim = new ArrayList<>(windowSimulations).get(idx);
                windowFitness.add(evaluateOneCaseWithElites(recSim, topIndividuals, parameters));
            }

            // 用窗口数据计算 spearman（注意：只基于 windowFitness）
            spearmanMatrix = calcSpearmanMatrix(windowFitness);

            // 代表性 case 选择（返回的是窗口内 index；如需映射回 recordedSimulations 的全局索引，请在 window 里存一个 wrapper 带原始全局 idx）
            finalSelectedCasesIndex = selectCases(spearmanMatrix, 5, 1);*/

            int n = windowSimulations.size();
            int fromIndex = Math.max(0, n - caseInOneInstance); // 确保不越界
            for (int i = fromIndex; i < n; i++) {
                finalSelectedCasesIndex.add(i);
            }

            //2.comparison: randomly select
/*            List<Integer> numbers = new ArrayList<>();
            for (int i = 0; i < windowSimulations.size(); i++) {
                numbers.add(i);
            }
            // 打乱
            Collections.shuffle(numbers, new Random(19980818));
            // 取前 count 个
            finalSelectedCasesIndex = new ArrayList<>(numbers.subList(0, 15));*/
        }else{ // in generation0
            int n = windowSimulations.size();
            int fromIndex = Math.max(0, n - caseInOneInstance); // 确保不越界
            for (int i = fromIndex; i < n; i++) {
                finalSelectedCasesIndex.add(i);
            }
        }
            // 将选中的窗口内案例映射回真正的 simulation 对象
            finalSelectedSimulations = finalSelectedCasesIndex.stream()
                    .map(i -> ((List<DynamicSimulation>) new ArrayList<>(windowSimulations)).get(i))
                    .collect(Collectors.toList());
            selectedCaseNum.add(finalSelectedCasesIndex.size());

            int[] caseIndex = new int[finalSelectedSimulations.size()];
            double[] caseStage = new double[finalSelectedSimulations.size()];
            for (int sim = 0; sim < finalSelectedSimulations.size(); sim++) {
                caseIndex[sim] = (int) ((DynamicSimulation) finalSelectedSimulations.get(sim)).getSeed() / 10000;
                double raw = ((double) ((BatchArrivalEvent) finalSelectedSimulations
                        .get(sim).getEventQueue().iterator().next()).getBatch().getId()) / wholeBatchNum;
                caseStage[sim] = Math.round(raw * 100.0) / 100.0;
            }
            selectedCasesIndex.add(caseIndex);
            selectedCasesStage.add(caseStage);


        finalSelectedCasesIndex.clear();

        // EVALUATION
        statistics.preEvaluationStatistics(this);


        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

//        if (generation >= switchGen) {
////            transferFitness2Rank();
//            transferFitnessToNormalizedSum();
//        }

        topIndividuals.clear();
//        PopulationUtils.sort(population.subpops[0].individuals);
        for (int i = 0; i < 10; i++) {
            topIndividuals.add(population.subpops[0].individuals[i]);
        }

        statistics.postEvaluationStatistics(this);

//After evaluate all individuals, we record feature information about 5 elites
// SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }

        long finish = System.currentTimeMillis();// yimei.util.Timer.getCpuTime();
        duration = (finish - start) / 1000;//000000;

        if (totalTime + duration >= endTime) {
            switchGen = 10000;
            triggerEnd++;
        }

        // SHOULD WE QUIT?
//        if (generation == numGenerations - 1) {
        if (totalTime + duration >= endTime && triggerEnd > 1) {

            writeAveGenRulesizeToFile(parameters.getInt(new Parameter("pop.subpop.0.species.ind.numtrees"), null, 2));
            writeDiversityToFile(entropyDiversity);
            writeSelectedCasesIndexToFile(selectedCasesIndex);
            writeSelectedCasesStageToFile(selectedCasesStage);
            writeSelectedCasesNumToFile(selectedCaseNum);

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

        if (triggerEnd < 1) {
            // BREEDING
            statistics.preBreedingStatistics(this);

            population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N


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

    private void transferFitness2Rank() {

        int numIndividuals = population.subpops[0].individuals.length;
        // trials 是 ArrayList<Double>，所以先拿第一个个体的长度
        int numInstances = ((GPIndividual) population.subpops[0].individuals[0]).fitness.trials.size();

        // 保存所有个体在每个 instance 的 fitness
        double[][] trialFitness = new double[numIndividuals][numInstances];
        for (int ind = 0; ind < numIndividuals; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            ArrayList<Double> trials = individual.fitness.trials;
            for (int inst = 0; inst < numInstances; inst++) {
                trialFitness[ind][inst] = trials.get(inst);
            }
        }

        // 保存 rank
        double[][] ranks = new double[numIndividuals][numInstances];

        for (int inst = 0; inst < numInstances; inst++) {
            Double[] values = new Double[numIndividuals];
            Integer[] indices = new Integer[numIndividuals];
            for (int ind = 0; ind < numIndividuals; ind++) {
                values[ind] = trialFitness[ind][inst];
                indices[ind] = ind;
            }

            // 假设 fitness 越小越好 → rank 越高
            Arrays.sort(indices, (a, b) -> Double.compare(values[a], values[b]));

            for (int rank = 0; rank < numIndividuals; rank++) {
                int ind = indices[rank];
                ranks[ind][inst] = rank + 1; // rank 从 1 开始
            }
        }

        // 计算平均 rank，更新 fitness
        for (int ind = 0; ind < numIndividuals; ind++) {
            double sumRank = 0.0;
            for (int inst = 0; inst < numInstances; inst++) {
                sumRank += ranks[ind][inst];
            }
            double meanRank = sumRank / numInstances;

            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            // 更新 fitness，用 meanRank 替代
            double[] fitnesses = new double[]{meanRank};
            ((MultiObjectiveFitness) individual.fitness).setObjectives(this, fitnesses);
        }
    }

    public void transferFitnessToNormalizedSum() {
        int numIndividuals = population.subpops[0].individuals.length;

        // 取实例数量（这里用 trials 的长度）
        GPIndividual first = (GPIndividual) population.subpops[0].individuals[0];
        int numInstances = ((MultiObjectiveFitness) first.fitness).trials.size();

        // 收集原始 fitness：shape = [individual][instance]
        double[][] trialFitness = new double[numIndividuals][numInstances];
        for (int ind = 0; ind < numIndividuals; ind++) {
            GPIndividual gi = (GPIndividual) population.subpops[0].individuals[ind];
            ArrayList<Double> trials = ((MultiObjectiveFitness) gi.fitness).trials;
            for (int inst = 0; inst < numInstances; inst++) {
                trialFitness[ind][inst] = trials.get(inst);
            }
        }

        // 归一化（按实例做 min-max）：f_norm = (f - min) / (max - min)
        // 注意：若 max==min（全相等），设为 0，避免除零
        double[][] norm = new double[numIndividuals][numInstances];

        for (int inst = 0; inst < numInstances; inst++) {
            double min = Double.POSITIVE_INFINITY;
            double max = Double.NEGATIVE_INFINITY;

            // 先找到该 instance 下有效的 min 和 max
            for (int ind = 0; ind < numIndividuals; ind++) {
                double v = trialFitness[ind][inst];
                if (v == 1.0E12) continue; // 跳过无效
                if (v < min) min = v;
                if (v > max) max = v;
            }

            double denom = max - min;
            for (int ind = 0; ind < numIndividuals; ind++) {
                double v = trialFitness[ind][inst];
                if (v == 1.0E12) {
                    // 如果你希望完全惩罚这种个体，可以给它一个极大值
                    norm[ind][inst] = 1.0;  // 最差
                } else if (denom == 0.0) {
                    norm[ind][inst] = 0.0;  // 所有人都一样，设为 0
                } else {
                    norm[ind][inst] = (v - min) / denom; // 正常归一化
                }
            }
        }

        // 对每个体在所有 instance 上求和（你也可以改成求平均）
        for (int ind = 0; ind < numIndividuals; ind++) {
            double sumNorm = 0.0;
            for (int inst = 0; inst < numInstances; inst++) {
                sumNorm += norm[ind][inst];
            }
            // 如果更想要“平均归一化 fitness”，把 sumNorm /= numInstances;

            GPIndividual gi = (GPIndividual) population.subpops[0].individuals[ind];
            MultiObjectiveFitness mf = (MultiObjectiveFitness) gi.fitness;

            // 写回为单目标（最小化）；若你有多目标，可把其他指标一并放进数组
            mf.setObjectives(state, new double[]{sumNorm});
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


    public ArrayList<Integer> selectCases(double[][] spearmanMatrix, int maxCount, double corrThreshold) {
        int n = spearmanMatrix.length;
        double[][] distanceMatrix = new double[n][n];

        // Spearman → 距离矩阵: 1 - |ρ|
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distanceMatrix[i][j] = 1 - Math.abs(spearmanMatrix[i][j]);
            }
        }

        ArrayList<Integer> selected = new ArrayList<>();

        // ① 随机选择第一个
//        int firstCase = this.random[0].nextInt(n);

        // ① 第一个不随机：选择与所有其它样本平均距离最大的样本（最“多样”）
        int firstCase = -1;
        double bestAvgDist = -1.0;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                sum += distanceMatrix[i][j];
            }
            double avg = sum / (n - 1);
            if (avg > bestAvgDist) {
                bestAvgDist = avg;
                firstCase = i;
            }
        }
        selected.add(firstCase);

//        int consecutiveFail = 0; // 连续“全部高相关”的次数

        // ② max–min 选择
        while (selected.size() < maxCount) {
            int bestCase = -1;
            double bestMinDist = -1;

            for (int candidate = 0; candidate < n; candidate++) {
                if (selected.contains(candidate)) continue;

                double minDist = Double.MAX_VALUE;
                for (int chosen : selected) {
                    minDist = Math.min(minDist, distanceMatrix[candidate][chosen]);
                }

                if (minDist > bestMinDist) {
                    bestMinDist = minDist;
                    bestCase = candidate;
                }
            }

            if (bestCase == -1) break;

            // 如果需要相关性阈值可以打开这一行
             if (bestMinDist < 1 - corrThreshold) break;
            System.out.println(1-bestMinDist);
            selected.add(bestCase);
        }

        return selected;
    }

    // 评估“一个 case”在当前精英集下的适应度向量
    public double[] evaluateOneCaseWithElites(DynamicSimulation recordedSim,
                                               List<Individual> elites,
                                               ParameterDatabase parameters) {
        double[] fitnessInOneCase = new double[elites.size()];

        // 这里只克隆一次“记录的 simulation”，每个精英都基于同一个“原始记录”做深拷贝运行
        // （或者：对每个精英都从 recordedSim 再 clone，避免相互干扰）
        for (int i = 0; i < elites.size(); i++) {
            GPIndividual indi = (GPIndividual) elites.get(i);
            GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
            GPRule routingRule = new GPRule(RuleType.ROUTING, indi.trees[1]);

            DynamicSimulation sim = recordedSim.deepCloneAllSame();
            sim.setBothRuleWithoutReset(sequencingRule, routingRule);
            sim.run();

            String objectiveName = parameters.getStringWithDefault(
                    new Parameter("eval.problem.eval-model.objectives.0"), null, "");
            Objective objective = Objective.get(objectiveName);
            double ObjValue = sim.objectiveValue(objective);

            // 双保险 / double check (通常 run() 内部已中止)
            for (WorkCenter w : sim.getSystemState().getWorkCenters()) {
                if (w.numOpsInQueue() > 100) {
                    if (objective.getName().endsWith("profit")) ObjValue = -Double.MAX_VALUE;
                    else ObjValue = Double.MAX_VALUE;
                    break;
                }
            }
            sim.reset();
            fitnessInOneCase[i] = ObjValue;
        }
        return fitnessInOneCase;
    }

    public void writeSelectedCasesNumToFile(ArrayList<Integer> casesNum) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".selectedCasesNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,casesNum");
            writer.newLine();
            for (int i = 0; i < casesNum.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + switchGen-10000 + ", " + casesNum.get(i) + "\n");
            }
            casesNum.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    

    public void writeSelectedCasesIndexToFile(ArrayList<int[]> caseIndex) {
        // 输出目录
        File outDir = new File(out_dir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            System.err.println("Failed to create output directory: " + out_dir);
            return;
        }
        File csvFile = new File(outDir, "job." + jobSeed + ".selectedCasesIndex.csv");

        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(csvFile, false), StandardCharsets.UTF_8))) {

            // 写表头
            writer.write("Gen");
            if (caseIndex != null && !caseIndex.isEmpty()) {
                int cols = caseIndex.get(0).length;
                for (int c = 0; c < cols; c++) {
                    writer.write(",CaseIndex" + c);
                }
            }
            writer.newLine();

            // 写内容
            if (caseIndex != null) {
                for (int gen = 0; gen < caseIndex.size(); gen++) {
                    int[] indices = caseIndex.get(gen);
                    writer.write(String.valueOf(gen + switchGen-10000)); // 第一列：Gen

                    if (indices != null) {
                        for (int idx : indices) {
                            writer.write("," + idx);
                        }
                    }
                    writer.newLine();
                }
            }

            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeRepresentativeIndividualsNumToFile(ArrayList<int[]> indsNum) {
        // 输出目录
        File outDir = new File(out_dir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            System.err.println("Failed to create output directory: " + out_dir);
            return;
        }
        File csvFile = new File(outDir, "job." + jobSeed + ".representativeIndividualsNum.csv");

        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(csvFile, false), StandardCharsets.UTF_8))) {

            // 写表头
            writer.write("Gen");
            if (indsNum != null && !indsNum.isEmpty()) {
                int cols = indsNum.get(0).length;
                for (int c = 0; c < cols; c++) {
                    writer.write(",indsNum" + c);
                }
            }
            writer.newLine();

            // 写内容
            if (indsNum != null) {
                for (int gen = 0; gen < indsNum.size(); gen++) {
                    int[] indices = indsNum.get(gen);
                    writer.write(String.valueOf(gen)); // 第一列：Gen

                    if (indices != null) {
                        for (int idx : indices) {
                            writer.write("," + idx);
                        }
                    }
                    writer.newLine();
                }
            }

            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void writeSelectedCasesStageToFile(ArrayList<double[]> caseStage) {
        // 输出目录
        File outDir = new File(out_dir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            System.err.println("Failed to create output directory: " + out_dir);
            return;
        }
        File csvFile = new File(outDir, "job." + jobSeed + ".selectedCasesStage.csv");

        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(csvFile, false), StandardCharsets.UTF_8))) {

            // 表头
            writer.write("Gen");
            if (caseStage != null && !caseStage.isEmpty()) {
                int cols = caseStage.get(0).length;
                for (int c = 0; c < cols; c++) {
                    writer.write(",Case" + c);
                }
            }
            writer.newLine();

            // 内容
            if (caseStage != null) {
                for (int gen = 0; gen < caseStage.size(); gen++) {
                    double[] values = caseStage.get(gen);
                    writer.write(String.valueOf(gen + switchGen-10000)); // 第一列：Gen

                    if (values != null) {
                        for (double v : values) {
                            writer.write("," + String.format(Locale.US, "%.2f", v)); // 每个值单独一格
                        }
                    }
                    writer.newLine();
                }
            }

            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSelectedCasesWorkloadToFile(ArrayList<double[]> caseWorkload) {
        // 输出目录
        File outDir = new File(out_dir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            System.err.println("Failed to create output directory: " + out_dir);
            return;
        }
        File csvFile = new File(outDir, "job." + jobSeed + ".selectedCaseWorkload.csv");

        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(csvFile, false), StandardCharsets.UTF_8))) {

            // 表头
            writer.write("Gen");
            if (caseWorkload != null && !caseWorkload.isEmpty()) {
                int cols = caseWorkload.get(0).length;
                for (int c = 0; c < cols; c++) {
                    writer.write(",Case" + c);
                }
            }
            writer.newLine();

            // 内容
            if (caseWorkload != null) {
                for (int gen = 0; gen < caseWorkload.size(); gen++) {
                    double[] values = caseWorkload.get(gen);
                    writer.write(String.valueOf(gen + switchGen-10000)); // 第一列：Gen

                    if (values != null) {
                        for (double v : values) {
                            writer.write("," + String.format(Locale.US, "%.2f", v)); // 每个值单独一格
                        }
                    }
                    writer.newLine();
                }
            }

            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSelectedCasesCorrelationToFile(ArrayList<double[]> caseCorrelation) {
        // 输出目录
        File outDir = new File(out_dir);
        if (!outDir.exists() && !outDir.mkdirs()) {
            System.err.println("Failed to create output directory: " + out_dir);
            return;
        }
        File csvFile = new File(outDir, "job." + jobSeed + ".selectedCaseCorrelation.csv");

        try (BufferedWriter writer = new BufferedWriter(
                new OutputStreamWriter(new FileOutputStream(csvFile, false), StandardCharsets.UTF_8))) {

            // 表头
            writer.write("Gen");
            if (caseCorrelation != null && !caseCorrelation.isEmpty()) {
                int cols = caseCorrelation.get(0).length;
                for (int c = 0; c < cols; c++) {
                    writer.write(",Case" + c);
                }
            }
            writer.newLine();

            // 内容
            if (caseCorrelation != null) {
                for (int gen = 0; gen < caseCorrelation.size(); gen++) {
                    double[] values = caseCorrelation.get(gen);
                    writer.write(String.valueOf(gen + switchGen-10000)); // 第一列：Gen

                    if (values != null) {
                        for (double v : values) {
                            writer.write("," + String.format(Locale.US, "%.2f", v)); // 每个值单独一格
                        }
                    }
                    writer.newLine();
                }
            }

            writer.flush();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}


