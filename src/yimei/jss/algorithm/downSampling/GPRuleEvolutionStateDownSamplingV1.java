package yimei.jss.algorithm.downSampling;

import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.ParameterDatabase;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Job;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.MultipleTreeRuleOptimizationProblem;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.Simulation;
import yimei.jss.simulation.state.SystemState;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static yimei.jss.algorithm.downSampling.KMedoidsPAM.clusterMembers;
import static yimei.jss.algorithm.downSampling.KMedoidsPAM.fit;
import static yimei.jss.algorithm.downSampling.SimpleStatisticsSaveRulesizeGen.writeAveGenRulesizeToFile;

public class GPRuleEvolutionStateDownSamplingV1 extends GPRuleEvolutionStateDownSampling {

    ArrayList<Individual> representativeIndividual = new ArrayList<>();
    List<List<Simulation>> recordedSimulationsClusters = new ArrayList<>();

    HashMap<List<Integer>, List<Simulation>> simulationsForIndividualsInOneCluster = new HashMap<>();

    Weights w;

    public int evolve() {

        long start = System.currentTimeMillis();

        if (generation > 0)
            output.message("Generation " + generation);

        //in each generation, calculate phenoCharacterisation
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;

        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        //use 10 representative individuals to collect individuals

        if (generation >= switchGen) {
            HashMap<List<Integer>,Double> realFitnessRepresentativeIndividuals = new HashMap<>();
            simulationsForIndividualsInOneCluster.clear();
            int k = 50;
            KMedoidsPAM.Result r = fit(indsCharListsMultiTree, k, 100, jobSeed);
//            System.out.println("Medoids (indices): " + Arrays.toString(r.medoids));
            List<List<Integer>> groups = clusterMembers(r.labels, r.medoids, k);
            for (int c = 0; c < k; c++) {

                List<Integer> indIndex; //collect all index of individual in the same cluster

//            System.out.println("Cluster " + c + " members " + groups.get(c));
//            System.out.println("Cluster " + c + " medoid instance index = " + r.medoids[c]);

                indIndex = groups.get(c);
                representativeIndividual.add(population.subpops[0].individuals[r.medoids[c]]);


                //-------1. use one reference rule to run the simulation and record the sub-simulation cases-----------
                DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
                recordedSimulations.clear();
                List<Simulation> simulationsInOneCluster = new ArrayList<>();

                AbstractRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) (representativeIndividual.get(c))).trees[0]);
                AbstractRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) (representativeIndividual.get(c))).trees[1]);
                simulation.setSequencingRule(sequencingRule);
                simulation.setRoutingRule(routingRule);
                simulation.runRecordSimulation();
                if (recordedSimulations.size() < caseInOneInstance) { //means this is a bad run using this rule pair
                    simulationsInOneCluster.clear();
                } else {
                    for (Simulation sim : recordedSimulations) {
                        simulationsInOneCluster.add(((DynamicSimulation) sim).deepCloneAllSame());
                    }

                    double[] caseFitness = multiCaseMeanFlowtime(caseInOneInstance,simulation.getSystemState());

                    windowFitness.add(caseFitness);
                }
                simulation.reset();
                (representativeIndividual.get(c)).evaluated = true;
                simulationsForIndividualsInOneCluster.put(indIndex, simulationsInOneCluster);
            }

/*            for (Map.Entry<List<Integer>, List<Simulation>> entry : simulationsForIndividualsInOneCluster.entrySet()) {
                if (entry.getValue().size() > 0) {
                    int indIndex = entry.getKey().get(0); //the first individual is the representative one
                    Individual representativeIndividual = population.subpops[0].individuals[indIndex];
                    double[] fitnessInCases = new double[entry.getValue().size()+1];
                    for (int sim = 0; sim < entry.getValue().size(); sim++) {
                        DynamicSimulation simulation = (DynamicSimulation) entry.getValue().get(sim);
                        fitnessInCases[sim] = evaluateOneCaseWithOneRepresentative(simulation, representativeIndividual, parameters);

                    }
                    fitnessInCases[fitnessInCases.length-1] = realFitnessRepresentativeIndividuals.get(entry.getKey());
                    windowFitness.add(fitnessInCases);
                }
            }*/
            // 用窗口数据计算 spearman（注意：只基于 windowFitness）
            spearmanMatrix = calcSpearmanMatrixFromIndividuals(windowFitness);

            for (int a=0; a<spearmanMatrix[spearmanMatrix.length-1].length; a++){
                System.out.println(spearmanMatrix[spearmanMatrix.length-1][a] + ",");
            }

            // 代表性 case 选择（返回的是窗口内 index；如需映射回 recordedSimulations 的全局索引，请在 window 里存一个 wrapper 带原始全局 idx）
/*            finalSelectedCasesIndex = selectCases(spearmanMatrix, caseInOneInstance, 1);
            Collections.sort(finalSelectedCasesIndex);
            //assign weight to different cases
            w = assignCoverageWeights(
                    spearmanMatrix,
                    finalSelectedCasesIndex,
                    caseInOneInstance,   // 例如 10
                    true                 // true=相似度，false=距离
            );*/

            int n = caseInOneInstance;
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

            // 将选中的窗口内案例映射回真正的 simulation 对象

            selectedCaseNum.add(finalSelectedCasesIndex.size());

            final Set<Integer> keep = new HashSet<>(finalSelectedCasesIndex);

            // 遍历并替换每个 entry 的 List<Simulation>
            simulationsForIndividualsInOneCluster.replaceAll((i, sims) ->
                    (sims == null || sims.isEmpty())
                            ? sims
                            : IntStream.range(0, sims.size())
                            .filter(keep::contains)   // 只保留 keep 集合中的下标
                            .mapToObj(sims::get)
                            .collect(Collectors.toList())
            );


            int[] caseIndex = new int[finalSelectedCasesIndex.size()];
            double[] caseStage = new double[finalSelectedCasesIndex.size()];
            for (int sim = 0; sim < finalSelectedCasesIndex.size(); sim++) {
                caseIndex[sim] = finalSelectedCasesIndex.get(sim);
                caseStage[sim] = ((double) finalSelectedCasesIndex.get(sim)) /caseInOneInstance;
            }
            selectedCasesIndex.add(caseIndex);
            selectedCasesStage.add(caseStage);



        }

        // EVALUATION
        statistics.preEvaluationStatistics(this);

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

//        if (generation >= switchGen) {
////            transferFitness2Rank();
////            transferFitnessToNormalizedSum();
//            assignWeightedFitness();
////            assignNormalisedFitness();
//        }

        finalSelectedCasesIndex.clear();
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

    private void assignWeightedFitness() {

        if(!((MultipleTreeRuleOptimizationProblem)state.evaluator.p_problem).getEvaluationModel().getObjectives().get(0).getName().startsWith("max")) {
            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                double weightedFitness = 0;
                for (int c = 0; c < individual.fitness.trials.size(); c++) {
                    weightedFitness += ((ArrayList<Double>) individual.fitness.trials).get(c) * w.weights[c];
                }
                // 更新 fitness，用 meanRank 替代
                double[] fitnesses = new double[]{weightedFitness};
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, fitnesses);
            }
        }
    }

    private void assignNormalisedFitness() {
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            double NormalisedFitness = 0;
            for (int c=0; c < individual.fitness.trials.size(); c++){
                double baseFitness = windowFitness.get(0)[finalSelectedCasesIndex.get(c)];
                NormalisedFitness += ((ArrayList<Double>)individual.fitness.trials).get(c) / baseFitness ;
            }
            // 更新 fitness，用 meanRank 替代
            double[] fitnesses = new double[]{NormalisedFitness};
            ((MultiObjectiveFitness) individual.fitness).setObjectives(this, fitnesses);
        }

    }


    public double evaluateOneCaseWithOneRepresentative(
            DynamicSimulation recordedSim,
            Individual representative,
            ParameterDatabase parameters) {

        GPIndividual indi = (GPIndividual) representative;
        GPRule seq = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
        GPRule rou = new GPRule(RuleType.ROUTING, indi.trees[1]);

        DynamicSimulation sim = recordedSim.deepCloneAllSame();
        sim.setBothRuleWithoutReset(seq, rou);
        sim.run();

        String objectiveName = parameters.getStringWithDefault(
                new Parameter("eval.problem.eval-model.objectives.0"), null, "");
        Objective objective = Objective.get(objectiveName);
        double v = sim.objectiveValue(objective);

        for (WorkCenter w : sim.getSystemState().getWorkCenters()) {
            if (w.numOpsInQueue() > 100) {
                v = objective.getName().endsWith("profit") ? -Double.MAX_VALUE : Double.MAX_VALUE;
                break;
            }
        }
        sim.reset();
        return v;
    }
    public static double[][] calcSpearmanMatrixFromIndividuals(List<double[]> indivVectors) {
        if (indivVectors == null || indivVectors.isEmpty()) return new double[0][0];

        int numIndividuals = indivVectors.size();       // 个体数
        int numCases = indivVectors.get(0).length;      // case 数（每个向量长度应一致）

        // 校验列数一致
        for (int r = 1; r < numIndividuals; r++) {
            if (indivVectors.get(r).length != numCases) {
                throw new IllegalArgumentException("每个个体的向量长度（case 数）必须一致。");
            }
        }

        // 转置：case x individual
        double[][] caseByInd = new double[numCases][numIndividuals];
        for (int r = 0; r < numIndividuals; r++) {
            double[] row = indivVectors.get(r); // 这个个体在所有 cases 的得分
            for (int c = 0; c < numCases; c++) {
                caseByInd[c][r] = row[c];       // 第 c 个 case 的第 r 个个体得分
            }
        }

        SpearmansCorrelation sc = new SpearmansCorrelation();
        double[][] result = new double[numCases][numCases];
        for (int i = 0; i < numCases; i++) {
            result[i][i] = 1.0;
            for (int j = i + 1; j < numCases; j++) {
                double corr = sc.correlation(caseByInd[i], caseByInd[j]);
                result[i][j] = corr;
                result[j][i] = corr;
            }
        }
        return result;
    }

    // 计算代表 case 的覆盖权重：每个已选 case 初始=1（覆盖自己），
// 每个未选 case 归到与之最相似的已选 case，那个代表的权重+1。
// 返回两个数组：counts（整数计数），weights（归一化，总和=1）
    static class Weights {
        int[] counts;
        double[] weights;
        Weights(int[] c, double[] w) { this.counts = c; this.weights = w; }
    }

    static Weights assignCoverageWeights(
            double[][] spearmanMatrix,           // 相似度矩阵：ρ(i,j)，建议取绝对值
            List<Integer> selected,              // 已选代表的索引（例如 5 个）
            int totalCases,                      // 本实例的总 case 数（例如 10）
            boolean isSimilarity                 // true: 相似度（大好）；false: 距离（小好）
    ) {
        int K = selected.size();
        // map: caseIdx -> selected内部的位置 [0..K-1]
        Map<Integer,Integer> posInSelected = new HashMap<>();
        for (int s = 0; s < K; s++) posInSelected.put(selected.get(s), s);

        // 计数：代表覆盖数，先给每个代表+1（覆盖自己）
        int[] counts = new int[K];
        Arrays.fill(counts, 0);
        for (int s = 0; s < K; s++) counts[s] += 1;

        // 处理未选中的 case
        for (int c = 0; c < totalCases; c++) {
            if (posInSelected.containsKey(c)) continue; // 代表自己已计入
            // 在所有代表中找与 c 最相近的那个
            int bestS = -1;
            double bestVal = isSimilarity ? -Double.MAX_VALUE : Double.MAX_VALUE;
            for (int s = 0; s < K; s++) {
                int repIdx = selected.get(s);
                double val = spearmanMatrix[c][repIdx];
                // 通常用 |ρ| 当相似度更稳
                val = Math.abs(val);
                if (isSimilarity) {
                    if (val > bestVal) { bestVal = val; bestS = s; }
                } else { // 距离：越小越好
                    if (val < bestVal) { bestVal = val; bestS = s; }
                }
            }
            counts[bestS] += 1; // 把未选 c 归到 bestS 代表上
        }

        // 归一化：总和应为 totalCases（例如 10）
        double[] weights = new double[K];
        for (int s = 0; s < K; s++) weights[s] = counts[s] / (double) totalCases;
        return new Weights(counts, weights);
    }

    public double[] multiCaseMeanFlowtime(int numCase, SystemState systemState) {
        List<Job> done = systemState.getJobsCompleted();
        int totalJobs = done.size();
        double[] caseFitness = new double[numCase];

        // 稳定顺序：按 Job ID 排序
        done.sort(Comparator.comparingInt(Job::getId));

        // === 第一步：把 jobs 切成 numCase 个“窗口” ===
        int base = (numCase == 0) ? 0 : totalJobs / numCase;   // 每个窗口的基础 job 数
        int remainder = (numCase == 0) ? 0 : totalJobs % numCase; // 余数分给靠前的窗口
        int[] winSizes = new int[numCase];
        for (int w = 0; w < numCase; w++) {
            winSizes[w] = base + (w < remainder ? 1 : 0);
        }

        // 每个窗口的平均 flow time
        double[] winMeans = new double[numCase];
        int idx = 0;
        for (int w = 0; w < numCase; w++) {
            int k = winSizes[w];
            if (k <= 0) {                     // 没有数据的窗口
                winMeans[w] = Double.POSITIVE_INFINITY;
                continue;
            }
            double sum = 0.0;
            for (int j = 0; j < k; j++) {
                if (idx >= totalJobs) {       // 防御：不够取
                    sum = Double.POSITIVE_INFINITY;
                    break;
                }
                Job job = done.get(idx++);
                if (job == null) {            // 防御：异常空指针
                    sum = Double.POSITIVE_INFINITY;
                    break;
                }
                sum += job.flowTime();
            }
            winMeans[w] = (sum == Double.POSITIVE_INFINITY) ? Double.POSITIVE_INFINITY : (sum / k);
        }

        // === 第二步：第 i 个 case = 前 i+1 个窗口的平均（最后一个就是所有窗口的平均）===
        double cum = 0.0;
        int count = 0;
        for (int i = 0; i < numCase; i++) {
            cum += winMeans[i];
            count++;
            caseFitness[i] = (Double.isInfinite(cum)) ? Double.POSITIVE_INFINITY : (cum / count);
        }
        return caseFitness;
    }




}
