package yimei.jss.algorithm.downSampling;

import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
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

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static yimei.jss.algorithm.downSampling.KMedoidsPAM.clusterMembers;
import static yimei.jss.algorithm.downSampling.KMedoidsPAM.fit;
import static yimei.jss.algorithm.downSampling.SimpleStatisticsSaveRulesizeGen.writeAveGenRulesizeToFile;
import static yimei.jss.algorithm.surrogateCaseLS.GPRuleEvolutionStateSavedSurrogateCasesV6.calcSpearmanMatrix;

public class GPRuleEvolutionStateDownSamplingV2 extends GPRuleEvolutionStateDownSamplingV1 {

    SparseCaseSelector.Result rLasso;

    int triggerEnd = 0; //record the time that bigger than end time

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
            ArrayList<double[]> individualsCaseFitness = new ArrayList<>();
            ArrayList<Double> individualsRealFitness = new ArrayList<>();
            KMedoidsPAM.Result r = fit(indsCharListsMultiTree, k, 100, jobSeed);
//            System.out.println("Medoids (indices): " + Arrays.toString(r.medoids));
            List<List<Integer>> groups = clusterMembers(r.labels, r.medoids, k);
            for (int c = 0; c < k; c++) {

                List<Integer> indIndex; //collect all index of individual in the same cluster
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
                if (simulation.getSystemState().getClockTime() == Double.MAX_VALUE) { //means this is a bad run using this rule pair
                    simulationsInOneCluster.clear();
                } else {
                    for (Simulation sim : recordedSimulations) {
                        simulationsInOneCluster.add(((DynamicSimulation) sim).deepCloneAllSame());
//                        outputSystemState(sim);
                    }
//                    System.out.println("------------------------------");
                    String objectiveName = parameters.getStringWithDefault(
                            new Parameter("eval.problem.eval-model.objectives.0"), null, "");
                    Objective objective = Objective.get(objectiveName);
                    double[] caseFitness = simulation.objectiveValueMultiCase(objective,caseInOneInstance);
                    windowFitness.add(caseFitness);
                }
                simulation.reset();
                simulationsForIndividualsInOneCluster.put(indIndex, simulationsInOneCluster);
            }

/*            double[][] spearmanMatrixTemp = calcSpearmanMatrixFromIndividuals(windowFitness);
            finalSelectedCasesIndex = selectCases(spearmanMatrixTemp, caseInOneInstance, 0.8);*/
            windowFitness.clear();
            int n = caseInOneInstance;
            int fromIndex = Math.max(0, n - caseInOneInstance); // 确保不越界
            for (int i = fromIndex; i < n; i++) {
                finalSelectedCasesIndex.add(i);
            }

//            linearRegression = SimpleLinearRegression.fit(individualsCaseFitness, individualsRealFitness);
//            System.out.println("Intercept: " + linearRegression.intercept);
//            System.out.println("Coefficients: " + Arrays.toString(linearRegression.coef));
//            RegressionMetrics.report(linearRegression,individualsCaseFitness,individualsRealFitness);
            // 代表性 case 选择（use regression)
            // 建模器
           /* SparseCaseSelector selector = new SparseCaseSelector(
                    individualsCaseFitness, individualsRealFitness);

            // 方案A：LASSO（只用 L1，稀疏更强）
            rLasso = selector.fitK(false, 1.0,caseInOneInstance/2);
            System.out.println(rLasso);
            finalSelectedCasesIndex = (ArrayList<Integer>) rLasso.selectedColumns;*/
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

        if(generation >= switchGen) {
            //now, we have fitness of all individual in all cases, calculate their spearman correlation
            for (int caseIndex = 0; caseIndex < this.population.subpops[0].individuals[0].fitness.trials.size(); caseIndex++) {
                ArrayList<Double> fitnessInOneCase = new ArrayList<>();
                for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
                    if (this.population.subpops[0].individuals[ind].fitness.fitness() < Double.MAX_VALUE) {
                        fitnessInOneCase.add((double) this.population.subpops[0].individuals[ind].fitness.trials.get(caseIndex));
                    }
                }
                double[] arr = fitnessInOneCase.stream()
                        .mapToDouble(Double::doubleValue)
                        .toArray();
                windowFitness.add(arr);
            }
            spearmanMatrix = calcSpearmanMatrix(windowFitness);

            windowFitness.clear();
        }

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

        // SHOULD WE QUIT?
//        if (generation == numGenerations - 1) {

        if (totalTime + duration >= endTime) {
            switchGen = 10000;
            triggerEnd++;
        }


        if (totalTime + duration >= endTime && triggerEnd > 1) {


            //realEvaluateFinalPopulation

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

        if(triggerEnd<1) {
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

    private void outputSystemState(Simulation sim) {
        int operationWaiting = 0;
        double workloadWaiting = 0;
        int operationRemaining = 0;
        double workloadRemaining = 0;

        for (WorkCenter workCenter:sim.getSystemState().getWorkCenters()) {
            workloadWaiting += workCenter.getWorkInQueue();
            operationWaiting += workCenter.getNumOpsInQueue();
        }

        for (Job job:sim.getSystemState().getJobsInSystem()){
            if(job.getId() < sim.getSystemState().getBatchesInSystem().get(sim.getSystemState().getBatchesInSystem().size()-1).jobs.get(0).getId()) {
                workloadRemaining = job.workLoadRemaining;
                operationRemaining = job.operationRemaining;
            }
        }

        System.out.println("The waiting and remaining workload are " + workloadWaiting + " , " + workloadRemaining + ". The sum is " + (workloadWaiting + workloadRemaining));
        System.out.println("The waiting and remaining operations are " + operationWaiting + " , " + operationRemaining + ". The sum is " + (operationRemaining + operationWaiting));
    }

    private void assignNewFitness() {

            if(!((MultipleTreeRuleOptimizationProblem)evaluator.p_problem).getEvaluationModel().getObjectives().get(0).getName().startsWith("max")) {
                for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                    double weightedFitness = rLasso.predict(individual.fitness.trials);
                    double[] fitnesses = new double[]{weightedFitness};
                    ((MultiObjectiveFitness) individual.fitness).setObjectives(this, fitnesses);
                }
            }

    }

    public static void selectColumnsInPlace(ArrayList<double[]> X,
                                            List<Integer> colIndex,
                                            boolean keepSelected) {
        if (X == null || X.isEmpty()) return;
        int p = X.get(0).length;

        // 1) 去重且保持原顺序（用 LinkedHashSet）
        LinkedHashSet<Integer> set = new LinkedHashSet<>();
        for (int idx : colIndex) {
            if (idx >= 0 && idx < p) set.add(idx);
        }
        if (set.isEmpty()) {
            // 全删：保留 0 列
            for (int i = 0; i < X.size(); i++) X.set(i, new double[0]);
            return;
        }

        // 2) 构造 keep 掩码
        boolean[] keep = new boolean[p];
        if (keepSelected) {
            for (int idx : set) keep[idx] = true;
        } else {
            Arrays.fill(keep, true);
            for (int idx : set) keep[idx] = false;
        }

        // 3) 统计保留列数并生成新列顺序（按原列顺序输出，便于稳定）
        int m = 0;
        for (int j = 0; j < p; j++) if (keep[j]) m++;
        int[] newCols = new int[m];
        for (int j = 0, k = 0; j < p; j++) if (keep[j]) newCols[k++] = j;

        // 4) 重建每一行
        for (int i = 0; i < X.size(); i++) {
            double[] row = X.get(i);
            if (row.length != p)
                throw new IllegalStateException("行列数不一致：第 " + i + " 行是 " + row.length + ", 期望 " + p);
            double[] nr = new double[m];
            for (int k = 0; k < m; k++) nr[k] = row[newCols[k]];
            X.set(i, nr);
        }

    }
}
