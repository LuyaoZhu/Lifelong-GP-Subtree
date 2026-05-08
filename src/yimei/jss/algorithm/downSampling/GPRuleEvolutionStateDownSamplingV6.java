package yimei.jss.algorithm.downSampling;

import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Job;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.Simulation;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static yimei.jss.algorithm.downSampling.KMedoidsPAM.fit;
import static yimei.jss.algorithm.downSampling.SimpleStatisticsSaveRulesizeGen.writeAveGenRulesizeToFile;

/**
 * 1. first use some individuals to select representative cases (10 -> 5)， and collect sub-simulations
 * 2. then determine a random/curriculum learning training order
 * 3. when selecting the startpoints, use the nearest individual
 * <p>
 * <p>
 *
 * @author Luyao Zhu
 */


public class GPRuleEvolutionStateDownSamplingV6 extends GPRuleEvolutionStateDownSamplingV2 {

    public int caseIndex;

    ArrayList<Integer> caseNumGen = new ArrayList<>();

    HashMap<int[], List<Simulation>> PCSimulations = new HashMap<>();

    ArrayList<int[]> representativeIndividualPC = new ArrayList<>();

    public int evolve() {

        long start = System.currentTimeMillis();

        if (generation > 0)
            output.message("Generation " + generation);

        //in each generation, calculate phenoCharacterisation
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);

        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        int k = 50;
        KMedoidsPAM.Result r = fit(indsCharListsMultiTree, k, 100, jobSeed);
        for (int c = 0; c < k; c++) {
            representativeIndividual.add(population.subpops[0].individuals[r.medoids[c]]);
            representativeIndividualPC.add(indsCharListsMultiTree[r.medoids[c]]);
        }

        if (generation >= switchGen) {

            PCSimulations.clear();

            for (int c = 0; c < representativeIndividual.size(); c++) {

                List<Simulation> simulations = new ArrayList<>();

                //-------1. use some rules to run the simulation and record the sub-simulation cases-----------
                recordedSimulations.clear();

                AbstractRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) (representativeIndividual.get(c))).trees[0]);
                AbstractRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) (representativeIndividual.get(c))).trees[1]);
                simulation.numBatchesRecorded = caseInOneInstance * batchInOneCase;
                simulation.setSequencingRule(sequencingRule);
                simulation.setRoutingRule(routingRule);
                simulation.runRecordSimulation();
                if (!(simulation.getSystemState().getClockTime() == Double.MAX_VALUE)) { //means this is a bad run using this rule pair

                    for (Simulation sim : recordedSimulations) {
                        simulations.add(((DynamicSimulation) sim).deepCloneAllSame());
//                        outputSystemState(sim);
                    }
//                    System.out.println("------------------------------");
                    String objectiveName = parameters.getStringWithDefault(
                            new Parameter("eval.problem.eval-model.objectives.0"), null, "");
                    Objective objective = Objective.get(objectiveName);
                    double[] caseFitness = simulation.objectiveValueMultiCase(objective, caseInOneInstance);
                    windowFitness.add(caseFitness);
                    PCSimulations.put(representativeIndividualPC.get(c), simulations);
                }
                simulation.reset();
            }

            double[][] spearmanMatrixTemp = calcSpearmanMatrixFromIndividuals(windowFitness);
            double ratio = parameters.getDouble(new Parameter("down-sampling-ratio"), null);
            finalSelectedCasesIndex = selectCases(spearmanMatrixTemp, (int) (caseInOneInstance*ratio), correlationThreshold);
            windowFitness.clear();

/*            for (int i = 0; i < caseInOneInstance; i++) {
                finalSelectedCasesIndex.add(i);
            }*/

            selectedCaseNum.add(finalSelectedCasesIndex.size());

            final Set<Integer> keep = new HashSet<>(finalSelectedCasesIndex);

            // 遍历并替换每个 entry 的 List<Simulation>
            PCSimulations.replaceAll((i, sims) ->
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
                caseStage[sim] = ((double) finalSelectedCasesIndex.get(sim)) / caseInOneInstance;
            }
            selectedCasesIndex.add(caseIndex);
            selectedCasesStage.add(caseStage);

        }

        representativeIndividual.clear();
        representativeIndividualPC.clear();

        // EVALUATION
        statistics.preEvaluationStatistics(this);


        if (generation >= switchGen) { //down-sampling
            downSampling(simulation);
        } else {
            evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population
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

        // SHOULD WE QUIT?
//        if (generation == numGenerations - 1) {
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

    private void downSampling(DynamicSimulation simulation) {

        String objectiveName = parameters.getStringWithDefault(
                new Parameter("eval.problem.eval-model.objectives.0"), null, "");
        Objective objective = Objective.get(objectiveName);

        //determine the order of cases
        // 1.randomly
        for (Map.Entry<int[], List<Simulation>> entry : PCSimulations.entrySet()) {
            Collections.shuffle(entry.getValue(), new Random(818));
        }

        //2. curriculum learning based on the remaining workload

/*        Random rand = new Random(818);
        List<Map.Entry<int[], List<Simulation>>> entries = new ArrayList<>(PCSimulations.entrySet());
        Map.Entry<int[], List<Simulation>> baseEntry = entries.get(rand.nextInt(entries.size()));
        List<Simulation> base = baseEntry.getValue();

        int n = base.size();
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) order[i] = i;

        Arrays.sort(order, Comparator.comparingDouble(i -> getWorkloadRemaining(base.get(i))));

        for (Map.Entry<int[], List<Simulation>> e : PCSimulations.entrySet()) {
            List<Simulation> sims = e.getValue();
            int m = Math.min(n, sims.size()); // 防止不同长度导致越界
            List<Simulation> reordered = new ArrayList<>(m);
            for (int k = 0; k < m; k++) {
                reordered.add(sims.get(order[k]));
            }
            sims.clear();
            sims.addAll(reordered);
        }

        for (Map.Entry<int[], List<Simulation>> e : PCSimulations.entrySet()) {
            for (Simulation sim : e.getValue()) {
                outputSystemState(sim);
            }
            System.out.println("-----------------------------------");
        }*/

        int selectedCaseIndex = 0;
        while (finalSelectedCasesIndex.size() > 0) {

            finalSelectedCasesIndex.remove(0);

            int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);

            //then evaluate and record simulation
            for (int i = 0; i < population.subpops[0].individuals.length; i++) {
                Individual ind = population.subpops[0].individuals[i];
                Simulation sim;

                int[] indPC = indsCharListsMultiTree[i];
                double minDistance = Double.MAX_VALUE;
                Simulation nearestSimulation = null;
                for (Map.Entry<int[], List<Simulation>> entry : PCSimulations.entrySet()) {
                    int[] indsCharList = entry.getKey();
                    double distance = PhenoCharacterisation.distance(indPC, indsCharList);
                    if (distance == 0) {
                        nearestSimulation = entry.getValue().get(selectedCaseIndex);
                        break;
                    } else if (distance < minDistance) {
                        minDistance = distance;
                        nearestSimulation = entry.getValue().get(selectedCaseIndex);
                    }
                }

                sim = ((DynamicSimulation) nearestSimulation).deepCloneAllSame();
                sim.setBothRuleWithoutReset(new GPRule(RuleType.SEQUENCING, ((GPIndividual) ind).trees[0]), new GPRule(RuleType.ROUTING, ((GPIndividual) ind).trees[1]));

                sim.run();

                double[] caseFitness;
                if (!(sim.getSystemState().getClockTime() == Double.MAX_VALUE)) {
                    caseFitness = new double[]{sim.objectiveValue(objective)};
                } else {
                    caseFitness = new double[]{Double.MAX_VALUE};
                }
                ((MultiObjectiveFitness) ind.fitness).setObjectives(this, caseFitness);
                ind.fitness.trials = new ArrayList<>();
                ind.fitness.trials.add(caseFitness[0]);
                ind.evaluated = true;
                sim.reset();
            }
            selectedCaseIndex ++;
            if (!finalSelectedCasesIndex.isEmpty()) {
                population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N
            }
        }
    }

    public double getWorkloadRemaining(Simulation sim) {
        double workloadWaiting = 0;
        double workloadRemaining = 0;

        for (WorkCenter workCenter : sim.getSystemState().getWorkCenters()) {
            workloadWaiting += workCenter.getWorkInQueue();
        }

        for (Job job : sim.getSystemState().getJobsInSystem()) {
            if (job.getId() < sim.getSystemState().getBatchesInSystem().get(sim.getSystemState().getBatchesInSystem().size() - 1).jobs.get(0).getId()) {
                workloadRemaining = job.workLoadRemaining;
            }
        }
        return workloadRemaining + workloadWaiting;
    }

    private void outputSystemState(Simulation sim) {
        int operationWaiting = 0;
        double workloadWaiting = 0;
        int operationRemaining = 0;
        double workloadRemaining = 0;

        for (WorkCenter workCenter : sim.getSystemState().getWorkCenters()) {
            workloadWaiting += workCenter.getWorkInQueue();
            operationWaiting += workCenter.getNumOpsInQueue();
        }

        for (Job job : sim.getSystemState().getJobsInSystem()) {
            if (job.getId() < sim.getSystemState().getBatchesInSystem().get(sim.getSystemState().getBatchesInSystem().size() - 1).jobs.get(0).getId()) {
                workloadRemaining = job.workLoadRemaining;
                operationRemaining = job.operationRemaining;
            }
        }

        System.out.println("The waiting and remaining workload are " + workloadWaiting + " , " + workloadRemaining + ". The sum is " + (workloadWaiting + workloadRemaining));
        System.out.println("The waiting and remaining operations are " + operationWaiting + " , " + operationRemaining + ". The sum is " + (operationRemaining + operationWaiting));
    }


    public static double[] spearmanWithLastCase(List<double[]> indivVectors) {
        if (indivVectors == null || indivVectors.isEmpty()) return new double[0];

        int numIndividuals = indivVectors.size();
        int numCases = indivVectors.get(0).length;

        // 校验每行长度一致
        for (int r = 1; r < numIndividuals; r++) {
            if (indivVectors.get(r).length != numCases) {
                throw new IllegalArgumentException("每个个体的向量长度（case 数）必须一致。");
            }
        }

        SpearmansCorrelation sc = new SpearmansCorrelation();

        // 提取“最后一个 case”的列向量
        double[] lastCol = new double[numIndividuals];
        for (int r = 0; r < numIndividuals; r++) {
            lastCol[r] = indivVectors.get(r)[numCases - 1];
        }

        double[] corr = new double[numCases];
        corr[numCases - 1] = 1.0; // 与自身相关为1

        // 复用缓冲，避免重复分配
        double[] tmp = new double[numIndividuals];
        for (int c = 0; c < numCases - 1; c++) {
            for (int r = 0; r < numIndividuals; r++) {
                tmp[r] = indivVectors.get(r)[c];
            }
            double v = sc.correlation(tmp, lastCol);
            corr[c] = Double.isNaN(v) ? 0.0 : v; // 如遇常数列导致NaN，可按需设为0或保留NaN
        }

        return corr;
    }


}
