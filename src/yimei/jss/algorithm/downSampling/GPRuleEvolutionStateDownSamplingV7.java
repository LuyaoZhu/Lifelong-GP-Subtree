package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
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

import static yimei.jss.algorithm.downSampling.SimpleStatisticsSaveRulesizeGen.writeAveGenRulesizeToFile;

/**
 * 1. first use some individuals to select representative cases (10 -> 5)， and collect sub-simulations
 * 2. then determine a random/curriculum learning training order
 * 3. when selecting the startpoints, use the nearest individual
 * <p>
 * a universal version
 *
 * <p>
 * <p>
 *
 * @author Luyao Zhu
 */


public class GPRuleEvolutionStateDownSamplingV7 extends GPRuleEvolutionStateDownSamplingV2 {

    public int caseIndex;

    ArrayList<double[]> selectedCaseWorkload = new ArrayList<>();

    ArrayList<double[]> selectedCaseCorrelation = new ArrayList<>();

    HashMap<int[], List<Simulation>> PCSimulations = new LinkedHashMap<>();;

    ArrayList<int[]> representativeIndsNumber = new ArrayList<>();

    Boolean useReferenceRule;
    Boolean useRepresentativeInds;
    Boolean loopInOneGen;

    public boolean useRandomDownSampling;
    public boolean useInformedDownSampling;

    public static final String RANDOM_DOWN_SAMPLING = "randomDownSampling";
    public static final String INFORMED_DOWN_SAMPLING = "informedDownSampling";

    int[][] indsCharListsMultiTree;


    public void setup(EvolutionState state, Parameter base) {
        super.setup(this, base);

        //1.use which rule to collect simulations
        useReferenceRule = parameters.getBoolean(new Parameter("useReferenceRule"), null, false);
        useRepresentativeInds = parameters.getBoolean(new Parameter("useRepresentativeInds"), null, false);
        //2. whether there is a loop in one generation
        loopInOneGen = parameters.getBoolean(new Parameter("loopInOneGeneration"), null, false);

        useRandomDownSampling = parameters.getBoolean(new Parameter(RANDOM_DOWN_SAMPLING), null, false);
        useInformedDownSampling = parameters.getBoolean(new Parameter(INFORMED_DOWN_SAMPLING), null, false);
    }


    public int evolve() {

        long start = System.currentTimeMillis();

        if (generation > 0)
            output.message("Generation " + generation);

        //get the problem and simulation to evaluate individuals outside the evaluation class
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
        String objectiveName = parameters.getStringWithDefault(
                new Parameter("eval.problem.eval-model.objectives.0"), null, "");
//        List<Objective> objectives = problem.getObjectives();
        Objective objective = Objective.get(objectiveName);

        indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);

        for(int ind=0; ind<population.subpops[0].individuals.length; ind++){
            ((GPIndividual)population.subpops[0].individuals[ind]).PC = indsCharListsMultiTree[ind];
        }

        if (generation >= switchGen) {

            PCSimulations.clear();
            recordedSimulations.clear();

            //determine which rule to run the simulation and collect the sub-simulations
            if (useReferenceRule) { //means we use SPT and WIQ to run the simulation, the evaluation of all individual will use the same sub-simulations
                List<Simulation> simulations = new ArrayList<>(); //to save collected simulations

                simulation.setRoutingRule(simulation.defaultRoutingRule);
                simulation.setSequencingRule(simulation.defaultSequencingRule);
                simulation.runRecordSimulation();
                simulation.reset();

                for (Simulation sim : recordedSimulations) {
                    simulations.add(((DynamicSimulation) sim).deepCloneAllSame());
                }

                int[] a = new int[40];
                java.util.Arrays.fill(a, 1);
                PCSimulations.put(a, simulations); //so for reference rule, there is only one sample

            }

            //so after that, we will calculate the fitness of some representative individuals on all sub-simulations

            for (int c = 0; c < representativeIndividual.size(); c++) {
                //to save the fitness information
                double[] ObjValue = new double[1];
                double[] caseObjective = new double[caseInOneInstance];

                List<Simulation> simulations = new ArrayList<>();
                //-------1. use some rules to run the simulation and record the sub-simulation cases-----------
                recordedSimulations.clear();
                AbstractRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) (representativeIndividual.get(c))).trees[0]);
                AbstractRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) (representativeIndividual.get(c))).trees[1]);

                simulation.setSequencingRule(sequencingRule);
                simulation.setRoutingRule(routingRule);
                simulation.runRecordSimulation();
                if (!(simulation.getSystemState().getClockTime() == Double.MAX_VALUE)) { //means this is a bad run using this rule pair
                    for (Simulation sim : recordedSimulations) {
                        simulations.add(((DynamicSimulation) sim).deepCloneAllSame());
//                        outputSystemState(sim);
                    }
//                    System.out.println("------------------------------");
                    ObjValue[0] = simulation.objectiveValue(objective);
                    caseObjective = simulation.objectiveValueMultiCase(objective, caseInOneInstance);
                    windowFitness.add(caseObjective);

                    if (useRepresentativeInds) {
                        //means we use representative individuals to run the simulation, the evaluation of all individual will use different sub-simulations
                        PCSimulations.put(((GPIndividual)representativeIndividual.get(c)).PC, simulations);
                    }
                }
                simulation.reset();
            }

            //then select cases
            double ratio = parameters.getDouble(new Parameter("down-sampling-ratio"), null);
            int numSelectedCases = (int) Math.ceil(ratio * caseInOneInstance); // 向上取整

            //2. with random down-sampling
            if (useRandomDownSampling) {
                while (finalSelectedCasesIndex.size() < numSelectedCases) {
                    int randomCaseIndex = random[0].nextInt(caseInOneInstance);
                    if (!finalSelectedCasesIndex.contains(randomCaseIndex)) {
                        finalSelectedCasesIndex.add(randomCaseIndex);
                    }
                }
            } else if (useInformedDownSampling) { //2. informed down sampling
                double[][] spearmanMatrixTemp = calcSpearmanMatrixFromIndividuals(windowFitness);
                finalSelectedCasesIndex = selectCases(spearmanMatrixTemp, numSelectedCases, correlationThreshold);
            } else {
                //1.without down-sampling
                for (int i = 0; i < caseInOneInstance; i++) {
                    finalSelectedCasesIndex.add(i);
                }
            }

            //after obtaining the fitness, to select cases
            windowFitness.clear();

            selectedCaseNum.add(finalSelectedCasesIndex.size());

            int[] caseIndex = new int[finalSelectedCasesIndex.size()];
            double[] caseStage = new double[finalSelectedCasesIndex.size()];
            double[] caseWorkload = new double[finalSelectedCasesIndex.size()];
            //random to select one map to report the workload
            List<Map.Entry<int[], List<Simulation>>> entries = new ArrayList<>(PCSimulations.entrySet());
            List<Simulation> valueSimulations = entries.get(0).getValue();

            for (int sim = 0; sim < finalSelectedCasesIndex.size(); sim++) {
                caseIndex[sim] = finalSelectedCasesIndex.get(sim);
                caseStage[sim] = ((double) finalSelectedCasesIndex.get(sim)) / caseInOneInstance;
                caseWorkload[sim] = getWorkloadRemaining(valueSimulations.get(finalSelectedCasesIndex.get(sim)));
            }
            selectedCasesIndex.add(caseIndex);
            selectedCasesStage.add(caseStage);
            selectedCaseWorkload.add(caseWorkload);

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

        }

        representativeIndividual.clear();


        // EVALUATION
        statistics.preEvaluationStatistics(this);

        if (loopInOneGen) {
            if (generation >= switchGen)  //down-sampling
                downSampling(simulation);
            else
                evaluator.evaluatePopulation(this);

        } else {
            evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population
        }

        PopulationUtils.sort(population);
/*        ArrayList<double[]> trainTestFitness = new ArrayList<>();
        if(generation >= switchGen) {
            //after evaluate all individuals, calculate the test fitness
            String batchSize = parameters.getString(new Parameter("eval.problem.eval-model.sim-models.0.batch-size"), null);
            SchedulingSet testSet = SchedulingSet.generateSet(968356, new String("dynamic-job-shop"), new String("missing-0.85-1.2-" + batchSize), objectives, 50);

            for (int i = 0; i < population.subpops[0].individuals.length*0.3; i++) {
                double[] TrainTestFitness = new double[2];
                GPIndividual ind = (GPIndividual) population.subpops[0].individuals[i];
                TrainTestFitness[0] = ind.fitness.fitness();
                AbstractRule sequencingRule = new GPRule(RuleType.SEQUENCING, ind.trees[0]);
                AbstractRule routingRule = new GPRule(RuleType.ROUTING, ind.trees[1]);
                Fitness fitness = new MultiObjectiveFitness();
                sequencingRule.calcFitness(  //in calcFitness(), it will check which one is routing/sequencing rule
                        fitness, null,
                        testSet, routingRule, objectives);
                TrainTestFitness[1] = fitness.fitness();
                trainTestFitness.add(TrainTestFitness);
            }
        }*/

        //after evaluate all individuals, we select representative individuals
//        PopulationUtils.sort(population);
        indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        for(int ind=0; ind<population.subpops[0].individuals.length; ind++){
            ((GPIndividual)population.subpops[0].individuals[ind]).PC = indsCharListsMultiTree[ind];
        }

        int[] PCNum = new int[2]; //the first one to save the whole unique PC number, the second is to save the unique PC number for the 50% individuals

        DedupResult result1 = dedupRowsWithIndices(indsCharListsMultiTree);
        PCNum[0] = result1.keptIndices.length;

        int[][] first50 = Arrays.copyOfRange(indsCharListsMultiTree, 0, (int) (indsCharListsMultiTree.length * 0.5));
        DedupResult result = dedupRowsWithIndices(first50);

        if(result.keptIndices.length < 5){  //if the diversity is too low, maybe there is only one individual
            result = result1;
        }

        int[][] deduped = result.distinct;     // 去重后的内容
        int[] idx = result.keptIndices;        // 原始下标（例如 [0, 3, 7, ...]）
        PCNum[1] = idx.length;
        for (int i = 0; i < idx.length; i++) {
            representativeIndividual.add((GPIndividual)population.subpops[0].individuals[idx[i]].clone());
        }
        System.out.println(idx.length);
        representativeIndsNumber.add(PCNum);

/*      int k = representativeIndsNum;
        KMedoidsPAM.Result r = fit(first50, k, 100, jobSeed);
        for (int c = 0; c < k; c++) {
            representativeIndividual.add((GPIndividual)population.subpops[0].individuals[r.medoids[c]].clone());
        }*/

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
        if (totalTime + duration >= endTime) { //use full-evaluation to output the best individual
            switchGen = 10000;
            triggerEnd++;
        }


        if (totalTime + duration >= endTime && triggerEnd > 1) {

            writeAveGenRulesizeToFile(parameters.getInt(new Parameter("pop.subpop.0.species.ind.numtrees"), null, 2));
            writeDiversityToFile(entropyDiversity);
            writeSelectedCasesIndexToFile(selectedCasesIndex);
            writeSelectedCasesStageToFile(selectedCasesStage);
            writeSelectedCasesNumToFile(selectedCaseNum);
            writeRepresentativeIndividualsNumToFile(representativeIndsNumber);
            writeSelectedCasesWorkloadToFile(selectedCaseWorkload);
            if(useInformedDownSampling){
                writeSelectedCasesCorrelationToFile(selectedCaseCorrelation);
            }

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

       /* Random rand = new Random(818);
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
        }*/

//        for (Map.Entry<int[], List<Simulation>> e : PCSimulations.entrySet()) {
//            for (Simulation sim : e.getValue()) {
//                outputSystemState(sim);
//            }
//            System.out.println("-----------------------------------");
//        }

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
            selectedCaseIndex++;
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


/*    public ArrayList<Integer> selectCases(double[][] spearmanMatrix, int maxCount, double corrThreshold) {
        int n = spearmanMatrix.length;
        double[][] distanceMatrix = new double[n][n];

        // Spearman → 距离矩阵: 1 - |ρ|
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distanceMatrix[i][j] = 1 - Math.abs(spearmanMatrix[i][j]);
            }
        }

        double[] correlation = new double[maxCount];

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
        correlation[selected.size()] = 1-bestAvgDist;
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
            correlation[selected.size()] = 1-bestMinDist;
            selected.add(bestCase);
        }

        selectedCaseCorrelation.add(correlation);
        return selected;
    }*/

    public ArrayList<Integer> selectCases(double[][] spearmanMatrix, int maxCount, double corrThreshold) {
        if (spearmanMatrix == null || spearmanMatrix.length == 0)
            throw new IllegalArgumentException("spearmanMatrix is null/empty");
        if (maxCount <= 0) throw new IllegalArgumentException("maxCount must be > 0");
        final int n = spearmanMatrix.length;

        // 校验方阵与ρ∈[-1,1]且有限
        for (int i = 0; i < n; i++) {
            if (spearmanMatrix[i] == null || spearmanMatrix[i].length != n)
                throw new IllegalArgumentException("spearmanMatrix must be square at row " + i);
            for (int j = 0; j < n; j++) {
                double rho = spearmanMatrix[i][j];
                if (!Double.isFinite(rho)) throw new IllegalStateException("rho non-finite at ["+i+"]["+j+"]: "+rho);
                if (Math.abs(rho) > 1.0000001) throw new IllegalStateException("rho out of [-1,1] at ["+i+"]["+j+"]: "+rho);
            }
        }

        // 距离矩阵 d=1-|ρ| 并夹到[0,1]
        double[][] distanceMatrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double d = 1.0 - Math.abs(spearmanMatrix[i][j]);
                if (d < 0) d = 0; else if (d > 1) d = 1;
                if (!Double.isFinite(d)) throw new IllegalStateException("distance non-finite at ["+i+"]["+j+"]: "+d);
                distanceMatrix[i][j] = d;
            }
        }

        double[] correlation = new double[maxCount];
        ArrayList<Integer> selected = new ArrayList<>();

        // ① 选第一个：平均距离最大
        int firstCase = -1; double bestAvgDist = -1.0;
        if (n == 1) { firstCase = 0; bestAvgDist = 0.0; }
        else {
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) if (i != j) sum += distanceMatrix[i][j];
                double avg = sum / (n - 1);
                if (!Double.isFinite(avg)) throw new IllegalStateException("avg non-finite for i="+i+": "+avg);
                if (avg > bestAvgDist) { bestAvgDist = avg; firstCase = i; }
            }
        }
        if (firstCase < 0) throw new IllegalStateException("firstCase stayed -1");
        correlation[0] = 1.0 - bestAvgDist;
        if (!Double.isFinite(correlation[0])) throw new IllegalStateException("correlation[0] non-finite");
        selected.add(firstCase);

        // ② max–min 选择
        while (selected.size() < maxCount) {
            int bestCase = -1; double bestMinDist = -1.0;
            for (int cand = 0; cand < n; cand++) {
                if (selected.contains(cand)) continue;
                double minDist = 1.0;
                for (int chosen : selected) {
                    if (chosen < 0 || chosen >= n) throw new IllegalStateException("chosen out of range: "+chosen);
                    double d = distanceMatrix[cand][chosen];
                    if (!Double.isFinite(d)) throw new IllegalStateException("distance non-finite at ["+cand+"]["+chosen+"]: "+d);
                    if (d < minDist) minDist = d;
                }
                if (!Double.isFinite(minDist)) throw new IllegalStateException("minDist non-finite for cand="+cand);
                if (minDist > bestMinDist) { bestMinDist = minDist; bestCase = cand; }
            }
            if (bestCase == -1) break;
            if (bestMinDist < 1.0 - corrThreshold) break; // 相似度过高则停止

            int idx = selected.size();
            if (idx >= maxCount) break;
            correlation[idx] = 1.0 - bestMinDist;
            if (!Double.isFinite(correlation[idx])) throw new IllegalStateException("correlation["+idx+"] non-finite");
            System.out.println("The max correlation is " + correlation[idx]);
            selected.add(bestCase);
        }

        // 补齐未写位置，确保无 NaN/Inf
        double fill = 1.0 - corrThreshold;
        if (!Double.isFinite(fill)) throw new IllegalStateException("fill non-finite from corrThreshold");
        for (int k = selected.size(); k < maxCount; k++) correlation[k] = fill;

        // 最终再验一次
        for (int k = 0; k < maxCount; k++)
            if (!Double.isFinite(correlation[k])) throw new IllegalStateException("final correlation["+k+"] non-finite");

        // 如果你维护了这个外部列表就保留，否则可删
        selectedCaseCorrelation.add(correlation);

        return selected;
    }

    static class DedupResult {
        final int[][] distinct;   // 去重后的行（保留第一次出现）
        final int[] keptIndices;  // 对应保留行在原数组中的下标
        DedupResult(int[][] d, int[] k) { this.distinct = d; this.keptIndices = k; }
    }

    public static DedupResult dedupRowsWithIndices(int[][] a) {
        if (a == null) return new DedupResult(new int[0][], new int[0]);
        Set<String> seen = new LinkedHashSet<>();  // 保持原顺序
        List<int[]> distinct = new ArrayList<>();
        List<Integer> kept = new ArrayList<>();

        for (int i = 0; i < a.length; i++) {
            int[] row = a[i];
            if (row == null) continue;                  // 如需保留null可自定义
            String key = Arrays.toString(row);          // 行内容签名
            if (seen.add(key)) {                        // 第一次出现
                distinct.add(row);
                kept.add(i);
            }
        }
        return new DedupResult(
                distinct.toArray(new int[0][]),
                kept.stream().mapToInt(Integer::intValue).toArray()
        );
    }


}
