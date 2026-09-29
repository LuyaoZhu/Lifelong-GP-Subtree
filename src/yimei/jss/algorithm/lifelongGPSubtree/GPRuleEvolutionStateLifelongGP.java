package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.algorithm.lifelongGP.GPRuleEvolutionStateLifelongGPV10N1;
import yimei.jss.algorithm.lifelongGP.SimpleKMedoids;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.AbstractRuleHelper;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.MultipleRuleEvaluationModel;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static yimei.jss.algorithm.lifelongGPSubtree.FeatureConstructionCrossoverPipelineSimplyProtect.LispContains;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * in each task, select individuals with good quality and diversity as the initial population of next task
 * in addition, the final population will be evaluated in 10 unseen instances to build surrogate
 * <p>
 * adaptively change the thresholds
 *
 * @author luyao
 */

public class GPRuleEvolutionStateLifelongGP extends GPRuleEvolutionStateLifelongGPV10N1 {

    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;

    public double thresholdSeq;
    public double thresholdRou;

    ArrayList<Integer> numSeqArchive = new ArrayList<>();
    ArrayList<Integer> numRouArchive = new ArrayList<>();

    ArrayList<int[]> seqArchiveDepth = new ArrayList<>();
    ArrayList<int[]> rouArchiveDepth = new ArrayList<>();

    public ArrayList<GPTree> seqFeatureArchive = new ArrayList<>();
    public ArrayList<GPTree> rouFeatureArchive = new ArrayList<>();      // Save selected high-level features ( will be updated in every generation)

    public ArrayList<Individual[][]> elitesEveryGen = new ArrayList<>();

    public int[] seqFeatureDepth;
    public int[] rouFeatureDepth;

    public int pcDistance;

    public static final String P_REPLICATIONS = "num-Rep";
    public int numRep;

    public static final String P_ELITES = "num-elites";
    public static int numElites;

    public double coefficientCurrentTask;

    public double seedingRatio;

    public static int generationPerTask;

    public ArrayList<ArrayList<ArrayList<GPTree>>> TaskSpecificBuildingBlocks = new ArrayList<>();
    public ArrayList<ArrayList<GPTree>> GeneralBuildingBlocks = new ArrayList<>();
    //define PC
    PhenoCharacterisation[] pc = new PhenoCharacterisation[2];

    RuleOptimizationProblem problem;
    public static ArrayList<Double> BlockOccurrenceRateOneGen = new ArrayList<>();

    ArrayList<Double> BlockOccurrenceRateEveryGen = new ArrayList<>();

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

        surrogateThreshold = parameters.getDoubleWithDefault(new Parameter("surrogateThreshold"), null, 5);

        seedingRatio = parameters.getDoubleWithDefault(new Parameter("seedingRatio"), null, 0.1);

        thresholdSeq = parameters.getDoubleWithDefault(new Parameter("thresholdSeq"), null, 0.6);
        thresholdRou = parameters.getDoubleWithDefault(new Parameter("thresholdRou"), null, 0.6);

        phenoCharacterisation[0] =
                SequencingPhenoCharacterisation.defaultPhenoCharacterisation();
        phenoCharacterisation[1] =
                RoutingPhenoCharacterisation.defaultPhenoCharacterisation();

    }

    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        problem = (RuleOptimizationProblem) evaluator.p_problem;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(generation / generationPerTask);

        //in each generation, calculate phenoCharacterisation
        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);

        //assign PC to individuals
        for (int s = 0; s < indsCharListsMultiTree.length; s++) {
            ((GPIndividual) population.subpops[0].individuals[s]).PC = indsCharListsMultiTree[s];
        }

        //record population diversity
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
//        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        // EVALUATION
        statistics.preEvaluationStatistics(this);

        refFit = calculateReferenceRuleFitness(simulation); //this is for normalising individuals' fitness

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

        double[] fitness = new double[population.subpops[0].individuals.length];
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            fitness[ind] = individual.fitness.fitness();
        }

        double minCurrentTask = Arrays.stream(fitness).min().getAsDouble();
        double judgeValue;
        if(minCurrentTask == 0){
            judgeValue = 1;
        } else {
            judgeValue = refFit;
        }

        //only preserve good individual (better than reference rule) with unique PC
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind].clone();
            if (individual.fitness.fitness() < judgeValue) {
                List<Integer> key = Arrays.stream(individual.PC)
                        .boxed()
                        .collect(Collectors.toList());
                PCIndividualMap.put(key, individual);
            }
        }
        archiveSampleNumber.add(PCIndividualMap.size());

        //then the fitness of one individual should be normalised raw fitness + normalised estimated fitness
        if (generation >= generationPerTask) {

            double[] thresholds = new double[surrogateSamples.size() + 1];

            double[][] estimatedFitness = new double[surrogateFitness.size() + 1][population.subpops[0].individuals.length];

            double[] minFitnessGap = new double[surrogateSamples.size() + 1];

            double[] meanPCDistance = calculateMeanPCDistance(surrogateSamples);
//            double[] mediumPCDistance = calculateMediumPCDistance(surrogateSamples);

            meanPCDistanceEveryGeneration.add(meanPCDistance);

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                estimatedFitness[estimatedFitness.length - 1][ind] = individual.fitness.fitness();
            }

            if(minCurrentTask == 0){ //means already normalised
                //based on the average fitness of top30% to determine whether we need to study this task only
                double meanFitnessCurrentTask = top30PercentMean(population.subpops[0].individuals);
//                if (generation%generationPerTask==1 && onlyCurrentTaskPhase) {
                if (checker.check(meanFitnessCurrentTask) && onlyCurrentTaskPhase) {
                    System.out.println("Converged at generation: " + generation);
                    switchGen.add(generation);
                    onlyCurrentTaskPhase = false;
                    ((surrogateClearingMultitreeEvaluatorV10N1)evaluator).nonIntermediatePop = false;
                }
            }

            if (onlyCurrentTaskPhase) {
                for (int a = 0; a < thresholds.length; a++) {
                    thresholds[a] = 0;
                    if (a == thresholds.length - 1) {
                        thresholds[a] = 1;
                    }
                }
            } else {

                if(thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1)[0] == 0){ //means this is the first time to calculate the thresholds
                    for (int t = 0; t < surrogateSamples.size(); t++) {
//                    estimatedFitness[t] = evaluatePopulation(this, surrogateSamples.get(t), surrogateFitness.get(t), 10);
                        estimatedFitness[t] = evaluatePopulationV1(this, surrogateSamples.get(t), surrogateFitness.get(t), surrogateThreshold);
//                System.out.println(Arrays.stream(estimatedFitness[t]).min().getAsDouble());
                        double minInPreviousTask = Arrays.stream(estimatedFitness[t]).min().getAsDouble();
                        for (int a = 0; a < estimatedFitness[t].length; a++) {
                            if (estimatedFitness[t][a] < Double.MAX_VALUE) {
                                estimatedFitness[t][a] = (estimatedFitness[t][a] - minInPreviousTask) / (1 - minInPreviousTask);
                            }
                        }
                        // this is to calculate the forgetting ratio and calculate the thresholds
                        if(estimatedFitness[t][estimatedFitness[t].length-1] < Double.MAX_VALUE){
                            double lowBound = Arrays.stream(surrogateFitness.get(t)).min().getAsDouble();
                            minFitnessGap[t] = (estimatedFitness[t][estimatedFitness[t].length-1] - lowBound) / (1 - lowBound);
                        } else {
                            minFitnessGap[t] = Double.MAX_VALUE;
                        }


                    }

                    //then based on fitness of top 50% individuals to calculate the thresholds
//                thresholds = calculateThresholdsV2(estimatedFitness,minFitness);
//                thresholds = calculateThresholdsV5(estimatedFitness, minFitnessGap);

//                thresholds = calculateThresholdsV3(estimatedFitness, minFitness);
//                thresholds = calculateThresholdsV4(minFitness);
                    thresholds = calculateThresholdsSimilarityV1(estimatedFitness);
//                    thresholds = calculateThresholdsSimilarityV2(estimatedFitness);

//                    double[] a = new double[thresholds.length];
//                    Arrays.fill(a, 1.0);
//                    thresholds = a;

                    ArrayList<Individual> elites = Arrays.stream(population.subpops[0].individuals)
                            .map(ind -> (GPIndividual) ind.clone())
                            .sorted((a, b) -> {
                                if (a.fitness.betterThan(b.fitness)) return -1;
                                else if (b.fitness.betterThan(a.fitness)) return 1;
                                else return 0;
                            })
                            .limit(Math.max(1, (int) Math.ceil(
                                    population.subpops[0].individuals.length * seedingRatio)))
                            .collect(Collectors.toCollection(ArrayList::new));

                    //task-specific building blocks collection
                    featureConstruction(elites);

                    //general building blocks collection
                    extractGeneralBlocks(TaskSpecificBuildingBlocks,generation/generationPerTask+1);

                } else {
                    thresholds = thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1);
                }

            }

            thresholdsEveryGeneration.add(thresholds);

            for (int t = 0; t < thresholds.length; t++) {
                System.out.println("The threshold of Task " + t + " is " + thresholds[t]);
            }

            double[] combinedFitness = new double[population.subpops[0].individuals.length];

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {

                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];

                double[] objective = new double[1];

                for (int t = 0; t < thresholds.length; t++) {
                    objective[0] += thresholds[t] * estimatedFitness[t][ind];
                    if (objective[0] >= Double.POSITIVE_INFINITY || objective[0] <= Double.NEGATIVE_INFINITY || Double.isNaN(objective[0])) {
//                        state.output.warning("Bad objective #" + ": " + objective[0] + ", setting to worst value for that objective.");
                        objective[0] = Double.MAX_VALUE;
                    }
                }
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                combinedFitness[ind] = objective[0];
            }
            //it's possible that all individuals' combined fitness in this population are Double.Max.
            // Thus, we still use the fitness in the current task to output the best individual and select parents, and preselection is based on fitness in all tasks.
            double min = Arrays.stream(combinedFitness).min().getAsDouble();
            if(min > 10) {
                for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                    double[] objective = new double[1];
                    objective[0] = estimatedFitness[estimatedFitness.length-1][ind];
                    ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                }
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

            if (PCs.length >= 500) {

                //then do k-means to select 500 individuals
                SimpleKMedoids.Result r = SimpleKMedoids.fit(PCs, population.subpops[0].individuals.length, 20250101L, 68);

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
            double[] refFitness = new double[simulationsPerTask];

            for (int s = 0; s < simulationsPerTask; s++) {

                simulation.reset(608 + s);
                refFitness[s] = calculateReferenceRuleFitness(simulation);

                // （可选但推荐）避免 refFit1[1] 为 0 导致除零
                if (refFitness[s] == 0.0) {
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
                    fitnessOneSurrogate[ind] += ObjValue / refFitness[s];
                }
            }

            for (int ind = 0; ind < fitnessOneSurrogate.length; ind++) {
                fitnessOneSurrogate[ind] /= simulationsPerTask;
            }

//            System.out.println(Arrays.stream(fitnessOneSurrogate).min().getAsDouble());

            surrogateFitness.add(fitnessOneSurrogate);
            surrogateSamples.add(PC);   // ✅ 这里 PCs[ind] 与 QDIndividuals.get(ind) 严格一一对应

            //investigate if the large distance, the bad fitness
//            investigateDistanceFitnessRelation(PC,fitnessOneSurrogate);
//            investigateTopKLocalRelation(PC,fitnessOneSurrogate,100);

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

            writeTaskSpecificBlockDepth(TaskSpecificBuildingBlocks);
            writeGeneralBlockDepth(GeneralBuildingBlocks);

            writeDiversityToFile(entropyDiversity);

            writeBlockOccurrenceRateCrossoverToFile(BlockOccurrenceRateEveryGen);

//            writeAccuracyToFile(MSEGen, SpearmanCorrelationGen, SamePCNum, AveragePCDistance);

            writeArchiveSampleNumToFile(archiveSampleNumber);

            writeThresholdsToFile(thresholdsEveryGeneration);

            writeMeanPCDistanceToFile(meanPCDistanceEveryGeneration);

            writeThresholdSwitchGenToFile(switchGen);

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


            savedTopIndividuals = new ArrayList<>();
            PopulationUtils.sort(population);

            for (int i = 0; i < population.subpops[0].individuals.length * seedingRatio; i++) {
                savedTopIndividuals.add(population.subpops[0].individuals[i]);
            }

            if(generation == generationPerTask - 1){  //only for the first task
//                TaskSpecific_BuildingBlocks_Collection(savedTopIndividuals,decisionSituations);
                featureConstruction(savedTopIndividuals);
            }

            population.clear();

            population = initializer.initialPopulation(this, 0);

            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                for (int replace = 0; replace < savedTopIndividuals.size(); replace++) {
                    population.subpops[sub].individuals[replace] = savedTopIndividuals.get(replace);
                }
            }

            savedIndividuals.clear();
            PCIndividualMap.clear();

            checker = new ConvergenceChecker();

            onlyCurrentTaskPhase = true;

        } else {

            if(onlyCurrentTaskPhase) {
                population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N
            } else {
                population = preselection(); //aims to select the individuals with fitness<Double.MaxValue
            }

        }

        double avg = BlockOccurrenceRateOneGen.stream()
                .mapToDouble(Double::doubleValue)
                .average()
                .orElse(0.0);

        BlockOccurrenceRateOneGen.clear();

        BlockOccurrenceRateEveryGen.add(avg);


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

    public void normalisePopulation(double min, double refFit) {

        ArrayList<Double> fitness = new ArrayList<>();

        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            if (individual.fitness.fitness() < Double.MAX_VALUE) {
                double[] objective = new double[1];
                objective[0] = (individual.fitness.fitness() - min) / (refFit - min);
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                fitness.add(individual.fitness.fitness());
            }
        }

    }

    public void featureConstruction(ArrayList<Individual> elitesEveryGen) {

        seqFeatureArchive.clear();
        rouFeatureArchive.clear();

        PhenoCharacterisation[] pcTask = new PhenoCharacterisation[2];

        DynamicSimulation simulation = (DynamicSimulation) ((MultipleRuleEvaluationModel)problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(generation/generationPerTask);

        // then calculate phenoCharacterisation and obtain some decision points
        pcTask[0] = SequencingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation);
        pcTask[1] = RoutingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation);


        //delete the features that are not related to the new elites in the archive
//        if (this.generation > 0) {
//            deleteUnimportantFeature(pc);
//        }
//        if(this.generation%5 == 0)
//            System.out.println("1");
        //select features into the archive

        for (int sub = 0; sub < 2; sub++) {
            for (int ind = 0; ind < elitesEveryGen.size(); ind++) {
                if (sub == 0) {
                    decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pcTask[sub]).decisionSituations;
                } else {
                    decisionSituationsRouting = ((RoutingPhenoCharacterisation) pcTask[sub]).decisionSituations;
                }
                AbstractRule rule = null;
                GPTree tree = ((GPIndividual) (elitesEveryGen.get(ind))).trees[sub];
                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                double[] score = new double[nonterminals];
                double[] frequency = new double[nonterminals];
                if (nonterminals > 1) {
                    if (sub == 0) {
                        double[][] score_matrix = new double[nonterminals][];
                        for (int i = 0; i < decisionSituationsSequencing.size(); i++) {

                            int decisionSize = decisionSituationsSequencing.get(i).getQueue().size();

                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                                score_matrix[numSubtree] = new double[decisionSize];

                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.SEQUENCING, nodeToTree);

                                SequencingDecisionSituation situation = decisionSituationsSequencing.get(i);
                                List<OperationOption> queue = situation.getQueue();
                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    queue.get(candiateNum).setPriority(rule.priority(queue.get(candiateNum), situation.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }

/*                                if(i == 0){ //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for(int b=0; b<elitesEveryGen.get(this.generation)[sub].length; b++){
                                        GPTree elite = ((GPIndividual)elitesEveryGen.get(this.generation)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num=0; num < subtreeNum; num++){
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if(eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
                                }*/

                            }

                            //2019.10.24 use the rank rather than the priority value directly
                            //========================================start==========================================
                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }
                            //========================================end============================================

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }

                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsSequencing.size();
                        }

                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {

                            if (Math.abs(score[numScore]) >= thresholdSeq) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(sequencingFeature.size() < 30){
//                                    sequencingFeature.add(nodeToTree);
//                                }
                                seqFeatureArchive.add(nodeToTree);
                            }
                        }

                    } else {

                        double[][] score_matrix = new double[nonterminals][];
                        //for (int i = 0; i < 1; i++) {
                        for (int i = 0; i < decisionSituationsRouting.size(); i++) {

                            int decisionSize = decisionSituationsRouting.get(i).getQueue().size();

                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {

                                score_matrix[numSubtree] = new double[decisionSize];

                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.ROUTING, nodeToTree);

                                RoutingDecisionSituation situation = decisionSituationsRouting.get(i);
                                List<OperationOption> queue = situation.getQueue();

                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    OperationOption operationOption = queue.get(candiateNum);
                                    queue.get(candiateNum).setPriority(rule.priority(operationOption, operationOption.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }

/*                                if (i == 0) { //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for (int b = 0; b < elitesEveryGen.get(this.generation)[sub].length; b++) {
                                        GPTree elite = ((GPIndividual) elitesEveryGen.get(this.generation)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num = 0; num < subtreeNum; num++) {
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if (eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
                                }*/


                            }

                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }
                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsRouting.size();
                        }
                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {

                            if (Math.abs(score[numScore]) >= thresholdRou) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(routingFeature.size() < 30){
//                                    routingFeature.add(nodeToTree);
//                                }
                                rouFeatureArchive.add(nodeToTree);
                            }
                        }

                    }
                }
            }
        }


        //delete duplicate features
        deleteDuplicateFeatures();

//        if (this.generation < 30)
        deleteBigRedundantFeatures();
//        else {
//            deleteDuplicateFeatures();
//        }

        // record the depth of the features
       /* seqFeatureDepth = new int[6];
        rouFeatureDepth = new int[6];
        int sumSeqFeature = 0;

        int gen = generation/generationPerTask;
        if (!(seqFeatureArchive.size() == 0)) {
//            System.out.println("the seqFeature number in generation " + (gen + 1) + " is " + seqFeatureArchive.get(gen).size());
            for (int i = 0; i < seqFeatureArchive.size(); i++) {
//                System.out.println(seqFeatureArchive.get(gen).get(i).child.makeLispTree());
                sumSeqFeature++;
                if (seqFeatureArchive.get(i).child.depth() == 2)
                    seqFeatureDepth[0]++;
                if (seqFeatureArchive.get(i).child.depth() == 3)
                    seqFeatureDepth[1]++;
                if (seqFeatureArchive.get(i).child.depth() == 4)
                    seqFeatureDepth[2]++;
                if (seqFeatureArchive.get(i).child.depth() == 5)
                    seqFeatureDepth[3]++;
                if (seqFeatureArchive.get(i).child.depth() == 6)
                    seqFeatureDepth[4]++;
                if (seqFeatureArchive.get(i).child.depth() == 7)
                    seqFeatureDepth[5]++;
            }
        }

        seqArchiveDepth.add(seqFeatureDepth);
        numSeqArchive.add(sumSeqFeature);
        System.out.println("The number of feature in sequencing Archive is  " + sumSeqFeature);

        int sumRouFeature = 0;

        if (!(rouFeatureArchive.size() == 0)) {
//            System.out.println("the rouFeature number in generation " + (gen + 1) + " is " + rouFeatureArchive.get(gen).size());
            for (int i = 0; i < rouFeatureArchive.size(); i++) {
//                System.out.println(rouFeatureArchive.get(gen).get(i).child.makeLispTree());
                sumRouFeature++;
                if (rouFeatureArchive.get(i).child.depth() == 2)
                    rouFeatureDepth[0]++;
                if (rouFeatureArchive.get(i).child.depth() == 3)
                    rouFeatureDepth[1]++;
                if (rouFeatureArchive.get(i).child.depth() == 4)
                    rouFeatureDepth[2]++;
                if (rouFeatureArchive.get(i).child.depth() == 5)
                    rouFeatureDepth[3]++;
                if (rouFeatureArchive.get(i).child.depth() == 6)
                    rouFeatureDepth[4]++;
                if (rouFeatureArchive.get(i).child.depth() == 7)
                    rouFeatureDepth[5]++;
            }
        }

        rouArchiveDepth.add(rouFeatureDepth);
        numRouArchive.add(sumRouFeature);
        System.out.println("The number of feature in routing Archive is  " + sumRouFeature);*/

        //now preserve them
        ArrayList<ArrayList<GPTree>> BlocksCurrentTask = new ArrayList<>();

        ArrayList<GPTree> seqBlocksThisTask = new ArrayList<>();
        for (GPTree tree : seqFeatureArchive) {
            seqBlocksThisTask.add((GPTree) tree.clone());
        }

        ArrayList<GPTree> rouBlocksThisTask = new ArrayList<>();
        for (GPTree tree : rouFeatureArchive) {
            rouBlocksThisTask.add((GPTree) tree.clone());
        }

        BlocksCurrentTask.add(seqBlocksThisTask);
        BlocksCurrentTask.add(rouBlocksThisTask);

        TaskSpecificBuildingBlocks.add(BlocksCurrentTask);
    }

    private void deleteBigRedundantFeatures() {

        // sequencing
        for (int a = 0; a < seqFeatureArchive.size(); a++) {

            GPTree base = seqFeatureArchive.get(a);

            for (int b = a + 1; b < seqFeatureArchive.size(); b++) {

                GPTree count = seqFeatureArchive.get(b);

                // base 包含 count → 删除 base（大的）
                if (LispContains(base.child, count.child)) {
                    seqFeatureArchive.remove(a);
                    a--;
                    break;
                }

                // count 包含 base → 删除 count（大的）
                else if (LispContains(count.child, base.child)) {
                    seqFeatureArchive.remove(b);
                    b--;
                }
            }
        }

        // routing
        for (int a = 0; a < rouFeatureArchive.size(); a++) {

            GPTree base = rouFeatureArchive.get(a);

            for (int b = a + 1; b < rouFeatureArchive.size(); b++) {

                GPTree count = rouFeatureArchive.get(b);

                if (LispContains(base.child, count.child)) {
                    rouFeatureArchive.remove(a);
                    a--;
                    break;
                }

                else if (LispContains(count.child, base.child)) {
                    rouFeatureArchive.remove(b);
                    b--;
                }
            }
        }
    }

    private void deleteDuplicateFeatures() {

        // sequencing
        for (int a = 0; a < seqFeatureArchive.size(); a++) {

            GPTree base = seqFeatureArchive.get(a);

            for (int b = a + 1; b < seqFeatureArchive.size(); b++) {

                GPTree count = seqFeatureArchive.get(b);

                if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                    seqFeatureArchive.remove(b);
                    b--;
                }
            }
        }

        // routing
        for (int a = 0; a < rouFeatureArchive.size(); a++) {

            GPTree base = rouFeatureArchive.get(a);

            for (int b = a + 1; b < rouFeatureArchive.size(); b++) {

                GPTree count = rouFeatureArchive.get(b);

                if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                    rouFeatureArchive.remove(b);
                    b--;
                }
            }
        }
    }

    public void extractGeneralBlocks(
            ArrayList<ArrayList<ArrayList<GPTree>>> TaskSpecificBuildingBlocks,
            int threshold   // 推荐 >= 2
    ) {

        // 初始化
        GeneralBuildingBlocks.clear();
        GeneralBuildingBlocks.add(new ArrayList<>()); // seq
        GeneralBuildingBlocks.add(new ArrayList<>()); // rou

        int numTrees = 2; // seq + rou

        for (int treeID = 0; treeID < numTrees; treeID++) {

            HashMap<String, Integer> freq = new HashMap<>();
            HashMap<String, GPTree> map = new HashMap<>();

            // 遍历所有 task
            for (int t = 0; t < TaskSpecificBuildingBlocks.size(); t++) {

                ArrayList<ArrayList<GPTree>> taskBlocks =
                        TaskSpecificBuildingBlocks.get(t);

                if (taskBlocks == null || taskBlocks.size() <= treeID) continue;

                // 防止同一个 task 重复计数
                HashSet<String> seenInTask = new HashSet<>();

                for (GPTree tree : taskBlocks.get(treeID)) {

                    String key = tree.child.makeLispTree();

                    if (!seenInTask.contains(key)) {

                        freq.put(key, freq.getOrDefault(key, 0) + 1);
                        map.put(key, tree);

                        seenInTask.add(key);
                    }
                }
            }

            // 选出 general blocks
            for (String key : freq.keySet()) {
                if (freq.get(key) >= threshold) {
                    GeneralBuildingBlocks.get(treeID).add(map.get(key));
                }
            }
        }
    }

    public void writeArchiveInformationToFile(ArrayList<Integer> numSeqArchive, ArrayList<Integer> numRouArchive, ArrayList<int[]> seqArchiveDepth, ArrayList<int[]> rouArchiveDepth) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File( out_dir + "/job." + jobSeed +  ".archiveInformation.csv");
//        File weightFile = new File( "/job." + jobSeed + ".archiveInformation.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,numSeqArchive,numRouArchive,seqDepth2,seqDepth3,seqDepth4,seqDepth5,seqDepth6,seqDepth7,rouDepth2,rouDepth3,rouDepth4,rouDepth5,rouDepth6,rouDepth7");
            writer.newLine();
            for (int i = 0; i < numSeqArchive.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + numSeqArchive.get(i) + ", " + numRouArchive.get(i) + ", "
                        + seqArchiveDepth.get(i)[0] + ", " + seqArchiveDepth.get(i)[1] + ", "+ seqArchiveDepth.get(i)[2] + ", "
                        + seqArchiveDepth.get(i)[3] + ", " + seqArchiveDepth.get(i)[4] + ", "+ seqArchiveDepth.get(i)[5] + ", "
                        + rouArchiveDepth.get(i)[0] + ", " + rouArchiveDepth.get(i)[1] + ", "+ rouArchiveDepth.get(i)[2] + ", "
                        + rouArchiveDepth.get(i)[3] + ", " + rouArchiveDepth.get(i)[4] + ", "+ rouArchiveDepth.get(i)[5] +"\n");
            }
            numSeqArchive.clear();
            numRouArchive.clear();
            seqArchiveDepth.clear();
            rouArchiveDepth.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeBlockOccurrenceRateCrossoverToFile(ArrayList<Double> BlockOccurrenceRateEveryGen) {

        File weightFile = new File(out_dir + "/job." + jobSeed + ".blockOccurrenceRateCrossover.csv");

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // header
            writer.write("Gen,OccurrenceRate");
            writer.newLine();

            // data
            for (int i = 0; i < BlockOccurrenceRateEveryGen.size(); i++) {

                double value = BlockOccurrenceRateEveryGen.get(i);

                String line = i + "," + String.format("%.4f", value);

                writer.write(line);
                writer.newLine();
            }

            BlockOccurrenceRateEveryGen.clear();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void writeTaskSpecificBlockDepth(
            ArrayList<ArrayList<ArrayList<GPTree>>> TaskSpecificBuildingBlocks) {

        File file = new File(out_dir + "/job." + jobSeed + ".taskSpecificBlockDepth.csv");

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {

            writer.write("Task,Type,Total,Depth2,Depth3,Depth4,Depth5,Depth6,Depth7");
            writer.newLine();

            for (int task = 0; task < TaskSpecificBuildingBlocks.size(); task++) {

                ArrayList<ArrayList<GPTree>> oneTask =
                        TaskSpecificBuildingBlocks.get(task);

                if (oneTask == null) continue;

                for (int treeID = 0; treeID < 2; treeID++) {

                    String type = (treeID == 0) ? "SEQ" : "ROU";

                    int total = 0;
                    int[] depthCount = new int[6];

                    if (oneTask.size() > treeID && oneTask.get(treeID) != null) {

                        ArrayList<GPTree> blocks = oneTask.get(treeID);
                        total = blocks.size();

                        for (GPTree tree : blocks) {

                            if (tree == null || tree.child == null) continue;

                            int depth = tree.child.depth();

                            if (depth >= 2 && depth <= 7) {
                                depthCount[depth - 2]++;
                            }
                        }
                    }

                    writer.write(task + ","
                            + type + ","
                            + total + ","
                            + depthCount[0] + ","
                            + depthCount[1] + ","
                            + depthCount[2] + ","
                            + depthCount[3] + ","
                            + depthCount[4] + ","
                            + depthCount[5]);

                    writer.newLine();
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeGeneralBlockDepth(
            ArrayList<ArrayList<GPTree>> GeneralBuildingBlocks) {

        File file = new File(out_dir + "/job." + jobSeed + ".generalBlockDepth.csv");

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {

            writer.write("Type,Total,Depth2,Depth3,Depth4,Depth5,Depth6,Depth7");
            writer.newLine();

            for (int treeID = 0; treeID < 2; treeID++) {

                String type = (treeID == 0) ? "SEQ" : "ROU";

                int total = 0;
                int[] depthCount = new int[6];

                if (GeneralBuildingBlocks != null
                        && GeneralBuildingBlocks.size() > treeID
                        && GeneralBuildingBlocks.get(treeID) != null) {

                    ArrayList<GPTree> blocks = GeneralBuildingBlocks.get(treeID);
                    total = blocks.size();

                    for (GPTree tree : blocks) {

                        if (tree == null || tree.child == null) continue;

                        int depth = tree.child.depth();

                        if (depth >= 2 && depth <= 7) {
                            depthCount[depth - 2]++;
                        }
                    }
                }

                writer.write(type + ","
                        + total + ","
                        + depthCount[0] + ","
                        + depthCount[1] + ","
                        + depthCount[2] + ","
                        + depthCount[3] + ","
                        + depthCount[4] + ","
                        + depthCount[5]);

                writer.newLine();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}



