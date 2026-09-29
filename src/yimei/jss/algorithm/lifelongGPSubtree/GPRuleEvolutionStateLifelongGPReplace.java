package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.SortComparatorL;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import smile.clustering.HierarchicalClustering;
import smile.clustering.linkage.CompleteLinkage;
import yimei.jss.algorithm.featureConstruction.FeatureConstructionCrossoverPipeline;
import yimei.jss.algorithm.lifelongGP.GPRuleEvolutionStateLifelongGPV10N1;
import yimei.jss.algorithm.lifelongGP.SimpleKMedoids;
import yimei.jss.algorithm.surrogateCaseLS.SilhouetteScore;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.*;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.AbstractRuleHelper;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.SPT;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
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

import static yimei.jss.algorithm.surrogateAccuracy.surrogateClearingMultitreeEvaluatorV1.SamePCNum;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * parent selection
 * crossover point
 *
 * @author luyao
 */

public class GPRuleEvolutionStateLifelongGPReplace extends GPRuleEvolutionStateLifelongGPV10N1 {

    ArrayList<GPIndividual> historicalGeneralists = new ArrayList<>(); //preserve some elites after training one task
    ArrayList<GPIndividual> currentSpecialists = new ArrayList<>();  //preserve some elites after the first stage of training
    ArrayList<GPIndividual> currentGeneralists = new ArrayList<>();
    public ArrayList<int[]> armSelectionHistory = new ArrayList<>();
    public ArrayList<double[]> armRewardHistory = new ArrayList<>();
    public ArrayList<double[]> armSurvivalHistory = new ArrayList<>();

    private double[] hgSum = new double[3]; private int hgCount = 0;
    private double[] csSum = new double[3]; private int csCount = 0;
    private double[] cgSum = new double[3]; private int cgCount = 0;

    public ArrayList<double[]> HGSubtreeImportance = new ArrayList<>();
    public ArrayList<double[]> CSSubtreeImportance = new ArrayList<>();
    public ArrayList<double[]> CGSubtreeImportance = new ArrayList<>();


    public boolean firstStageSwitchGeneration = false;

    public static double UCBExplorationConstant;
    public static double currentRewardWeight;
    public double individualGroupRatio;
    // =====================================
    // importance statistics
    // [tree0 donor,
    //  tree0 replace,
    //  tree1 donor,
    //  tree1 replace]
    // =====================================

    public double[] importanceSum = new double[4];

    public int[] importanceCount = new int[4];

    // =====================================
    // depth statistics
    // =====================================

    public double[] depthSum = new double[4];

    public int[] depthCount = new int[4];

    // =====================================
    // history
    // =====================================

    public ArrayList<double[]> importanceHistory = new ArrayList<>();

    public ArrayList<double[]> depthHistory = new ArrayList<>();

    public ArrayList<Integer> invalidIndsNumber = new ArrayList<>();

    public ArrayList<Double> surrogateThresholds = new ArrayList<>();

    public void resetGenerationStatistics() {

        importanceSum = new double[4];
        importanceCount = new int[4];

        depthSum = new double[4];
        depthCount = new int[4];
    }

    public void setup(EvolutionState state, Parameter base) {
        super.setup(state, base);
        fixedArm = parameters.getIntWithDefault(new Parameter("fixedArm"), null, 1);
        UCBExplorationConstant = parameters.getDoubleWithDefault(new Parameter("UCBExplorationConstant"), null, 0.5);
        currentRewardWeight = parameters.getDoubleWithDefault(new Parameter("currentRewardWeight"), null, 0.8);
        individualGroupRatio = parameters.getDoubleWithDefault(new Parameter("individualGroupRatio"), null, 0.1);
    }
    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        Arrays.fill(AdaptiveParentSelection.generationSelectionCounts, 0);
        resetGenerationStatistics();

        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(generation / generationPerTask);

        //collect different decision points for different tasks
        if (generation % generationPerTask == 0) {
            phenoCharacterisation = new PhenoCharacterisation[2];
            phenoCharacterisation[0] = SequencingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation);
            phenoCharacterisation[1] = RoutingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation, ((SequencingPhenoCharacterisation) phenoCharacterisation[0]).decisionSituations.size());
            decisionSituationsEachTask.add(phenoCharacterisation);
        }

        //record population diversity
        double[] diversityValue = new double[decisionSituationsEachTask.size()];

        for (int s = 0; s < population.subpops[0].individuals.length; s++) {
            ((GPIndividual) population.subpops[0].individuals[s]).PCs = new ArrayList<>();
        }

        for (int i = 0; i < decisionSituationsEachTask.size(); i++) {

            PhenoCharacterisation[] pc = decisionSituationsEachTask.get(i);

            //in each generation, calculate phenoCharacterisation
            int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, pc);

            //assign PC to individuals
            for (int s = 0; s < population.subpops[0].individuals.length; s++) {
                ((GPIndividual) population.subpops[0].individuals[s]).PCs.add(indsCharListsMultiTree[s]);
            }

            diversityValue[i] = PopulationUtils.entropy(indsCharListsMultiTree);
//        System.out.println(diversityValue[0]);
        }
        double[] diversity = new double[population.subpops.length];
        diversity[0] = Arrays.stream(diversityValue).average().getAsDouble();
        entropyDiversity.add(diversity);

        // EVALUATION
        statistics.preEvaluationStatistics(this);

        refFit = calculateReferenceRuleFitness(simulation); //this is for normalising individuals' fitness

        evaluator.evaluatePopulation(this);  // here, after this we evaluate the population

        //only preserve good individual (better than reference rule) with unique PC
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind].clone();
            if (individual.fitness.fitness() < refFit) {
                List<Integer> key = Arrays.stream(individual.PCs.get(individual.PCs.size() - 1))
                        .boxed()
                        .collect(Collectors.toList());
                PCIndividualMap.put(key, individual);
            }
        }
        archiveSampleNumber.add(PCIndividualMap.size());

        //find the min fitness value of the population and normalise the population
        double[] fitness = new double[population.subpops[0].individuals.length];
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            fitness[ind] = individual.fitness.fitness();
        }

        double minCurrentTask = Arrays.stream(fitness).min().getAsDouble();
        if(minCurrentTask < refFit) {
            normalisePopulation(minCurrentTask, refFit);
        }

        //then the fitness of one individual should be normalised raw fitness + normalised estimated fitness
        if (generation >= generationPerTask) {

            double[] thresholds = new double[surrogateSamples.size() + 1];

            double[][] estimatedFitness = new double[surrogateFitness.size() + 1][population.subpops[0].individuals.length];

            double[] meanPCDistance = calculateMeanPCDistance(surrogateSamples);

            meanPCDistanceEveryGeneration.add(meanPCDistance);

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                estimatedFitness[estimatedFitness.length - 1][ind] = individual.fitness.fitness();
            }

            //based on the average fitness of top30% to determine whether we need to study this task only
            double meanFitnessCurrentTask = top30PercentMean(population.subpops[0].individuals);
//                if (generation%generationPerTask==1 && onlyCurrentTaskPhase) {
            if (checker.check(meanFitnessCurrentTask) && onlyCurrentTaskPhase) {
                System.out.println("Converged at generation: " + generation);
                switchGen.add(generation);
                onlyCurrentTaskPhase = false;
                ((surrogateClearingMultitreeEvaluatorV10N1) evaluator).nonIntermediatePop = false;
                firstStageSwitchGeneration = true;
            }


            if (onlyCurrentTaskPhase) {
                for (int a = 0; a < thresholds.length; a++) {
                    thresholds[a] = 0;
                    if (a == thresholds.length - 1) {
                        thresholds[a] = 1;
                    }
                }
            } else {

                if (thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1)[0] == 0) { //means this is the first time to calculate the thresholds

                    for (int t = 0; t < surrogateSamples.size(); t++) {

//                        estimatedFitness[t] = evaluatePopulationV1(this, surrogateSamples.get(t), surrogateFitness.get(t), surrogateThresholds.get(t), t);
                        estimatedFitness[t] = evaluatePopulationV1(this, surrogateSamples.get(t), surrogateFitness.get(t), surrogateThreshold, t);

                        double minInPreviousTask = Arrays.stream(estimatedFitness[t]).min().getAsDouble();
                        for (int a = 0; a < estimatedFitness[t].length; a++) {
                            if (estimatedFitness[t][a] < Double.MAX_VALUE) {
                                estimatedFitness[t][a] = ( estimatedFitness[t][a] - minInPreviousTask ) / (1- minInPreviousTask);
                            }
                        }
                    }

                    thresholds = calculateThresholdsSimilarityV1(estimatedFitness);
                    currentSpecialists = selectQDIndividuals(population, individualGroupRatio);

                } else {
                    thresholds = thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1);
                }
            }

            thresholdsEveryGeneration.add(thresholds);

            for (int t = 0; t < thresholds.length; t++) {
                System.out.println("The threshold of Task " + t + " is " + thresholds[t]);
            }

/*            double[] combinedFitness = new double[population.subpops[0].individuals.length];

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
            if (min > 10) {
                for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                    double[] objective = new double[1];
                    objective[0] = estimatedFitness[estimatedFitness.length - 1][ind];
                    ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                }
            }*/

            if (!onlyCurrentTaskPhase && !firstStageSwitchGeneration) {
                currentGeneralists = new ArrayList<>();
                for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                    currentGeneralists.add(individual);
                }
            }

        }

        statistics.postEvaluationStatistics(this);

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

            if (PCs.length >= population.subpops[0].individuals.length) {

                //then do k-means to select 500 individuals
                SimpleKMedoids.Result r = SimpleKMedoids.fit(PCs, population.subpops[0].individuals.length, 20250101L, 68);
                surrogateThresholds.add( Math.ceil(SimpleKMedoids.printAndGetAverageDistance(PCs, r)) );

                PC = new int[population.subpops[0].individuals.length][PCs[0].length];

                for (int a = 0; a < r.medoids.length; a++) {

                    int m = r.medoids[a];
                    savedIndividuals.add(QDIndividuals.get(m));   // 每个簇中心个体

                    for (int i = 0; i < PCs[m].length; i++) {
                        PC[a][i] = (int) PCs[m][i];
                    }
                }
            } else {
                surrogateThresholds.add(2.0);

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
                    fitnessOneSurrogate[ind] += ObjValue/refFitness[s];
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

        //After evaluate all individuals, we record feature information about 5 elites
        // SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }


        // SHOULD WE QUIT?
        if (generation == numGenerations - 1) {

            writeDiversityToFile(entropyDiversity);

            writeArmSelectionToFile(armSelectionHistory);

            writeArmRewardToFile(armRewardHistory);

            writeArmSurvivalToFile(armSurvivalHistory);

            writeImportanceToFile(importanceHistory);

            writeDepthToFile(depthHistory);

//            writeAccuracyToFile(MSEGen, SpearmanCorrelationGen, SamePCNum, AveragePCDistance);

            writeArchiveSampleNumToFile(archiveSampleNumber);

            writeThresholdsToFile(thresholdsEveryGeneration);

            writeMeanPCDistanceToFile(meanPCDistanceEveryGeneration);

            writeThresholdSwitchGenToFile(switchGen);

            writeInvalidIndsNumGenToFile(invalidIndsNumber);

            writeSurrogateThresholdsToFile(surrogateThresholds);

            writeSubtreeImportanceToFile(HGSubtreeImportance, CSSubtreeImportance, CGSubtreeImportance);

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

//            PopulationUtils.sort(population);
//
//            for (int i = 0; i < population.subpops[0].individuals.length * seedingRatio; i++) {
//                historicalGeneralists.add(population.subpops[0].individuals[i]);
//            }

            historicalGeneralists = selectQDIndividuals(population, individualGroupRatio);

            population.clear();
            population = initializer.initialPopulation(this, 0);
//            population = breeder.breedPopulation(this);  //only use mutation

            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                for (int replace = 0; replace < historicalGeneralists.size(); replace++) {
                    population.subpops[sub].individuals[replace] = historicalGeneralists.get(replace);
                }
            }

            PCIndividualMap.clear();

            checker = new ConvergenceChecker();

            onlyCurrentTaskPhase = true;

            AdaptiveParentSelection.resetMAB();

            TaskwiseImportanceNodeSelector.importanceCache.clear();

            invalidIndsNumber.add(-100);

        } else {

            if (onlyCurrentTaskPhase) {
                population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N

                // =====================================
                // fill zero for inactive generations
                // =====================================

                armSelectionHistory.add(new int[]{0, 0, 0, 0});

                armRewardHistory.add(new double[]{0, 0, 0, 0});

                armSurvivalHistory.add(new double[]{0, 0, 0, 0});

                importanceHistory.add(new double[]{0, 0, 0, 0});

                depthHistory.add(new double[]{0, 0, 0, 0});

                invalidIndsNumber.add(-100);

            } else {
                population = preselection(); //aims to select the individuals with fitness<Double.MaxValue

                // =========================================
                // MAB reward:
                // survived lifelong preselection
                // =========================================
                updateMABRewardBySurvivors(population);

                if(firstStageSwitchGeneration == true)
                    AdaptiveParentSelection.resetMAB();

                firstStageSwitchGeneration = false;

                AdaptiveParentSelection.printStatistics();

                armSelectionHistory.add(AdaptiveParentSelection.generationSelectionCounts.clone());
                armRewardHistory.add(AdaptiveParentSelection.getRewards());
                armSurvivalHistory.add(AdaptiveParentSelection.getCurrentGenerationSurvival());

                // =====================================
                // importance averages
                // =====================================

                double[] importanceStats = new double[4];

                for (int i = 0; i < 4; i++) {

                    importanceStats[i] =
                            importanceCount[i] == 0 ? 0 : importanceSum[i] / importanceCount[i];
                }

                importanceHistory.add(importanceStats);

                // =====================================
                // depth averages
                // =====================================

                double[] depthStats = new double[4];

                for (int i = 0; i < 4; i++) {

                    depthStats[i] = depthCount[i] == 0 ? 0 : depthSum[i] / depthCount[i];
                }

                depthHistory.add(depthStats);

            }

        }

        // 每代结束时，无条件调用——不管这代有没有触发过 adaptive crossover
        finalizeSubtreeImportanceStats();

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

    public void addSubtreeImportanceSample(String groupLabel, double previous, double current, double all) {

        double[] sum;
        int count;

        switch (groupLabel) {
            case "HG": sum = hgSum; break;
            case "CS": sum = csSum; break;
            case "CG": sum = cgSum; break;
            default: return;
        }

        sum[0] += previous;
        sum[1] += current;
        sum[2] += all;

        if (groupLabel.equals("HG")) hgCount++;
        else if (groupLabel.equals("CS")) csCount++;
        else cgCount++;
    }

    public void finalizeSubtreeImportanceStats() {

        HGSubtreeImportance.add(finalizeGroup(hgSum, hgCount));
        CSSubtreeImportance.add(finalizeGroup(csSum, csCount));
        CGSubtreeImportance.add(finalizeGroup(cgSum, cgCount));

        hgSum = new double[3]; hgCount = 0;
        csSum = new double[3]; csCount = 0;
        cgSum = new double[3]; cgCount = 0;
    }

    private double[] finalizeGroup(double[] sum, int count) {
        return (count == 0) ? new double[3] : new double[]{ sum[0]/count, sum[1]/count, sum[2]/count };
    }

    public ArrayList<GPIndividual> selectQDIndividuals(
            Population population,
            double ratio) {

        ArrayList<GPIndividual> qdIndividuals = new ArrayList<>();

        // sort by fitness (best first)
        PopulationUtils.sort(population.subpops[0].individuals);

        HashSet<String> behaviourSet = new HashSet<>();

        int targetSize = (int) (
                population.subpops[0].individuals.length * ratio);

        for (int i = 0;
             i < population.subpops[0].individuals.length
                     && qdIndividuals.size() < targetSize;
             i++) {

            GPIndividual ind =
                    (GPIndividual) population.subpops[0].individuals[i];

            // build behavioural signature
            StringBuilder sb = new StringBuilder();

            for (double v : ind.PCs.get(ind.PCs.size() - 1)) {
                sb.append(String.format("%.3f", v)).append(",");
            }

            String signature = sb.toString();

            // preserve only behaviourally unique individuals
            if (!behaviourSet.contains(signature)) {

                behaviourSet.add(signature);

                qdIndividuals.add(ind);
            }
        }

        return qdIndividuals;
    }

    public void updateMABRewardBySurvivors(
            Population population) {

        int armNum =
                AdaptiveParentSelection.ARM_NUM;

        int[] generated =
                new int[armNum];

        int[] survived =
                new int[armNum];

        // =====================================
        // count generated offspring
        // =====================================

        for (int s = 0;
             s < population.subpops.length;
             s++) {

            for (int i = 0;
                 i < population.subpops[s].individuals.length;
                 i++) {

                GPIndividual ind =
                        (GPIndividual)
                                population.subpops[s]
                                        .individuals[i];

                int arm =
                        ind.producedByArm;

                if (ind.creationType
                        != CreationType.CROSSOVER) {

                    continue;
                }

                if (ind.birthGeneration
                        != this.generation) {

                    continue;
                }

                if (arm >= 0
                        && arm < armNum) {

                    survived[arm]++;
                }
            }
        }

        // =====================================
        // use global generated count
        // =====================================

        for (int arm = 0;
             arm < armNum;
             arm++) {

            generated[arm] =
                    AdaptiveParentSelection
                            .generatedByArm[arm];
        }

        // =====================================
        // reward = survival rate
        // =====================================

        for (int arm = 0;
             arm < armNum;
             arm++) {

            if (generated[arm] <= 0) {
                continue;
            }

            double reward =
                    (double) survived[arm]
                            / generated[arm];

            AdaptiveParentSelection
                    .currentGenerationSurvival[arm]
                    = reward;

            AdaptiveParentSelection.updateReward(
                    arm,
                    reward);
        }
    }

    public void writeArmSelectionToFile(
            ArrayList<int[]> armSelectionHistory) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + ".armSelection.csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,Arm0,Arm1,Arm2,Arm3");

            writer.newLine();

            for (int i = 0;
                 i < armSelectionHistory.size();
                 i++) {

                int[] values =
                        armSelectionHistory.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (int v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(line.toString());

                writer.newLine();
            }

            armSelectionHistory.clear();

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    public void writeArmRewardToFile(
            ArrayList<double[]> armRewardHistory) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + ".armReward.csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,Arm0,Arm1,Arm2,Arm3");

            writer.newLine();

            for (int i = 0;
                 i < armRewardHistory.size();
                 i++) {

                double[] values =
                        armRewardHistory.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (double v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(line.toString());

                writer.newLine();
            }

            armRewardHistory.clear();

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }



    public void writeArmSurvivalToFile(
            ArrayList<double[]> armSurvivalHistory) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + ".armSurvival.csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,Arm0,Arm1,Arm2,Arm3");

            writer.newLine();

            for (int i = 0;
                 i < armSurvivalHistory.size();
                 i++) {

                double[] values =
                        armSurvivalHistory.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (double v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(line.toString());

                writer.newLine();
            }

            armSurvivalHistory.clear();

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    public void writeImportanceToFile(
            ArrayList<double[]> history) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + ".importance.csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,"
                            + "Tree1Donor,"
                            + "Tree1Replace,"
                            + "Tree2Donor,"
                            + "Tree2Replace");

            writer.newLine();

            for (int i = 0; i < history.size(); i++) {

                double[] values =
                        history.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (double v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(
                        line.toString());

                writer.newLine();
            }

            writer.close();

        } catch (Exception e) {

            e.printStackTrace();
        }
    }

    public void writeDepthToFile(
            ArrayList<double[]> history) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + ".depth.csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,"
                            + "Tree1Donor,"
                            + "Tree1Replace,"
                            + "Tree2Donor,"
                            + "Tree2Replace");

            writer.newLine();

            for (int i = 0; i < history.size(); i++) {

                double[] values =
                        history.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (double v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(
                        line.toString());

                writer.newLine();
            }

            writer.close();

        } catch (Exception e) {

            e.printStackTrace();
        }
    }

    public void normaliseBestPopulation(double min) {

        ArrayList<Double> fitness = new ArrayList<>();

        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            if (individual.fitness.fitness() < Double.MAX_VALUE) {
                double[] objective = new double[1];
                objective[0] = individual.fitness.fitness() / min;
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                fitness.add(individual.fitness.fitness());
            }
        }
    }

    public void writeInvalidIndsNumGenToFile(ArrayList<Integer> invalidIndsNum) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".invalidIndsNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!invalidIndsNum.isEmpty()) {
                header.append(",invalidIndsNum");

            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < invalidIndsNum.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                line.append(", ").append(invalidIndsNum.get(i));

                writer.write(line.toString());
                writer.newLine();
            }

            invalidIndsNum.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSurrogateThresholdsToFile(ArrayList<Double> surrogateThresholds) {
        // 假设 out_dir 和 jobSeed 已经在类中定义
        File weightFile = new File(out_dir + "/job." + jobSeed + ".surrogateThresholds.csv");

        // 使用 try-with-resources 自动管理资源关闭，防止内存/句柄泄漏
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // 创建表头
            StringBuilder header = new StringBuilder("TaskID");
            if (surrogateThresholds != null && !surrogateThresholds.isEmpty()) {
                header.append(",surrogateThresholds");
            }
            writer.write(header.toString());
            writer.newLine();

            // 写入数据（去掉了逗号后面的空格，符合严格的 CSV 规范）
            if (surrogateThresholds != null) {
                for (int i = 0; i < surrogateThresholds.size(); i++) {
                    StringBuilder line = new StringBuilder();
                    line.append(i).append(",").append(surrogateThresholds.get(i));
                    writer.write(line.toString());
                    writer.newLine();
                }
            }

            // 建议删掉 surrogateThresholds.clear(); 避免产生意料之外的副作用

        } catch (IOException e) {
            System.err.println("写入 CSV 文件失败: " + weightFile.getAbsolutePath());
            e.printStackTrace();
        }
    }

    public Population preselection() {

        java.util.Arrays.fill(AdaptiveParentSelection.generatedByArm, 0);

        //ensure the elites can be saved in next generation
        ArrayList<Individual[]> elites = new ArrayList<Individual[]>(this.population.subpops.length);
        PopulationUtils.sort(this.population);

        numElites = this.parameters.getIntWithDefault(new Parameter(P_ELITES), null, 10);

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
        numRep = this.parameters.getIntWithDefault(new Parameter(P_REPLICATIONS), null, 3);
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

        PopulationUtils.sort(population); //sort the inds in the intermediate pop by estimated fitness

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            population.subpops[sub].resize(population.subpops[sub].individuals.length / numRep);
        }

        int invalidIndNum = 0;

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            for (int ind = 0; ind < this.population.subpops[sub].individuals.length; ind++) {
                Individual individual = this.population.subpops[sub].individuals[ind];
                if(individual.fitness.fitness() == Double.MAX_VALUE) {
                    invalidIndNum = this.population.subpops[sub].individuals.length - ind;
                    break;
                }
            }
        }

        System.out.println("invalidIndNum = " + invalidIndNum);

        invalidIndsNumber.add(invalidIndNum);

        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            int e = 0;
            for (int replace = population.subpops[sub].individuals.length - 1; replace >= population.subpops[sub].individuals.length - numElites; replace--) {
                population.subpops[sub].individuals[replace] = elites.get(sub)[e];
                e++;
            }
        }


        return population;
    }

    public void writeSubtreeImportanceToFile(
            ArrayList<double[]> hgHistory,
            ArrayList<double[]> csHistory,
            ArrayList<double[]> cgHistory) {

        writeSubtreeImportanceGroupToFile(hgHistory, "HGSubtreeImportance");
        writeSubtreeImportanceGroupToFile(csHistory, "CSSubtreeImportance");
        writeSubtreeImportanceGroupToFile(cgHistory, "CGSubtreeImportance");
    }

    private void writeSubtreeImportanceGroupToFile(
            ArrayList<double[]> history,
            String fileTag) {

        File file =
                new File(out_dir
                        + "/job."
                        + jobSeed
                        + "."
                        + fileTag
                        + ".csv");

        try {

            BufferedWriter writer =
                    new BufferedWriter(
                            new FileWriter(file));

            writer.write(
                    "Gen,Previous,Current,All");

            writer.newLine();

            for (int i = 0;
                 i < history.size();
                 i++) {

                double[] values =
                        history.get(i);

                StringBuilder line =
                        new StringBuilder(i + "");

                for (double v : values) {

                    line.append(",")
                            .append(v);
                }

                writer.write(line.toString());

                writer.newLine();
            }

            history.clear();

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

}


