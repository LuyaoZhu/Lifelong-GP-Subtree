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
import yimei.jss.algorithm.surrogateCaseLS.SilhouetteScore;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRule;
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
import java.util.stream.IntStream;

import static yimei.jss.algorithm.surrogateAccuracy.surrogateClearingMultitreeEvaluatorV1.SamePCNum;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * in each task, select individuals with good quality and diversity as the initial population of next task
 * in addition, the final population will be evaluated in 10 unseen instances to build surrogate
 * <p>
 * adaptively change the thresholds
 *
 * @author luyao
 */

public class GPRuleEvolutionStateLifelongGPV69 extends GPRuleEvolutionStateLifelongGPV68 {

    ConvergenceChecker checker;

    List<GPRSurrogateModel> GPRSurrogateModels = new ArrayList<>();


    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
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
        if (minCurrentTask == 0) {
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
            double[][] estimatedFitness1 = new double[surrogateFitness.size() + 1][population.subpops[0].individuals.length];

            double[] minFitnessGap = new double[surrogateSamples.size() + 1];

            double[] minInPreviousTask = new double[surrogateSamples.size()];

            double[] meanPCDistance = calculateMeanPCDistance(surrogateSamples);

            meanPCDistanceEveryGeneration.add(meanPCDistance);

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                estimatedFitness[estimatedFitness.length - 1][ind] = individual.fitness.fitness();
            }

            if (minCurrentTask == 0) { //means already normalised
                //based on the average fitness of top30% to determine whether we need to study this task only
                double meanFitnessCurrentTask = top30PercentMean(population.subpops[0].individuals);
//                if (onlyCurrentTaskPhase) {
                if (checker.check(meanFitnessCurrentTask) && onlyCurrentTaskPhase) {
                    System.out.println("Converged at generation: " + generation);
                    switchGen.add(generation);
                    onlyCurrentTaskPhase = false;
                    ((surrogateClearingMultitreeEvaluatorV69) evaluator).nonIntermediatePop = false;
                }
            }
            if(!onlyCurrentTaskPhase) {
                for (int t = 0; t < GPRSurrogateModels.size(); t++) {
                    for (int ind=0; ind < population.subpops[0].individuals.length; ind++) {
                        GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                        double[] doublePC =  Arrays.stream(individual.PC).asDoubleStream().toArray();
                        estimatedFitness1[t][ind] = GPRSurrogateModels.get(t).estimate(doublePC);
                    }
                    estimatedFitness[t] = evaluatePopulationV1(this, surrogateSamples.get(t), surrogateFitness.get(t), surrogateThreshold);

                    minInPreviousTask[t] = Arrays.stream(estimatedFitness[t]).min().getAsDouble();
                    if (minInPreviousTask[t] < Double.MAX_VALUE) {
                        for (int a = 0; a < estimatedFitness[t].length; a++) {
                            if (estimatedFitness[t][a] < Double.MAX_VALUE) {
                                estimatedFitness[t][a] = (estimatedFitness[t][a] - minInPreviousTask[t]) / (1 - minInPreviousTask[t]);
                            }
                        }
                    }
                }
            }
            if(generation%generationPerTask==0){
                //set thresholds
/*                for (int a = 0; a < thresholds.length; a++) {
                    thresholds[a] = 1.0/(thresholds.length-1);
                    if (a == thresholds.length - 1) {
                        thresholds[a] = 1;
                    }
                }*/
                for (int a = 0; a < thresholds.length; a++) {
                    thresholds[a] = 0;
                    if (a == thresholds.length - 1) {
                        thresholds[a] = 1;
                    }
                }
            } else {
                thresholds = thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1);
            }

            for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {

                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];

                individual.fitness.trials = new ArrayList<>();

                for (int task=0; task<estimatedFitness.length; task++) {
                    individual.fitness.trials.add(estimatedFitness[task][ind]);
                }

                double[] objective = new double[1];

                for (int t = 0; t < thresholds.length; t++) {
                    objective[0] += thresholds[t] * estimatedFitness[t][ind];
                    if (objective[0] >= Double.POSITIVE_INFINITY || objective[0] <= Double.NEGATIVE_INFINITY || Double.isNaN(objective[0])) {
                        objective[0] = Double.MAX_VALUE;
                    }
                }
                ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
            }

            //it's possible that all individuals' combined fitness in this population are Double.Max.
            // Thus, we still use the fitness in the current task to output the best individual and select parents, and preselection is based on fitness in all tasks.
            PopulationUtils.sort(population.subpops[0].individuals);
            if (population.subpops[0].individuals[0].fitness.fitness() > 10) {
                for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                    double[] objective = new double[1];
                    objective[0] = estimatedFitness[estimatedFitness.length - 1][ind];
                    ((MultiObjectiveFitness) individual.fitness).setObjectives(this, objective);
                }
            }

            //now we get the best individual in this generation,use it to calculate the threshold
            if (onlyCurrentTaskPhase) {
                for (int a = 0; a < thresholds.length; a++) {
                    thresholds[a] = 0;
                    if (a == thresholds.length - 1) {
                        thresholds[a] = 1;
                    }
                }
            } else {

//                if(thresholdsEveryGeneration.get(thresholdsEveryGeneration.size() - 1)[0] == 0){ //means this is the first time to calculate the thresholds
                for (int t = 0; t < minInPreviousTask.length; t++) {

                    // this is to calculate the forgetting ratio and calculate the thresholds
                    if (((Double)population.subpops[0].individuals[0].fitness.trials.get(t)) < Double.MAX_VALUE) {
                        double lowBound = Arrays.stream(surrogateFitness.get(t)).min().getAsDouble();
                        minFitnessGap[t] = ((Double)population.subpops[0].individuals[0].fitness.trials.get(t)) - (lowBound-minInPreviousTask[t]) / (1 - minInPreviousTask[t]);
                    } else {
                        minFitnessGap[t] = Double.MAX_VALUE;
                    }
                }
                thresholds = calculateThresholdsV5(estimatedFitness, minFitnessGap);

            }

            thresholdsEveryGeneration.add(thresholds);

            for (int t = 0; t < thresholds.length; t++) {
                System.out.println("The threshold of Task " + t + " is " + thresholds[t]);
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

            surrogateFitness.add(fitnessOneSurrogate);
            surrogateSamples.add(PC);   // ✅ 这里 PCs[ind] 与 QDIndividuals.get(ind) 严格一一对应

//            System.out.println(Arrays.stream(fitnessOneSurrogate).min().getAsDouble());
            double[][] newPC = new double[PC.length][PC[0].length];
            for (int a=0; a<PC.length; a++) {
                for (int b=0; b<PC[0].length; b++) {
                    newPC[a][b] = (double)PC[a][b];
                }
            }

            //use Gaussian Process
            GPRSurrogateModel surrogate = new GPRSurrogateModel();
            surrogate.train(newPC, fitnessOneSurrogate);

            GPRSurrogateModels.add(surrogate);


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

            if (onlyCurrentTaskPhase) {
                population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N
            } else {
                population = preselection(); //aims to select the individuals with fitness<Double.MaxValue
            }

        }


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

    private double[] calculateMeanPCDistance(ArrayList<int[][]> surrogateSamples) {
        double[] meanPCDistance = new double[surrogateSamples.size()];

        //calculate the fitness based on surrogate model
        for (int t = 0; t < surrogateSamples.size(); t++) {
            int[][] indsCharListsMultiTree = surrogateSamples.get(t);
            ArrayList<Double> PCDistances = new ArrayList<>();

            for (int i = 0; i < population.subpops[0].individuals.length; i++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[i];
                double dMin = Double.MAX_VALUE;
                int[] pcIntermediate = individual.PC;

                for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                    int[] pcModel = indsCharListsMultiTree[pc];
                    double d = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                    if (d == 0) {
                        dMin = d;
                        break;
                    }

                    if (d < dMin) {
                        dMin = d;
                    }

                }
                PCDistances.add(dMin);
                //record the medium

            }

            meanPCDistance[t] = PCDistances.stream()
                    .mapToDouble(Double::doubleValue)
                    .average()
                    .orElse(Double.NaN);
        }
        return meanPCDistance;
    }

    private double[] calculateMediumPCDistance(ArrayList<int[][]> surrogateSamples) {
        double[] mediumPCDistance = new double[surrogateSamples.size()];

        //calculate the fitness based on surrogate model
        for (int t = 0; t < surrogateSamples.size(); t++) {
            int[][] indsCharListsMultiTree = surrogateSamples.get(t);
            ArrayList<Double> PCDistances = new ArrayList<>();

            for (int i = 0; i < population.subpops[0].individuals.length; i++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[i];
                double dMin = Double.MAX_VALUE;
                int[] pcIntermediate = individual.PC;

                for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                    int[] pcModel = indsCharListsMultiTree[pc];
                    double d = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                    if (d == 0) {
                        dMin = d;
                        break;
                    }

                    if (d < dMin) {
                        dMin = d;
                    }

                }
                PCDistances.add(dMin);
                //record the medium

            }

            // 拷贝一份，避免改原数据
            ArrayList<Double> sorted = new ArrayList<>(PCDistances);
            Collections.sort(sorted);

            int n = sorted.size();
            if (n % 2 == 1) {
                // 奇数个
                mediumPCDistance[t] = sorted.get(n / 2);
            } else {
                // 偶数个
                mediumPCDistance[t] = (sorted.get(n / 2 - 1) + sorted.get(n / 2)) / 2.0;
            }

        }
        return mediumPCDistance;
    }

    private double[] calculateThresholds(double[][] estimatedFitness) {

        double[] thresholds = new double[estimatedFitness.length];

        double sumThreshold = 0;

        for (int i = 0; i < estimatedFitness.length; i++) {
            ArrayList<Double> fitnessInOneTask = new ArrayList<>();
            for (int j = 0; j < estimatedFitness[i].length; j++) {
                if (estimatedFitness[i][j] < Double.MAX_VALUE) {
                    fitnessInOneTask.add(estimatedFitness[i][j]);
                }
            }

            // 从小到大排序（fitness 越小越好）
            Collections.sort(fitnessInOneTask);

            // 取前 50%
            int k = fitnessInOneTask.size() / 2;
            List<Double> topHalf = fitnessInOneTask.subList(0, k);

            // 计算 mean 和 min
            double sum = 0.0;
            double min = Double.MAX_VALUE;

            for (double v : topHalf) {
                sum += v;
                min = Math.min(min, v);
            }

            double mean = topHalf.isEmpty() ? Double.NaN : sum / topHalf.size();

            if (i < estimatedFitness.length - 1) {
                min = Arrays.stream(surrogateFitness.get(i)).min().getAsDouble();
            }

            thresholds[i] = Math.max(0, (mean - min) / min);

            sumThreshold += thresholds[i];

        }

        for (int i = 0; i < thresholds.length; i++) {
            thresholds[i] /= sumThreshold;
        }

        return thresholds;

    }

    private double[] calculateThresholdsV1(double[][] estimatedFitness) {
        int n = estimatedFitness.length;
        double[] r = new double[n]; // 存储每个 task 的原始信号强度
        double[] finalWeights = new double[n];
        double sumThreshold = 0;

        for (int i = 0; i < estimatedFitness.length; i++) {
            ArrayList<Double> fitnessInOneTask = new ArrayList<>();
            for (int j = 0; j < estimatedFitness[i].length; j++) {
                if (estimatedFitness[i][j] < Double.MAX_VALUE) {
                    fitnessInOneTask.add(estimatedFitness[i][j]);
                }
            }

            // 从小到大排序（fitness 越小越好）
            Collections.sort(fitnessInOneTask);

            // 取前 50%
            int k = fitnessInOneTask.size() / 2;
            List<Double> topHalf = fitnessInOneTask.subList(0, k);

            // 计算 mean 和 min
            double sum = 0.0;
            double min = Double.MAX_VALUE;

            for (double v : topHalf) {
                sum += v;
                min = Math.min(min, v);
            }

            double mean = topHalf.isEmpty() ? Double.NaN : sum / topHalf.size();

            if (i < estimatedFitness.length - 1) {
                min = Arrays.stream(surrogateFitness.get(i)).min().getAsDouble();
            }

            r[i] = Math.max(0, (mean - min) / min);

            sumThreshold += r[i];

        }

        // 分配最终权重
        for (int i = 0; i < n; i++) {
            if (i == n - 1) {
                // 当前任务：0.5 固定 + 0.5 竞争
                finalWeights[i] = 0.5 + 0.5 * (r[i] / sumThreshold);
            } else {
                // 历史任务：共享剩下的 0.5 竞争额度
                finalWeights[i] = 0.5 * (r[i] / sumThreshold);
            }
        }

        return finalWeights;
    }

    private double[] calculateThresholdsV2(double[][] estimatedFitness, double[] minFitness) {

        double[] thresholds = new double[estimatedFitness.length];
        double[] variance = new double[estimatedFitness.length];

        double sumThreshold = 0;

        for (int i = 0; i < estimatedFitness.length; i++) {
            ArrayList<Double> fitnessInOneTask = new ArrayList<>();
            for (int j = 0; j < estimatedFitness[i].length; j++) {
                if (estimatedFitness[i][j] < Double.MAX_VALUE) {
                    fitnessInOneTask.add(estimatedFitness[i][j]);
                }
            }

//            variance[i] = varianceStream(fitnessInOneTask,false);
//            System.out.println(variance);

            // 计算 mean 和 min
            double sum = 0.0;

            for (double v : fitnessInOneTask) {
                sum += v;
            }

            double mean = fitnessInOneTask.isEmpty() ? Double.NaN : sum / fitnessInOneTask.size();

            thresholds[i] = Math.max(0, mean - minFitness[i]);

            sumThreshold += thresholds[i];

        }

        for (int i = 0; i < thresholds.length; i++) {
            if (i == thresholds.length - 1) {
                thresholds[i] = thresholds[i] / sumThreshold * 0.5 + 0.5;
            } else {
                thresholds[i] = thresholds[i] / sumThreshold * 0.5;
            }
        }

//        for (int i = 0; i < thresholds.length; i++) {
//            thresholds[i] = thresholds[i]/sumThreshold;
////            thresholds[i] = thresholds[i]/variance[i];
//        }

        return thresholds;
    }

    private double[] calculateThresholdsV3(
            double[][] estimatedFitness,
            double[] minFitness
    ) {

        int T = estimatedFitness.length;
        int N = estimatedFitness[0].length;

        // 存储符合条件的个体索引
        ArrayList<Integer> selectedIndividuals = new ArrayList<>();

        // =================================================
        // Step 0: 筛选在所有 task 上都有效的个体
        // =================================================
        for (int j = 0; j < N; j++) {
            boolean isAllValid = true;

            // 检查当前个体 j 在所有任务 i 下的 fitness
            for (int i = 0; i < T; i++) {
                if (estimatedFitness[i][j] >= Double.MAX_VALUE) {
                    isAllValid = false;
                    break; // 只要有一个无效，跳出内层循环
                }
            }

            // 只有全部有效才加入
            if (isAllValid) {
                selectedIndividuals.add(j);
            }
        }

        double[] thresholds = new double[T];
        double sumThreshold = 0.0;

        // =================================================
        // Step 1: 计算均值与阈值
        // =================================================
        if (selectedIndividuals.size() > 0) { // 如果有符合条件的个体
            for (int i = 0; i < T - 1; i++) {
                double sum = 0.0;
                for (int j : selectedIndividuals) {
                    sum += estimatedFitness[i][j];
                }

                double mean = sum / selectedIndividuals.size();

                // 计算与最小值的偏差作为原始阈值
                thresholds[i] = Math.max(0, mean - minFitness[i]);
                sumThreshold += thresholds[i];
            }

            // =================================================
            // Step 2: 归一化逻辑 (保持原样)
            // =================================================
            for (int i = 0; i < thresholds.length; i++) {
                if (i == thresholds.length - 1) {
                    thresholds[i] = 1.0;
                } else {
                    // 防止 sumThreshold 为 0 导致除以零（NaN）
                    thresholds[i] = (sumThreshold > 0) ? (thresholds[i] / sumThreshold) : 0;
                }
            }
        } else {
            // 如果没有个体满足“全task有效”的条件，进入保底分配
            for (int i = 0; i < thresholds.length; i++) {
                if (i == thresholds.length - 1) {
                    thresholds[i] = 1.0;
                } else {
                    thresholds[i] = 1.0 / (T - 1);
                }
            }
        }
        return thresholds;
    }

    private double[] calculateThresholdsV4(
            double[] minFitness
    ) {

        int T = minFitness.length;

        double[] thresholds = new double[T];
        double sumThreshold = 0.0;

        // =================================================
        // Step 1: 计算均值与阈值
        // =================================================

        for (int i = 0; i < T - 1; i++) {

            if (minFitness[i] == Double.MAX_VALUE) { // means that all fitness is Double.MaxValue in this task
                thresholds[i] = 1;
            } else {
                thresholds[i] = 0 - minFitness[i];
            }
            sumThreshold += thresholds[i];
        }

        // =================================================
        // Step 2: 归一化逻辑 (保持原样)
        // =================================================
        for (int i = 0; i < thresholds.length; i++) {
            if (i == thresholds.length - 1) {
                thresholds[i] = 1.0;
            } else {
                // 防止 sumThreshold 为 0 导致除以零（NaN）
                thresholds[i] = (sumThreshold > 0) ? (thresholds[i] / sumThreshold) : 0;
            }
        }

        return thresholds;
    }

    private double[] calculateThresholdsV5(
            double[][] estimatedFitness,   // [task][individual]
            double[] minFitness         // task 权重

    ) {

        int Task = estimatedFitness.length;
        double[] thresholds = new double[Task];
        double sumThreshold = 0.0;

        for (int t = 0; t < Task - 1; t++) {

            if (minFitness[t] == Double.MAX_VALUE) {
                thresholds[t] = 1.0;
            } else {
                thresholds[t] = minFitness[t];  // 仍按你原逻辑
            }

            sumThreshold += thresholds[t];
        }

        // =================================================
        // Step 5: 归一化
        // =================================================
        for (int t = 0; t < Task; t++) {

            if (t == Task - 1) {
                thresholds[t] = 1.0;  // 当前任务固定为1
            } else {
                thresholds[t] = (sumThreshold > 0)
                        ? thresholds[t] / sumThreshold
                        : 0.0;
            }
        }

        return thresholds;
    }

    private double[] calculateThresholdsSimilarity(
            double[][] estimatedFitness
    ) {

        double[] thresholds = new double[estimatedFitness.length];

        int numTasks = estimatedFitness.length;
        int numIndividuals = estimatedFitness[0].length;

        int[][] PC = new int[numIndividuals][];

        //first remove some individuals with the same PC
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind].clone();
            PC[ind] = individual.PC;
        }


        // ========= 1️⃣ 标记需要删除的个体 =========
        boolean[] validIndividual = new boolean[numIndividuals];
        int validCount = 0;

        Set<String> seenPCs = new HashSet<>(); // 用字符串表示 PC，方便判断重复

        for (int j = 0; j < numIndividuals; j++) {

            boolean valid = true;

            // 1. 删除 fitness = Double.MAX_VALUE
            for (int i = 0; i < numTasks; i++) {
                if (estimatedFitness[i][j] == Double.MAX_VALUE) {
                    valid = false;
                    break;
                }
            }

            // 2. 删除重复 PC
            String pcKey = Arrays.toString(PC[j]);
            if (seenPCs.contains(pcKey)) {
                valid = false;
            } else if (valid) {
                seenPCs.add(pcKey);
            }

            validIndividual[j] = valid;
            if (valid) validCount++;
        }

        // ========= 2️⃣ 构造过滤后的 fitness =========
        double[][] filtered = new double[numTasks][validCount];

        int newCol = 0;
        for (int j = 0; j < numIndividuals; j++) {
            if (validIndividual[j]) {
                for (int i = 0; i < numTasks; i++) {
                    filtered[i][newCol] = estimatedFitness[i][j];
                }
                newCol++;
            }
        }

        // ========= 3️⃣ 计算 Spearman =========

        SpearmansCorrelation spearmansCorrelation = new SpearmansCorrelation();

        double[] lastTask = filtered[numTasks - 1];

        for (int i = 0; i < numTasks; i++) {
            if (i == numTasks - 1 || lastTask.length <= 5) {
                thresholds[i] = 1;
            } else {
                double corr = spearmansCorrelation.correlation(filtered[i], lastTask);
                if (Double.isNaN(corr) || Double.isInfinite(corr)) {
                    thresholds[i] = 1;   // treat as no correlation
                } else {
                    thresholds[i] = 1 - Math.abs(corr);
                }
            }
        }

        return thresholds;
    }

    private double[] calculateThresholdsSimilarityV1(
            double[][] estimatedFitness
    ) {

        double[] thresholds = new double[estimatedFitness.length];

        int numTasks = estimatedFitness.length;
        int numIndividuals = estimatedFitness[0].length;

        int[][] PC = new int[numIndividuals][];

        //first remove some individuals with the same PC
        for (int ind = 0; ind < population.subpops[0].individuals.length; ind++) {
            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind].clone();
            PC[ind] = individual.PC;
        }


        // ========= 1️⃣ 标记需要删除的个体 =========
        boolean[] validIndividual = new boolean[numIndividuals];
        int validCount = 0;

        Set<String> seenPCs = new HashSet<>(); // 用字符串表示 PC，方便判断重复

        for (int j = 0; j < numIndividuals; j++) {

            boolean valid = true;

            // 1. 删除 fitness = Double.MAX_VALUE
            for (int i = 0; i < numTasks; i++) {
                if (estimatedFitness[i][j] == Double.MAX_VALUE) {
                    valid = false;
                    break;
                }
            }

            // 2. 删除重复 PC
            String pcKey = Arrays.toString(PC[j]);
            if (seenPCs.contains(pcKey)) {
                valid = false;
            } else if (valid) {
                seenPCs.add(pcKey);
            }

            validIndividual[j] = valid;
            if (valid) validCount++;
        }

        // ========= 2️⃣ 构造过滤后的 fitness =========
        double[][] filtered = new double[numTasks][validCount];

        int newCol = 0;
        for (int j = 0; j < numIndividuals; j++) {
            if (validIndividual[j]) {
                for (int i = 0; i < numTasks; i++) {
                    filtered[i][newCol] = estimatedFitness[i][j];
                }
                newCol++;
            }
        }

        // ========= 3️⃣ 计算 Spearman =========

        SpearmansCorrelation spearmansCorrelation = new SpearmansCorrelation();

        double[] lastTask = filtered[numTasks - 1];
        double sum = 0.0;

        for (int i = 0; i < numTasks; i++) {
            if (i == numTasks - 1 || lastTask.length <= 5) {
                thresholds[i] = 1;
            } else {
                double corr = spearmansCorrelation.correlation(filtered[i], lastTask);
                if (Double.isNaN(corr) || Double.isInfinite(corr)) {
                    thresholds[i] = 1;   // treat as no correlation
                } else {
                    thresholds[i] = 1 - Math.abs(corr);
                }
            }
            sum += thresholds[i];
        }

        for (int j = 0; j < thresholds.length - 1; j++) {
            thresholds[j] = thresholds[j] / (sum - 1);
        }


        return thresholds;
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

    public void writeThresholdSwitchGenToFile(ArrayList<Integer> switchGen) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".thresholdSwitchGen.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!switchGen.isEmpty()) {
                header.append(",switchGen");

            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < switchGen.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                line.append(", ").append(switchGen.get(i));

                writer.write(line.toString());
                writer.newLine();
            }

            switchGen.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeThresholdsToFile(ArrayList<double[]> thresholdsEveryGen) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".thresholdsEveryGen.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of thresholdsEveryGen elements
            StringBuilder header = new StringBuilder("Gen");
            if (!thresholdsEveryGen.isEmpty()) {
                for (int j = 0; j < thresholdsEveryGen.get(0).length; j++) {
                    header.append(",threshold").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < thresholdsEveryGen.size(); i++) {
                StringBuilder line = new StringBuilder(i + 100 + "");
                for (double value : thresholdsEveryGen.get(i)) {
                    line.append(", ").append(value);
                }
                writer.write(line.toString());
                writer.newLine();
            }

            thresholdsEveryGen.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeMeanPCDistanceToFile(ArrayList<double[]> mediumPCDistance) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".meanPCDistance.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of mediumPCDistance elements
            StringBuilder header = new StringBuilder("Gen");
            if (!mediumPCDistance.isEmpty()) {
                for (int j = 0; j < mediumPCDistance.get(0).length; j++) {
                    header.append(",meanPCDistance").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < mediumPCDistance.size(); i++) {
                StringBuilder line = new StringBuilder(i + 100 + "");
                for (double value : mediumPCDistance.get(i)) {
                    line.append(", ").append(value);
                }
                writer.write(line.toString());
                writer.newLine();
            }

            mediumPCDistance.clear();
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
        int countInfinity = 0;
        for (int i=0; i<population.subpops[0].individuals.length; i++) {
            Individual individual = population.subpops[0].individuals[i];
            if(individual.fitness.fitness() > 10){
                countInfinity++;
            }
        }

        System.out.println("countInfinity = " + countInfinity);

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
        double[] realFitness = new double[population.subpops[0].individuals.length];

        int[] referencePC = new int[40];
        Arrays.fill(referencePC, 1);

        for (int sub = 0; sub < state.population.subpops.length; sub++) {

            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                double rawEstimated = 0;
                double confidenceEstimatedFitness = 0;
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

                    if (dMin < threshold) {
                        rawEstimated = fitnessesForModel[index];
                        estimatedFitness[i] = rawEstimated + (1 - rawEstimated) / threshold * dMin;
//
                    } else {
                        estimatedFitness[i] = Double.MAX_VALUE;
                        rawEstimated = Double.MAX_VALUE;
                        confidenceEstimatedFitness = Double.MAX_VALUE;
                    }
                } else {
                    estimatedFitness[i] = Double.MAX_VALUE;
                }

                RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
                DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
                double[] refFitness = new double[simulationsPerTask];

                //real fitness
                for (int s = 0; s < simulationsPerTask; s++) {

                    simulation.reset(608 + s);
                    refFitness[s] = calculateReferenceRuleFitness(simulation);

                    // （可选但推荐）避免 refFit1[1] 为 0 导致除零
                    if (refFitness[s] == 0.0) {
                        throw new IllegalStateException("refFit1[1] is 0, cannot normalise ObjValue / refFit1[1].");
                    }

                    GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) individual).trees[0]);
                    GPRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) individual).trees[1]);
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
                    realFitness[i] += (ObjValue / refFitness[s]) / simulationsPerTask;
                }
                System.out.println("estimated Fitness: " + rawEstimated + " , confidence Fitness: " + estimatedFitness[i] + " , realFitness: " + realFitness[i]);
            }

        }

        return estimatedFitness;

    }

    public double[] evaluatePopulationV1(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel, double threshold) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        double[] estimatedFitness = new double[population.subpops[0].individuals.length];
        double[] realFitness = new double[population.subpops[0].individuals.length];

        for (int sub = 0; sub < state.population.subpops.length; sub++) {

            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                double rawEstimated = 0;

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
                    rawEstimated = fitnessesForModel[index];
                    estimatedFitness[i] = fitnessesForModel[index];
                /*} else if (threshold / 2 < dMin && dMin < threshold) {
                    rawEstimated = fitnessesForModel[index];
                    //find the neighbourhood of this sample, and based on the distance to determine the final fitness, x=distance, y=fitness
                    ArrayList<Double> x_distance = new ArrayList<>();
                    ArrayList<Double> y_fitness = new ArrayList<>();

                    for (int pc = 0; pc < indsCharListsMultiTree.length; pc++) {
                        int[] pcModel = indsCharListsMultiTree[pc];
                        double d = PhenoCharacterisation.distance(pcModel, indsCharListsMultiTree[index]);
                        //determine the neighbour (distance<10)
                        if (d <= threshold && fitnessesForModel[pc] < Double.MAX_VALUE) {
                            x_distance.add(d);
                            y_fitness.add(fitnessesForModel[pc]);
                        }
                    }
                    //then build regression model and set estimatde fitness
                    LinearRegression lr = new LinearRegression();
                    lr.fit(x_distance, y_fitness);
                    estimatedFitness[i] = lr.predict(dMin);*/
//
                } else {
                    estimatedFitness[i] = Double.MAX_VALUE;
                    rawEstimated = Double.MAX_VALUE;
                }


                /*RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
                DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
                double[] refFitness = new double[simulationsPerTask];

                //real fitness
                for (int s = 0; s < simulationsPerTask; s++) {

                    simulation.reset(608 + s);
                    refFitness[s] = calculateReferenceRuleFitness(simulation);

                    // （可选但推荐）避免 refFit1[1] 为 0 导致除零
                    if (refFitness[s] == 0.0) {
                        throw new IllegalStateException("refFit1[1] is 0, cannot normalise ObjValue / refFit1[1].");
                    }

                    GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) individual).trees[0]);
                    GPRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) individual).trees[1]);
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
                    realFitness[i] += ObjValue / refFitness[s];
                }

                realFitness[i] = realFitness[i] / simulationsPerTask;

                System.out.println("raw Fitness: " + rawEstimated + " , estimated Fitness: " + estimatedFitness[i] + " , realFitness: " + realFitness[i] + " , distance: " + dMin);*/
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


    public double median(ArrayList<Double> list) {
        if (list == null || list.isEmpty()) {
            throw new IllegalArgumentException("List is null or empty");
        }

        // 拷贝一份，避免改变原列表顺序
        ArrayList<Double> sorted = new ArrayList<>(list);
        Collections.sort(sorted);

        int n = sorted.size();
        int mid = n / 2;

        if (n % 2 == 1) {
            // 奇数个
            return sorted.get(mid);
        } else {
            // 偶数个，取中间两个的平均
            return (sorted.get(mid - 1) + sorted.get(mid)) / 2.0;
        }
    }

    public static double varianceStream(List<Double> list, boolean sample) {
        int n = list.size();
        if (n == 0 || (sample && n == 1)) return 0.0;

        double mean = list.stream()
                .mapToDouble(Double::doubleValue)
                .average()
                .orElse(0.0);

        double sum = list.stream()
                .mapToDouble(x -> (x - mean) * (x - mean))
                .sum();

        return sample ? sum / (n - 1) : sum / n;
    }

    public static void investigateDistanceFitnessRelation(
            int[][] PC,
            double[] fitnessOneSurrogate
    ) {
        int n = PC.length;

        // =========================
        // 1. 找 baseline（最小 fitness）
        // =========================
        int bestIdx = 0;
        double bestFitness = Double.MAX_VALUE;
        for (int i = 0; i < n; i++) {
            if (fitnessOneSurrogate[i] < bestFitness) {
                bestFitness = fitnessOneSurrogate[i];
                bestIdx = i;
            }
        }

        int[] bestPC = PC[bestIdx];

        // =========================
        // 2. 计算 distance 和 Δfitness
        // =========================
        List<Double> distances = new ArrayList<>();
        List<Double> deltaFitness = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            if (i == bestIdx) continue;

            double d = PhenoCharacterisation.distance(PC[i], bestPC);
            double df = fitnessOneSurrogate[i] - bestFitness;

            distances.add(d);
            deltaFitness.add(df);
        }

        // =========================
        // 3. 相关性分析
        // =========================

        double pearson = pearsonCorrelation(distances, deltaFitness);
        double spearman = spearmanCorrelation(distances, deltaFitness);

        // =========================
        // 4. 输出
        // =========================
        System.out.println("Baseline index: " + bestIdx);
        System.out.println("Baseline fitness: " + bestFitness);
        System.out.println("Pearson(distance, Δfitness)   = " + pearson);
        System.out.println("Spearman(distance, Δfitness)  = " + spearman);
    }

    public static void investigateTopKLocalRelation(int[][] PC, double[] fitness, int Kbest) {

        int N = fitness.length;

        // 1️⃣ 找 fitness 最小的 Kbest 个体
        Integer[] order = IntStream.range(0, N).boxed().toArray(Integer[]::new);
        Arrays.sort(order, Comparator.comparingDouble(i -> fitness[i]));

        int[] bestIdx = new int[Kbest];
        for (int i = 0; i < Kbest; i++) {
            bestIdx[i] = order[i];
        }

        // 2️⃣ 遍历所有个体，找最近 best
        List<Double> distances = new ArrayList<>();
        List<Double> deltaFitness = new ArrayList<>();

        for (int i = 0; i < N; i++) {

            double minDist = Double.MAX_VALUE;
            int nearestBest = -1;

            for (int b : bestIdx) {
                double d = PhenoCharacterisation.distance(PC[i], PC[b]);
                if (d < minDist) {
                    minDist = d;
                    nearestBest = b;
                }
            }

            if (nearestBest >= 0 &&
                    fitness[i] != Double.MAX_VALUE &&
                    fitness[nearestBest] != Double.MAX_VALUE) {

                distances.add(minDist);
                deltaFitness.add(fitness[i] - fitness[nearestBest]);
            }
        }

        double pearson = pearsonCorrelation(distances, deltaFitness);
        double spearman = spearmanCorrelation(distances, deltaFitness);

        // =========================
        // 4. 输出
        // =========================
        Arrays.sort(fitness);
        System.out.println("Pearson(distance, Δfitness)   = " + pearson);
        System.out.println("Spearman(distance, Δfitness)  = " + spearman);


    }

    private static double pearsonCorrelation(List<Double> x, List<Double> y) {
        int n = x.size();
        double meanX = x.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double meanY = y.stream().mapToDouble(Double::doubleValue).average().orElse(0);

        double num = 0, denX = 0, denY = 0;
        for (int i = 0; i < n; i++) {
            double dx = x.get(i) - meanX;
            double dy = y.get(i) - meanY;
            num += dx * dy;
            denX += dx * dx;
            denY += dy * dy;
        }
        return num / Math.sqrt(denX * denY + 1e-12);
    }

    private static double spearmanCorrelation(List<Double> x, List<Double> y) {
        List<Double> rx = rank(x);
        List<Double> ry = rank(y);
        return pearsonCorrelation(rx, ry);
    }

    private static List<Double> rank(List<Double> v) {
        int n = v.size();
        Integer[] idx = IntStream.range(0, n).boxed().toArray(Integer[]::new);

        Arrays.sort(idx, Comparator.comparingDouble(v::get));

        double[] ranks = new double[n];
        for (int i = 0; i < n; i++) {
            ranks[idx[i]] = i + 1;
        }

        List<Double> res = new ArrayList<>();
        for (double r : ranks) res.add(r);
        return res;
    }

    public class LinearRegression {
        private double slope;      // k
        private double intercept;  // b
        private boolean fitted = false;

        // 拟合
        public void fit(ArrayList<Double> x, ArrayList<Double> y) {
            if (x.size() != y.size()) {
                throw new IllegalArgumentException("x and y must have same size");
            }
            int n = x.size();
            double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

            for (int i = 0; i < n; i++) {
                double xi = x.get(i);
                double yi = y.get(i);
                sumX += xi;
                sumY += yi;
                sumXY += xi * yi;
                sumXX += xi * xi;
            }

            slope = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX);
            intercept = (sumY - slope * sumX) / n;
            fitted = true;
        }

        // 预测
        public double predict(double x) {
            if (!fitted) throw new IllegalStateException("Model not fitted yet");
            return slope * x + intercept;
        }

        // 获取参数
        public double getSlope() {
            return slope;
        }

        public double getIntercept() {
            return intercept;
        }
    }

    public static double top30PercentMean(Individual[] GPIndividuals) {

        double[] a = new double[GPIndividuals.length];
        for (int i = 0; i < GPIndividuals.length; i++) {
            a[i] = GPIndividuals[i].fitness.fitness();
        }

        double[] copy = Arrays.copyOf(a, a.length);
        Arrays.sort(copy); // 从小到大排序

        int count = (int) Math.ceil(copy.length * 0.3);

        double sum = 0.0;
        for (int i = 0; i < count; i++) {
            sum += copy[i];
        }

        return sum / count;
    }

    public class ConvergenceChecker {

        private LinkedList<Double> window = new LinkedList<>();

        private int windowSize = 5;
        private double improvementThreshold = 0.001; // 0.1%

        public boolean check(double currentMean) {

            window.add(currentMean);

            if (window.size() < windowSize) {
                return false;
            }

            if (window.size() > windowSize) {
                window.removeFirst();
            }

            double totalImprovement = 0.0;

            for (int i = 1; i < window.size(); i++) {
                double prev = window.get(i - 1);
                double curr = window.get(i);

                double improvement = (prev - curr) / Math.abs(prev);
                totalImprovement += improvement;
            }

            double avgImprovement = totalImprovement / (windowSize - 1);

            return avgImprovement < improvementThreshold;
        }
    }
}



