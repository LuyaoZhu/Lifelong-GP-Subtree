package yimei.jss.algorithm.lifelongGP;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleEvaluator;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.Clearable;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import static yimei.jss.algorithm.surrogateLSnew.GPRuleEvolutionStateSavedSurrogateCases.realEvaluationIntermediate;
import static yimei.jss.gp.GPRuleEvolutionState.phenoCharacterisation;


/**
 * Created by fzhang on 2019.9.24.
 */
//a class can not extend from more than one class
public class surrogateClearingMultitreeEvaluatorV68 extends SimpleEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

//    public static final String P_PERCENTAGESURRACC = "percentage-surrogateAcc";
//    protected double percentageSurrogateAcc;

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    //ArrayList<Double[]> fitnessesForModel = new ArrayList<>();
    double[] tempfitnessesForModel; //for transfer the fitness into the surrogate model to estimate the fitness
    int[][] tempindsCharListsMultiTree;

    protected boolean clear = true;
    protected boolean nonIntermediatePop = true;

    protected double radius;
    protected int capacity;

    public static ArrayList<Double> correlationSurrogateAccuracy = new ArrayList<>();
    public static ArrayList<Double> correlationSurrogateAccuracyWhole = new ArrayList<>();
    public static ArrayList<Double> distanceRankSurrogateAccuracy = new ArrayList<>();
    public static ArrayList<Double> distanceRankSurrogateAccuracyWhole = new ArrayList<>();

    public ArrayList<Double> estimatedFitness;
    public ArrayList<Double> realFitness;

    public static ArrayList<Double> MSEGen = new ArrayList<>();
    public static ArrayList<Double> SpearmanCorrelationGen = new ArrayList<>();
    public static ArrayList<Integer> SamePCNum = new ArrayList<>();
    public static ArrayList<Double> AveragePCDistance = new ArrayList<>();

    HashMap<List<Integer>, List<Double>> PCFitnessMap = new HashMap<>();


//    protected PhenoCharacterisation[] phenoCharacterisation;


    public static ArrayList<double[]> entropyDiversity = new ArrayList<>();

    public double getRadius() {
        return radius;
    }

    public int getCapacity() {
        return capacity;
    }

    public PhenoCharacterisation[] getPhenoCharacterisation() {
        return phenoCharacterisation;
    }

    public PhenoCharacterisation getPhenoCharacterisation(int index) {
        return phenoCharacterisation[index];
    }

    public static ArrayList<Double> pcDistance = new ArrayList<>();

    public static ArrayList<Double> averageFitnessGens = new ArrayList<>();



    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        radius = state.parameters.getDoubleWithDefault(
                base.push(P_RADIUS), null, 0.0);
        capacity = state.parameters.getIntWithDefault(
                base.push(P_CAPACITY), null, 1);
//        percentageSurrogateAcc = state.parameters.getDoubleWithDefault(
//                base.push(P_PERCENTAGESURRACC), null, 0.0);

        String filePath = state.parameters.getString(new Parameter("filePath"), null);
        //It's a little tricky to know whether we have 1 or 2 populations here, so we will assume
        //2 for the purpose of the phenoCharacterisation, and ignore the second object if only
        //1 is used
    }

    @Override
    public void evaluatePopulation(final EvolutionState state) {

        //if it is a normal evaluation
        if (realEvaluationIntermediate == true) {
//            multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);
            super.evaluatePopulation(state);
        } else {

            if (nonIntermediatePop) {
                super.evaluatePopulation(state);

                double[] fitness = new double[state.population.subpops[0].individuals.length];
                for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {
                    GPIndividual individual = (GPIndividual) state.population.subpops[0].individuals[ind];
                    fitness[ind] = individual.fitness.fitness();
                }

                double min = Arrays.stream(fitness).min().getAsDouble();

                if( min < ((GPRuleEvolutionStateLifelongGPV68) state).refFit) {
                    //then normalise
                    ((GPRuleEvolutionStateLifelongGPV68) state).normalisePopulation(min, ((GPRuleEvolutionStateLifelongGPV68) state).refFit);
                }

                int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

                //after normalisation, we build surrogate
                double[] fitnessesForModel = new double[state.population.subpops[0].individuals.length];
                int[][] indsCharListsMultiTreeCopy = indsCharListsMultiTree.clone();
                for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {
                    fitnessesForModel[ind] = state.population.subpops[0].individuals[ind].fitness.fitness();
                }

                // int removeIdx = 0;
                for (int i = 0; i < fitnessesForModel.length; i++) {
                    if (fitnessesForModel[i] == Double.MAX_VALUE) {
                   /* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*/
                        indsCharListsMultiTreeCopy = ArrayUtils.remove(indsCharListsMultiTreeCopy, i);
                        fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i);
                        i--;
                        // removeIdx++;
                    }
                }

                tempfitnessesForModel = fitnessesForModel;
                tempindsCharListsMultiTree = indsCharListsMultiTreeCopy;


                if(!((GPRuleEvolutionStateLifelongGPV68) state).onlyCurrentTaskPhase) {
                    nonIntermediatePop = false;
                }

                if(state.generation == ((GPRuleEvolutionStateLifelongGPV68) state).generationPerTask - 1 ||
                        state.generation == ((GPRuleEvolutionStateLifelongGPV68) state).generationPerTask * 2 - 1 ||
                        state.generation == ((GPRuleEvolutionStateLifelongGPV68) state).generationPerTask * 3 - 1 ||
                        state.generation % ((GPRuleEvolutionStateLifelongGPV68) state).generationPerTask == 0 ) {
                    nonIntermediatePop = true;
                }

//            surrogateFitnessEveryGen.add(tempfitnessesForModel);
//            surrogatePCEveryGen.add(tempindsCharListsMultiTree);

            } else {
//                multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);//the individuals do not have fitness or they are evaluated yet
//            multitreeClearingFirst.clearPopulationRanks(state, radius, capacity, phenoCharacterisation);
                realFitness = new ArrayList<>();
                estimatedFitness = new ArrayList<>();
//                this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel);

                //we need to assign estimated fitness of all tasks based on weights
                double[][] estimatedFitness = new double[state.generation/((GPRuleEvolutionStateLifelongGPV68) state).generationPerTask + 1][state.population.subpops[0].individuals.length];
                for (int task=0; task<estimatedFitness.length; task++) {
                    if (task == estimatedFitness.length - 1) {
                        estimatedFitness[task] = this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel, ((GPRuleEvolutionStateLifelongGPV68) state).surrogateThreshold);
                    } else {
                        estimatedFitness[task] = this.evaluatePopulation(state, ((GPRuleEvolutionStateLifelongGPV68) state).surrogateSamples.get(task), ((GPRuleEvolutionStateLifelongGPV68) state).surrogateFitness.get(task), ((GPRuleEvolutionStateLifelongGPV68) state).surrogateThreshold);
                        //normalise
                        double minInPreviousTask = Arrays.stream(estimatedFitness[task]).min().getAsDouble();
                        for (int a = 0; a < estimatedFitness[task].length; a++) {
                            if (estimatedFitness[task][a] < Double.MAX_VALUE) {
                                estimatedFitness[task][a] = (estimatedFitness[task][a] - minInPreviousTask) / (1 - minInPreviousTask);
                            }
                        }
                    }
                }

                double[] thresholds = ((GPRuleEvolutionStateLifelongGPV68) state).thresholdsEveryGeneration.get(((GPRuleEvolutionStateLifelongGPV68) state).thresholdsEveryGeneration.size()-1);

                //then assign combined fitness
                for (int ind = 0; ind < state.population.subpops[0].individuals.length; ind++) {

                    GPIndividual individual = (GPIndividual) state.population.subpops[0].individuals[ind];

                    double[] objective = new double[1];

                    for (int t = 0; t < thresholds.length; t++) {
                        objective[0] += thresholds[t] * estimatedFitness[t][ind];
                        if (objective[0] >= Double.POSITIVE_INFINITY || objective[0] <= Double.NEGATIVE_INFINITY || Double.isNaN(objective[0])) {
//                        state.output.warning("Bad objective #" + ": " + objective[0] + ", setting to worst value for that objective.");
                            objective[0] = Double.MAX_VALUE;
                        }
                    }
                    ((MultiObjectiveFitness) individual.fitness).setObjectives(state, objective);
                }


//            calculateSurrogateAccuracy();
                nonIntermediatePop = true;
            }
        }

    }

    private void calculateSurrogateAccuracy() {

        double MSE = 0.0;
        double infinityNum = 0;
        for (int i = 0; i < estimatedFitness.size(); i++) {
            if (estimatedFitness.get(i).equals(Double.MAX_VALUE) || realFitness.get(i).equals(Double.MAX_VALUE)) {
                infinityNum++;
            } else {
                MSE += Math.pow(estimatedFitness.get(i) - realFitness.get(i), 2);
            }
        }
        MSE /= (estimatedFitness.size() - infinityNum);
        double RMSE = Math.sqrt(MSE);

        MSEGen.add(RMSE);

        SpearmansCorrelation spearmansCorrelation = new SpearmansCorrelation();
        double[] realFit = realFitness.stream().mapToDouble(Double::doubleValue).toArray();
        double[] estimatedFit = estimatedFitness.stream().mapToDouble(Double::doubleValue).toArray();
        double correlation = spearmansCorrelation.correlation(estimatedFit, realFit);

        SpearmanCorrelationGen.add(correlation);

    }

    //==========================================KNN========================================
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int i = 0; i < state.population.subpops[0].individuals.length; i++) //
        {
            GPIndividual individual = (GPIndividual) state.population.subpops[0].individuals[i];
            individual.PC = indsCharListsIntermediatePop[i];
        }

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            double PCDistance = 0;
            int samePCNum = 0;
            int numEstimatedInd = 0;
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
                    if (individual.evaluated == false) {
                        numEstimatedInd++;
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

                            if (dMin == 0)
                                samePCNum++;

                        }
                        PCDistance += dMin;
                        ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index]);
                        //individual.evaluated = true;
                        //individual.fitness.trials = new ArrayList();
                        //individual.fitness.trials.add(individual.fitness.fitness());

                        estimatedFitness.add(individual.fitness.fitness());
//                        realFitness.add(realEvaluation((GPRuleEvolutionState)state,(GPIndividual)individual));

                    }
                    if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                            || state.generation == 40 || state.generation == 49) {
                        pcDistance.add(dMin);
                    }
                } else {
                    ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
            SamePCNum.add(samePCNum);
            AveragePCDistance.add(PCDistance / numEstimatedInd);
        }
    }

    public double[] evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel, double threshold) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        double[] estimatedFitnesses = new double[state.population.subpops[0].individuals.length];

        for (int i = 0; i < state.population.subpops[0].individuals.length; i++) //
        {
            GPIndividual individual = (GPIndividual) state.population.subpops[0].individuals[i];
            individual.PC = indsCharListsIntermediatePop[i];
        }

        for (int sub = 0; sub < state.population.subpops.length; sub++) {

            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
                    if (individual.evaluated == false) {
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

                        if(dMin <= threshold) {
                            estimatedFitnesses[i] = fitnessesForModel[index];
                        } else {
                            estimatedFitnesses[i] = Double.MAX_VALUE;
                        }
                    } else {
                        estimatedFitnesses[i] = Double.MAX_VALUE;
                    }
                } else {
                    estimatedFitnesses[i] = Double.MAX_VALUE;
                }
            }
        }
        return estimatedFitnesses;
    }

    public void evaluatePopulation(final EvolutionState state, int[][][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int[][][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulationRanks(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            double PCDistance = 0;
            int samePCNum = 0;
            int numEstimatedInd = 0;
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if (indsCharListsMultiTree.length != 0) {
                    //KNN
                    //===============================start==============================================
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
                    if (individual.evaluated == false) {
                        numEstimatedInd++;
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
                            if (dMin == 0)
                                samePCNum++;

                        }
                        PCDistance += dMin;
                        ((Clearable) individual.fitness).surrogateFitness(fitnessesForModel[index]);
                        //individual.evaluated = true;
                        //individual.fitness.trials = new ArrayList();
                        //individual.fitness.trials.add(individual.fitness.fitness());

                        estimatedFitness.add(individual.fitness.fitness());
                        realFitness.add(realEvaluation((GPRuleEvolutionState) state, (GPIndividual) individual));

                    }
                    if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                            || state.generation == 40 || state.generation == 49) {
                        pcDistance.add(dMin);
                    }
                } else {
                    ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
            SamePCNum.add(samePCNum);
            AveragePCDistance.add(PCDistance / numEstimatedInd);
        }
    }

    private double realEvaluation(GPRuleEvolutionState state, GPIndividual indi) {

        //then calculate the real fitness
        GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
        GPRule routingRule = new GPRule(RuleType.ROUTING, indi.trees[1]);

        RuleOptimizationProblem problem = (RuleOptimizationProblem) state.evaluator.p_problem;

        DynamicSimulation simulation = (DynamicSimulation) (((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0));
        simulation.setSequencingRule(sequencingRule); //indicate different individuals
        simulation.setRoutingRule(routingRule);

        simulation.run();
        String objectiveName = state.parameters.getStringWithDefault(new Parameter("eval.problem.eval-model.objectives.0"), null, "");
        Objective objective = Objective.get(objectiveName);

        double ObjValue = simulation.objectiveValue(objective); // this line: the value of makespan

        //in essence, here is useless. because if w.numOpsInQueue() > 100, the simulation has been canceled in run(). here is a double check
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


        return ObjValue;

    }

    //===================SVM=======================
/*    public void evaluatePopulation(final EvolutionState state, double[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        double[][] indsCharListsIntermediatePop = phenotypicForSurrogate.phenotypicPopulation(state, phenoCharacterisation, nonIntermediatePop); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsMultiTree.length != 0) {
                //===============================start==============================================
                GaussianKernel mercerKernel = new GaussianKernel(100);
                SVR.Trainer<double[]> trainer = new SVR.Trainer<double[]>(mercerKernel, 0.001, 1); //change the C from 1 to 10, no effect on fitness
                SVR<double[]> network = trainer.train(indsCharListsMultiTree, fitnessesForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {// do not use the surrogate, in this way, equals to use the original way to get offsprings
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/


    //=====================================RBF==========================================2019.9.15
/*  public void evaluatePopulation(final EvolutionState state, double[][] indsCharListsMultiTree, double[] fitnessesForModel) {

      double[][] indsCharListsIntermediatePop = phenotypicForSurrogate.phenotypicPopulation(state, phenoCharacterisation, nonIntermediatePop); //3. calculate the phenotypic characteristic

      for (int sub = 0; sub < state.population.subpops.length; sub++) {
          //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
          int maxNeighbour = (int) Math.round(0.01 * 10 * (indsCharListsMultiTree.length - 10)); //insure the RBF can work
          if (indsCharListsMultiTree.length != 0 && maxNeighbour >= 1) {
              //===============================start==============================================
              Metric<double[]> metric = new EuclideanDistance();
              RBFNetwork.Trainer<double[]> trainer = new RBFNetwork.Trainer<double[]>(metric);
              RBFNetwork<double[]> network = trainer.train(indsCharListsMultiTree, fitnessesForModel);

              for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
              {
                  Individual individual = state.population.subpops[sub].individuals[i];
                  if (individual.evaluated == false) {
                      double[] pcIntermediate = indsCharListsIntermediatePop[i];
                      //calculate the fitness based on surrogate model
                      double estimatedFitness = network.predict(pcIntermediate);
                      ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);

                        if(estimatedFitness <= 0){
                            ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                        }
                        else{
                            ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                        }
                  }
              }
          } else {
              for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
              {
                  Individual individual = state.population.subpops[sub].individuals[i];
                  if (individual.evaluated == false) {
                      ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                  }
              }
          }
      }
  }*/

//2019.9.16
    //=======================================Gaussian Process Regression==============================================
/*    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================
                MercerKernel<double[]> mercerKernel = new MercerKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };
                GaussianProcessRegression.Trainer<double[]> trainer = new GaussianProcessRegression.Trainer<double[]>(mercerKernel,0.1);
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                GaussianProcessRegression<double[]> network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/


    //LASSO
   /* public void evaluatePopulation(final EvolutionState state, double[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        double[][] indsCharListsIntermediatePop = phenotypicForSurrogate.phenotypicPopulation(state, phenoCharacterisation, nonIntermediatePop); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsMultiTree.length != 0) {
                //===============================start==============================================
                LASSO.Trainer trainer = new LASSO.Trainer(1);
                LASSO network = trainer.train(indsCharListsMultiTree, fitnessesForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {//large uncertainty
                        double[] pcIntermediate = indsCharListsIntermediatePop[i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                        if (estimatedFitness <= 0) {
                            ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE);
                        } else {
                            ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                        }
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just make it different with the duplicated ones
                    }
                }
            }
        }
    }*/

    //NN
   /* public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 ArrayList<double[][]> indsCharListsForModel) {

        double[][][] indsCharListsIntermediatePop = surrogateClearing.phenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            //if there is training data, train the model. Otherwise, do not need to train model, just assign them the same fitness (Double.MAX_VALUE)
            if (indsCharListsForModel.get(sub).length != 0) {
                //===============================start==============================================

                int numInput = indsCharListsForModel.get(sub)[0].length;
                NeuralNetwork.Trainer trainer = new NeuralNetwork.Trainer(numInput,200, 200, 200,1);//the number of input, ...the number of nodes in each layer..., the nunmber of output
                double[] fitnessForModel = ArrayUtils.toPrimitive(fitnessesForModel.get(sub));
                NeuralNetwork network = trainer.train(indsCharListsForModel.get(sub), fitnessForModel);

                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if (individual.evaluated == false) {
                        double[] pcIntermediate = indsCharListsIntermediatePop[sub][i];
                        //calculate the fitness based on surrogate model
                        double estimatedFitness = network.predict(pcIntermediate);
                        ((Clearable) individual.fitness).surrogateFitness(estimatedFitness);
                    }
                }
            } else {
                for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //512
                {
                    Individual individual = state.population.subpops[sub].individuals[i];
                    if(individual.evaluated == false){
                        ((Clearable) individual.fitness).surrogateFitness(Double.MAX_VALUE - 1000); //just remove the duplicated ones
                    }
                }
            }
        }
    }*/

    public void setClear(boolean clear) {
        this.clear = clear;
    }

}
