package yimei.jss.algorithm.surrogateLSnew;

import ec.EvolutionState;
import ec.Individual;
import ec.simple.SimpleEvaluator;
import ec.util.Parameter;
import org.apache.commons.lang3.ArrayUtils;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.niching.Clearable;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.multitreeClearingFirst;
import yimei.jss.niching.phenotypicForSurrogate;

import java.util.ArrayList;

import static yimei.jss.algorithm.surrogateLSnew.GPRuleEvolutionStateSavedSurrogateCases.*;
import static yimei.jss.gp.GPRuleEvolutionState.phenoCharacterisation;


/**
 * Created by fzhang on 2019.9.24.
 */
//a class can not extend from more than one class
public class surrogateClearingMultitreeEvaluatorV1 extends SimpleEvaluator {
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
        }
        else{
                if (nonIntermediatePop) {
            super.evaluatePopulation(state);

            //2021.8.3 record the average fitness in the pop
                    //2021.8.3 average fitness value
                   /* ArrayList<Double> averageFitnessCurrentGen = new ArrayList<>();
                    for (Individual ind : state.population.subpops[0].individuals) {
                        if (ind.fitness.fitness() != Double.MAX_VALUE) {
                            averageFitnessCurrentGen.add(ind.fitness.fitness());
                        }
                    }
                    averageFitnessGens.add(averageFitnessCurrentGen.stream().mapToDouble(i -> i).average().orElse(0));*/

            //2021.7.7 after evaluate the individuals in the new population
/*            if(state.generation != 0){
                ArrayList<Pair<Integer, Double>> indexIndFitnessRank = new ArrayList<>();
                for (int ind = 0; ind < state.population.subpops[0].individuals.length - numElites; ind++){
                    indexIndFitnessRank.add(new ImmutablePair<>(ind,state.population.subpops[0].individuals[ind].fitness.fitness()));
                }
                indexIndFitnessRank.sort(Comparator.comparingDouble(Pair::getValue)); // in a increase order, the smaller the better

                //2021.8.3 calculate the accuracy of surrogate based on all inds in the population
                //=======================================begin=========================================
                double[] evaluatedRankWhole = new double[state.population.subpops[0].individuals.length - numElites];
                for(int i = 0; i < indexIndFitnessRank.size(); i++){
                    evaluatedRankWhole[i] = indexIndFitnessRank.get(i).getKey();
                }
                double [] preSelectedRankWhole = new double[state.population.subpops[0].individuals.length - numElites];
                for(int j = 0; j < preSelectedRankWhole.length; j++){
                    preSelectedRankWhole[j] = j;
                }

                Double correlationWhole = new SpearmansCorrelation().correlation(preSelectedRankWhole, evaluatedRankWhole);
//                System.out.println(correlation);
                correlationSurrogateAccuracyWhole.add(correlationWhole);
                //==============================end===========================

                //2021.7.28 calculate the surrogate accuracy based on top 0.3*popsize inds
                //========================================begin================================================
                double[] evaluatedRank = new double[(int)(state.population.subpops[0].individuals.length * 0.3)];
                for(int i = 0; i < evaluatedRank.length; i++){
                    evaluatedRank[i] = indexIndFitnessRank.get(i).getKey();
                }
                double [] preSelectedRank = new double[(int)(state.population.subpops[0].individuals.length * 0.3)];
                for(int j = 0; j < preSelectedRank.length; j++){
                    preSelectedRank[j] = j;
                }
                Double correlation = new SpearmansCorrelation().correlation(preSelectedRank, evaluatedRank);
//                System.out.println(correlation);
                correlationSurrogateAccuracy.add(correlation);
                //======================================end=================================================

                //2021.7.8 also define another way to measure the accuracy of the surrogate
                Double distanceRank = 0.0;
                for(int i = 0; i < preSelectedRank.length; i++){
                    distanceRank += Math.abs(preSelectedRank[i] - evaluatedRank[i]);
                }
                distanceRank = distanceRank / preSelectedRank.length;
//                System.out.println(distanceRank);
                distanceRankSurrogateAccuracy.add(distanceRank);

                //2021.8.3 calculate the surrogate accuracy based on the whole pop inds
                Double distanceRankWhole = 0.0;
                for(int i = 0; i < preSelectedRankWhole.length; i++){
                    distanceRankWhole += Math.abs(preSelectedRankWhole[i] - evaluatedRankWhole[i]);
                }
                distanceRankWhole = distanceRankWhole / preSelectedRankWhole.length;
//                System.out.println(distanceRank);
                distanceRankSurrogateAccuracyWhole.add(distanceRankWhole);
            }*/

            int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic
            double[] diversityValue = new double[state.population.subpops.length];
            diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
            entropyDiversity.add(diversityValue);
            System.out.println("Entropy diversity: " + diversityValue[0]);

//            for (int i = 0; i < state.population.subpops[0].individuals.length; i++) {
//                GPIndividual ind = (GPIndividual)state.population.subpops[0].individuals[i];
//                ind.PC = indsCharListsMultiTree[i];
//            }

            //get the fitness for training model
            double[] fitnessesForModel = new double[state.population.subpops[0].individuals.length];
            for(int ind = 0; ind < state.population.subpops[0].individuals.length; ind++){
                fitnessesForModel[ind] = state.population.subpops[0].individuals[ind].fitness.fitness();
            }

           // int removeIdx = 0;
            for (int i = 0; i < fitnessesForModel.length; i++) {
                if (fitnessesForModel[i] == Double.MAX_VALUE) {
                   /* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*/
                    indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i);
                    i--;
                   // removeIdx++;
                }
            }
            tempfitnessesForModel = fitnessesForModel;
            tempindsCharListsMultiTree = indsCharListsMultiTree;
            nonIntermediatePop = false;

            surrogateFitnessEveryGen.add(tempfitnessesForModel);
            surrogatePCEveryGen.add(tempindsCharListsMultiTree);

        } else {
            multitreeClearingFirst.clearPopulation(state, radius, capacity, phenoCharacterisation);//the individuals do not have fitness or they are evaluated yet
                                                                                                   //but has a fitness there
            this.evaluatePopulation(state, tempindsCharListsMultiTree, tempfitnessesForModel);
            nonIntermediatePop = true;
            }
        }

        }

    //==========================================KNN========================================
    public void evaluatePopulation(final EvolutionState state, int[][] indsCharListsMultiTree, double[] fitnessesForModel) {

        int[][] indsCharListsIntermediatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulation(state, phenoCharacterisation); //3. calculate the phenotypic characteristic

        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int i = 0; i < state.population.subpops[sub].individuals.length; i++) //
            {
                Individual individual = state.population.subpops[sub].individuals[i];
                if(indsCharListsMultiTree.length != 0){
                    //KNN
                    //===============================start==============================================
                    double dMin  = Double.MAX_VALUE;
                    int index = 0;
                    if(individual.evaluated == false){
                        int[] pcIntermediate = indsCharListsIntermediatePop[i];
                        //calculate the fitness based on surrogate model
                        for(int pc = 0; pc < indsCharListsMultiTree.length; pc++){
                            int[] pcModel = indsCharListsMultiTree[pc];
                            double d  = PhenoCharacterisation.distance(pcIntermediate, pcModel);
                            if(d == 0){
                                dMin = d;
                                index = pc;
                                break;
                            }

                            if (d < dMin){
                                dMin = d;
                                index = pc;
                            }
                        }
                        ((Clearable)individual.fitness).surrogateFitness(fitnessesForModel[index]);
                        //individual.evaluated = true;
                        //individual.fitness.trials = new ArrayList();
                        //individual.fitness.trials.add(individual.fitness.fitness());
                    }
                    if (state.generation == 0 || state.generation == 10 || state.generation == 20 || state.generation == 30
                            || state.generation == 40 || state.generation == 49) {
                        pcDistance.add(dMin);
                    }
                }
                else{
                    ((Clearable)individual.fitness).surrogateFitness(Double.MAX_VALUE);
                }
            }
        }
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
