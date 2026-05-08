package yimei.jss.niching;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.coevolve.GroupedProblemForm;
import ec.coevolve.MultiPopCoevolutionaryEvaluator;
import ec.multiobjective.MultiObjectiveFitness;
import ec.util.Parameter;
import smile.math.kernel.MercerKernel;
import smile.regression.GaussianProcessRegression;

import java.util.ArrayList;

import static yimei.jss.niching.surrogateClearing.indsCharLists;

/**
 * Created by dyska on 26/09/17.
 */
//a class can not extend from more than one class
public class surrogateClearingMultiPopCoevolutionaryEvaluator extends MultiPopCoevolutionaryEvaluator {
    public static final String P_RADIUS = "radius";
    public static final String P_CAPACITY = "capacity";

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";
    private int preGenerations;

    protected boolean clear = true;

    protected double radius;
    protected int capacity;

    protected PhenoCharacterisation[] phenoCharacterisation;

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

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        radius = state.parameters.getDoubleWithDefault(
                base.push(P_RADIUS), null, 0.0);
        capacity = state.parameters.getIntWithDefault(
                base.push(P_CAPACITY), null, 1);
        String filePath = state.parameters.getString(new Parameter("filePath"), null);
        //It's a little tricky to know whether we have 1 or 2 populations here, so we will assume
        //2 for the purpose of the phenoCharacterisation, and ignore the second object if only
        //1 is used
        phenoCharacterisation = new PhenoCharacterisation[2];
        if (filePath == null) {
            //dynamic simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation();
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation();
        } else {
            //static simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
        }
    }

    @Override
    public void evaluatePopulation(final EvolutionState state) {

        if (clear) {
            //  System.out.println("MultiPopCoevolution");
            double[][][] indsCharLists = surrogateClearing.clearPopulation(state,phenoCharacterisation);
        }

         super.evaluatePopulation(state, indsCharLists); //fzhang  all the evolution process is the same. for one individual, get a fitness value
    }

    // individuals to evaluate together
    Individual[] inds = null;
    // which individual should have its fitness updated as a result
    boolean[] updates = null;
    public static final String POP_FULL_EVALUATED_INDS = "frac-full-evaluated-inds";
    private double fracInds;
    double[] tempFitnessesForModel = null;
    double[] fitnessesForModel = null;

    public void performCoevolutionaryEvaluation( final EvolutionState state,
                                                 final Population population,
                                                 final GroupedProblemForm prob,
                                                 final double[][][] indsCharLists) {
        int evaluations = 0;

        inds = new Individual[population.subpops.length];
        updates = new boolean[population.subpops.length];//check if we need to evaluate the individuals, because it is so expensive, should be careful.
        //numCurrent = 0
        if (numCurrent > 0) {
            for (int i = 0; i < selectionMethodCurrent.length; i++) {
                selectionMethodCurrent[i].prepareToProduce(state, i, 0);  //A default version of prepareToProduce which does nothing.
            }
        }

        //==========================here skip: step 2: load the selectionMethod for previous generation==================================
        if (numPrev > 0) {
            for (int i = 0; i < selectionMethodPrev.length; i++) {
                // do a hack here
                Population currentPopulation = state.population; //state.population is current population.
                state.population = previousPopulation;
                selectionMethodPrev[i].prepareToProduce(state, i, 0);  //A default version of prepareToProduce which does nothing.  here state is important, it is the state
                //for previous population
                state.population = currentPopulation; //currentPopulaiton is a temp parameter to save current population
            }
        }

        //step 3: build subPopulaiton, subpops[0] and subpops[1]
        // build subpopulation array to pass in each time
        int[] subpops = new int[state.population.subpops.length];
        //System.out.println(subpops.length);
        for(int j = 0; j < subpops.length; j++) {
            subpops[j] = j;
        }

        //System.out.println(prob);  //yimei.jss.ruleoptimisation.RuleCoevolutionProblem@5ae63ade
        //here skip: step 3: setup the shuffle: here num-shuffled = 0  here, we do not use it.
        if (numShuffled > 0) {
            int[/*numShuffled*/][/*subpop*/][/*shuffledIndividualIndexes*/] ordering = null;
            // build shuffled orderings
            ordering = new int[numShuffled][state.population.subpops.length][state.population.subpops[0].individuals.length]; // if num-shuffled =1  [1][2][512]
            for(int c = 0; c < numShuffled; c++)
                for(int m = 0; m < state.population.subpops.length; m++)
                {
                    for(int i = 0; i < state.population.subpops[0].individuals.length; i++)
                        ordering[c][m][i] = i; //ordering[0][0][0] = 0, ordering[0][0][1] = 1, ordering[0][0][2] = 2, ordering[0][0][3] = 3
                    if (m != 0)
                        shuffle(state, ordering[c][m]); //ordering = new int[numShuffled][state.population.subpops.length]
                }

            // for each individual
            for (int i = 0; i < state.population.subpops[0].individuals.length; i++) {
                for (int k = 0; k < numShuffled; k++) {
                    for (int ind = 0; ind < inds.length; ind++) {
                        inds[ind] = state.population.subpops[ind].individuals[ordering[k][ind][i]];
                        updates[ind] = true;
                    }
                    prob.evaluate(state, inds, updates, false, subpops, 0);  ////yimei.jss.ruleoptimisation.RuleCoevolutionProblem@5ae63ade
                    evaluations++;
                }
            }
        }

        fracInds = state.parameters.getDoubleWithDefault(new Parameter(POP_FULL_EVALUATED_INDS),null,0.0);
        for (int j = 0; j < state.population.subpops.length; j++) //0,1
        {
            if (!shouldEvaluateSubpop(state, j, 0)) continue;  // don't evaluate this subpopulation

            // for each individual
            tempFitnessesForModel = new double[(int)(state.population.subpops[j].individuals.length*fracInds)];
            for (int i = 0; i < (int)(state.population.subpops[j].individuals.length*fracInds); i++) //512
            {
                Individual individual = state.population.subpops[j].individuals[i];
                // Test against all the elites
                for (int k = 0; k < eliteIndividuals[j].length; k++) { //2
                    for (int ind = 0; ind < inds.length; ind++) { //2
                        if (ind == j) {   //j = 0, 1  (ind j) ---> (0 0) or (1 1) that is to say, this is the subpopulation1
                            inds[ind] = individual; //inds[0] = individual = state.population.subpops[0].individuals[0];
                            updates[ind] = true;   // updates[0] = true    updates[1] = true   evaluate
                        }
                        else {  // this is subpopulation2
                            inds[ind] = eliteIndividuals[ind][k];   // (ind j) ---> (0 1) or (1 0)
                            //inds[1] = eliteIndividuals[1][*]   inds[0] = eliteIndividuals[0][*]
                            updates[ind] = false;  // do not evaluate
                        }
                    }

                    prob.evaluate(state,inds,updates, false, subpops, 0);
                    tempFitnessesForModel[i] = individual.fitness.fitness();
                    evaluations++;
                }
                //System.out.println(evaluations);  //4  8  12 16 20 24 28 32 ... 4096   2*512*4 = 4096  inds[] is used to save the individuals we want to evaluated.

                //here, skip this part: test against random selected individuals of the current population
                for(int k = 0; k < numCurrent; k++) //0  skip this part
                {
                    for(int ind = 0; ind < inds.length; ind++) //2
                    {
                        if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                        else { inds[ind] = produceCurrent(ind, state, 0); updates[ind] = true; }
                    }
                    prob.evaluate(state,inds,updates, false, subpops, 0);
                    evaluations++;
                }

                // here, skip this part. Test against random individuals of previous population
                for(int k = 0; k < numPrev; k++)  // 0  skip this part
                {
                    for(int ind = 0; ind < inds.length; ind++)
                    {
                        if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                        else { inds[ind] = producePrevious(ind, state, 0); updates[ind] = false; }
                    }
                    prob.evaluate(state,inds,updates, false, subpops, 0);
                    evaluations++;
                }
            }

            //fzhang 2019.8.25 choose the individuals from the first half population to predict the fitness of other individuals
            int countNonDoubleMax = 0;
            for(int count = 0; count < tempFitnessesForModel.length; count ++){
                if(tempFitnessesForModel[count] != Double.MAX_VALUE){
                    countNonDoubleMax++;
                }
            }

//            int maxNeighbour = (int) Math.round(0.01 * 10 * (countNonDoubleMax - 10));
//            if (maxNeighbour >= 1){
            if (countNonDoubleMax >= 1){
                double[][][] indsCharListsForModel = new double[state.population.subpops.length][countNonDoubleMax][];
                fitnessesForModel = new double[countNonDoubleMax];
                int countGoodInds = 0;
                for (int idxGoodInds = 0; idxGoodInds < (int)(state.population.subpops[j].individuals.length*fracInds); idxGoodInds++){
                    if(tempFitnessesForModel[idxGoodInds] != Double.MAX_VALUE){
                        indsCharListsForModel[j][countGoodInds] = indsCharLists[j][idxGoodInds].clone();
                        fitnessesForModel[countGoodInds] = tempFitnessesForModel[idxGoodInds];
                        countGoodInds++;
                    }
                }

                //evolve model
                //RBF
//                Metric<double[]> metric = new EuclideanDistance();
//                RBFNetwork.Trainer<double[]> trainer = new RBFNetwork.Trainer<double[]>(metric);
//                RBFNetwork<double[]> network = trainer.train(indsCharListsForModel[j], fitnessesForModel);

                //SVM
             /*   MercerKernel<double[]> mercerKernel = new MercerKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };
                SVR.Trainer<double[]> trainer = new SVR.Trainer<double[]>(mercerKernel,1,1.0);
                SVR<double[]> network = trainer.train(indsCharListsForModel[j], fitnessesForModel);*/

             //Gaussian Process Regression
                MercerKernel<double[]> mercerKernel = new MercerKernel<double[]>() {
                    @Override
                    public double k(double[] x, double[] y) {
                        return 0;
                    }
                };
                GaussianProcessRegression.Trainer<double[]> trainer = new GaussianProcessRegression.Trainer<double[]>(mercerKernel,0.1);
                GaussianProcessRegression<double[]> network = trainer.train(indsCharListsForModel[j], fitnessesForModel);


                //===============================it is the same for every model builder========================================
                for (int remainInds = (int)(state.population.subpops[j].individuals.length*fracInds); remainInds < state.population.subpops[j].individuals.length; remainInds++) {
                    Individual estimatedInd = state.population.subpops[j].individuals[remainInds];
                    double[] remainInd = indsCharLists[j][remainInds];
                    double estimatedFitness = network.predict(remainInd);
                    ((Clearable)state.population.subpops[j].individuals[remainInds].fitness).surrogateFitness(estimatedFitness);
                    ((MultiObjectiveFitness)(estimatedInd.fitness)).trials = new ArrayList();
                    estimatedInd.fitness.trials.add(estimatedInd.fitness.fitness());
                }
            }
            else// in case there is no training data
                {
                for (int i = (int)(state.population.subpops[j].individuals.length*fracInds); i < (int)(state.population.subpops[j].individuals.length); i++) //512
                {
                    Individual individual = state.population.subpops[j].individuals[i];
                    // Test against all the elites
                    for (int k = 0; k < eliteIndividuals[j].length; k++) { //2
                        for (int ind = 0; ind < inds.length; ind++) { //2
                            if (ind == j) {   //j = 0, 1  (ind j) ---> (0 0) or (1 1) that is to say, this is the subpopulation1
                                inds[ind] = individual; //inds[0] = individual = state.population.subpops[0].individuals[0];
                                updates[ind] = true;   // updates[0] = true    updates[1] = true   evaluate
                            }
                            else {  // this is subpopulation2
                                inds[ind] = eliteIndividuals[ind][k];   // (ind j) ---> (0 1) or (1 0)
                                //inds[1] = eliteIndividuals[1][*]   inds[0] = eliteIndividuals[0][*]
                                updates[ind] = false;  // do not evaluate
                            }
                        }

                        prob.evaluate(state,inds,updates, false, subpops, 0);
                        evaluations++;
                    }
                    //System.out.println(evaluations);  //4  8  12 16 20 24 28 32 ... 4096   2*512*4 = 4096  inds[] is used to save the individuals we want to evaluated.

                    //here, skip this part: test against random selected individuals of the current population
                    for(int k = 0; k < numCurrent; k++) //0  skip this part
                    {
                        for(int ind = 0; ind < inds.length; ind++) //2
                        {
                            if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                            else { inds[ind] = produceCurrent(ind, state, 0); updates[ind] = true; }
                        }
                        prob.evaluate(state,inds,updates, false, subpops, 0);
                        evaluations++;
                    }

                    // here, skip this part. Test against random individuals of previous population
                    for(int k = 0; k < numPrev; k++)  // 0  skip this part
                    {
                        for(int ind = 0; ind < inds.length; ind++)
                        {
                            if (ind == j) { inds[ind] = individual; updates[ind] = true; }
                            else { inds[ind] = producePrevious(ind, state, 0); updates[ind] = false; }
                        }
                        prob.evaluate(state,inds,updates, false, subpops, 0);
                        evaluations++;
                    }
                }
            }
        }
        //============================================================================================================================

        //here, skip this part
        // now shut down the selection methods
        if (numCurrent > 0)
            for( int i = 0 ; i < selectionMethodCurrent.length; i++)
                selectionMethodCurrent[i].finishProducing( state, i, 0 );  //A default version of finishProducing, which does nothing.

        if (numPrev > 0)
            for( int i = 0 ; i < selectionMethodPrev.length ; i++ )
            {
                // do a hack here
                Population currentPopulation = state.population;
                state.population = previousPopulation;
                selectionMethodPrev[i].finishProducing( state, i, 0 );
                state.population = currentPopulation;
            }

        state.output.message("Evaluations: " + evaluations);
    }

    public void setClear(boolean clear) {
        this.clear = clear;
    }
}
