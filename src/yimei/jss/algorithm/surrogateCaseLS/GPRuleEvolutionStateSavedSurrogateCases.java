package yimei.jss.algorithm.surrogateCaseLS;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.SortComparatorL;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.util.Pair;
import smile.clustering.HierarchicalClustering;
import smile.clustering.linkage.CompleteLinkage;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.phenotypicForSurrogate;
import yimei.jss.rule.AbstractRuleHelper;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static yimei.jss.gp.GPRun.out_dir;

/**
 * The evolution state of evolving dispatching rules with GP.
 *
 * @author yimei
 */

public class GPRuleEvolutionStateSavedSurrogateCases extends GPRuleEvolutionState {

    public int parentSelectionNum;
    public ArrayList<Integer> firstSelectedCasesIndex = new ArrayList<>();
    public double sumDistanceInReal;
    public ArrayList<Integer> finalSelectedCasesIndex = new ArrayList<>();
    ArrayList<double[]> surrogateFitnessEveryGen = new ArrayList<>();
    ArrayList<int[][]> surrogatePCEveryGen = new ArrayList<>();
    ArrayList<Integer> uniquePCNum = new ArrayList<>();

    ArrayList<Double> averageDistance = new ArrayList<>();

    ArrayList<Integer> caseNum = new ArrayList<>();

    public ArrayList<Integer> offspringLS = new ArrayList<>();
    public ArrayList<Integer> offspringTS = new ArrayList<>();

    public ArrayList<Integer> offspringMixed = new ArrayList<>();

    Map<List<Integer>, List<Integer>> pcToIndexMap = new HashMap<>();
    public int switchGen;

    double QTS;
    double QLS;

    double useLS;

    ArrayList<double[]> spearmanRealFitnessCurrentGen = new ArrayList<>();
    ArrayList<double[]> spearmanEstimatedFitnessCurrentGen = new ArrayList<>();
    ArrayList<double[]> spearmanEstimatedRealFitness = new ArrayList<>();

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

        String filePath = this.parameters.getString(new Parameter("filePath"), null);//有无filepath决定PhenoCharacterisation（静态 or 动态）
        phenoCharacterisation = new PhenoCharacterisation[2];
        //this is the baseline PhenoCharacterisation with baseline rule "SPT" "WIQ", it will be set again by the best rule, so it is useful here
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

        switchGen = parameters.getIntWithDefault(new Parameter("switchGen"), null, 5);

    }

    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        //record the diversity
        int[][] indsCharListsMultiTree = phenotypicForSurrogate.muchBetterPhenotypicPopulation(this, phenoCharacterisation);
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree);
        System.out.println(diversityValue[0]);
        entropyDiversity.add(diversityValue);

        //System.out.println("generation "+generation);
        // EVALUATION

        statistics.preEvaluationStatistics(this);

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

        //adaptively change the ratio


        // 记录唯一PC及对应fitness值列表
        Map<List<Integer>, List<Double>> pcToFitnessMap = new HashMap<>();

        for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {

            GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
            individual.PC = indsCharListsMultiTree[ind].clone();  //give PC to individuals
            List<Integer> pcKey = arrayToList(individual.PC);
//            pcToIndexMap.putIfAbsent(pcKey, new ArrayList<>());
//            pcToIndexMap.get(pcKey).add(ind);
//            pcToFitnessMap.putIfAbsent(pcKey, new ArrayList<>()); // 初始化列表
//            pcToFitnessMap.get(pcKey).add(individual.fitness.fitness()); // 添加fitness
            if (!(individual.fitness.fitness() == Double.MAX_VALUE)) {
                pcToIndexMap.putIfAbsent(pcKey, new ArrayList<>());
                pcToIndexMap.get(pcKey).add(ind);
                pcToFitnessMap.putIfAbsent(pcKey, new ArrayList<>()); // 初始化列表
                pcToFitnessMap.get(pcKey).add(individual.fitness.fitness()); // 添加fitness
            }
        }

        //-----------------------begin--------------------------
        //after evaluate all individuals, we record the information to build surrogate

        int[][] indsCharListsMultiTreeCopy = new int[pcToFitnessMap.size()][indsCharListsMultiTree[0].length];
        double[] fitnessesForModelCopy = new double[pcToFitnessMap.size()];
        int[][] indexArray = new int[pcToIndexMap.size()][];

        int PCIndex = 0;
        for (List<Integer> key : pcToFitnessMap.keySet()) {
            indsCharListsMultiTreeCopy[PCIndex] = key.stream().mapToInt(i -> i).toArray();
            List<Double> fitnessList = pcToFitnessMap.get(key);
//            double fitness = fitnessList.get(0);       // choose the first
//            double fitness = fitnessList.stream().mapToDouble(d -> d).average().orElse(0.0);    // choose average
            double fitness = fitnessList.stream().mapToDouble(d -> d).min().getAsDouble();
            fitnessesForModelCopy[PCIndex] = fitness;
            indexArray[PCIndex] = pcToIndexMap.get(key).stream().mapToInt(i -> i).toArray();
            PCIndex++;
        }

/*        for (int i = 0; i < fitnessesForModelCopy.length; i++) {
            if (fitnessesForModelCopy[i] == Double.MAX_VALUE) {
                   *//* indsCharListsMultiTree = ArrayUtils.remove(indsCharListsMultiTree, i - removeIdx);
                    fitnessesForModel = ArrayUtils.remove(fitnessesForModel, i - removeIdx);*//*
                indsCharListsMultiTreeCopy = ArrayUtils.remove(indsCharListsMultiTreeCopy, i);
                fitnessesForModelCopy = ArrayUtils.remove(fitnessesForModelCopy, i);
                i--;
                // removeIdx++;
            }
        }*/

//        System.out.println(Arrays.stream(fitnessesForModelCopy).average().getAsDouble());

        surrogateFitnessEveryGen.add(fitnessesForModelCopy);
        surrogatePCEveryGen.add(indsCharListsMultiTreeCopy);
        uniquePCNum.add(fitnessesForModelCopy.length);
        System.out.println("The number of unique PC is " + fitnessesForModelCopy.length);

        //-----------------------end--------------------------

        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;

        //------------------------begin to select surrogate cases----------------------
        if (generation >= switchGen ) {

            // for each individual, use KNN to estimate their fitness
            for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                individual.caseFitness = new ArrayList<>(); // to save the estimated fitness in all cases
                individual.PCdistance = new ArrayList<>(); // to save the estimated fitness in all cases
                individual.realCaseFitness = new ArrayList<>(Collections.nCopies(generation, 0.0));
            }

            //pcDistance

            double[][] PCdistance = new double[surrogatePCEveryGen.size()-1][pcToIndexMap.size()];

            //we calculate the average distance between candidates and the nearest sample (measure the accuracy of surrogate)

            List<Pair<ArrayList<Double>, ArrayList<Double>>> fitnessPair = new ArrayList<>();

            for (int instance = 0; instance < surrogatePCEveryGen.size() - 1; instance++) {

                ArrayList<Double> estimatedFitness = new ArrayList<>();
                ArrayList<Double> realFitness = new ArrayList<>();

                double distanceInOneCase = 0;
                int num = 0;
                for (Map.Entry<List<Integer>, List<Integer>> pcToIndex : pcToIndexMap.entrySet()) {
                    int[] estimatedPC = pcToIndex.getKey().stream().mapToInt(i -> i).toArray();
                    double dMin = Double.MAX_VALUE;
                    int index = 0;
                    for (int pc = 0; pc < surrogatePCEveryGen.get(instance).length; pc++) {
                        int[] pcModel = surrogatePCEveryGen.get(instance)[pc];
//                        double d = PhenoCharacterisation.hammingDistance(estimatedPC, pcModel);
                        double d = PhenoCharacterisation.distance(estimatedPC, pcModel);
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

                    if(dMin == 0){
                        estimatedFitness.add(surrogateFitnessEveryGen.get(instance)[index]);
                        realFitness.add(population.subpops[0].individuals[pcToIndex.getValue().get(0).intValue()].fitness.fitness());

                        //----------------------

                        /*int indsIndex = pcToIndex.getValue().get(0);

                        GPIndividual indi = (GPIndividual) population.subpops[0].individuals[indsIndex];
                        GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
                        GPRule routingRule = new GPRule(RuleType.ROUTING, indi.trees[1]);

                        long seed = 10000 * instance;
                        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);

                        simulation.reset(seed);
                        simulation.setSequencingRule(sequencingRule); //indicate different individuals
                        simulation.setRoutingRule(routingRule);

                        simulation.run();
                        String objectiveName = parameters.getStringWithDefault(new Parameter("eval.problem.eval-model.objectives.0"), null, "");
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

                        simulation.reset();*/
                        //-----------------------

                    }

                    PCdistance[instance][num] = dMin;

                        distanceInOneCase += dMin;
//                    individual.caseFitness.add(surrogateFitnessEveryGen.get(instance)[index]);
                        for (int indsIndex : pcToIndex.getValue()) {
                            ((GPIndividual) population.subpops[0].individuals[indsIndex]).caseFitness.add(surrogateFitnessEveryGen.get(instance)[index]);
                            ((GPIndividual) population.subpops[0].individuals[indsIndex]).PCdistance.add(dMin);
                        }
                    num++;
                }
                averageDistance.add(distanceInOneCase / pcToIndexMap.size());
//                estimatedRankInAllCases[instance] = computeRanks(estimatedFitnessInAllCases[instance]);

                Pair<ArrayList<Double>, ArrayList<Double>> fitnessPairCurrentGen = new Pair<>(estimatedFitness, realFitness);

                fitnessPair.add(fitnessPairCurrentGen);

            }

            double meanSameSize = 0;
            for(int a=0; a<fitnessPair.size(); a++){
                meanSameSize += fitnessPair.get(a).getFirst().size();
            }
            meanSameSize /= fitnessPair.size();

            for (int i=0; i<fitnessPair.size(); i++) {
                if(fitnessPair.get(i).getFirst().size() > meanSameSize){ //do regression
//                    if(i >= generation-15){

                    finalSelectedCasesIndex.add(i);

//                    double[] coefLinear = LinearRegression.fit(fitnessPair.get(i).getSecond(),fitnessPair.get(i).getFirst());
                    double[] coef = QuadraticRegression.fit(fitnessPair.get(i).getSecond(),fitnessPair.get(i).getFirst());

                    System.out.printf("拟合公式: y = %.3f * x^2 + %.3f * x + %.3f\n", coef[0], coef[1], coef[2]);

                    KernelRegressor reg = new KernelRegressor(fitnessPair.get(i).getSecond(),fitnessPair.get(i).getFirst(), 0.8);

                    //after doing regression, we need to estimate some inds

                    for (Map.Entry<List<Integer>, List<Integer>> pcToIndex : pcToIndexMap.entrySet()) {
                        int indsIndex = pcToIndex.getValue().get(0);
                        double distance = ((GPIndividual)population.subpops[0].individuals[indsIndex]).PCdistance.get(i);

                        //-----------------------begin: investigate the gap between real fitness and estimated fitness
/*                            GPIndividual indi = (GPIndividual) population.subpops[0].individuals[indsIndex];
                            GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, indi.trees[0]);
                            GPRule routingRule = new GPRule(RuleType.ROUTING, indi.trees[1]);

                            long seed = 10000 * i;
                            DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);

                            simulation.reset(seed);
                            simulation.setSequencingRule(sequencingRule); //indicate different individuals
                            simulation.setRoutingRule(routingRule);

                            simulation.run();
                            String objectiveName = parameters.getStringWithDefault(new Parameter("eval.problem.eval-model.objectives.0"), null, "");
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

                            for (int index : pcToIndex.getValue()) {
                                ((GPIndividual)population.subpops[0].individuals[index]).caseFitness.set(i,ObjValue);
                            }*/

                        //------------------------------end-------------------------------------------
                            if(!(distance == 0)){
                                double newX  = population.subpops[0].individuals[indsIndex].fitness.fitness();
//                                double predictedYLinear = LinearRegression.predict(coefLinear[0],coefLinear[1], newX);
                                double predictedY = QuadraticRegression.predict(coef, newX);
//                                double predictedY = reg.predict(newX);
                                for (int index : pcToIndex.getValue()) {
                                    ((GPIndividual)population.subpops[0].individuals[index]).caseFitness.set(i,predictedY);
                                }
                            }
                    }
                }
            }

            //assign the real fitness in the current gen to caseFitness
            for (int ind = 0; ind < this.population.subpops[0].individuals.length; ind++) {
                GPIndividual individual = (GPIndividual) population.subpops[0].individuals[ind];
                individual.caseFitness.add(individual.fitness.fitness());
            }

            finalSelectedCasesIndex.add(generation);


            //after get estimated fitness, we need to calculate the spearman correlation
            double[][] estimatedFitness = new double[finalSelectedCasesIndex.size()][pcToIndexMap.size()];
            for (int a=0; a<finalSelectedCasesIndex.size(); a++){
                int index = 0;
                for (Map.Entry<List<Integer>, List<Integer>> pcToIndex : pcToIndexMap.entrySet()) {
                    int indsIndex = pcToIndex.getValue().get(0);
                    estimatedFitness[a][index] = ((GPIndividual)population.subpops[0].individuals[indsIndex]).caseFitness.get(finalSelectedCasesIndex.get(a));
                    index++;
                }
            }

            double[][] matrix = computeSpearmanMatrix(estimatedFitness);
            for (double[] row : matrix) {
                System.out.println(Arrays.toString(row));
            }

        }

        System.out.println("The number of cases is " + finalSelectedCasesIndex.size());


//---------------------end------------------------------
//use different cases to generate more individuals, and

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

//            writeCaseNumToFile(caseNum);

            writeUniquePCNumToFile(uniquePCNum);

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

        population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N

        long seed = 10000 * generation;
        DynamicSimulation simulation = (DynamicSimulation) ((MultipleTreeMultipleRuleEvaluationModel) problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(0);
        simulation.reset(seed);

        averageDistance.clear();
        firstSelectedCasesIndex.clear();
        finalSelectedCasesIndex.clear();
        pcToIndexMap.clear();

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
            statistics. preCheckpointStatistics(this);
            Checkpoint.setCheckpoint(this);
            statistics.postCheckpointStatistics(this);
        }

        return R_NOTDONE;
    }

    private double calculateReward(ArrayList<Integer> offspring) {

        List<Double> allFitness = new ArrayList<>();    // the smaller, the better
        List<Double> allDiversity = new ArrayList<>();  // the smaller, the better

        for (int i = 0; i < offspring.size(); i++) {
            double fitness = population.subpops[0].individuals[offspring.get(i)].fitness.fitness();
            if (fitness < Double.MAX_VALUE) {
                allFitness.add(fitness);

                //calculate inds' diversity
                int[] basePC = ((GPIndividual) population.subpops[0].individuals[offspring.get(i)]).PC;
                double diversity = 0;
                for (Individual individual : population.subpops[0].individuals) {
                    int[] otherPC = ((GPIndividual) individual).PC;
                    double distance = PhenoCharacterisation.hammingDistance(basePC, otherPC);
                    diversity += (40 - distance);

                }

                allDiversity.add(diversity);

            }
        }
//        List<Double> normFitness = zScoreNormalize(allFitness);
//        List<Double> normDiversity = zScoreNormalize(allDiversity);

        System.out.println(" Average fitness is " + allFitness.stream().mapToDouble(Double::doubleValue).average());
        System.out.println(" Average Diversity is " + allDiversity.stream().mapToDouble(Double::doubleValue).average());

        List<Double> logFitness = logTransform(allFitness);
        List<Double> normFitness = minMaxNormalize(logFitness);
        List<Double> normDiversity = minMaxNormalize(allDiversity);

        double weight = ((double) generation / numGenerations);
        double meanNormFitness = 5 * normFitness.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double meanNormDiversity = normDiversity.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);

        double reward = weight * meanNormFitness + (1 - weight) * meanNormDiversity;

        return reward;
    }

    private double[][] filterArray(double[][] A, List<Integer> indicesToKeep) {
        double[][] B = new double[indicesToKeep.size()][];
        for (int i = 0; i < indicesToKeep.size(); i++) {
            B[i] = A[indicesToKeep.get(i)]; // 直接引用 A 中的行
        }
        return B;
    }

    public String findTreeString(String tree, int f) {
        String start = "T" + String.valueOf(f) + ":";
        int startIndex = tree.indexOf(start);

        int endIndex;

        if (f == 2) {
            endIndex = tree.length();
        } else {
            String end = "T" + String.valueOf(f + 1);
            endIndex = tree.indexOf(end);
            if (endIndex == -1)
                endIndex = tree.length();
        }

        if (startIndex != -1 && endIndex != -1) {
            return tree.substring(startIndex + start.length(), endIndex);
        } else {
            return "";
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

    public void writeUniquePCNumToFile(ArrayList<Integer> uniquePCNum) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".uniquePCNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!uniquePCNum.isEmpty()) {
                for (int j = 0; j < 1; j++) {
                    header.append(",uniquePCNum").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < uniquePCNum.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                line.append(", ").append(uniquePCNum.get(i));
                writer.write(line.toString());
                writer.newLine();
            }

            caseNum.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeCaseNumToFile(ArrayList<Integer> caseNum) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".caseNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!caseNum.isEmpty()) {
                for (int j = 0; j < 1; j++) {
                    header.append(",caseNum").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < caseNum.size(); i++) {
                StringBuilder line = new StringBuilder(i + switchGen + "");
                line.append(", ").append(caseNum.get(i));
                writer.write(line.toString());
                writer.newLine();
            }

            caseNum.clear();
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

    // 进行层次聚类
    public static int[] hierarchicalClustering(double[][] distanceMatrix, int numClusters) {
        // 使用 Complete Linkage 进行层次聚类
        HierarchicalClustering hc = new HierarchicalClustering(new CompleteLinkage(distanceMatrix));

        // 分割成 numClusters 个簇
        return hc.partition(numClusters);
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

}


