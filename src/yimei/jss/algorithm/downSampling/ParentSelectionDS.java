package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Individual;
import ec.SelectionMethod;
import ec.select.SelectDefaults;
import ec.steadystate.SteadyStateBSourceForm;
import ec.steadystate.SteadyStateEvolutionState;
import ec.util.MersenneTwisterFast;
import ec.util.Parameter;
import yimei.jss.ruleoptimisation.MultipleTreeRuleOptimizationProblem;

import java.util.ArrayList;
import java.util.List;

public class ParentSelectionDS extends SelectionMethod implements SteadyStateBSourceForm {

    /**
     * default base
     */
    public static final String P_LEXICASE = "lexicase";

    //    public static final String P_NUMCASE = "caseInOneInstance";
    public int numCase;
    public static final String P_PRO_LEXICASE = "pro-lexicase";
    public double proLexicase;
    public static final String P_RANDOM_PRO = "random-pro";
    public String randomPro;

    public static final String P_ELITIST_SELECTION_RATE = "elitist-selection-rate";
    public double elitistSR;

//    public int genForUseLS = 5;

    /** tournament */
    /**
     * default base
     */
    public static final String P_TOURNAMENT_SIZE = "select.tournament.size";

    public static final String P_CANDIDATE_SIZE = "candidate-size";

    public static final String P_PICKWORST = "pick-worst";

    public static final String USE_TS = "useTS";

    public static final String USE_ParetoTS = "useParetoTS";

    /** size parameter */
//    public static final String P_SIZE = "size";

    /**
     * Base size of the tournament; this may change.
     */
    int tournamentSize;

    int candidateSize;

    /**
     * Probablity of picking the size plus one
     */
    public double probabilityOfPickingSizePlusOne;

    /**
     * Do we pick the worst instead of the best?
     */
    public boolean pickWorst;

    public boolean useTS;
    public boolean useParetoTS;
    public boolean useRandomDownSampling;
    public boolean useInformedDownSampling;

    public static final String RANDOM_DOWN_SAMPLING = "randomDownSampling";
    public static final String INFORMED_DOWN_SAMPLING = "informedDownSampling";
//    public Parameter defaultBase()
//    {
//        return SelectDefaults.base().push(P_TOURNAMENT);
//    }

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

//        int numIns;
        Parameter def = defaultBase();
//        numIns = state.parameters.getInt(new Parameter(P_NUMCASE), def.push(P_NUMCASE), 0);
//        if (numIns <= 0)
//            state.output.fatal("The number of cases must be an integer >= 1.", base.push(P_NUMCASE), def.push(P_NUMCASE));

//        this.numCase = state.parameters.getInt(new Parameter(P_NUMCASE), null, 1);

        proLexicase = state.parameters.getDouble(new Parameter(P_PRO_LEXICASE), null, 0);
        randomPro = state.parameters.getStringWithDefault(new Parameter(P_RANDOM_PRO), null, "no");


        //tournament
        double val = state.parameters.getDouble(new Parameter(P_TOURNAMENT_SIZE), null, 1.0);
        if (val < 1.0)
            state.output.fatal("Tournament size must be >= 1.");
        else if (val == (int) val)  // easy, it's just an integer
        {
            tournamentSize = (int) val;
            probabilityOfPickingSizePlusOne = 0.0;
        } else {
            tournamentSize = (int) Math.floor(val);
            probabilityOfPickingSizePlusOne = val - tournamentSize;  // for example, if we have 5.4, then the probability of picking *6* is 0.4
        }

        //candidate for LS
        candidateSize = state.parameters.getInt(new Parameter(P_CANDIDATE_SIZE), null, 1);

        pickWorst = state.parameters.getBoolean(base.push(P_PICKWORST), def.push(P_PICKWORST), false);

        elitistSR = state.parameters.getDoubleWithDefault(new Parameter(P_ELITIST_SELECTION_RATE), null, 1.0);

        useTS = state.parameters.getBoolean(new Parameter(USE_TS), new Parameter(USE_TS), false);
        useParetoTS = state.parameters.getBoolean(new Parameter(USE_ParetoTS), new Parameter(USE_ParetoTS), false);

        useRandomDownSampling = state.parameters.getBoolean(base.push(RANDOM_DOWN_SAMPLING), new Parameter(RANDOM_DOWN_SAMPLING), false);
        useInformedDownSampling = state.parameters.getBoolean(new Parameter(INFORMED_DOWN_SAMPLING), new Parameter(INFORMED_DOWN_SAMPLING), false);

    }

    //todo: maybe can be modified!
    public List<Integer> getRandomSequenceInstance(EvolutionState state, int thread) {
        List<Integer> Instance = new ArrayList<>();
        for (int i = 0; i < this.numCase; i++) {
            Instance.add(i);
        }
        List<Integer> randomInstance = new ArrayList<>();
        while (!Instance.isEmpty()) {
            int index = state.random[thread].nextInt(Instance.size());
            randomInstance.add(Instance.get(index));
            Instance.remove(index);
        }

////        test
//        for(int i=0; i<randomInstance.size();i++){
//            System.out.print(randomInstance.get(i)+", ");
//        }
//        System.out.println();

        return randomInstance;
    }

    @Override
    public int produce(int subpopulation, EvolutionState state, int thread) { //epsilon lexicase selection

        // pick size random individuals, then pick the best.
        int index;
        numCase = state.population.subpops[0].individuals[0].fitness.trials.size();

        if (state.generation < ((GPRuleEvolutionStateDownSampling) state).switchGen) {
            index = produceTournament(subpopulation, state, thread);
        } else {
            if (useTS) {
                index = tournamentSelection(subpopulation, state, thread);
            } else if (useParetoTS) {
                index = paretoTournamentSelection(subpopulation, state, thread);
            } else {
                index = produceLS(subpopulation, state, thread);
            }
        }

        return index;

    }


    //add 2021.11.19
    int bestIndexOfIndividual(List<Individual> individuals, List<Integer> individualIndex) {
        int index = individualIndex.get(0);
        Individual bestIndividual = individuals.get(0);
        for (int i = 1; i < individuals.size(); i++) {
            if (individuals.get(i).fitness.betterThan(bestIndividual.fitness)) {
                index = individualIndex.get(i);
                //modified by mengxu 2021.11.29
                bestIndividual = individuals.get(i);
            }
        }
        return index;
    }


    public int produceLS(int subpopulation, EvolutionState state, int thread) { //epsilon lexicase selection
        Individual[] oldinds = state.population.subpops[subpopulation].individuals.clone();
        List<Individual> allIndividuals = new ArrayList<>();
        List<Integer> allIndividualsIndex = getCandidates(subpopulation, state, thread);
        for (int i = 0; i < allIndividualsIndex.size(); i++) {
            allIndividuals.add(oldinds[allIndividualsIndex.get(i)]);
        }

        //determine which cases are used to select parents
//        int numSelectedCases = (int) Math.round(state.parameters.getDouble(new Parameter("down-sampling-ratio"),null)*numCase);

        double ratio = state.parameters.getDouble(new Parameter("down-sampling-ratio"), null);
        int numSelectedCases = (int) Math.ceil(ratio * numCase); // 向上取整

        List<Integer> finalSelectedCases = new ArrayList<>();

//        List<Integer> finalSelectedCases =  ((GPRuleEvolutionStateDownSampling)state).finalSelectedCasesIndex;

        //2. with random down-sampling
/*        if(useRandomDownSampling) {
            while (finalSelectedCases.size() < numSelectedCases) {
                int randomCaseIndex = state.random[thread].nextInt(numCase);
                if (!finalSelectedCases.contains(randomCaseIndex)) {
                    finalSelectedCases.add(randomCaseIndex);
                }
            }
        } else if(useInformedDownSampling) { //2. informed down sampling
            finalSelectedCases = selectCases(((GPRuleEvolutionStateDownSampling) state).spearmanMatrix, numSelectedCases, 1, state);
        } else {
            //1.without down-sampling
            finalSelectedCases = getRandomSequenceInstance(state,0);
        }*/

        finalSelectedCases = getRandomSequenceInstance(state, 0);

        for (int i = 0; i < finalSelectedCases.size(); i++) {
            if (allIndividuals.size() == 1) {
//                System.out.println("epsilon selection index: " + allIndividualsIndex.get(0));
                return allIndividualsIndex.get(0);
            }

            int instanceIndex = finalSelectedCases.get(i);

            int indexBest = getIndexOfbestIndsInstance(allIndividuals, instanceIndex);
            double bestFitness = (double) allIndividuals.get(indexBest).fitness.trials.get(instanceIndex);

            //calculate the epsilon------------------------
            //todo: need to check if the location to calculate epsilon is right? 2021.10.14
            double epsilon = calculateEpsilon(allIndividuals, instanceIndex);
//            System.out.println("Epsilon: " + epsilon);

            int j = 0;
            while (j < allIndividuals.size()) {
                if (allIndividuals.size() == 1) {
//                    System.out.println("epsilon selection index: " + allIndividualsIndex.get(0));

                    return allIndividualsIndex.get(0);
                }
                if ((double) allIndividuals.get(j).fitness.trials.get(instanceIndex) > bestFitness + epsilon) {
                    allIndividualsIndex.remove(j);
                    allIndividuals.remove(j);
                } else {
                    j++;
                }
            }
        }

        //modified by mengxu 2021.11.19
        //modify this part by adding a strategy to select individual based on the number of best case of fitness.
//        int index = bestIndexOfIndividual(allIndividuals, allIndividualsIndex);
//        System.out.println("size: " + allIndividuals.size());
//        if(state instanceof GPRuleEvolutionStateMV0){
//            ((GPRuleEvolutionStateMV0)state).parentIndex.add(index);
//        }
//        return index;

        int index = state.random[thread].nextInt(allIndividuals.size());
//        System.out.println("size: " + allIndividuals.size());
//        System.out.println("epsilon selection index by random: " + allIndividualsIndex.get(index));
        return allIndividualsIndex.get(index);
    }

    public int produceTournament(final int subpopulation,
                                 final EvolutionState state,
                                 final int thread) {
        // pick size random individuals, then pick the best.
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
        int best = getRandomIndividual(0, subpopulation, state, thread);

        int s = getTournamentSizeToUse(state.random[thread]);

        if (pickWorst)
            for (int x = 1; x < s; x++) {
                int j = getRandomIndividual(x, subpopulation, state, thread);
                if (!betterThan(oldinds[j], oldinds[best], subpopulation, state, thread))  // j is at least as bad as best
                    best = j;
            }
        else
            for (int x = 1; x < s; x++) {
                int j = getRandomIndividual(x, subpopulation, state, thread);
                if (betterThan(oldinds[j], oldinds[best], subpopulation, state, thread))  // j is better than best
                    best = j;
            }

        return best;
    }

    public double calculateEpsilon(List<Individual> allIndividuals, int i) {
        List<Double> allFitness = new ArrayList<>();
        for (int ref = 0; ref < allIndividuals.size(); ref++) {
            double fit = (double) allIndividuals.get(ref).fitness.trials.get(i);
            allFitness.add(fit);
        }
        allFitness.sort(Double::compareTo);

//        //test
//        System.out.println("Fit sort: ");
//        for(double fit:allFitness){
//            System.out.print(fit + ", ");
//        }
//        System.out.println();

        double median;
        if (allFitness.size() % 2 == 0) {
            int indexMedian = allFitness.size() / 2;
            median = (allFitness.get(indexMedian - 1) + allFitness.get(indexMedian)) / 2;
        } else {
            int indexMedian = allFitness.size() / 2;
            median = allFitness.get(indexMedian);
        }
        for (int ref2 = 0; ref2 < allFitness.size(); ref2++) {
            allFitness.set(ref2, Math.abs(allFitness.get(ref2) - median));
        }
        allFitness.sort(Double::compareTo);
        if (allFitness.size() % 2 == 0) {
            int indexMedian = allFitness.size() / 2;
            median = (allFitness.get(indexMedian - 1) + allFitness.get(indexMedian)) / 2;
        } else {
            int indexMedian = allFitness.size() / 2;
            median = allFitness.get(indexMedian);
        }
        return median;
    }


    public double calculateEpsilonNormalisation(List<double[][]> allIndividualsNormalisation, int i) {
        List<Double> allFitness = new ArrayList<>();
        for (int ref = 0; ref < allIndividualsNormalisation.size(); ref++) {
            double fit = allIndividualsNormalisation.get(ref)[i][0];
            allFitness.add(fit);
        }
        allFitness.sort(Double::compareTo);

//        //test
//        System.out.println("Fit sort: ");
//        for(double fit:allFitness){
//            System.out.print(fit + ", ");
//        }
//        System.out.println();

        double median;
        if (allFitness.size() % 2 == 0) {
            int indexMedian = allFitness.size() / 2;
            median = (allFitness.get(indexMedian - 1) + allFitness.get(indexMedian)) / 2;
        } else {
            int indexMedian = allFitness.size() / 2;
            median = allFitness.get(indexMedian);
        }
        for (int ref2 = 0; ref2 < allFitness.size(); ref2++) {
            allFitness.set(ref2, Math.abs(allFitness.get(ref2) - median));
        }
        allFitness.sort(Double::compareTo);
        if (allFitness.size() % 2 == 0) {
            int indexMedian = allFitness.size() / 2;
            median = (allFitness.get(indexMedian - 1) + allFitness.get(indexMedian)) / 2;
        } else {
            int indexMedian = allFitness.size() / 2;
            median = allFitness.get(indexMedian);
        }
        return median;
    }


    /**
     * Returns a tournament size to use, at random, based on base size and probability of picking the size plus one.
     */
    public int getTournamentSizeToUse(MersenneTwisterFast random) {
        double p = probabilityOfPickingSizePlusOne;   // pulls us to under 35 bytes
        if (p == 0.0) return tournamentSize;
        return tournamentSize + (random.nextBoolean(p) ? 1 : 0);
    }

    public int getCandidateSizeToUse(MersenneTwisterFast random) {
        return candidateSize;
    }


    /**
     * Produces the index of a (typically uniformly distributed) randomly chosen individual
     * to fill the tournament.  <i>number</> is the position of the individual in the tournament.
     */
    public int getRandomIndividual(int number, int subpopulation, EvolutionState state, int thread) {
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
        return state.random[thread].nextInt(oldinds.length);
    }

    /**
     * Returns true if *first* is a better (fitter, whatever) individual than *second*.
     */
    public boolean betterThan(Individual first, Individual second, int subpopulation, EvolutionState state, int thread) {
        return first.fitness.betterThan(second.fitness);
    }

    public static void randomPermutation(EvolutionState state, int thread, int[] perm, int size) {
//    JMetalRandom randomGenerator = JMetalRandom.getInstance();
        int[] index = new int[size];
        boolean[] flag = new boolean[size];

        for (int n = 0; n < size; n++) {
            index[n] = n;
            flag[n] = true;
        }

        int num = 0;
        while (num < size) {
            int start = state.random[thread].nextInt(size);//todo: need to check
//      int start = randomGenerator.nextInt(0, size - 1);
            while (true) {
                if (flag[start]) {
                    perm[num] = index[start];
                    flag[start] = false;
                    num++;
                    break;
                }
                if (start == (size - 1)) {
                    start = 0;
                } else {
                    start++;
                }
            }
        }
    }

    public List<Integer> getCandidates(final int subpopulation,
                                       final EvolutionState state,
                                       final int thread) {
        int populationSize = state.population.subpops[0].individuals.length;

        int permSize = (int) (populationSize * elitistSR);
        int[] permutation = new int[permSize];
        randomPermutation(state, 0, permutation, permSize);

        // pick size random individuals, then pick the best.
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
        int best = getRandomIndividual(0, subpopulation, state, thread);

        int s = getCandidateSizeToUse(state.random[thread]);
//        System.out.println("tournament size for ELSMO:" + s);

        List<Integer> tourIndividuals = new ArrayList<>();

        for (int x = 0; x < permSize; x++) {
            //modified by mengxu 2021.09.24
            //need to exclude the bad fitness, or will make error when calculate the epsilon.
            int j = permutation[x];
            Individual ind = state.population.subpops[0].individuals[j];
            double[][] multiCaseFitness = new double[ind.fitness.trials.size()][1];
            for (int k = 0; k < ind.fitness.trials.size(); k++) {
                multiCaseFitness[k][0] = (double) ind.fitness.trials.get(k);
            }
            if (multiCaseFitness[0][0] >= Double.POSITIVE_INFINITY || multiCaseFitness[0][0] >= Double.MAX_VALUE) {
                continue;
            } else {
                tourIndividuals.add(j);
            }
            if (tourIndividuals.size() == s) {
                return tourIndividuals;
            }
        }

//        System.out.println("TourIndividuals' size has not arrive at the set size in advance.");
        return tourIndividuals;
    }

    public int getIndexOfbestIndsInstance(List<Individual> individuals, int indexInstance) {

        int best_index = 0;
        double best_fitness = (double) individuals.get(0).fitness.trials.get(indexInstance);
        for (int ind = 0; ind < individuals.size(); ind++) {
            if ((double) individuals.get(ind).fitness.trials.get(indexInstance) < best_fitness) {
                best_fitness = (double) individuals.get(ind).fitness.trials.get(indexInstance);
                best_index = ind;
            }
        }
        return best_index;

    }

    public ArrayList<Integer> selectCases(double[][] spearmanMatrix, int maxCount, double corrThreshold, final EvolutionState state) {
        int n = spearmanMatrix.length;
        double[][] distanceMatrix = new double[n][n];

        // Spearman → 距离矩阵: 1 - |ρ|
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distanceMatrix[i][j] = 1 - Math.abs(spearmanMatrix[i][j]);
            }
        }

        ArrayList<Integer> selected = new ArrayList<>();

        // ① 随机选择第一个

        int firstCase = state.random[0].nextInt(n);
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

            selected.add(bestCase);
        }

        return selected;
    }

    public int tournamentSelection(final int subpopulation, final EvolutionState state, final int thread) {
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
        int best = getRandomIndividual(0, subpopulation, state, thread);

        int s = getTournamentSizeToUse(state.random[thread]);

        //determine which cases are used to select parents
//        int numSelectedCases = (int) Math.round(state.parameters.getDouble(new Parameter("down-sampling-ratio"),null)*numCase);

        double ratio = state.parameters.getDouble(new Parameter("down-sampling-ratio"), null);
        int numSelectedCases = (int) Math.ceil(ratio * numCase); // 向上取整

        ArrayList<Integer> finalSelectedCases = new ArrayList<>();
//        ArrayList<Integer> finalSelectedCases = ((GPRuleEvolutionStateDownSampling)state).finalSelectedCasesIndex;

        //2. with random down-sampling
/*        if(useRandomDownSampling) {
            while (finalSelectedCases.size() < numSelectedCases) {
                int randomCaseIndex = state.random[thread].nextInt(numCase);
                if (!finalSelectedCases.contains(randomCaseIndex)) {
                    finalSelectedCases.add(randomCaseIndex);
                }
            }
        } else if(useInformedDownSampling) { //2. informed down sampling
            finalSelectedCases = selectCases(((GPRuleEvolutionStateDownSampling) state).spearmanMatrix, numSelectedCases, 1, state);
        } else {
            //1.without down-sampling
            finalSelectedCases = new ArrayList<>(numCase); // 预分配容量
            for (int i = 0; i < numCase; i++) {
                finalSelectedCases.add(i);
            }
        }*/

        finalSelectedCases = new ArrayList<>(numCase); // 预分配容量
        for (int i = 0; i < numCase; i++) {
            finalSelectedCases.add(i);
        }

        if (pickWorst)
            for (int x = 1; x < s; x++) {
                int j = getRandomIndividual(x, subpopulation, state, thread);
                if (!betterThan(oldinds[j], oldinds[best], subpopulation, state, thread, finalSelectedCases))  // j is at least as bad as best
                    best = j;
            }
        else
            for (int x = 1; x < s; x++) {
                int j = getRandomIndividual(x, subpopulation, state, thread);
                if (betterThan(oldinds[j], oldinds[best], subpopulation, state, thread, finalSelectedCases))  // j is better than best
                    best = j;
            }

        return best;
    }

    public int paretoTournamentSelection(final int subpopulation, final EvolutionState state, final int thread) {
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;

        int s = getTournamentSizeToUse(state.random[thread]);

        int[] candidates = new int[s];
        for (int i = 0; i < s; i++) {
            candidates[i] = getRandomIndividual(i, subpopulation, state, thread);
        }
        //determine which cases are used to select parents
//        int numSelectedCases = (int) Math.round(state.parameters.getDouble(new Parameter("down-sampling-ratio"),null)*numCase);

        double ratio = state.parameters.getDouble(new Parameter("down-sampling-ratio"), null);
        int numSelectedCases = (int) Math.ceil(ratio * numCase); // 向上取整

        ArrayList<Integer> finalSelectedCases = new ArrayList<>();
//        ArrayList<Integer> finalSelectedCases = ((GPRuleEvolutionStateDownSampling)state).finalSelectedCasesIndex;

        //2. with random down-sampling
/*        if(useRandomDownSampling) {
            while (finalSelectedCases.size() < numSelectedCases) {
                int randomCaseIndex = state.random[thread].nextInt(numCase);
                if (!finalSelectedCases.contains(randomCaseIndex)) {
                    finalSelectedCases.add(randomCaseIndex);
                }
            }
        } else if(useInformedDownSampling) { //2. informed down sampling
            finalSelectedCases = selectCases(((GPRuleEvolutionStateDownSampling) state).spearmanMatrix, numSelectedCases, 1, state);
        } else {
            //1.without down-sampling
            finalSelectedCases = new ArrayList<>(numCase); // 预分配容量
            for (int i = 0; i < numCase; i++) {
                finalSelectedCases.add(i);
            }
        }*/

        finalSelectedCases = new ArrayList<>(numCase); // 预分配容量
        for (int i = 0; i < numCase; i++) {
            finalSelectedCases.add(i);
        }

        // 取出候选的 trials（值越小越好）
        double[][] vals = new double[s][numCase];
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < oldinds[candidates[i]].fitness.trials.size(); j++) {
                vals[i][j] = (double) oldinds[candidates[i]].fitness.trials.get(j); // 若需强制转换，请按你的Fitness类型改写
            }
            if (vals[i] == null) vals[i] = new double[]{Double.POSITIVE_INFINITY};
        }

        // 朴素O(s^2) 判定：被任何候选支配则标记为 dominated
        boolean[] dominated = new boolean[s];
        for (int i = 0; i < s; i++) {
            if (dominated[i]) continue;
            for (int j = 0; j < s; j++) {
                if (i == j) continue;
                // 判断 j 是否支配 i（全维不差且至少一维严格更好；越小越好）
                double[] a = vals[j];
                double[] b = vals[i];
                int m = Math.min(a.length, b.length);
                boolean notWorseAll = true;
                boolean strictlyBetter = false;
                for (int k = 0; k < m; k++) {
                    if (a[k] > b[k]) {
                        notWorseAll = false;
                        break;
                    }
                    if (a[k] < b[k]) strictlyBetter = true;
                }
                if (notWorseAll && strictlyBetter) {
                    dominated[i] = true;
                    break;
                }
            }
        }

        // 收集 Pareto front 在 candidates 中的位置
        java.util.ArrayList<Integer> front = new java.util.ArrayList<>();
        for (int i = 0; i < s; i++) if (!dominated[i]) front.add(i);

        // 从前沿中随机选一个（若意外为空，则在所有候选中随机）
        int winnerPos = front.isEmpty() ? state.random[thread].nextInt(s)
                : front.get(state.random[thread].nextInt(front.size()));

        // 最终选择的父代（在 subpopulation 中的全局索引）
        int chosen = candidates[winnerPos];

        return chosen;
    }

    public boolean betterThan(Individual first,
                              Individual second,
                              int subpopulation,
                              EvolutionState state,
                              int thread,
                              ArrayList<Integer> selectedCases) {
        // 1) 判别目标方向（示例：名字以 "max" 开头视为最大化）
        boolean isMaximizing =
                ((MultipleTreeRuleOptimizationProblem) state.evaluator.p_problem)
                        .getEvaluationModel()
                        .getObjectives()               // 注意：一般应为 getObjectives()
                        .get(0)
                        .getName()
                        .toLowerCase()
                        .startsWith("max");

        // 2) 聚合器初始化
        double aggFirst = 0.0;
        double aggSecond = 0.0;

        if (selectedCases == null || selectedCases.isEmpty()) {
            // 无可比对样本时，保守返回 false（认为 first 不优于 second）
            return false;
        }

        // 3) 遍历所选 case，聚合两者的表现
        for (int c = 0; c < selectedCases.size(); c++) {
            int caseIndex = selectedCases.get(c);

            // 安全取值：trials 里常为 Number 或 Double
            double f1 = ((Number) first.fitness.trials.get(caseIndex)).doubleValue();
            double f2 = ((Number) second.fitness.trials.get(caseIndex)).doubleValue();

/*            if (isMaximizing) {
                // 最大化：取各自的最大值作为代表
                if (f1 > aggFirst)  aggFirst  = f1;
                if (f2 > aggSecond) aggSecond = f2;
            } else {
                // 最小化：累计求和（如需平均可在循环后除以样本数）
                aggFirst  += f1;
                aggSecond += f2;
            }*/

            aggFirst += f1;
            aggSecond += f2;

        }

        // 如果你希望最小化时用“平均”而不是“和”，取消下面两行注释：
        // if (!isMaximizing) {
        //     double n = selectedCases.size();
        //     aggFirst /= n; aggSecond /= n;
        // }

        // 4) 比较与返回（相等视为“不更好”）
        if (aggFirst < aggSecond) {
            return true;
        } else
            return false;
    }

    @Override
    public Parameter defaultBase() {
        return SelectDefaults.base().push(P_LEXICASE);
    }

    @Override
    public void individualReplaced(SteadyStateEvolutionState state, int subpopulation, int thread, int individual) {

    }

    @Override
    public void sourcesAreProperForm(SteadyStateEvolutionState state) {

    }
}
