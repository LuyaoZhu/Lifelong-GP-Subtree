package yimei.jss.algorithm.lifelongGP;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.select.TournamentSelection;
import ec.util.Parameter;

import java.util.*;

/**
 * @author Luyao
 * the current task is always the first.
 */
public class EpsilonLexicaseSelection extends TournamentSelection {
    private static final long serialVersionUID = 1;

    public double sumDistance;

//    int casesNum;

    @Override
    public int produce(final int subpopulation, final EvolutionState state, final int thread) {

//        if (state.generation < 200) {
        if (state.generation < 100) {
//        if (state.generation < ((GPRuleEvolutionStateSavedSurrogateCases) state).switchGen) {
            // pick size random individuals, then pick the best.
            int index = tournamentSelection(subpopulation, state, thread);

//            ((GPIndividual)state.population.subpops[0].individuals[index]).fromTS = true;

            return index;

        } else {

            Individual[] oldinds = state.population.subpops[subpopulation].individuals.clone();
            List<Individual> allIndividuals = new ArrayList<>();

            List<Integer> allIndividualsIndex = getCandidates(subpopulation, state, thread);
            for (int i = 0; i < allIndividualsIndex.size(); i++) {
                allIndividuals.add(oldinds[allIndividualsIndex.get(i)]);
            }

            //rotate the cases to increase the diversity
            List<Integer> randomInstanceIndex = this.getRandomSequenceInstance(state, thread);

            //always let the current task as the first case
//            randomInstanceIndex.remove(state.population.subpops[0].individuals[0].fitness.trials.size()-1);
//            randomInstanceIndex.add(0, state.population.subpops[0].individuals[0].fitness.trials.size()-1);

            //-----------------------------------------------------before----------------------------------------------------------
            for (int i = 0; i < randomInstanceIndex.size(); i++) {

                if (allIndividuals.size() == 1) {
//                System.out.println("epsilon selection index: " + allIndividualsIndex.get(0));
//                    System.out.println( i + " out of " + randomInstanceIndex.size() + " cases are all used." );
                    int index = allIndividualsIndex.get(0);
                    return index;
                }

                int instanceIndex = randomInstanceIndex.get(i);

                // after choose one instance, we calculate the estimated fitness of all candidates.
                ArrayList<Double> estimatedFitness = new ArrayList<>();
                for (int ind = 0; ind < allIndividualsIndex.size(); ind++) {
                    Individual candidateInd = state.population.subpops[0].individuals[allIndividualsIndex.get(ind)];
                    estimatedFitness.add((Double) candidateInd.fitness.trials.get(instanceIndex));
                }

                double bestFitness = Collections.min(estimatedFitness);

                //calculate the epsilon------------------------
                double epsilon = calculateEpsilon(estimatedFitness);
//                if(i == randomInstanceIndex.size()-1) { //means this is the last case, we just select the min!!
//                    epsilon = 0;
//                }
//                System.out.println("Epsilon: " + epsilon);
//
                int j = 0;
                while (j < allIndividuals.size()) {
                    if (allIndividuals.size() == 1) {
//                        System.out.println("current case index is "+ i + ", the whole number of cases is " + randomInstanceIndex.size() );
//                            System.out.println("epsilon selection index: " + allIndividualsIndex.get(0));
                        return allIndividualsIndex.get(0);
                    }
                    if (estimatedFitness.get(j) > bestFitness + epsilon) {
                        estimatedFitness.remove(j);
                        allIndividualsIndex.remove(j);
                        allIndividuals.remove(j);
                    } else {
                        j++;
                    }
                }
            }

            //if there are multiple remaining individuals, we will select the individual with the best current task fitness
//            int finalIndex = 0;
//            double bestFitness = Double.MAX_VALUE;
//            for (int i = 0; i < allIndividualsIndex.size(); i++) {
//                int index = allIndividualsIndex.get(i);
//                double fitness = state.population.subpops[0].individuals[index].fitness.fitness();
//                if (fitness <= bestFitness) {
//                    bestFitness = fitness;
//                    finalIndex = index;
//                }
//            }

            int finalIndex = allIndividualsIndex.get(state.random[thread].nextInt(allIndividuals.size()));

//            System.out.println(finalIndex);
//        System.out.println("size: " + allIndividuals.size());
//        System.out.println("epsilon selection index by random: " + index);
//            System.out.println(randomInstanceIndex.size() + " cases are all used." );

            //if multiple individuals remain, we use real evaluation to select the best

//            return index;
            return finalIndex;
        }
        //--------------------------------------------------------------end--------------------------------------------------------------------

    }

    private void reverseCases(final EvolutionState state, List<Integer> randomInstanceIndex, List<Individual> allIndividuals, List<Integer> allIndividualsIndex) {

        for (int i = 0; i < randomInstanceIndex.size(); i++) {
            if (allIndividuals.size() == 1) {
//                System.out.println("epsilon selection index: " + allIndividualsIndex.get(0));
//                    System.out.println( i + " out of " + randomInstanceIndex.size() + " cases are all used." );
//                System.out.println(allIndividualsIndex.get(0));
            }

            int instanceIndex = randomInstanceIndex.get(i);

            // after choose one instance, we calculate the estimated fitness of all candidates.
            ArrayList<Double> estimatedFitness = new ArrayList<>();
            for (int ind = 0; ind < allIndividualsIndex.size(); ind++) {
                Individual candidateInd = state.population.subpops[0].individuals[allIndividualsIndex.get(ind)];

                estimatedFitness.add(((GPIndividual) candidateInd).caseFitness.get(instanceIndex));
            }

            double bestFitness = Collections.min(estimatedFitness);

            //calculate the epsilon------------------------
            double epsilon = calculateEpsilon(estimatedFitness);
//            System.out.println("Epsilon: " + epsilon);

            int j = 0;
            while (j < allIndividuals.size()) {
                if (allIndividuals.size() == 1) {
//                        System.out.println("current case index is "+ i + ", the whole number of cases is " + randomInstanceIndex.size() );
                    System.out.println(allIndividualsIndex.get(0));
                }
                if (estimatedFitness.get(j) > bestFitness + epsilon) {
                    estimatedFitness.remove(j);
                    allIndividualsIndex.remove(j);
                    allIndividuals.remove(j);
                } else {
                    j++;
                }
            }
        }
        int index = state.random[0].nextInt(allIndividuals.size());
//        System.out.println("size: " + allIndividuals.size());
//        System.out.println("epsilon selection index by random: " + allIndividualsIndex.get(index));
        System.out.println(randomInstanceIndex.size() + " cases are all used.");
//        System.out.println(allIndividualsIndex.get(index));

    }

    public static Set<Integer> generateUniqueRandomNumbers(int M, int N, long seed) {
        if (M >= N + 1) {
            throw new IllegalArgumentException("M must be less than N+1 to ensure uniqueness.");
        }

        Set<Integer> uniqueNumbers = new HashSet<>();
        Random random = new Random(seed); // 固定种子

        while (uniqueNumbers.size() < M) {
            int num = random.nextInt(N); // 生成 [0, N) 之间的随机数
            uniqueNumbers.add(num);
        }

        return uniqueNumbers;
    }

    public double calculateEpsilon(ArrayList<Double> estimatedFitness) {
        List<Double> allFitness = new ArrayList<>();
        for (int ref = 0; ref < estimatedFitness.size(); ref++) {
            allFitness.add(estimatedFitness.get(ref));
        }
        allFitness.sort(Double::compareTo);

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

    public List<Integer> getCandidates(final int subpopulation,
                                       final EvolutionState state,
                                       final int thread) {
        int populationSize = state.population.subpops[0].individuals.length;

        int permSize = (int) (populationSize * 1);
        int[] permutation = new int[permSize];
        randomPermutation(state, 0, permutation, permSize);

        // pick size random individuals, then pick the best.
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
//        int best = getRandomIndividual(0, subpopulation, state, thread);

//        int s = state.parameters.getIntWithDefault(new Parameter("candidate"), null, 800); // the number of candidates
        int s = state.parameters.getIntWithDefault(new Parameter("pop.subpop.0.size"), null, 250); // the number of candidates
//        System.out.println("tournament size for ELSMO:" + s);

        List<Integer> tourIndividuals = new ArrayList<>();

        for (int x = 0; x < permSize; x++) {
            //modified by mengxu 2021.09.24
            //need to exclude the bad fitness, or will make error when calculate the epsilon.
            int j = permutation[x];
            Individual ind = state.population.subpops[0].individuals[j];

            if (ind.fitness.fitness() >= Double.POSITIVE_INFINITY || ind.fitness.fitness() >= Double.MAX_VALUE) {
                continue;
            } else {
                tourIndividuals.add(j);
            }
            if (tourIndividuals.size() == s*0.5) {
                return tourIndividuals;
            }
        }

//        System.out.println("TourIndividuals' size has not arrive at the set size in advance.");
        return tourIndividuals;
    }

    public List<Integer> getRandomSequenceInstance(EvolutionState state, int thread) {
        List<Integer> Instance = new ArrayList<>();

        for (int i = 0; i < state.population.subpops[0].individuals[0].fitness.trials.size(); i++) {
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

    public int tournamentSelection(final int subpopulation, final EvolutionState state, final int thread) {
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


}
