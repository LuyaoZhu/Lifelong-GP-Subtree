package yimei.jss.algorithm.lifelongGP;

import ec.EvolutionState;
import ec.Individual;
import ec.select.TournamentSelection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * @author Luyao
 * the current task is always the first.
 */
public class ParetoTournamentSelection extends TournamentSelection {
    private static final long serialVersionUID = 1;

    public double sumDistance;

//    int casesNum;

    @Override
    public int produce(final int subpopulation, final EvolutionState state, final int thread) {

        int index;

        if(state.generation<100){
            index = tournamentSelection(subpopulation, state, thread);
        } else {
            index = paretoTournamentSelection(subpopulation, state, thread);
        }

        return index;

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

    public int paretoTournamentSelection(final int subpopulation, final EvolutionState state, final int thread) {
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;

        List<Integer> candidateIndices = new ArrayList<>();
        List<Individual> candidates = new ArrayList<>();
        int s = getTournamentSizeToUse(state.random[thread]);

        //randomly get s candidates
        while (candidates.size() < s) {
            int j = getRandomIndividual(candidates.size(), subpopulation, state, thread);
            if(oldinds[j].fitness.fitness() < Double.MAX_VALUE){
                candidateIndices.add(j);
                candidates.add(oldinds[j]);
            }
        }

        // 找非支配解
        List<Integer> frontLocalIdx = getNonDominatedIndices(candidates);
        List<Individual> front = new ArrayList<>();
        for (int idx : frontLocalIdx) front.add(candidates.get(idx));

        if (front.size() == 1) {
            return candidateIndices.get(frontLocalIdx.get(0)); // 只有一个直接返回
        } else if (front.size() == 2) {
            if(oldinds[candidateIndices.get(frontLocalIdx.get(0))].fitness.fitness() < oldinds[candidateIndices.get(frontLocalIdx.get(1))].fitness.fitness()){
                return candidateIndices.get(frontLocalIdx.get(0));
            } else {
                return candidateIndices.get(frontLocalIdx.get(1));
            }
        } else {
            // 计算拥挤度，选拥挤度最xiao
            double[] crowding = computeCrowdingDistance(front);
            double minDistance = Double.MAX_VALUE;
            int selectedIdxInFront = 0;
            for (int i = 0; i < front.size(); i++) {
                if (crowding[i] < minDistance) {
                    minDistance = crowding[i];
                    selectedIdxInFront = i;
                }
            }
            return candidateIndices.get(frontLocalIdx.get(selectedIdxInFront));
        }

    }


    public static boolean dominates(Individual a, Individual b) {
        boolean betterInAny = false;
        for (int i = 0; i < a.fitness.trials.size(); i++) {
            if (((double)a.fitness.trials.get(i)) > ((double)b.fitness.trials.get(i))) return false; // a 在某目标更差
            if (((double)a.fitness.trials.get(i)) < ((double)b.fitness.trials.get(i))) betterInAny = true;
        }
        return betterInAny;
    }

    // 找局部 Front 1
    public static List<Integer> getNonDominatedIndices(List<Individual> candidates) {
        List<Integer> front = new ArrayList<>();
        for (int i = 0; i < candidates.size(); i++) {
            Individual ind = candidates.get(i);
            boolean dominated = false;
            for (int j = 0; j < candidates.size(); j++) {
                if (i == j) continue;
                if (dominates(candidates.get(j), ind)) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) front.add(i);
        }
        return front;
    }

    // 计算拥挤度，返回与 front 对应的 double[]，不修改 individual
    public static double[] computeCrowdingDistance(List<Individual> frontIndividuals) {
        int size = frontIndividuals.size();
        int numObj = frontIndividuals.get(0).fitness.trials.size();
        double[] distance = new double[size];
        Arrays.fill(distance, 0);

        for (int m = 0; m < numObj; m++) {
            final int obj = m;
            List<Integer> sortedIdx = new ArrayList<>();
            for (int i = 0; i < size; i++) sortedIdx.add(i);

            // 按当前目标排序
            sortedIdx.sort(Comparator.comparingDouble(i -> (double) frontIndividuals.get(i).fitness.trials.get(obj)));

            distance[sortedIdx.get(0)] = Double.POSITIVE_INFINITY;
            distance[sortedIdx.get(size - 1)] = Double.POSITIVE_INFINITY;

            double fMin = (double) frontIndividuals.get(sortedIdx.get(0)).fitness.trials.get(obj);
            double fMax = (double) frontIndividuals.get(sortedIdx.get(size - 1)).fitness.trials.get(obj);
            if (fMax - fMin == 0) continue;

            for (int i = 1; i < size - 1; i++) {
                double nextVal = (double) frontIndividuals.get(sortedIdx.get(i + 1)).fitness.trials.get(obj);
                double prevVal = (double) frontIndividuals.get(sortedIdx.get(i - 1)).fitness.trials.get(obj);
                distance[sortedIdx.get(i)] += (nextVal - prevVal) / (fMax - fMin);
            }
        }
        return distance;
    }


}
