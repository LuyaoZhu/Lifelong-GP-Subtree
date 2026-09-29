package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.select.TournamentSelection;
import yimei.jss.algorithm.lifelongGP.GPRuleEvolutionStateLifelongGPV10N1;

import java.util.ArrayList;
import java.util.Arrays;

public class AdaptiveParentSelection extends TournamentSelection {

    public enum ParentRole {
        PARENT1,
        PARENT2
    }

    public static final int ARM_HG_CS = 0;
    public static final int ARM_HG_CG = 1;
    public static final int ARM_CS_CG = 2;
    public static final int ARM_CG_CG = 3;

    public static final int ARM_NUM = 4;

    private static double[] rewards = new double[ARM_NUM];
    private static int[] counts = new int[ARM_NUM];
    private static int totalSelections = 0;

    public static int[] generationSelectionCounts = new int[ARM_NUM];

    private static double explorationConstant;

    private static int currentArm = ARM_HG_CS;

    public static int[] generatedByArm =
            new int[ARM_NUM];

    public static double[] currentGenerationSurvival =
            new double[ARM_NUM];

    public static double[]
            currentSelectedImportance =
            new double[2];

    public static void reset() {
        rewards = new double[ARM_NUM];
        counts = new int[ARM_NUM];
        totalSelections = 0;
        currentArm = ARM_HG_CS;
    }

    public static void resetMAB() {

        Arrays.fill(rewards, 0.0);

        Arrays.fill(counts, 0);

        Arrays.fill(
                generationSelectionCounts,
                0);

        totalSelections = 0;
    }

    public static double[] getRewards() {
        return rewards.clone();
    }

    public static void setCurrentSelectedImportance(
            double parent1,
            double parent2) {

        currentSelectedImportance[0] = parent1;

        currentSelectedImportance[1] = parent2;
    }

    public static double[] getCurrentSelectedImportance() {

        return currentSelectedImportance.clone();
    }

    public static double[] getCurrentGenerationSurvival() {

        return currentGenerationSurvival.clone();
    }

    public static void setExplorationConstant(double c) {
        explorationConstant = c;
    }

    public static int selectArm(EvolutionState state, int thread) {

        GPRuleEvolutionStateLifelongGPReplace s = (GPRuleEvolutionStateLifelongGPReplace) state;

        // =====================================
        // first lifelong generation:
        // only HG + CS
        // =====================================

        if (s.firstStageSwitchGeneration) {
            currentArm = ARM_HG_CS;
        } else {
            currentArm = selectArmUCB(state, thread);
//            currentArm = state.random[thread].nextInt(ARM_NUM);
        }

//        currentArm = ((GPRuleEvolutionStateLifelongGPV10N1)state).fixedArm;

        counts[currentArm]++;

        totalSelections++;

        generationSelectionCounts[currentArm]++;

        return currentArm;
    }

    private static int selectArmUCB(EvolutionState state, int thread) {

        for (int i = 0; i < ARM_NUM; i++) {
            if (counts[i] < state.population.subpops[thread].individuals.length * 0.1 ) {
                return i;
            }
        }

        double[] ucbValues = new double[ARM_NUM];
        explorationConstant = ((GPRuleEvolutionStateLifelongGPReplace)state).UCBExplorationConstant;

        for (int i = 0; i < ARM_NUM; i++) {
            double avgReward = rewards[i];
            double exploration =
                    explorationConstant
                            * Math.sqrt(Math.log(Math.max(totalSelections, 1)) / counts[i]);

            ucbValues[i] = avgReward + exploration;
        }

        return sampleArmBySoftmax(
                ucbValues,
                state,
                thread);

        // 选取 ucbValues 中最大值的 Arm 索引 (Argmax)
//        int bestArm = 0;
//        double maxUCB = ucbValues[0];
//        for (int i = 1; i < ARM_NUM; i++) {
//            if (ucbValues[i] > maxUCB) {
//                maxUCB = ucbValues[i];
//                bestArm = i;
//            }
//        }
//        System.out.println("bestArm = " + bestArm);
//
//        return bestArm;

    }

    private static int sampleArmBySoftmax(
            double[] ucb,
            EvolutionState state,
            int thread) {

        double tau = 0.2;

        double[] probs =
                new double[ucb.length];

        double sum = 0.0;

        for (int i = 0; i < ucb.length; i++) {

            probs[i] =
                    Math.exp(ucb[i] / tau);

            sum += probs[i];
        }

        double r =
                state.random[thread]
                        .nextDouble() * sum;

        double cumulative = 0.0;

        for (int i = 0; i < probs.length; i++) {

            cumulative += probs[i];

            if (r <= cumulative) {
                return i;
            }
        }

        return probs.length - 1;
    }

    public static GPIndividual selectParentFromArm(
            int arm,
            ParentRole role,
            ArrayList<GPIndividual> historicalGeneralists,
            ArrayList<GPIndividual> currentSpecialists,
            ArrayList<GPIndividual> currentGeneralists,
            EvolutionState state,
            int thread,
            int tournamentSize) {

        ArrayList<GPIndividual> source =
                getSourceByArm(
                        arm,
                        role,
                        historicalGeneralists,
                        currentSpecialists,
                        currentGeneralists);

        if (source == null || source.isEmpty()) {
            return null;
        }


        // =====================================
        // HG / CS:
        // random selection
        // =====================================

        if (source == historicalGeneralists
                || source == currentSpecialists) {

            int idx =
                    state.random[thread]
                            .nextInt(source.size());

            return source.get(idx);
        }

        // =====================================
        // CG:
        // tournament selection
        // =====================================

        return (GPIndividual) tournamentSelect(
                source,
                state,
                thread,
                tournamentSize);
    }

    private static ArrayList<GPIndividual> getSourceByArm(
            int arm,
            ParentRole role,
            ArrayList<GPIndividual> historicalGeneralists,
            ArrayList<GPIndividual> currentSpecialists,
            ArrayList<GPIndividual> currentGeneralists) {

        switch (arm) {
            case ARM_HG_CS:
                return role == ParentRole.PARENT1
                        ? historicalGeneralists
                        : currentSpecialists;

            case ARM_HG_CG:
                return role == ParentRole.PARENT1
                        ? historicalGeneralists
                        : currentGeneralists;

            case ARM_CS_CG:
                return role == ParentRole.PARENT1
                        ? currentSpecialists
                        : currentGeneralists;

            case ARM_CG_CG:
                return currentGeneralists;

            default:
                return null;
        }
    }

    private static Individual tournamentSelect(
            ArrayList<? extends Individual> source,
            EvolutionState state,
            int thread,
            int tournamentSize) {

        int size = Math.max(1, tournamentSize);

        int best = state.random[thread].nextInt(source.size());

        for (int i = 1; i < size; i++) {
            int challenger = state.random[thread].nextInt(source.size());

            Individual bestInd = source.get(best);
            Individual challengerInd = source.get(challenger);

            if (challengerInd.fitness.betterThan(bestInd.fitness)) {
                best = challenger;
            }
        }

        return source.get(best);
    }

    public static void updateReward(int arm, double reward) {

        if (arm < 0 || arm >= ARM_NUM) {
            return;
        }

        double alpha = GPRuleEvolutionStateLifelongGPReplace.currentRewardWeight;

        rewards[arm] =
                (1-alpha)*rewards[arm]
                        + alpha*reward;

    }

    public static int getCurrentArm() {
        return currentArm;
    }

    public static double getAverageReward(int arm) {
        if (arm < 0 || arm >= ARM_NUM) {
            return 0.0;
        }

        return rewards[arm];
    }

    public static int getCount(int arm) {
        if (arm < 0 || arm >= ARM_NUM) {
            return 0;
        }

        return counts[arm];
    }

    public static void printStatistics() {

        System.out.println(
                "========== MAB Statistics ==========");

        for (int i = 0;
             i < ARM_NUM;
             i++) {

            System.out.println(
                    "Arm "
                            + i
                            + " count="
                            + counts[i]
                            + " reward="
                            + rewards[i]);
        }
    }

}