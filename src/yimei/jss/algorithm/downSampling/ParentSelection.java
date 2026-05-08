package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Individual;
import ec.select.TournamentSelection;
import ec.util.Parameter;
import yimei.jss.ruleoptimisation.MultipleTreeRuleOptimizationProblem;

import java.util.ArrayList;

public class ParentSelection extends TournamentSelection {

    public int produce(final int subpopulation, final EvolutionState state, final int thread) {

        // pick size random individuals, then pick the best.
        int index;

        if(state.generation < ((GPRuleEvolutionStateDownSampling)state).switchGen){
            index = super.produce(subpopulation, state, thread);
        } else {
            index = tournamentSelection(subpopulation, state, thread);
        }

            return index;
        }

    public int tournamentSelection(final int subpopulation, final EvolutionState state, final int thread) {
        Individual[] oldinds = state.population.subpops[subpopulation].individuals;
        int best = getRandomIndividual(0, subpopulation, state, thread);

        int s = getTournamentSizeToUse(state.random[thread]);

        //determine which cases are used to select parents
        int numSelectedCases = (int) Math.round(state.parameters.getDouble(new Parameter("down-sampling-ratio"),null)*((GPRuleEvolutionStateDownSampling)state).caseInOneInstance);
        ArrayList<Integer> finalSelectedCases = new ArrayList<>();

        //1. random down sampling
        while (finalSelectedCases.size() < numSelectedCases) {
            int randomCaseIndex = state.random[thread].nextInt(((GPRuleEvolutionStateDownSampling)state).caseInOneInstance);
            if(!finalSelectedCases.contains(randomCaseIndex)) {
                finalSelectedCases.add(randomCaseIndex);
            }
        }
        //2. informed down sampling
        finalSelectedCases = selectCases(((GPRuleEvolutionStateDownSampling)state).spearmanMatrix, numSelectedCases, 1, state);

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

    public boolean betterThan(Individual first,
                              Individual second,
                              int subpopulation,
                              EvolutionState state,
                              int thread,
                              ArrayList<Integer> selectedCases)
    {
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
        double aggFirst  = isMaximizing ? Double.NEGATIVE_INFINITY : 0.0;
        double aggSecond = isMaximizing ? Double.NEGATIVE_INFINITY : 0.0;

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

            if (isMaximizing) {
                // 最大化：取各自的最大值作为代表
                if (f1 > aggFirst)  aggFirst  = f1;
                if (f2 > aggSecond) aggSecond = f2;
            } else {
                // 最小化：累计求和（如需平均可在循环后除以样本数）
                aggFirst  += f1;
                aggSecond += f2;
            }
        }

        // 如果你希望最小化时用“平均”而不是“和”，取消下面两行注释：
        // if (!isMaximizing) {
        //     double n = selectedCases.size();
        //     aggFirst /= n; aggSecond /= n;
        // }

        // 4) 比较与返回（相等视为“不更好”）
        if(aggFirst < aggSecond) {
            return true;
        } else
            return false;
    }

}
