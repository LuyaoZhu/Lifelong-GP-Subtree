package yimei.jss.algorithm.surrogateCaseLS;

import java.util.*;

public class SpearmanCorrelationBatch {

    // 单数组秩转化
    public static double[] rank(double[] data) {
        int n = data.length;
        Double[] boxed = new Double[n];
        for (int i = 0; i < n; i++) {
            boxed[i] = data[i];
        }

        Integer[] indices = new Integer[n];
        for (int i = 0; i < n; i++) indices[i] = i;

        Arrays.sort(indices, Comparator.comparingDouble(i -> boxed[i]));

        double[] ranks = new double[n];
        int i = 0;
        while (i < n) {
            int j = i;
            while (j < n - 1 && boxed[indices[j]].equals(boxed[indices[j + 1]])) {
                j++;
            }
            double avgRank = (i + j + 2) / 2.0;
            for (int k = i; k <= j; k++) {
                ranks[indices[k]] = avgRank;
            }
            i = j + 1;
        }

        return ranks;
    }

    // 两个数组之间的 Spearman 相关系数
    public static double spearman(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Arrays must have same length.");
        }

        double[] rankX = rank(x);
        double[] rankY = rank(y);

        int n = x.length;
        double sumDiff2 = 0;
        for (int i = 0; i < n; i++) {
            double d = rankX[i] - rankY[i];
            sumDiff2 += d * d;
        }

        return 1.0 - (6 * sumDiff2) / (n * (n * n - 1));
    }

    // 输入多个surrogate case的估计fitness，输出相关系数矩阵
    public static double[][] computeSpearmanMatrix(double[][] estimatedFitness) {
        int m = estimatedFitness.length;
        double[][] corrMatrix = new double[m][m];

        for (int i = 0; i < m; i++) {
            corrMatrix[i][i] = 1.0;
            for (int j = i + 1; j < m; j++) {
                double rho = spearman(estimatedFitness[i], estimatedFitness[j]);
                corrMatrix[i][j] = rho;
                corrMatrix[j][i] = rho;
            }
        }

        return corrMatrix;
    }

    // 示例
    public static void main(String[] args) {
        double[][] estimatedFitness = {
                {1.2, 2.3, 3.1, 4.0},
                {1.1, 2.5, 3.3, 3.8},
                {4.0, 3.2, 2.0, 1.0}
        };

        double[][] matrix = computeSpearmanMatrix(estimatedFitness);

        for (double[] row : matrix) {
            System.out.println(Arrays.toString(row));
        }
    }
}

