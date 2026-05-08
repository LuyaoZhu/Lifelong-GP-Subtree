package yimei.jss.algorithm.surrogateCaseLS;

public class SilhouetteScore {
    public static double computeSilhouette(double[][] distanceMatrix, int[] labels) {
        int n = distanceMatrix.length;
        double totalScore = 0.0;

        for (int i = 0; i < n; i++) {
            double a = 0.0, b = Double.MAX_VALUE;
            int clusterI = labels[i];

            // 计算同簇内的平均距离 (a)
            int countA = 0;
            for (int j = 0; j < n; j++) {
                if (labels[j] == clusterI && i != j) {
                    a += distanceMatrix[i][j];
                    countA++;
                }
            }
            if (countA > 0) a /= countA;

            // 计算最近其他簇的平均距离 (b)
            for (int c = 0; c < max(labels) + 1; c++) {
                if (c == clusterI) continue;
                double sumB = 0.0;
                int countB = 0;
                for (int j = 0; j < n; j++) {
                    if (labels[j] == c) {
                        sumB += distanceMatrix[i][j];
                        countB++;
                    }
                }
                if (countB > 0) sumB /= countB;
                b = Math.min(b, sumB);
            }
            totalScore += (b - a) / Math.max(a, b);
        }
        return totalScore / n;
    }

    public static int max(int[] arr) {
        int max = arr[0];
        for (int n : arr) {
            if (n > max) max = n;
        }
        return max;
    }
}
