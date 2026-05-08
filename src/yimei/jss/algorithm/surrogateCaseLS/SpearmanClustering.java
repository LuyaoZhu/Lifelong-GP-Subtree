package yimei.jss.algorithm.surrogateCaseLS;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.util.*;

public class SpearmanClustering {

    // 计算 Spearman 相关系数距离矩阵（1 - 相关系数）
    public static double[][] computeSpearmanDistanceMatrix(double[][] fitnessMatrix) {

        // 👉 在方法内部进行转置（列变成行，行变成列）
        int rows = fitnessMatrix.length;        // the number of instances
        int cols = fitnessMatrix[0].length;     // the number of individuals
        double[][] transposed = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = fitnessMatrix[i][j];
            }
        }

        // 👉 现在 transposed 是 10×19，列为 instance，Spearman 将比较列之间的相关性
        RealMatrix dataMatrix = new BlockRealMatrix(transposed);  // 10×19
        SpearmansCorrelation spearman = new SpearmansCorrelation();
        RealMatrix corrMatrix = spearman.computeCorrelationMatrix(dataMatrix);  // 得到 19×19

        int n = corrMatrix.getRowDimension();  // ✅ 应该是 19
        double[][] distMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distMatrix[i][j] = 1.0 - corrMatrix.getEntry(i, j);
            }
        }
        return distMatrix;
    }

    // 简单的基于距离矩阵的层次聚类（最小距离合并法）
    public static Map<Integer, List<Integer>> agglomerativeClustering(double[][] distMatrix, int k) {
        int n = distMatrix.length;
        List<Set<Integer>> clusters = new ArrayList<>();
        for (int i = 0; i < n; i++) clusters.add(new HashSet<>(Collections.singleton(i)));

        while (clusters.size() > k) {
            double minDist = Double.MAX_VALUE;
            int c1 = -1, c2 = -1;
            for (int i = 0; i < clusters.size(); i++) {
                for (int j = i + 1; j < clusters.size(); j++) {
                    double dist = averageLinkage(clusters.get(i), clusters.get(j), distMatrix);
                    if (dist < minDist) {
                        minDist = dist;
                        c1 = i;
                        c2 = j;
                    }
                }
            }
            clusters.get(c1).addAll(clusters.get(c2));
            clusters.remove(c2);
        }

        Map<Integer, List<Integer>> result = new HashMap<>();
        for (int i = 0; i < clusters.size(); i++) {
            result.put(i, new ArrayList<>(clusters.get(i)));
        }
        return result;
    }

    // 计算两个 cluster 之间的平均距离（average linkage）
    private static double averageLinkage(Set<Integer> cluster1, Set<Integer> cluster2, double[][] distMatrix) {
        double total = 0;
        int count = 0;
        for (int i : cluster1) {
            for (int j : cluster2) {
                total += distMatrix[i][j];
                count++;
            }
        }
        return total / count;
    }

    public static void main(String[] args) {
        // 示例：26 个 instance，每个有 100 个规则的 fitness 表现
        int instances = 26;
        int rules = 100;
        double[][] fitnessMatrix = new double[instances][rules];
        Random rand = new Random(42);
        for (int i = 0; i < instances; i++) {
            for (int j = 0; j < rules; j++) {
                fitnessMatrix[i][j] = rand.nextDouble();
            }
        }

        double[][] distanceMatrix = computeSpearmanDistanceMatrix(fitnessMatrix);
        Map<Integer, List<Integer>> clusters = agglomerativeClustering(distanceMatrix, 4);

        for (Map.Entry<Integer, List<Integer>> entry : clusters.entrySet()) {
            System.out.println("Cluster " + entry.getKey() + ": " + entry.getValue());
        }
    }
}

