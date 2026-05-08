package yimei.jss.algorithm.surrogateCaseLS;

import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import java.util.Arrays;
import java.util.Random;

public class Kmeans {


    /**
     * K-Means 聚类核心函数
     *
     * @param data   原始数据，形状为 [nSamples, nFeatures]
     * @param k      聚类簇数
     * @param maxIter 最大迭代次数
     * @param tol    收敛阈值，当质心移动距离小于该阈值时停止
     * @return       每个样本对应的聚类标签 (0~k-1)
     */
    public static int[] kmeans(double[][] data, int k, int maxIter, double tol) {
        int nSamples = data.length;
        int nFeatures = data[0].length;

        // 随机初始化 k 个质心（可替换为更好的初始化方法，如kmeans++）
        double[][] centroids = new double[k][nFeatures];
        Random rand = new Random(1); // 固定随机种子

        // 随机挑选 k 个不同的样本点作为初始质心
        int[] initialIdx = randomUniqueIndices(nSamples, k, rand);
        for (int i = 0; i < k; i++) {
            centroids[i] = Arrays.copyOf(data[initialIdx[i]], nFeatures);
        }

        // 每个样本的当前聚类标签
        int[] labels = new int[nSamples];
        // 用来标记每次迭代后质心的移动量
        double centroidShift = Double.MAX_VALUE;

        // 迭代
        for (int iter = 0; iter < maxIter && centroidShift > tol; iter++) {
            // 1) 将每个样本归到与其最近的质心簇
            for (int i = 0; i < nSamples; i++) {
                labels[i] = findNearestCentroid(data[i], centroids);
            }

            // 2) 更新质心（计算每个簇的平均向量）
            double[][] newCentroids = new double[k][nFeatures];
            int[] counts = new int[k]; // 记录每个簇中样本数

            for (int i = 0; i < nSamples; i++) {
                int clusterId = labels[i];
                counts[clusterId]++;
                for (int f = 0; f < nFeatures; f++) {
                    newCentroids[clusterId][f] += data[i][f];
                }
            }

            for (int c = 0; c < k; c++) {
                if (counts[c] == 0) {
                    // 若某簇无分配样本，则可随机重置，或继续沿用旧质心
                    // 这里简单演示为：随机重置
                    newCentroids[c] = Arrays.copyOf(data[rand.nextInt(nSamples)], nFeatures);
                } else {
                    for (int f = 0; f < nFeatures; f++) {
                        newCentroids[c][f] /= counts[c];
                    }
                }
            }

            // 3) 计算本次迭代前后所有质心的移动距离
            centroidShift = 0.0;
            for (int c = 0; c < k; c++) {
                double dist = 1 - new SpearmansCorrelation().correlation(centroids[c], newCentroids[c]);
                centroidShift = Math.max(centroidShift, dist);
            }

            // 更新质心
            centroids = newCentroids;
        }

        return labels;
    }

    /**
     * 计算向量 x 到 k 个质心中距离最小的质心索引
     */
    private static int findNearestCentroid(double[] x, double[][] centroids) {
        int nearestIndex = 0;
        double minDist = Double.MAX_VALUE;

        for (int i = 0; i < centroids.length; i++) {
            double dist = 1 - new SpearmansCorrelation().correlation(x, centroids[i]);
            if (dist < minDist) {
                minDist = dist;
                nearestIndex = i;
            }
        }
        return nearestIndex;
    }

    /**
     * 计算向量 x 和 y 的欧几里得距离
     */
    private static double euclideanDistance(double[] x, double[] y) {
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            double diff = x[i] - y[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    /**
     * 在 [0, n) 区间内随机选择 k 个不重复的索引
     */
    private static int[] randomUniqueIndices(int n, int k, Random rand) {
        // 简易方法：先生成 [0,1,2,...,n-1]，再随机打乱，取前k个
        int[] indices = new int[n];
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }
        for (int i = n - 1; i >= 1; i--) {
            int j = rand.nextInt(i + 1);
            // swap
            int tmp = indices[i];
            indices[i] = indices[j];
            indices[j] = tmp;
        }
        return Arrays.copyOf(indices, k);
    }

}
