package yimei.jss.algorithm.lifelongGP;

import java.util.*;

public class SimpleKMedoids {

    public static class Result {
        public final int[] medoids; // length k, each is index in [0, n)
        public final int[] labels;  // length n, cluster id in [0, k)
        public Result(int[] medoids, int[] labels) {
            this.medoids = medoids;
            this.labels = labels;
        }
    }

    /** Squared Euclidean distance (faster, monotonic to Euclidean). */
    private static double dist2(double[] a, double[] b) {
        double s = 0.0;
        for (int i = 0; i < a.length; i++) {
            double d = a[i] - b[i];
            s += d * d;
        }
        return s;
    }

    /** k-medoids with reproducible initialization + iterative refine. */
    public static Result fit(double[][] X, int k, long seed, int maxIters) {
        int n = X.length;
        if (n == 0) throw new IllegalArgumentException("X is empty");
        if (k <= 0 || k > n) throw new IllegalArgumentException("k must be in [1, n]");

        Random rnd = new Random(seed);

        // 1) Initialize medoids (reproducible): sample k unique indices
        int[] medoids = initUniqueIndices(n, k, rnd);

        int[] labels = new int[n];
        Arrays.fill(labels, -1);

        // 2) Iterate: assign -> update medoids
        boolean changed = true;
        for (int iter = 0; iter < maxIters && changed; iter++) {
            changed = false;

            // 2.1 Assign each point to nearest medoid
            for (int i = 0; i < n; i++) {
                int bestC = -1;
                double bestD = Double.POSITIVE_INFINITY;
                for (int c = 0; c < k; c++) {
                    double d = dist2(X[i], X[medoids[c]]);
                    if (d < bestD) {
                        bestD = d;
                        bestC = c;
                    }
                }
                labels[i] = bestC;
            }

            // 2.2 Build cluster members list
            List<int[]> clusters = buildClusters(labels, k);

            // 2.3 Update medoid of each cluster: choose point minimizing sum distances within cluster
            for (int c = 0; c < k; c++) {
                int[] members = clusters.get(c);

                // Edge case: empty cluster -> re-seed medoid to a random point (reproducible)
                if (members.length == 0) {
                    int newMed = rnd.nextInt(n);
                    if (medoids[c] != newMed) {
                        medoids[c] = newMed;
                        changed = true;
                    }
                    continue;
                }

                int currentMed = medoids[c];
                int bestMed = currentMed;
                double bestScore = Double.POSITIVE_INFINITY;

                // Evaluate each member as candidate medoid
                for (int cand : members) {
                    double sum = 0.0;
                    // sum distance to all others in this cluster
                    for (int other : members) {
                        sum += dist2(X[cand], X[other]);
                        // small pruning
                        if (sum >= bestScore) break;
                    }
                    if (sum < bestScore) {
                        bestScore = sum;
                        bestMed = cand;
                    }
                }

                if (bestMed != currentMed) {
                    medoids[c] = bestMed;
                    changed = true;
                }
            }

            // Optional: you can print progress
            // System.out.println("iter=" + iter + " changed=" + changed);
        }

        // 3) Final assignment (ensure labels match final medoids)
        for (int i = 0; i < n; i++) {
            int bestC = -1;
            double bestD = Double.POSITIVE_INFINITY;
            for (int c = 0; c < k; c++) {
                double d = dist2(X[i], X[medoids[c]]);
                if (d < bestD) {
                    bestD = d;
                    bestC = c;
                }
            }
            labels[i] = bestC;
        }

        return new Result(medoids, labels);
    }

    private static int[] initUniqueIndices(int n, int k, Random rnd) {
        // Fisher-Yates on first k positions (memory-friendly)
        int[] idx = new int[n];
        for (int i = 0; i < n; i++) idx[i] = i;
        for (int i = 0; i < k; i++) {
            int j = i + rnd.nextInt(n - i);
            int tmp = idx[i]; idx[i] = idx[j]; idx[j] = tmp;
        }
        return Arrays.copyOf(idx, k);
    }

    private static List<int[]> buildClusters(int[] labels, int k) {
        int n = labels.length;
        int[] counts = new int[k];
        for (int i = 0; i < n; i++) counts[labels[i]]++;

        int[][] buckets = new int[k][];
        for (int c = 0; c < k; c++) buckets[c] = new int[counts[c]];

        int[] cursor = new int[k];
        for (int i = 0; i < n; i++) {
            int c = labels[i];
            buckets[c][cursor[c]++] = i;
        }

        List<int[]> clusters = new ArrayList<>(k);
        clusters.addAll(Arrays.asList(buckets));
        return clusters;
    }

    // Example usage
    public static void main(String[] args) {
        double[][] X = new double[5000][40];
        // TODO fill X

        int k = 500;
        long seed = 20250101L;
        int maxIters = 15;

        Result r = fit(X, k, seed, maxIters);
        System.out.println("medoids = " + r.medoids.length);
        System.out.println("labels  = " + r.labels.length);
    }
}

