package yimei.jss.algorithm.downSampling;

import java.util.*;
import java.util.stream.IntStream;

/** K-Medoids (PAM) clustering on int[][] (Euclidean distance). */
public class KMedoidsPAM {

    public static class Result {
        public final int[] medoids;  // indices of medoids (each is an instance index)
        public final int[] labels;   // cluster assignment for each instance
        public Result(int[] medoids, int[] labels) { this.medoids = medoids; this.labels = labels; }
    }

    /** Main API: fit k-medoids on int[][] data. */
    public static Result fit(int[][] data, int k, int maxIter, long seed) {
        if (k <= 0 || k > data.length) throw new IllegalArgumentException("Invalid k");
        double[][] X = toDouble(data);

        // --- init medoids with k-medoids++ style (better than random) ---
        int[] medoids = initKMedoidsPlusPlus(X, k, seed);

        int n = X.length;
        int[] labels = new int[n];
        boolean changed = true;
        int iter = 0;

        while (changed && iter++ < maxIter) {
            // 1) assign
            changed = assign(X, medoids, labels);

            // 2) update medoids (for each cluster, pick the point minimizing total distance)
            for (int c = 0; c < k; c++) {
                int bestMedoid = medoids[c];
                double bestCost = clusterCost(X, labels, c, bestMedoid);

                // consider every point in cluster as candidate medoid
                for (int i = 0; i < n; i++) {
                    if (labels[i] != c) continue;
                    double cost = clusterCost(X, labels, c, i);
                    if (cost < bestCost - 1e-12) {
                        bestCost = cost;
                        bestMedoid = i;
                    }
                }
                medoids[c] = bestMedoid;
            }
        }
        // final assign to be safe
        assign(X, medoids, labels);
        return new Result(medoids, labels);
    }

    // ---------- helpers ----------

    private static boolean assign(double[][] X, int[] medoids, int[] labels) {
        boolean changed = false;
        for (int i = 0; i < X.length; i++) {
            int best = -1;
            double bestD = Double.POSITIVE_INFINITY;
            for (int m = 0; m < medoids.length; m++) {
                double d = dist2(X[i], X[medoids[m]]);
                if (d < bestD) { bestD = d; best = m; }
            }
            if (labels[i] != best) { labels[i] = best; changed = true; }
        }
        return changed;
    }

    private static double clusterCost(double[][] X, int[] labels, int c, int medoidIndex) {
        double sum = 0.0;
        for (int i = 0; i < X.length; i++) {
            if (labels[i] == c) sum += dist2(X[i], X[medoidIndex]);
        }
        return sum;
    }

    private static double dist2(double[] a, double[] b) { // Euclidean (no sqrt for speed)
        double s = 0.0;
        for (int j = 0; j < a.length; j++) {
            double d = a[j] - b[j];
            s += d * d;
        }
        return s;
    }

    private static double[][] toDouble(int[][] arr) {
        double[][] x = new double[arr.length][arr[0].length];
        for (int i = 0; i < arr.length; i++)
            for (int j = 0; j < arr[i].length; j++)
                x[i][j] = arr[i][j];
        return x;
    }

    /** k-medoids++ init: pick first random, then probability ∝ minDist^2 to current set. */
    private static int[] initKMedoidsPlusPlus(double[][] X, int k, long seed) {
        Random rnd = new Random(seed);
        int n = X.length;
        int[] medoids = new int[k];
        boolean[] chosen = new boolean[n];

        // first medoid
        medoids[0] = rnd.nextInt(n);
        chosen[medoids[0]] = true;

        double[] minD2 = new double[n];
        Arrays.fill(minD2, Double.POSITIVE_INFINITY);
        for (int i = 0; i < n; i++) minD2[i] = dist2(X[i], X[medoids[0]]);

        for (int m = 1; m < k; m++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++) if (!chosen[i]) sum += minD2[i];

            double r = rnd.nextDouble() * sum;
            int next = -1;
            double acc = 0.0;
            for (int i = 0; i < n; i++) {
                if (chosen[i]) continue;
                acc += minD2[i];
                if (acc >= r) { next = i; break; }
            }
            if (next < 0) { // fallback
                do { next = rnd.nextInt(n); } while (chosen[next]);
            }
            medoids[m] = next;
            chosen[next] = true;

            for (int i = 0; i < n; i++) {
                double d2 = dist2(X[i], X[next]);
                if (d2 < minD2[i]) minD2[i] = d2;
            }
        }
        return medoids;
    }

    // utility: cluster -> member indices
    public static List<List<Integer>> clusterMembers(int[] labels, int[] medoids, int k) {
        List<List<Integer>> groups = new ArrayList<>();
        for (int c = 0; c < k; c++) groups.add(new ArrayList<>());

        for (int i = 0; i < labels.length; i++) {
            groups.get(labels[i]).add(i);
        }

        // 把 medoid 放到对应簇的首位
        for (int c = 0; c < k; c++) {
            int med = medoids[c];
            List<Integer> cluster = groups.get(c);
            if (cluster.remove((Integer) med)) {
                cluster.add(0, med);
            }
        }
        return groups;
    }

    // quick demo
    public static void main(String[] args) {
        int[][] indsCharListsMultiTree = {
                {1,2,3},
                {2,2,4},
                {10,11,12},
                {11,12,13},
                {9,10,11},
                {50,50,50}
        };
        int k = 2;
        Result r = fit(indsCharListsMultiTree, k, 100, 2025L);
        System.out.println("Medoids (indices): " + Arrays.toString(r.medoids));
        List<List<Integer>> groups = clusterMembers(r.labels, r.medoids, k);
        for (int c = 0; c < k; c++) {
            System.out.println("Cluster " + c + " members " + groups.get(c));
            System.out.println("Cluster " + c + " medoid instance index = " + r.medoids[c]);
        }
    }
}
