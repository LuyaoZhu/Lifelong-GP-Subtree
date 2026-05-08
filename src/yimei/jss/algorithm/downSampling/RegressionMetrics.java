package yimei.jss.algorithm.downSampling;

import java.util.ArrayList;
import java.util.Arrays;

public final class RegressionMetrics {

    private RegressionMetrics() {}

    /** 一键评估并打印：Spearman、Pearson、MSE、RMSE、MAE、R² */
    public static void report(SimpleLinearRegression.Model model,
                              ArrayList<double[]> X, ArrayList<Double> yList) {
        double[] y = yList.stream().mapToDouble(Double::doubleValue).toArray();
        double[] yhat = predict(model, X);

        double rho = spearman(yhat, y);
        double pr  = pearson(yhat, y);
        double mse = mse(yhat, y);
        double rmse = Math.sqrt(mse);
        double mae = mae(yhat, y);
        double r2  = r2(yhat, y);

        System.out.printf("Spearman ρ = %.6f%n", rho);
        System.out.printf("Pearson   r = %.6f%n", pr);
        System.out.printf("MSE         = %.6f%n", mse);
        System.out.printf("RMSE        = %.6f%n", rmse);
        System.out.printf("MAE         = %.6f%n", mae);
        System.out.printf("R²          = %.6f%n", r2);
    }

    /** 批量预测 */
    public static double[] predict(SimpleLinearRegression.Model model,
                                   ArrayList<double[]> X) {
        int n = X.size();
        double[] yhat = new double[n];
        for (int i = 0; i < n; i++) {
            yhat[i] = model.predict(X.get(i));
        }
        return yhat;
    }

    // ---------- 相关性/排名 ----------
    /** Spearman 等级相关（ties 取平均名次） */
    public static double spearman(double[] a, double[] b) {
        return pearson(rank(a), rank(b));
    }

    /** 皮尔逊相关 */
    public static double pearson(double[] x, double[] y) {
        int n = x.length;
        double mx = 0, my = 0;
        for (double v : x) mx += v; mx /= n;
        for (double v : y) my += v; my /= n;
        double num = 0, dx = 0, dy = 0;
        for (int i = 0; i < n; i++) {
            double ux = x[i] - mx, uy = y[i] - my;
            num += ux * uy; dx += ux * ux; dy += uy * uy;
        }
        double den = Math.sqrt(dx) * Math.sqrt(dy);
        return den > 0 ? num / den : 0.0;
    }

    /** 排名（1-based），相同值给平均名次 */
    public static double[] rank(double[] v) {
        int n = v.length;
        Integer[] idx = new Integer[n];
        for (int i = 0; i < n; i++) idx[i] = i;
        Arrays.sort(idx, (i, j) -> Double.compare(v[i], v[j]));
        double[] r = new double[n];
        for (int i = 0; i < n; ) {
            int j = i;
            double vi = v[idx[i]];
            while (j + 1 < n && v[idx[j + 1]] == vi) j++;
            double avg = (i + j + 2) / 2.0; // 1-based
            for (int k = i; k <= j; k++) r[idx[k]] = avg;
            i = j + 1;
        }
        return r;
    }

    // ---------- 误差指标 ----------
    public static double mse(double[] yhat, double[] y) {
        int n = y.length; double s = 0;
        for (int i = 0; i < n; i++) { double d = yhat[i] - y[i]; s += d * d; }
        return s / n;
    }

    public static double mae(double[] yhat, double[] y) {
        int n = y.length; double s = 0;
        for (int i = 0; i < n; i++) s += Math.abs(yhat[i] - y[i]);
        return s / n;
    }

    /** 决定系数 R² = 1 - SSE/SST */
    public static double r2(double[] yhat, double[] y) {
        int n = y.length;
        double my = 0; for (double v : y) my += v; my /= n;
        double sse = 0, sst = 0;
        for (int i = 0; i < n; i++) {
            double e = yhat[i] - y[i];
            sse += e * e;
            double d = y[i] - my;
            sst += d * d;
        }
        return sst > 0 ? 1.0 - sse / sst : 0.0;
    }
}
