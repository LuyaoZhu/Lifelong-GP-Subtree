package yimei.jss.algorithm.lifelongGP;

import smile.math.kernel.GaussianKernel;
import smile.regression.GaussianProcessRegression;
import smile.math.matrix.Matrix;
import smile.math.matrix.DenseMatrix;
import smile.math.matrix.Cholesky;

public class GPRSurrogateModel {

    private GaussianProcessRegression<double[]> model;
    private double[][] trainX;
    private GaussianKernel kernel;
    private Cholesky cholesky; // 我们需要自己存一份分解结果

    private final double lambda = 0.01;
    private final double THRESHOLD = 0.35;

    public void train(double[][] x, double[] y) {
        this.trainX = x;
        this.kernel = new GaussianKernel(1.0);

        // 1. 调用你源码中的 Trainer 进行训练
        GaussianProcessRegression.Trainer<double[]> trainer =
                new GaussianProcessRegression.Trainer<>(kernel, lambda);
        this.model = trainer.train(x, y);

        // 2. 为了计算 SD，必须手动复现源码中 138-146 行的矩阵逻辑
        int n = x.length;
        DenseMatrix K = Matrix.zeros(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                double kVal = kernel.k(x[i], x[j]);
                K.set(i, j, kVal);
                K.set(j, i, kVal);
            }
            K.add(i, i, lambda);
        }
        // 存下 Cholesky 分解结果用于 estimate
        this.cholesky = K.cholesky();
    }

    public double estimate(double[] x) {
        if (model == null) return Double.MAX_VALUE;

        // 预测均值 (使用 Smile 模型)
        double mu = model.predict(x);

        // 计算标准差 (手动实现)
        double sd = calculateSD(x);

        // 惩罚项逻辑
        double penalty = 0.0;
        if (sd > THRESHOLD) {
            penalty = 10.0 * (sd - THRESHOLD);
        }

        return mu + penalty;
    }

    private double calculateSD(double[] x) {
        int n = trainX.length;
        double[] kStar = new double[n];
        for (int i = 0; i < n; i++) {
            kStar[i] = kernel.k(x, trainX[i]);
        }

        // 解 L * v = kStar (对应 GP 回归中的方差计算)
        double[] v = kStar.clone();
        cholesky.solve(v); // 2.x 的 Cholesky.solve 直接解全系统

        double kSS = kernel.k(x, x);
        double dotProduct = 0;
        for (int i = 0; i < n; i++) {
            dotProduct += kStar[i] * v[i];
        }

        return Math.sqrt(Math.max(0, kSS - dotProduct));
    }
}
