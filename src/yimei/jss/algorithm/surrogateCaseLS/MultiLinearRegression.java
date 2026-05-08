package yimei.jss.algorithm.surrogateCaseLS;

import org.apache.commons.math3.linear.*;

public class MultiLinearRegression {

    // 训练模型：输入 int[][] X, double[] y；输出 beta（含截距）
    public static RealVector fit(int[][] X, double[] y) {
        int n = X.length;
        int dim = X[0].length;

        // 构建 X 矩阵并加一列 1（用于截距）
        double[][] Xd = new double[n][dim + 1];
        for (int i = 0; i < n; i++) {
            Xd[i][0] = 1.0; // intercept
            for (int j = 0; j < dim; j++) {
                Xd[i][j + 1] = X[i][j];
            }
        }

        RealMatrix XMatrix = MatrixUtils.createRealMatrix(Xd);
        RealVector yVector = MatrixUtils.createRealVector(y);

        // 正规方程：beta = (X^T X)^-1 X^T y
        RealMatrix Xt = XMatrix.transpose();
        RealMatrix XtX = Xt.multiply(XMatrix);
        RealMatrix XtXInv = new LUDecomposition(XtX).getSolver().getInverse();
        RealVector beta = XtXInv.multiply(Xt).operate(yVector);

        return beta;
    }

    // 预测新样本：输入 beta 和一个 int[] x（不含截距），输出预测值
    public static double predict(RealVector beta, int[] x) {
        double y = beta.getEntry(0); // 截距项
        for (int i = 0; i < x.length; i++) {
            y += beta.getEntry(i + 1) * x[i];
        }
        return y;
    }
}

