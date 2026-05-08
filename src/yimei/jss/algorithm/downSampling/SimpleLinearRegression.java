package yimei.jss.algorithm.downSampling;

import java.util.ArrayList;
import java.util.Arrays;

public class SimpleLinearRegression {

    public static class Model {
        public final double intercept;     // 截距
        public final double[] coef;        // 系数（不含截距，对应原始列）
        public Model(double intercept, double[] coef) {
            this.intercept = intercept; this.coef = coef;
        }
        public double predict(double[] x) {
            double s = intercept;
            for (int j = 0; j < coef.length; j++) s += coef[j] * x[j];
            return s;
        }

        public double predict(ArrayList<Double> x) {
            if (x.size() != coef.length) {
                throw new IllegalArgumentException("输入长度 " + x.size() + " 与系数长度 " + coef.length + " 不一致");
            }
            double s = intercept;
            for (int j = 0; j < coef.length; j++) {
                s += coef[j] * x.get(j);
            }
            return s;
        }

    }

    /** 用普通最小二乘拟合（带截距），对奇异/病态情形加极小ridge稳定项 */
    public static Model fit(ArrayList<double[]> Xlist, ArrayList<Double> ylist) {
        int n = Xlist.size();
        if (n == 0) throw new IllegalArgumentException("Empty X");
        int p = Xlist.get(0).length;

        // 1) 计算 XtX 与 Xty （带截距 => 维度 (p+1)）
        int d = p + 1;
        double[][] XtX = new double[d][d];
        double[] Xty = new double[d];

        for (int i = 0; i < n; i++) {
            double[] x = Xlist.get(i);
            double yi = ylist.get(i);
            // 扩展特征：z = [1, x0, x1, ..., xp-1]
            double[] z = new double[d];
            z[0] = 1.0;
            System.arraycopy(x, 0, z, 1, p);

            // 累加 XtX 与 Xty
            for (int a = 0; a < d; a++) {
                Xty[a] += z[a] * yi;
                double za = z[a];
                for (int b = 0; b < d; b++) {
                    XtX[a][b] += za * z[b];
                }
            }
        }

        // 2) 数值稳定：加入极小 ridge（不影响结果但防止奇异）
        double ridge = 1e-8;
        for (int k = 0; k < d; k++) XtX[k][k] += ridge;

        // 3) 解线性方程 XtX * beta = Xty  （用高斯消元 + 简单部分选主元）
        double[] beta = solveSymmetric(XtX, Xty); // 长度 d

        // 4) 封装：beta[0] 是截距，后面是各列系数
        double intercept = beta[0];
        double[] coef = Arrays.copyOfRange(beta, 1, d);
        return new Model(intercept, coef);
    }

    /** 简单的线性方程求解器（带部分选主元的高斯消元），适合小矩阵 */
    private static double[] solveSymmetric(double[][] A, double[] b) {
        int n = b.length;
        double[][] M = new double[n][n];
        for (int i = 0; i < n; i++) M[i] = Arrays.copyOf(A[i], n);
        double[] y = Arrays.copyOf(b, n);

        // 消元
        for (int k = 0; k < n; k++) {
            // 选主元
            int piv = k;
            double max = Math.abs(M[k][k]);
            for (int i = k + 1; i < n; i++) {
                double v = Math.abs(M[i][k]);
                if (v > max) { max = v; piv = i; }
            }
            if (max == 0.0) throw new RuntimeException("Singular matrix.");
            if (piv != k) {
                double[] tmpR = M[k]; M[k] = M[piv]; M[piv] = tmpR;
                double tmpV = y[k]; y[k] = y[piv]; y[piv] = tmpV;
            }

            // 归一化主元行
            double diag = M[k][k];
            for (int j = k; j < n; j++) M[k][j] /= diag;
            y[k] /= diag;

            // 消掉其他行的第k列
            for (int i = 0; i < n; i++) {
                if (i == k) continue;
                double factor = M[i][k];
                if (factor == 0) continue;
                for (int j = k; j < n; j++) M[i][j] -= factor * M[k][j];
                y[i] -= factor * y[k];
            }
        }
        return y; // 现在 y 就是解
    }




    // 简单演示
    public static void main(String[] args) {
        ArrayList<double[]> X = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        X.add(new double[]{1, 2, 3}); y.add(10.0);
        X.add(new double[]{2, 1, 0}); y.add( 8.0);
        X.add(new double[]{3, 3, 1}); y.add(15.0);

        Model m = fit(X, y);
        System.out.println("Intercept: " + m.intercept);
        System.out.println("Coefficients: " + Arrays.toString(m.coef));
        System.out.println("Pred[2,2,2] = " + m.predict(new double[]{2,2,2}));
    }
}