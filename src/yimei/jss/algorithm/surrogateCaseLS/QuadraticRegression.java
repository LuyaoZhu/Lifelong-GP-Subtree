package yimei.jss.algorithm.surrogateCaseLS;

import org.apache.commons.math3.linear.*;

import java.util.ArrayList;
import java.util.Arrays;

public class QuadraticRegression {

    // 拟合函数，返回系数 [a, b, c]
    public static double[] fit(ArrayList<Double> xList, ArrayList<Double> yList) {
        if (xList.size() != yList.size() || xList.isEmpty()) {
            throw new IllegalArgumentException("xList and yList must be same size and not empty");
        }

        int n = xList.size();
        double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
        double sumY = 0, sumXY = 0, sumX2Y = 0;

        for (int i = 0; i < n; i++) {
            double x = xList.get(i);
            double y = yList.get(i);
            double x2 = x * x;

            sumX += x;
            sumX2 += x2;
            sumX3 += x2 * x;
            sumX4 += x2 * x2;

            sumY += y;
            sumXY += x * y;
            sumX2Y += x2 * y;
        }

        // 构建系数矩阵 A 和右边常数向量 B
        double[][] A = {
                {sumX4, sumX3, sumX2},
                {sumX3, sumX2, sumX},
                {sumX2, sumX,  n}
        };

        double[] B = {sumX2Y, sumXY, sumY};

        // 解线性方程组 A * coef = B
        RealMatrix matrix = MatrixUtils.createRealMatrix(A);
        DecompositionSolver solver = new LUDecomposition(matrix).getSolver();
        RealVector constants = MatrixUtils.createRealVector(B);
        RealVector solution = solver.solve(constants);

        return solution.toArray(); // 返回 [a, b, c]
    }

    // 预测函数：输入 x 和系数，输出 y = ax^2 + bx + c
    public static double predict(double[] coef, double x) {
        return coef[0]*x*x + coef[1]*x + coef[2];
    }

    // 示例测试
    public static void main(String[] args) {
        ArrayList<Double> x = new ArrayList<>(Arrays.asList(1.0, 2.0, 3.0, 4.0, 5.0));
        ArrayList<Double> y = new ArrayList<>(Arrays.asList(2.0, 5.0, 10.0, 17.0, 26.0)); // y ≈ x^2 + 1

        double[] coef = fit(x, y);
        System.out.printf("拟合公式: y = %.3f * x^2 + %.3f * x + %.3f\n", coef[0], coef[1], coef[2]);

        double testX = 6.0;
        double predictedY = predict(coef, testX);
        System.out.printf("预测 y(%.1f) = %.3f\n", testX, predictedY);
    }
}

