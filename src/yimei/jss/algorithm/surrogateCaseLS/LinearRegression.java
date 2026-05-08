package yimei.jss.algorithm.surrogateCaseLS;

import java.util.ArrayList;

public class LinearRegression {

    // 拟合函数：返回 [a, b]
    public static double[] fit(ArrayList<Double> x, ArrayList<Double> y) {
        if (x.size() != y.size() || x.isEmpty()) {
            throw new IllegalArgumentException("x and y must have the same non-zero size");
        }

        int n = x.size();
        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

        for (int i = 0; i < n; i++) {
            double xi = x.get(i);
            double yi = y.get(i);
            sumX += xi;
            sumY += yi;
            sumXY += xi * yi;
            sumX2 += xi * xi;
        }

        double xMean = sumX / n;
        double yMean = sumY / n;

        double a = (sumXY - n * xMean * yMean) / (sumX2 - n * xMean * xMean);
        double b = yMean - a * xMean;

        return new double[] { a, b };
    }

    // 预测函数：输入 a、b、x，输出预测 y
    public static double predict(double a, double b, double x) {
        return a * x + b;
    }


}
