package yimei.jss.algorithm.surrogateCaseLS;

import java.util.ArrayList;

public class KernelRegressor {

    private ArrayList<Double> xList;
    private ArrayList<Double> yList;
    private double h; // 核带宽（smoothing bandwidth）

    public KernelRegressor(ArrayList<Double> xList, ArrayList<Double> yList, double h) {
        if (xList.size() != yList.size() || xList.isEmpty()) {
            throw new IllegalArgumentException("xList and yList must be same size and not empty.");
        }
        this.xList = xList;
        this.yList = yList;
        this.h = h;
    }

    // Gaussian 核函数
    private double gaussianKernel(double u) {
        return Math.exp(-0.5 * u * u) / Math.sqrt(2 * Math.PI);
    }

    // 预测单个 x 值对应的 y 值
    public double predict(double x) {
        double numerator = 0.0;
        double denominator = 0.0;

        for (int i = 0; i < xList.size(); i++) {
            double xi = xList.get(i);
            double yi = yList.get(i);
            double u = (x - xi) / h;
            double k = gaussianKernel(u);
            numerator += k * yi;
            denominator += k;
        }

        if (denominator == 0) return 0.0; // 避免除 0

        return numerator / denominator;
    }

    // 批量预测
    public ArrayList<Double> predictAll(ArrayList<Double> inputList) {
        ArrayList<Double> results = new ArrayList<>();
        for (double x : inputList) {
            results.add(predict(x));
        }
        return results;
    }

    // 均方误差
    public double computeMSE(ArrayList<Double> xTest, ArrayList<Double> yTest) {
        double sum = 0.0;
        for (int i = 0; i < xTest.size(); i++) {
            double pred = predict(xTest.get(i));
            double diff = pred - yTest.get(i);
            sum += diff * diff;
        }
        return sum / xTest.size();
    }
}
