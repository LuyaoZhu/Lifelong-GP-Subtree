package yimei.jss.algorithm.surrogateCaseLS;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.ArrayList;
import java.util.Collections;

public class SplineRegressor {

    private PolynomialSplineFunction spline;

    public SplineRegressor(ArrayList<Double> xList, ArrayList<Double> yList) {
        if (xList.size() != yList.size() || xList.size() < 3) {
            throw new IllegalArgumentException("xList and yList must match in size and have at least 3 points.");
        }

        // 排序 x 和 y
        ArrayList<Double> xSorted = new ArrayList<>(xList);
        ArrayList<Double> ySorted = new ArrayList<>(yList);
        sortXY(xSorted, ySorted);

        // 转为 double[]
        double[] x = toArray(xSorted);
        double[] y = toArray(ySorted);

        SplineInterpolator interpolator = new SplineInterpolator();
        this.spline = interpolator.interpolate(x, y);
    }

    public double predict(double x) {
        return spline.value(x);
    }

    public ArrayList<Double> predictAll(ArrayList<Double> testX) {
        ArrayList<Double> result = new ArrayList<>();
        for (double val : testX) {
            result.add(predict(val));
        }
        return result;
    }

    public double computeMSE(ArrayList<Double> xList, ArrayList<Double> yList) {
        double sum = 0;
        for (int i = 0; i < xList.size(); i++) {
            double diff = predict(xList.get(i)) - yList.get(i);
            sum += diff * diff;
        }
        return sum / xList.size();
    }

    private void sortXY(ArrayList<Double> x, ArrayList<Double> y) {
        int n = x.size();
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (x.get(i) > x.get(j)) {
                    Collections.swap(x, i, j);
                    Collections.swap(y, i, j);
                }
            }
        }
    }

    private double[] toArray(ArrayList<Double> list) {
        double[] arr = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            arr[i] = list.get(i);
        }
        return arr;
    }
}
