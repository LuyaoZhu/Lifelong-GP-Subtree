package yimei.jss.algorithm.featureweighted;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.regression.SimpleRegression;

public class PearsonsCorrelationV1 extends PearsonsCorrelation {

    public double correlation(double[] xArray, double[] yArray) {
        SimpleRegression regression = new SimpleRegression();
        if (xArray.length != yArray.length) {
            throw new DimensionMismatchException(xArray.length, yArray.length);
        } else if (xArray.length < 2) {
            throw new MathIllegalArgumentException(LocalizedFormats.INSUFFICIENT_DIMENSION, new Object[]{xArray.length, 2});
        } else {
            for(int i = 0; i < xArray.length; ++i) {
                regression.addData(xArray[i], yArray[i]);
            }

            Double result = regression.getR();
            if(result.isNaN()){
                return 0;
            }
            else{
                return regression.getR();
            }
        }
    }
}
