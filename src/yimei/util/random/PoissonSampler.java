package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

public class PoissonSampler extends AbstractRealSampler {

    private double mean;

    public PoissonSampler() {
        super();
    }

    public PoissonSampler(double mean) {
        super();
        this.mean = mean;
    }

    public void setMean(double mean) {
        this.mean = mean;
    }

    public double getMean() {
        return mean;
    }

    public double next(RandomDataGenerator rdg) {
        return rdg.nextPoisson(mean);
    }

    @Override
    public void setLower(double lower) {
        // do nothing.
    }

    @Override
    public void setUpper(double upper) {
        // do nothing.

    }

    public double getLower() {
        return 0;
    }

    @Override
    public double getUpper() {
        return  0;
    }

    @Override
    public AbstractRealSampler clone() {
        return new PoissonSampler(mean);
    }
}

