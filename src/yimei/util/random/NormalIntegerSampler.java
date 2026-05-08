package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

/**
 * Integer sampler from a (truncated) normal distribution.
 * Uses rejection sampling on half-integer bins so that rounding-to-int
 * stays unbiased near the bounds. Bounds are inclusive: [lower, upper].
 */
public class NormalIntegerSampler extends AbstractIntegerSampler {

    private int lower;
    private int upper;
    private double mean;
    private double sd;

    public NormalIntegerSampler() {
        super();
        this.lower = Integer.MIN_VALUE;
        this.upper = Integer.MAX_VALUE;
    }

    public NormalIntegerSampler(double mean, double sd) {
        super();
        this.mean = mean;
        this.sd = sd;
        this.lower = Integer.MIN_VALUE;
        this.upper = Integer.MAX_VALUE;
        validateParams();
    }

    public NormalIntegerSampler(double mean, double sd, int lower, int upper) {
        super();
        this.mean = mean;
        this.sd = sd;
        this.lower = lower;
        this.upper = upper;
        validateParams();
    }

    /** Set mean and standard deviation. */
    public void set(double mean, double sd) {
        this.mean = mean;
        this.sd = sd;
        validateParams();
    }

    /** Set both bounds (inclusive). */
    public void set(int lower, int upper) {
        this.lower = lower;
        this.upper = upper;
        validateParams();
    }

    public void setMean(double mean) { this.mean = mean; validateParams(); }
    public void setSd(double sd)     { this.sd = sd;     validateParams(); }

    public double getMean() { return mean; }
    public double getSd()   { return sd;   }
    public int getLower()   { return lower; }
    public int getUpper()   { return upper; }

    @Override
    public void setLower(int lower) {
        this.lower = lower;
        validateParams();
    }

    @Override
    public void setUpper(int upper) {
        this.upper = upper;
        validateParams();
    }

    private void validateParams() throws IllegalArgumentException {
        if (!(sd > 0.0)) {
            throw new IllegalArgumentException("sd must be > 0");
        }
        if (upper < lower) {
            throw new IllegalArgumentException("upper must be >= lower");
        }
    }

    @Override
    public int next(RandomDataGenerator rdg) {
        // Rejection on half-integer bins:
        // accept only if the real sample lies in [lower-0.5, upper+0.5),
        // then round to nearest integer -> guarantees result ∈ [lower, upper].
        final double lo = (double) lower - 0.5;
        final double hi = (double) upper + 0.5;

        int candidate;
        while (true) {
            double v = rdg.nextGaussian(mean, sd);
            if (v < lo || v >= hi) continue;   // outside truncation window
            candidate = (int) Math.round(v);   // safe: now within [lower, upper]
            break;
        }
        return candidate;
    }

    @Override
    public AbstractIntegerSampler clone() {
        return new NormalIntegerSampler(mean, sd, lower, upper);
    }
}
