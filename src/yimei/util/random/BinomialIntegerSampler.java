package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.distribution.BinomialDistribution;

/**
 * A binomial integer distribution sampler on [lower, upper] (endpoints included).
 * Internally: X = lower + Binomial(n = upper - lower, p).
 *
 * Examples:
 *   // Ops in {1..10} with mean 5.5: lower=1, upper=10 => n=9, p=0.5
 *   BinomialIntegerSampler s = new BinomialIntegerSampler(1, 10, 0.5);
 *
 * Notes:
 *   - Requires 0 <= p <= 1 and upper >= lower.
 *   - If n == 0, distribution is degenerate at 'lower'.
 *
 * @author yimei
 */
public class BinomialIntegerSampler extends AbstractIntegerSampler {

    private int lower;     // inclusive
    private int upper;     // inclusive
    private double p;      // success probability in Binomial(n, p), n = upper - lower

    public BinomialIntegerSampler() {
        super();
        this.lower = 0;
        this.upper = 0;
        this.p = 0.5;
    }

    public BinomialIntegerSampler(int lower, int upper, double p) {
        super();
        if (upper < lower) {
            throw new IllegalArgumentException("upper must be >= lower");
        }
        this.lower = lower;
        this.upper = upper;
        setP(p);
    }

    /** Convenience ctor: support {0..n}, i.e., lower=0, upper=n */
    public BinomialIntegerSampler(int n, double p) {
        this(0, n, p);
    }

    /** Set lower/upper, keep current p */
    public void set(int lower, int upper) {
        if (upper < lower) {
            throw new IllegalArgumentException("upper must be >= lower");
        }
        this.lower = lower;
        this.upper = upper;
    }

    /** Set lower/upper/p all at once */
    public void set(int lower, int upper, double p) {
        set(lower, upper);
        setP(p);
    }

    public void setLower(int lower) {
        // keep current upper, just shift the support
        this.lower = lower;
        if (this.upper < this.lower) {
            this.upper = this.lower; // collapse if needed
        }
    }

    public void setUpper(int upper) {
        if (upper < this.lower) {
            throw new IllegalArgumentException("upper must be >= lower");
        }
        this.upper = upper;
    }

    public void setP(double p) {
        if (p < 0.0) p = 0.0;
        if (p > 1.0) p = 1.0;
        this.p = p;
    }

    public int getLower() {
        return lower;
    }

    public int getUpper() {
        return upper;
    }

    public double getP() {
        return p;
    }

    private int n() {
        return upper - lower; // Binomial number of trials
    }

    @Override
    public int next(RandomDataGenerator rdg) {
        int n = n();
        if (n <= 0) {
            return lower; // degenerate at lower
        }
        // Use rdg's underlying RandomGenerator for consistent seeding with the framework
        BinomialDistribution bd = new BinomialDistribution(rdg.getRandomGenerator(), n, p);
        return lower + bd.sample();
    }

    @Override
    public double getMean() {
        return lower + n() * p;
    }

    @Override
    public void setMean(double mean) {
        // mean = lower + n * p  =>  p = (mean - lower)/n
        int n = n();
        if (n > 0) {
            double newP = (mean - lower) / (double) n;
            setP(newP); // clamps into [0,1]
        }
        // if n == 0, mean is fixed at 'lower'; do nothing.
    }

    @Override
    public AbstractIntegerSampler clone() {
        return new BinomialIntegerSampler(lower, upper, p);
    }

    @Override
    public String toString() {
        return String.format("BinomialIntegerSampler([%d,%d], p=%.6f)", lower, upper, p);
    }
}
