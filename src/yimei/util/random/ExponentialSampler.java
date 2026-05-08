package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

public class ExponentialSampler extends AbstractRealSampler {

	private double mean;
	private double lower = -Double.MAX_VALUE;
	private double upper = Double.MAX_VALUE;

	public ExponentialSampler() {
		super();
	}

	public ExponentialSampler(double mean) {
		super();
		this.mean = mean;
	}

	public ExponentialSampler(double mean, double lower, double upper) {
		super();
		this.mean = mean;
		this.lower = lower;
		this.upper = upper;
	}

	public void setMean(double mean) {
		this.mean = mean;
	}

	public double getMean() {
		return mean;
	}

	@Override
	public double getLower() {
		return lower;
	}

	@Override
	public double getUpper() {
		return  upper;
	}

	public double next(RandomDataGenerator rdg) {
		double value = rdg.nextExponential(mean);
		while (value > upper || value < lower) {
			value = rdg.nextExponential(mean);
		}
		return value;
	}

	@Override
	public void setLower(double lower) {
		// do nothing.
	}

	@Override
	public void setUpper(double upper) {
		// do nothing.

	}

	@Override
	public AbstractRealSampler clone() {
		return new ExponentialSampler(mean);
	}
}
