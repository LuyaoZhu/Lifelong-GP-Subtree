package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

public class NormalSampler extends AbstractRealSampler {

	private double mean;
	private double sd;
	private double lower;
	private double upper;

	public NormalSampler() {
		super();
	}

	public NormalSampler(double mean, double sd) {
		super();
		this.mean = mean;
		this.sd = sd;
	}

	public NormalSampler(double mean, double sd, double lower, double upper) {
		super();
		this.mean = mean;
		this.sd = sd;
		this.lower = lower;
		this.upper = upper;
	}

	public void set(double mean, double sd) {
		this.mean = mean;
		this.sd = sd;
	}

	public double getLower() {
		return lower;
	}

	@Override
	public double getUpper() {
		return  upper;
	}

	public void setMean(double mean) {
		this.mean = mean;
	}

	public void setSd(double sd) {
		this.sd = sd;
	}

	public double getMean() {
		return mean;
	}

	public double getSd() {
		return sd;
	}

	@Override
//	public double next(RandomDataGenerator rdg) {
//		return rdg.nextGaussian(mean, sd);
//	}

	public double next(RandomDataGenerator rdg) {
		double value = rdg.nextGaussian(mean, sd);
		while (value > upper || value < lower) {
			value = rdg.nextGaussian(mean, sd);
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
		return new NormalSampler(mean, sd);
	}
}
