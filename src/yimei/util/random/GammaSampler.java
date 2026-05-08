package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

public class GammaSampler extends AbstractRealSampler {

	private double shape;
	private double scale;
	private double lower;
	private double upper;

	public GammaSampler() {
		super();
	}

	public GammaSampler(double shape, double scale, double lower, double upper) {
		super();
		this.shape = shape;
		this.scale = scale;
		this.lower = lower;
		this.upper = upper;
	}

	public void setShape(double shape) {
		this.shape = shape;
	}

	public double getShape() {
		return shape;
	}

	public void setScale(double scale) {
		this.scale = scale;
	}

	public double getScale() {
		return scale;
	}

	@Override
	public double next(RandomDataGenerator rdg) {
		double value = rdg.nextGamma(shape, scale);
		while (value > upper || value < lower) {
			value = rdg.nextGamma(shape, scale);
		}
		return value;
	}

	@Override
	public void setLower(double lower) {

	}

	@Override
	public void setUpper(double upper) {

	}

	public double getLower() {
		return  lower;
	}

	public double getUpper() {
		return  upper;
	}

	@Override
	public void setMean(double mean) {

	}

	@Override
	public double getMean() {
		return shape * scale;
	}

	@Override
	public AbstractRealSampler clone() {
		return new GammaSampler(shape, scale,upper,lower);
	}
}
