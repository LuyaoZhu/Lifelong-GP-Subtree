package yimei.util.random;

import org.apache.commons.math3.random.RandomDataGenerator;

public class LogNormalSampler extends AbstractRealSampler {

    // 对数正态的底层参数：log(X) ~ N(logMean, logSd^2)
    private double logMean;
    private double logSd;

    // 截断区间（与参考 NormalSampler 风格一致）
    private double lower = -Double.MAX_VALUE;
    private double upper = Double.MAX_VALUE;

    public LogNormalSampler() {
        super();
    }

    public LogNormalSampler(double logMean, double logSd) {
        super();
        this.logMean = logMean;
        this.logSd = logSd;
    }

    public LogNormalSampler(double logMean, double logSd, double lower, double upper) {
        super();
        this.logMean = logMean;
        this.logSd = logSd;
        this.lower = lower;
        this.upper = upper;
    }

    /** 直接设置底层对数空间参数 */
    public void set(double logMean, double logSd) {
        this.logMean = logMean;
        this.logSd = logSd;
    }

    public double getLower() {
        return lower;
    }

    @Override
    public double getUpper() {
        return  upper;
    }

    public void setLogMean(double logMean) {
        this.logMean = logMean;
    }

    public void setLogSd(double logSd) {
        this.logSd = logSd;
    }

    public double getLogMean() {
        return logMean;
    }

    public double getLogSd() {
        return logSd;
    }

    /** AbstractRealSampler 要求实现：按算术均值设置（需要已有 logSd） */
    public void setMean(double mean) {
        if (mean <= 0) {
            throw new IllegalArgumentException("LogNormal arithmetic mean must be > 0");
        }
        // E[X] = exp(logMean + 0.5*logSd^2)  =>  logMean = ln(mean) - 0.5*logSd^2
        this.logMean = Math.log(mean) - 0.5 * (logSd * logSd);
    }

    /** 如果你的基类也会调用 setSd(double)，这里将其解释为“底层对数标准差” */
    public void setSd(double sd) {
        if (sd <= 0) {
            throw new IllegalArgumentException("LogNormal logSd must be > 0");
        }
        this.logSd = sd;
    }

    /** 可选：返回算术均值（便于调试/一致性） */
    public double getMean() {
        return Math.exp(logMean + 0.5 * logSd * logSd);
    }

    public double getSd() {
        double mu = logMean, s2 = logSd * logSd;
        // Var[X] = (e^{s^2}-1) e^{2mu + s^2}
        double var = (Math.expm1(s2)) * Math.exp(2 * mu + s2);
        return Math.sqrt(var);
    }

    @Override
    public double next(RandomDataGenerator rdg) {
        double value = Math.exp(rdg.nextGaussian(logMean, logSd));
        while (value > upper || value < lower) {
            value = Math.exp(rdg.nextGaussian(logMean, logSd));
        }
        return value;
    }

    @Override
    public void setLower(double lower) {
        // do nothing. （保持与你的 NormalSampler 一致）
    }

    @Override
    public void setUpper(double upper) {
        // do nothing. （保持与你的 NormalSampler 一致）
    }

    @Override
    public AbstractRealSampler clone() {
        // 与参考一致：不复制截断边界
        return new LogNormalSampler(logMean, logSd);
    }
}
