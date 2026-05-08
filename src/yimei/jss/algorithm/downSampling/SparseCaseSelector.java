package yimei.jss.algorithm.downSampling;

import smile.regression.LASSO;
import smile.regression.ElasticNet;

import java.util.*;

/**
 * SparseCaseSelector (K-sparse)
 * 目标：从 n 个特征里强制只选 K 个，其他系数精确为 0；以 Spearman ρ 最大化挑参数。
 * 实现：
 *  - 预筛去常数列
 *  - 标准化（均值-方差）
 *  - 在 λ 网格上用 LASSO/ENet 得到重要性
 *  - 将解“投影”为恰好 K 非零（>K 取Top-K；<K 用|corr|补足）
 *  - 在已选 K 特征上小Ridge重拟合（稳定），反标准化到原尺度
 *  - 以 LOOCV 的 Spearman ρ 选择最佳 λ
 *
 * 依赖：SMILE 4.3.0（LASSO/ElasticNet 构造器）
 */
public class SparseCaseSelector {

    private final ArrayList<double[]> Xlist; // 每行=一个体，列=case特征
    private final ArrayList<Double> ylist;   // 真实 y

    public SparseCaseSelector(ArrayList<double[]> Xlist, ArrayList<Double> ylist) {
        if (Xlist == null || ylist == null) throw new IllegalArgumentException("Xlist/ylist is null");
        if (Xlist.size() != ylist.size()) throw new IllegalArgumentException("X and y size mismatch");
        if (Xlist.isEmpty()) throw new IllegalArgumentException("Empty X");
        int p0 = Xlist.get(0).length;
        for (int i = 1; i < Xlist.size(); i++) {
            if (Xlist.get(i).length != p0) throw new IllegalArgumentException("Row " + i + " has length " + Xlist.get(i).length + ", expected " + p0);
        }
        this.Xlist = Xlist;
        this.ylist = ylist;
    }

    /** 公开：从 n 选 K，useElasticNet=true 用 ENet(α=alphaForENet)，否则 LASSO */
    public Result fitK(boolean useElasticNet, double alphaForENet, int K) {
        return fitK(useElasticNet, alphaForENet, K, /*useLOOCV=*/true);
    }

    /** 主入口：是否用LOOCV；如果 false 则用训练集 Spearman 直接挑 λ */
    public Result fitK(boolean useElasticNet, double alphaForENet, int K, boolean useLOOCV) {
        final int n = Xlist.size();
        final int p = Xlist.get(0).length;

        // —— 列表->数组
        final double[][] X = new double[n][p];
        final double[] y = new double[n];
        for (int i = 0; i < n; i++) { X[i] = Xlist.get(i); y[i] = ylist.get(i); }

        // —— 预筛常数/全NaN 列
        List<Integer> keptCols = new ArrayList<>();
        for (int j = 0; j < p; j++) {
            double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
            boolean allNaN = true;
            for (int i = 0; i < n; i++) {
                double v = X[i][j];
                if (!Double.isNaN(v)) {
                    allNaN = false;
                    if (v < min) min = v;
                    if (v > max) max = v;
                }
            }
            if (!allNaN && max > min) keptCols.add(j);
        }
        if (keptCols.isEmpty()) throw new IllegalStateException("All columns are constant/NaN.");
        final int pk = keptCols.size();

        final double[][] Xk = new double[n][pk];
        for (int i = 0; i < n; i++) {
            int c = 0;
            for (int j : keptCols) Xk[i][c++] = X[i][j];
        }

        // —— 标准化
        final double[] mean = colMean(Xk);
        final double[] sd   = colStd (Xk, mean);
        final double[][] Xs = standardize(Xk, mean, sd);

        // —— λ 网格
        final double[] lambdas = logspace(-2, 4, 60); // 1e-2..1e4

        // —— 预备：特征与 y 的 |corr| 排序（当 <K 时用于补足）
        final double[] absCorr = new double[pk];
        for (int j = 0; j < pk; j++) absCorr[j] = Math.abs(pearson(col(Xs, j), y));
        final Integer[] byAbsCorr = indexOrder(absCorr, /*desc=*/true);

        double bestRho = -2.0;
        double bestLambda = lambdas[0];
        double[] bestCoefStd = null;   // 标准化空间 K 稀疏系数（长度=pk，未选=0）
        double bestInterceptStd = 0.0;

        for (double lambda : lambdas) {
            // 1) 初始解：在标准化空间训练
            double[] wStd;
            if (useElasticNet) {
                ElasticNet en = new ElasticNet(Xs, y, lambda, alphaForENet);
                wStd = en.coefficients();
            } else {
                LASSO lasso = new LASSO(Xs, y, lambda);
                wStd = lasso.coefficients();
            }

            // 2) 选出恰好 K 个索引
            int[] picked = pickKIndicesByMagnitudeOrCorr(wStd, absCorr, K);

            // 3) 在 picked 上做一次小 ridge 重拟合（标准化空间），得到 K 稀疏系数
            double[] wKStd = refitRidgeOnSubset(Xs, y, picked, /*ridge=*/1e-6);

            // 4) LOOCV 或训练集 Spearman 评估（在标准化空间评估就行）
            double rho;
            if (useLOOCV && n >= 3) {
                // LOOCV：每次留一，子集固定为 picked，不重复选特征（更快更稳）
                double[] yhat = new double[n];
                for (int i = 0; i < n; i++) {
                    // 构造留一训练
                    int m = n - 1;
                    double[][] Xtr = new double[m][picked.length];
                    double[] ytr = new double[m];
                    for (int r = 0, t = 0; r < n; r++) {
                        if (r == i) continue;
                        for (int c = 0; c < picked.length; c++) Xtr[t][c] = Xs[r][picked[c]];
                        ytr[t] = y[r];
                        t++;
                    }
                    double[] wLOO = ridgeSolve(Xtr, ytr, 1e-6); // 仅系数，不含截距，因为我们构建了无截距形式
                    // 有截距的做法：我们这里 ridgeSolve 是带截距的形式（见实现），会返回完整解
                    double b0 = wLOO[0];
                    double[] b = Arrays.copyOfRange(wLOO, 1, wLOO.length);
                    // 预测
                    double s = b0;
                    for (int c = 0; c < picked.length; c++) s += b[c] * Xs[i][picked[c]];
                    yhat[i] = s;
                }
                rho = spearman(yhat, y);
            } else {
                // 训练集估计
                double[] yhat = new double[n];
                // 需要有截距：用全量再拟一次截距+系数
                double[] betaFull = ridgeSolve(buildSubset(Xs, picked), y, 1e-6);
                double b0 = betaFull[0];
                double[] b = Arrays.copyOfRange(betaFull, 1, betaFull.length);
                for (int i = 0; i < n; i++) {
                    double s = b0;
                    for (int c = 0; c < picked.length; c++) s += b[c] * Xs[i][picked[c]];
                    yhat[i] = s;
                }
                rho = spearman(yhat, y);
            }

            if (rho > bestRho) {
                bestRho = rho;
                bestLambda = lambda;
                // 保存标准化空间下的 K 稀疏解（对齐 pk 维，未选=0）
                bestCoefStd = new double[pk];
                for (int c = 0; c < picked.length; c++) bestCoefStd[picked[c]] = wKStd[c];
                // 记录对应的标准化空间截距（用全量再算，保证一致）
                double[] betaFull = ridgeSolve(buildSubset(Xs, picked), y, 1e-6);
                bestInterceptStd = betaFull[0];
            }
        }

        // —— 反标准化到原尺度（仅 keptCols，对应 pk 维；再补齐到 p 维）
        // y = b0_std + Σ b_std_j * (xj - meanj)/sdj
        //   = (b0_std - Σ b_std_j * meanj/sdj) + Σ (b_std_j / sdj) * xj
        final double[] coefKept = new double[pk];
        for (int j = 0; j < pk; j++) coefKept[j] = bestCoefStd[j] / sd[j];
        double interceptOrig = bestInterceptStd;
        for (int j = 0; j < pk; j++) interceptOrig -= bestCoefStd[j] * mean[j] / sd[j];

        final double[] coefFull = new double[p];
        Arrays.fill(coefFull, 0.0);
        for (int k = 0; k < pk; k++) coefFull[ keptCols.get(k) ] = coefKept[k];

        // —— 统计真正非零的原列索引
        List<Integer> selectedCols = new ArrayList<>();
        for (int j = 0; j < p; j++) if (Math.abs(coefFull[j]) > 0) selectedCols.add(j);

        Result out = new Result();
        out.bestLambda = bestLambda;
        out.bestSpearman = bestRho;
        out.intercept = interceptOrig;
        out.coef = coefFull;
        out.keptColumns = keptCols;
        out.selectedColumns = selectedCols;
        out.k = K;
        return out;
    }

    // ===================== 工具/内部实现 =====================

    private static double[] colMean(double[][] X) {
        int n = X.length, p = X[0].length;
        double[] m = new double[p];
        for (int j = 0; j < p; j++) {
            double s = 0; for (int i = 0; i < n; i++) s += X[i][j];
            m[j] = s / n;
        }
        return m;
    }
    private static double[] colStd(double[][] X, double[] mean) {
        int n = X.length, p = X[0].length;
        double[] sd = new double[p];
        for (int j = 0; j < p; j++) {
            double s2 = 0; for (int i = 0; i < n; i++) { double d = X[i][j]-mean[j]; s2 += d*d; }
            sd[j] = Math.sqrt(s2 / n);
            if (sd[j] == 0) sd[j] = 1.0;
        }
        return sd;
    }
    private static double[][] standardize(double[][] X, double[] mean, double[] sd) {
        int n = X.length, p = X[0].length;
        double[][] Z = new double[n][p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
                Z[i][j] = (X[i][j] - mean[j]) / sd[j];
        return Z;
    }
    private static double[] col(double[][] X, int j) {
        double[] out = new double[X.length];
        for (int i = 0; i < X.length; i++) out[i] = X[i][j];
        return out;
    }

    private static double pearson(double[] x, double[] y) {
        int n = x.length;
        double mx=0, my=0;
        for (double v : x) mx += v; mx /= n;
        for (double v : y) my += v; my /= n;
        double num=0, dx=0, dy=0;
        for (int i=0;i<n;i++){ double ux=x[i]-mx, uy=y[i]-my; num += ux*uy; dx += ux*ux; dy += uy*uy; }
        double den = Math.sqrt(dx) * Math.sqrt(dy);
        return den>0 ? num/den : 0.0;
    }
    private static double[] rank(double[] v) {
        int n=v.length; Integer[] idx=new Integer[n]; for(int i=0;i<n;i++) idx[i]=i;
        Arrays.sort(idx, java.util.Comparator.comparingDouble(i->v[i]));
        double[] r=new double[n];
        for (int i=0;i<n;){
            int j=i; double vi=v[idx[i]];
            while(j+1<n && v[idx[j+1]]==vi) j++;
            double avg=(i+j+2)/2.0;
            for (int k=i;k<=j;k++) r[idx[k]]=avg;
            i=j+1;
        }
        return r;
    }
    private static double spearman(double[] a, double[] b) {
        return pearson(rank(a), rank(b));
    }
    private static double[] logspace(double s, double e, int n) {
        double[] out = new double[n];
        double a = Math.pow(10, s), b = Math.pow(10, e);
        double r = Math.pow(b/a, 1.0/(n-1));
        out[0]=a; for(int i=1;i<n;i++) out[i]=out[i-1]*r;
        return out;
    }
    private static Integer[] indexOrder(double[] score, boolean desc) {
        Integer[] idx = new Integer[score.length];
        for (int i=0;i<score.length;i++) idx[i]=i;
        Arrays.sort(idx, (i,j)->desc?Double.compare(score[j],score[i]):Double.compare(score[i],score[j]));
        return idx;
    }

    /** 由初始 wStd 和 |corr| 产生恰好 K 个索引（>K 取Top-K；<K 用|corr|补足） */
    private static int[] pickKIndicesByMagnitudeOrCorr(double[] wStd, double[] absCorr, int K) {
        final int p = wStd.length;
        // 先按 |w| 排序
        Integer[] byW = indexOrder(Arrays.stream(wStd).map(Math::abs).toArray(), true);
        ArrayList<Integer> sel = new ArrayList<>();
        for (int j : byW) {
            if (Math.abs(wStd[j]) > 0) { sel.add(j); if (sel.size()==K) break; }
        }
        if (sel.size() == K) return sel.stream().mapToInt(Integer::intValue).toArray();

        // 如果不足K：先收集已选集合，按 |corr| 补足
        boolean[] used = new boolean[p];
        for (int j : sel) used[j]=true;
        Integer[] byCorr = indexOrder(absCorr, true);
        for (int j : byCorr) {
            if (!used[j]) { sel.add(j); used[j]=true; if (sel.size()==K) break; }
        }
        // 若仍不足（极端情况），补低序
        for (int j = 0; sel.size() < K && j < p; j++) if (!used[j]) sel.add(j);
        return sel.stream().mapToInt(Integer::intValue).toArray();
    }

    /** 在子集上做 ridge：返回 [b0, b1..bK]（带截距；输入是标准化后的 X） */
    private static double[] ridgeSolve(double[][] Xsub, double[] y, double ridge) {
        int n = Xsub.length, k = Xsub[0].length;
        int d = k + 1; // 截距 + k
        double[][] XtX = new double[d][d];
        double[] Xty = new double[d];

        // 先算 sums
        double sumY = 0; double[] sumX = new double[k];
        for (int i = 0; i < n; i++) {
            sumY += y[i];
            for (int j = 0; j < k; j++) sumX[j] += Xsub[i][j];
        }
        XtX[0][0] = n; Xty[0] = sumY;
        for (int j = 0; j < k; j++) {
            XtX[0][j+1] = sumX[j];
            XtX[j+1][0] = sumX[j];
            double sxy=0;
            for (int i=0;i<n;i++) sxy += Xsub[i][j]*y[i];
            Xty[j+1] = sxy;
        }
        for (int j = 0; j < k; j++) {
            for (int t = j; t < k; t++) {
                double sxx=0;
                for (int i=0;i<n;i++) sxx += Xsub[i][j]*Xsub[i][t];
                XtX[j+1][t+1]=sxx; XtX[t+1][j+1]=sxx;
            }
        }
        for (int j=1;j<d;j++) XtX[j][j] += ridge;
        for (int j=0;j<d;j++) XtX[j][j] += 1e-12;

        return solveSymmetric(XtX, Xty);
    }

    /** 从标准化 Xs 按 picked 列构建子矩阵 */
    private static double[][] buildSubset(double[][] Xs, int[] picked) {
        int n = Xs.length, k = picked.length;
        double[][] sub = new double[n][k];
        for (int i=0;i<n;i++) for (int c=0;c<k;c++) sub[i][c] = Xs[i][picked[c]];
        return sub;
    }

    /** 仅返回系数（不含截距），输入为标准化 Xs 与 picked 索引 */
    private static double[] refitRidgeOnSubset(double[][] Xs, double[] y, int[] picked, double ridge) {
        double[][] sub = buildSubset(Xs, picked);
        double[] beta = ridgeSolve(sub, y, ridge);
        return Arrays.copyOfRange(beta, 1, beta.length);
    }

    /** 高斯消元（部分选主元） */
    private static double[] solveSymmetric(double[][] A, double[] b) {
        int n = b.length;
        double[][] M = new double[n][n];
        for (int i=0;i<n;i++) M[i]=Arrays.copyOf(A[i], n);
        double[] y = Arrays.copyOf(b, n);

        for (int k=0;k<n;k++){
            int piv=k; double max=Math.abs(M[k][k]);
            for (int i=k+1;i<n;i++){ double v=Math.abs(M[i][k]); if (v>max){max=v;piv=i;} }
            if (max==0.0) throw new RuntimeException("Singular matrix.");
            if (piv!=k){ double[] tr=M[k]; M[k]=M[piv]; M[piv]=tr; double tv=y[k]; y[k]=y[piv]; y[piv]=tv; }

            double diag=M[k][k];
            for (int j=k;j<n;j++) M[k][j]/=diag; y[k]/=diag;

            for (int i=0;i<n;i++){
                if (i==k) continue;
                double f=M[i][k]; if (f==0) continue;
                for (int j=k;j<n;j++) M[i][j]-=f*M[k][j];
                y[i]-=f*y[k];
            }
        }
        return y;
    }

    // ===================== 结果对象 =====================

    public static class Result {
        public double bestLambda;
        public double bestSpearman;
        public double intercept;
        public double[] coef;                 // 长度 = 原始列数 p（未选=0）
        public List<Integer> keptColumns;     // 预筛保留列（常数列已剔除）
        public List<Integer> selectedColumns; // 最终选择的 K 个原始列索引
        public int k;

        @Override public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("bestLambda=").append(bestLambda)
                    .append(", bestSpearman=").append(bestSpearman)
                    .append(", intercept=").append(intercept)
                    .append(", K=").append(k)
                    .append(", selectedCols=").append(selectedColumns)
                    .append("\ncoef=");
            sb.append(Arrays.toString(coef));
            return sb.toString();
        }

        public double predict(ArrayList<Double> x) {
            if (x.size() == coef.length) {
                // 输入是完整长度
                double s = intercept;
                for (int j = 0; j < coef.length; j++) {
                    s += coef[j] * x.get(j);
                }
                return s;
            } else if (x.size() == selectedColumns.size()) {
                // 输入只有被选中的列
                double s = intercept;
                for (int i = 0; i < selectedColumns.size(); i++) {
                    int col = selectedColumns.get(i);
                    s += coef[col] * x.get(i);
                }
                return s;
            } else {
                throw new IllegalArgumentException(
                        "输入维度 " + x.size() + " 既不是完整长度 " + coef.length +
                                " 也不是选中长度 " + selectedColumns.size()
                );
            }
        }

    }



    // ===================== 示例 =====================
    public static void main(String[] args) {
        int n = 50, p = 25, K = 6;
        Random rnd = new Random(42);
        ArrayList<double[]> X = new ArrayList<>();
        ArrayList<Double> y = new ArrayList<>();
        for (int i=0;i<n;i++){
            double[] row=new double[p];
            for (int j=0;j<p;j++) row[j]=rnd.nextGaussian();
            X.add(row);
            // y = 0.8*x7 + 0.5*x12 - 0.3*x18 + 噪声
            double yy = 0.8*row[7]+0.5*row[12]-0.3*row[18]+rnd.nextGaussian()*0.3;
            y.add(yy);
        }
        SparseCaseSelector selector = new SparseCaseSelector(X, y);
        Result r = selector.fitK(/*useElasticNet=*/false, /*alpha=*/1.0, K, /*useLOOCV=*/true);
        System.out.println(r);
    }
}
