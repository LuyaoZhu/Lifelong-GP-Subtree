package yimei.jss.algorithm.surrogateCaseLS;

import java.util.*;

public class KneePointFinder {

    // 结果类：包含 knee 值（为 double）和原始 index
    public static class KneeResult {
        public double kneeValue;
        public int originalIndex;

        public KneeResult(double value, int index) {
            this.kneeValue = value;
            this.originalIndex = index;
        }

        @Override
        public String toString() {
            return "Knee value = " + kneeValue + " (original index = " + originalIndex + ")";
        }
    }

    // 找拐点（已排序值，double 版本）
    private static int findKneeIndexSorted(ArrayList<Double> sorted) {
        int n = sorted.size();
        if (n < 3) throw new IllegalArgumentException("List too short");

        double x0 = 0, y0 = sorted.get(0);
        double x1 = n - 1, y1 = sorted.get(n - 1);

        double maxDist = -1;
        int kneeIndex = -1;

        for (int i = 1; i < n - 1; i++) {
            double xi = i, yi = sorted.get(i);
            double num = Math.abs((y1 - y0) * xi - (x1 - x0) * yi + x1 * y0 - y1 * x0);
            double den = Math.sqrt(Math.pow(y1 - y0, 2) + Math.pow(x1 - x0, 2));
            double dist = num / den;

            if (dist > maxDist) {
                maxDist = dist;
                kneeIndex = i;
            }
        }

        return kneeIndex;
    }

    // 泛型方法，支持 Integer 或 Double 等 Number 类型
    public static <T extends Number & Comparable<T>> KneeResult findKneePoint(ArrayList<T> yList, int offset) {
        int n = yList.size();
        if (n < offset + 3) throw new IllegalArgumentException("List too short after offset");

        // 去除前 offset 个元素
        List<T> trimmed = yList.subList(offset, n);

        // 创建 (value, originalIndex) 对
        List<double[]> paired = new ArrayList<>();
        for (int i = 0; i < trimmed.size(); i++) {
            paired.add(new double[]{trimmed.get(i).doubleValue(), i + offset});
        }

        // 按 value 升序排序
        paired.sort(Comparator.comparingDouble(a -> a[0]));

        // 提取排序值
        ArrayList<Double> sortedValues = new ArrayList<>();
        for (double[] p : paired) {
            sortedValues.add(p[0]);
        }

        // 找拐点
        int kneeIndexInSorted = findKneeIndexSorted(sortedValues);

        // 映射回原始位置和值
        double kneeValue = sortedValues.get(kneeIndexInSorted);
        int originalIndex = (int) paired.get(kneeIndexInSorted)[1];

        return new KneeResult(kneeValue, originalIndex);
    }

    // 找拐点（已排序值，double 版本）
    private static int findKneeIndexSorted(double[] sorted) {
        int n = sorted.length;
        if (n < 3) throw new IllegalArgumentException("List too short");

        double x0 = 0, y0 = sorted[0];
        double x1 = n - 1, y1 = sorted[n - 1];

        double maxDist = -1;
        int kneeIndex = -1;

        for (int i = 1; i < n - 1; i++) {
            double xi = i, yi = sorted[i];
            double num = Math.abs((y1 - y0) * xi - (x1 - x0) * yi + x1 * y0 - y1 * x0);
            double den = Math.sqrt(Math.pow(y1 - y0, 2) + Math.pow(x1 - x0, 2));
            double dist = num / den;

            if (dist > maxDist) {
                maxDist = dist;
                kneeIndex = i;
            }
        }
        return kneeIndex;
    }

    // 核心方法：输入 double[]，输出 KneeResult
    public static KneeResult findKneePoint(double[] yArray, int offset) {
        int n = yArray.length;
        if (n < offset + 3) throw new IllegalArgumentException("List too short after offset");

        // 创建 (value, originalIndex) 对
        double[][] paired = new double[n - offset][2];
        for (int i = offset; i < n; i++) {
            paired[i - offset][0] = yArray[i];
            paired[i - offset][1] = i; // 原始索引
        }

        // 按 value 升序排序
        Arrays.sort(paired, Comparator.comparingDouble(a -> a[0]));

        // 提取排序值
        double[] sortedValues = new double[paired.length];
        for (int i = 0; i < paired.length; i++) {
            sortedValues[i] = paired[i][0];
        }

        // 找拐点
        int kneeIndexInSorted = findKneeIndexSorted(sortedValues);

        // 映射回原始位置和值
        double kneeValue = sortedValues[kneeIndexInSorted];
        int originalIndex = (int) paired[kneeIndexInSorted][1];

        return new KneeResult(kneeValue, originalIndex);
    }

    // 测试入口
    public static void main(String[] args) {
        ArrayList<Integer> intList = new ArrayList<>(Arrays.asList(
                0, 3, 1, 0, 337, 2, 21, 40, 14, 494,
                547, 119, 131, 506, 271, 504, 252, 180, 472, 550, 476
        ));

        ArrayList<Double> doubleList = new ArrayList<>();
        for (Integer i : intList) {
            doubleList.add(i * 1.0); // 构造 double 版
        }

        int offset = 5;

        // Integer 版本调用
        KneeResult resultInt = findKneePoint(intList, offset);
        System.out.println("Integer list: " + resultInt);

        // Double 版本调用
        KneeResult resultDouble = findKneePoint(doubleList, offset);
        System.out.println("Double list: " + resultDouble);

        // 额外：用拐点值筛选高值个体
        System.out.println("Values above knee point:");
        for (int i = 0; i < intList.size(); i++) {
            if (intList.get(i) > resultInt.kneeValue) {
                System.out.println("  index = " + i + ", value = " + intList.get(i));
            }
        }
    }
}
