package yimei.jss.ruleevaluation;

import yimei.jss.simulation.DynamicSimulation;

import java.util.HashMap;
import java.util.Map;

public class DFJSSProblem {

    // ---------- Step 1: instance 类型 ----------
    private enum InstanceType {
        E, G, U, N
    }

    // ---------- Step 2: instance 配置模板 ----------
    protected static class SimConfig {
        int numMachines;
        int numJobs;
        String batchSize;
        double utilLevel;
        double dueDateFactor;
        DynamicSimulation.ShopType distribution;

        SimConfig(int numMachines, int numJobs, String batchSize,
                  double utilLevel, double dueDateFactor,
                  DynamicSimulation.ShopType distribution) {
            this.numMachines = numMachines;
            this.numJobs = numJobs;
            this.batchSize = batchSize;
            this.utilLevel = utilLevel;
            this.dueDateFactor = dueDateFactor;
            this.distribution = distribution;
        }
    }

    protected static final Map<Character, SimConfig> TEMPLATE = new HashMap<>();

    static {
        TEMPLATE.put('G', new SimConfig(2, 500, "small", 0.60, 1.2,
                DynamicSimulation.ShopType.GAMMA));

        TEMPLATE.put('U', new SimConfig(10, 3000, "single", 0.85, 1.5,
                DynamicSimulation.ShopType.UNIFORM));

        TEMPLATE.put('E', new SimConfig(5, 3000, "medium", 0.85, 2.0,
                DynamicSimulation.ShopType.EXPONENTIAL));

        TEMPLATE.put('N', new SimConfig(15, 3000, "large", 0.95, 4.0,
                DynamicSimulation.ShopType.NORMAL));
    }



}

