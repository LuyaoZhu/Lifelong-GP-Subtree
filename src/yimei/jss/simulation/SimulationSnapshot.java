package yimei.jss.simulation;

import java.util.*;

public final class SimulationSnapshot {

    // —— 顶层标量（便于复现/审计） ——
    public final long seed;
    public final double clockTime;
    public final int numWorkCenters;
    public final int throughput;
    public final int warmupBatches;
    public final int numBatchesRecorded;
    public final double utilLevel;
    public final double dueDateFactor;

    // —— 核心状态 ——
    public final List<WorkCenterSnapshot> workCenters;      // 所有工作中心快照
    public final List<JobSnapshot> jobsInSystem;            // 当前系统内的作业
    public final List<EventSnapshot> futureEvents;          // 事件队列（按当前顺序）
    public final Map<Integer, Integer> workCenterCandidateOps; // WC id -> 候选操作数（=队列长度）

    public SimulationSnapshot(
            long seed, double clockTime, int numWorkCenters,
            int throughput, int warmupBatches, int numBatchesRecorded,
            double utilLevel, double dueDateFactor,
            List<WorkCenterSnapshot> workCenters,
            List<JobSnapshot> jobsInSystem,
            List<EventSnapshot> futureEvents,
            Map<Integer, Integer> workCenterCandidateOps) {
        this.seed = seed;
        this.clockTime = clockTime;
        this.numWorkCenters = numWorkCenters;
        this.throughput = throughput;
        this.warmupBatches = warmupBatches;
        this.numBatchesRecorded = numBatchesRecorded;
        this.utilLevel = utilLevel;
        this.dueDateFactor = dueDateFactor;
        this.workCenters = Collections.unmodifiableList(workCenters);
        this.jobsInSystem = Collections.unmodifiableList(jobsInSystem);
        this.futureEvents = Collections.unmodifiableList(futureEvents);
        this.workCenterCandidateOps = Collections.unmodifiableMap(workCenterCandidateOps);
    }

    // ===== 子快照结构 =====

    public static final class WorkCenterSnapshot {
        public final int id;
        public final int numMachines;
        public final List<Double> machineReadyTimes; // 与 WorkCenter 字段一致
        public final double readyTime;               // = min(machineReadyTimes)
        public final double workInQueue;             // 与 WorkCenter 字段一致
        public final double busyTime;                // 与 WorkCenter 字段一致
        public final List<OptionRef> queue;          // 队列里的 OperationOption 的关键信息
        public final int numOpsInQueue;              // 候选操作数（= queue.size()）

        public WorkCenterSnapshot(int id,
                                  int numMachines,
                                  List<Double> machineReadyTimes,
                                  double readyTime,
                                  double workInQueue,
                                  double busyTime,
                                  List<OptionRef> queue) {
            this.id = id;
            this.numMachines = numMachines;
            this.machineReadyTimes = Collections.unmodifiableList(machineReadyTimes);
            this.readyTime = readyTime;
            this.workInQueue = workInQueue;
            this.busyTime = busyTime;
            this.queue = Collections.unmodifiableList(queue);
            this.numOpsInQueue = queue.size();
        }
    }

    // 从 OperationOption 提取的最小必要引用
    public static final class OptionRef {
        public final int jobId;
        public final int opId;        // Operation.getId()
        public final int optionId;    // OperationOption.getOptionId()
        public final double procTime; // OperationOption.getProcTime()
        public final double readyTime;// OperationOption.getReadyTime()
        public OptionRef(int jobId, int opId, int optionId, double procTime, double readyTime) {
            this.jobId = jobId;
            this.opId = opId;
            this.optionId = optionId;
            this.procTime = procTime;
            this.readyTime = readyTime;
        }
    }

    public static final class JobSnapshot {
        public final int id;
        public final double arrivalTime;
        public final double releaseTime;
        public final double dueDate;
        public final double weight;
        public final double revenue;
        public final double totalProcTime;   // Job.getTotalProcTime()
        public final double avgProcTime;     // Job.getAvgProcTime()
        public final double completionTime;  // Job.getCompletionTime()
        public final double priority;        // Job.getPriority()
        public final double workLoadRemaining; // Job.workLoadRemaining
        public final int operationRemaining;   // Job.operationRemaining
        public final Integer batchId;          // 若有批次
        public final List<OperationSnapshot> operations; // 每道工序与其选项

        public JobSnapshot(int id, double arrivalTime, double releaseTime, double dueDate,
                           double weight, double revenue, double totalProcTime, double avgProcTime,
                           double completionTime, double priority, double workLoadRemaining,
                           int operationRemaining, Integer batchId, List<OperationSnapshot> operations) {
            this.id = id;
            this.arrivalTime = arrivalTime;
            this.releaseTime = releaseTime;
            this.dueDate = dueDate;
            this.weight = weight;
            this.revenue = revenue;
            this.totalProcTime = totalProcTime;
            this.avgProcTime = avgProcTime;
            this.completionTime = completionTime;
            this.priority = priority;
            this.workLoadRemaining = workLoadRemaining;
            this.operationRemaining = operationRemaining;
            this.batchId = batchId;
            this.operations = Collections.unmodifiableList(operations);
        }
    }

    public static final class OperationSnapshot {
        public final int opId; // Operation.getId()
        public final List<OperationOptionSnapshot> options;
        public OperationSnapshot(int opId, List<OperationOptionSnapshot> options) {
            this.opId = opId;
            this.options = Collections.unmodifiableList(options);
        }
    }

    public static final class OperationOptionSnapshot {
        public final int optionId;
        public final int workCenterId;
        public final double procTime;
        public OperationOptionSnapshot(int optionId, int workCenterId, double procTime) {
            this.optionId = optionId;
            this.workCenterId = workCenterId;
            this.procTime = procTime;
        }
    }

    public static final class EventSnapshot {
        public final String type;     // 事件类型名
        public final double time;     // 触发时间
        public final Integer entityId;// 可关联 batch/job id（能取则取）
        public EventSnapshot(String type, double time, Integer entityId) {
            this.type = type;
            this.time = time;
            this.entityId = entityId;
        }
    }
}