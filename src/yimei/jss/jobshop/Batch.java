package yimei.jss.jobshop;

import java.util.ArrayList;
import java.util.List;


/**
 * A batch consists of some jobs.
 *
 * Created by luyao on 20/12/2024.
 */

public class Batch {

    private int id;

    public int numOfJobs;
    public List<Job> jobs = new ArrayList<>();
    public double arrivalTime; //change it from private final to public
    public double releaseTime;
    public double dueDate;

    private double weight;

    public double totalProcTime;
    private double avgProcTime;

    private double completionTime;

    private double priority;

    public int completedJobs;

    public double revenue;

    public double estimatedProcTime;

    public Batch(double avgProcTime) {
        this.avgProcTime = avgProcTime;
    }

    public double getPriority() {
        return priority;
    }

    public void setPriority(double priority) {
        this.priority = priority;
    }

    public Batch(int id,
               int numberOfJobs,
               double totalProcTime,
               double arrivalTime,
               double releaseTime,
               double dueDate,
               double weight) {
        this.id = id;
        this.numOfJobs = numberOfJobs;
        this.totalProcTime = totalProcTime;
        this.arrivalTime = arrivalTime;
        this.releaseTime = releaseTime;
        this.dueDate = dueDate;
        this.weight = weight;
    }


    public int getId() {
        return id;
    }

    public Job getJob(int idx) {
        return jobs.get(idx);
    }

    public List<Job> getJobs() {
        return jobs;
    }

    public double getArrivalTime() {
        return arrivalTime;
    }

    public double getReleaseTime() {
        return releaseTime;
    }

    public double getDueDate() {
        return dueDate;
    }

    public double getWeight() {
        return weight;
    }

    public double getTotalProcTime() {
        return totalProcTime; //the time that really eas used to process the jobs
    }

    public int getNumOfJobs() {return numOfJobs;}

    public double getAvgProcTime() {
        return avgProcTime;
    }

    //fzhang 29.8.2018 get the unprocessingtime (waiting time)
    public double getWaitingTime() {
        return this.flowTime() - totalProcTime;
    }

    public double getCompletionTime() {
        return completionTime; //completionTime: the finish time(a time points)
    }

    public double getWorkLoadRemaining(){
        double value = 0;

        for (Job job : jobs) {
            value += job.workLoadRemaining;
        }
        return value;
    }

    public int getOperationRemaining(){
        int value = 0;

        for (Job job : jobs) {
            value += job.operationRemaining;
        }
        return value;
    }

    public void setDueDate(double dueDate) {
        this.dueDate = dueDate;
    }

    public void setCompletionTime(double completionTime) {
        this.completionTime = completionTime;
    }

    public double flowTime() {
        return completionTime - arrivalTime; // the time period between the job arrives and the job is finished. Including the waiting time
    }

    public double weightedFlowTime() {
        return weight * flowTime();
    }

    public double tardiness() {
        double tardiness = completionTime - dueDate;
        if (tardiness < 0)
            tardiness = 0;

        return tardiness;
    }

    public double weightedTardiness() {
        return weight * tardiness();
    }

    public void addJob(Job job) {
        jobs.add(job);
    }

    public String toString() {
        String string = String.format("Batch %d, arrives at %.1f, due at %.1f, weight is %.1f. It has %d jobs:\n",
                getId(), arrivalTime, dueDate, weight, jobs.size());

        return string;
    }

}
