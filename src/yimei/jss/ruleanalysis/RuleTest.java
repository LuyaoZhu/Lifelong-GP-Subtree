package yimei.jss.ruleanalysis;

import ec.gp.GPNode;
import ec.multiobjective.MultiObjectiveFitness;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.SchedulingSet;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.simulation.Simulation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class RuleTest {

    public static final long simSeed = 968356;

    protected String trainPath; //the directory of training things
    protected RuleType ruleType;
    protected int numRuns;
    protected String testScenario;
    protected String testSetName;
    protected List<Objective> objectives; // The objectives to test.
    protected int numPopulations;

    public RuleTest(String trainPath, RuleType ruleType, int numRuns,
                    String testScenario, String testSetName,
                    List<Objective> objectives, int numPopulations) {
        this.trainPath = trainPath;
        this.ruleType = ruleType;
        this.numRuns = numRuns;
        this.testScenario = testScenario;
        this.testSetName = testSetName;
        this.objectives = objectives;
        this.numPopulations = numPopulations;
    }

    public RuleTest(String trainPath, RuleType ruleType, int numRuns,
                    String testScenario, String testSetName, int numPopulations) {
        this(trainPath, ruleType, numRuns, testScenario, testSetName, new ArrayList<>(), numPopulations);
    }

    public String getTrainPath() {
        return trainPath;
    }

    public RuleType getRuleType() {
        return ruleType;
    }

    public int getNumRuns() {
        return numRuns;
    }

    public int getNumPopulations() {
        return numPopulations;
    }

    public String getTestScenario() {
        return testScenario;
    }

    public List<Objective> getObjectives() {
        return objectives;
    }

    public void setObjectives(List<Objective> objectives) {
        this.objectives = objectives;
    }

    public void addObjective(Objective objective) {
        this.objectives.add(objective);
    }

    public void addObjective(String objective) {
        addObjective(Objective.get(objective));
    }

    //generate testset using simseed, replications
    public SchedulingSet generateTestSet() {
        return SchedulingSet.generateSet(simSeed, testScenario,
                testSetName, objectives, 50);
    }

    public void writeToCSV() {
        SchedulingSet testSet = generateTestSet();
        Simulation.optimisationObjective = objectives.get(0).getName();

        File targetPath = new File(trainPath + "test"); //create a folder named "test" in trainPath
        if (!targetPath.exists()) {
            targetPath.mkdirs();
        }

        File csvFile = new File(targetPath + "/" + testSetName + ".csv"); //create a .csv to save the test result

        List<TestResult> testResults = new ArrayList<>();

        //for test: which machines are chosen by routing rule, CCGP. Scenario: run1, rule in generation 51
        for (int i = 0; i < numRuns; i++) {
            System.out.println("Run " + i);
        /*for (int i = 0; i < 3; i++) {
        	System.out.println("Run "+ i);*/

            File sourceFile = new File(trainPath + "job." + i + ".out.stat");  //this file keeps the rule
            TestResult result = TestResult.readFromFile(sourceFile, ruleType, numPopulations);

            //Didn't bother saving time files
            File timeFile = new File(trainPath + "job." + i + ".time.csv");
            result.setGenerationalTimeStat(ResultFileReader.readTimeFromFile(timeFile));

            //read average rule size in to .csv
            File aveRuleSizeFile = new File(trainPath + "job." + i + ".aveGenRulesize.csv");
            result.setGenerationalAveSequencingRuleSizeStat(ResultFileReader.readAveRuleSizeFromFile(aveRuleSizeFile, 1));
            result.setGenerationalAveRoutingRuleSizeStat(ResultFileReader.readAveRuleSizeFromFile(aveRuleSizeFile, 2));
//            if (numPopulations == 2) {
//                result.setGenerationalAveRuleSizeStat(ResultFileReader.readAveRuleSizeFromFile(aveRuleSizeFile, 3));
//            } else if (numPopulations == 3) {
                result.setGenerationalAveOrderingRuleSizeStat(ResultFileReader.readAveRuleSizeFromFile(aveRuleSizeFile, 3));
                result.setGenerationalAveRuleSizeStat(ResultFileReader.readAveRuleSizeFromFile(aveRuleSizeFile, 4));
//            }

            //24.8.2018 fzhang read badrun in CSV
//            File badrunsFile = new File(trainPath + "job." + i + ".BadRun.csv");
//            result.setGenerationalBadRunStat(ResultFileReader.readBadRunFromFile(badrunsFile));

            //2019.8.24 read the number of cleared inds
            //File clearedInds = new File(trainPath + "job." + i + ".clearedInds.csv");
//            File clearedInds = new File(trainPath + "job." + i + ".badInds.csv");
//            result.setGenerationalSeqClearedIndsStat(ResultFileReader.readClearedIndsFromFile(clearedInds, 1));
//            result.setGenerationalRouClearedIndsStat(ResultFileReader.readClearedIndsFromFile(clearedInds, 2));

            long start = System.currentTimeMillis();

//            result.validate(objectives);

            //for (int j = 90; j < result.getGenerationalRules().size(); j++) {
            for (int j = 0; j < result.getGenerationalRules().size(); j++) {
                AbstractRule[] generationalRules = result.getGenerationalRules(j);
                if (numPopulations == 2) {
                    generationalRules[0].calcFitness(  //in calcFitness(), it will check which one is routing/sequencing rule
                            result.getGenerationalTestFitness(j), null,
                            testSet, generationalRules[1], objectives);
                } else if (numPopulations == 1) { //if only one rule, we will use WIQ as fixed routing rule
                    AbstractRule routingRule = new WIQ(yimei.jss.rule.RuleType.ROUTING);
                    generationalRules[0].calcFitness(result.getGenerationalTestFitness(j), null,
                            testSet, routingRule, objectives);
                } else if (numPopulations == 3) {
                    generationalRules[0].calcFitness(  //in calcFitness(), it will check which one is routing/sequencing rule
                            result.getGenerationalTestFitness(j), null,
                            testSet, generationalRules[1], generationalRules[2], objectives);
                }

              /*  System.out.println("Generation " + j + ": test fitness = " +
                        result.getGenerationalTestFitness(j).fitness());*/

                //only print out the last 51 generation, but from 0 to 50
                System.out.println("Generation " + j + ": test fitness = " +
                        result.getGenerationalTestFitness(j).fitness());
            }

            long finish = System.currentTimeMillis();
            long duration = (finish - start) / 1000;
            System.out.println("Duration = " + duration + " s.");

            testResults.add(result);
        }

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(csvFile.getAbsoluteFile()));
            //with badrun
        /*    writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
                    "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, BadRun");*/

            //without badrun
            /*writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
                    "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime");*/

            //with average rule size
//            writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
//                    "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, AveSeqRulesize, AveRouRuleSize, AveRuleSize");

            //with the third rule including average rule size
//            if (numPopulations == 2) {
//                writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
//                        "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, AveSeqRulesize, AveRouRuleSize, AveRuleSize");
//            } else if (numPopulations == 3) {
                writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
                        "RoutRuleUniqueTerminals,OrdRuleSize,OrdRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, AveSeqRulesize, AveRouRuleSize, AveOrdRuleSize, AveRuleSize");
//            }

            //with the number of cleared inds
//            writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
//                    "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, AveSeqRulesize, AveRouRuleSize, AveRuleSize, SeqClearedInds, RouClearedInds");

            /*writer.write("Run,Generation,SeqRuleSize,SeqRuleUniqueTerminals,RoutRuleSize," +
                    "RoutRuleUniqueTerminals,Obj,TrainFitness,TestFitness,TrainTime, AveSeqRulesize, AveRouRuleSize, AveRuleSize, SeqBadInds, RouBadInds");*/

            writer.newLine();

            //for (int i = 0; i < 1; i++) {
            for (int i = 0; i < numRuns; i++) {
                TestResult result = testResults.get(i);

                //for (int j = 50; j < result.getGenerationalRules().size(); j++) { //use rules in each generation for testing
                for (int j = 0; j < result.getGenerationalRules().size(); j++) { //use rules in each generation for testing

                    MultiObjectiveFitness trainFit =
                            (MultiObjectiveFitness) result.getGenerationalTrainFitness(j);
                    MultiObjectiveFitness testFit =
                            (MultiObjectiveFitness) result.getGenerationalTestFitness(j);
                    GPRule[] rules = (GPRule[]) result.getGenerationalRules(j);
                    GPRule seqRule = null;
                    GPRule routRule = null;
                    GPRule ordRule = null;
                    if (numPopulations == 2) {
                        if (rules[0].getType() == yimei.jss.rule.RuleType.SEQUENCING) {
                            seqRule = rules[0];
                            routRule = rules[1];
                        } else {
                            seqRule = rules[1];
                            routRule = rules[0];
                        }
                    } else if (numPopulations == 1) {
                        seqRule = rules[0];
                    } else if (numPopulations == 3) {
                        if (rules[0].getType() == yimei.jss.rule.RuleType.SEQUENCING) {
                            seqRule = rules[0];
                            if (rules[1].getType() == yimei.jss.rule.RuleType.ROUTING) {
                                routRule = rules[1];
                                ordRule = rules[2];
                            } else {
                                routRule = rules[2];
                                ordRule = rules[1];
                            }
                        } else if (rules[0].getType() == yimei.jss.rule.RuleType.ROUTING) {
                            routRule = rules[0];
                            if (rules[1].getType() == yimei.jss.rule.RuleType.SEQUENCING) {
                                seqRule = rules[1];
                                ordRule = rules[2];
                            } else {
                                seqRule = rules[2];
                                ordRule = rules[1];
                            }
                        } else if (rules[0].getType() == yimei.jss.rule.RuleType.ORDERING) {
                            ordRule = rules[0];
                            if (rules[1].getType() == yimei.jss.rule.RuleType.SEQUENCING) {
                                seqRule = rules[1];
                                routRule = rules[2];
                            } else {
                                seqRule = rules[2];
                                routRule = rules[1];
                            }
                        }
                    }


                    UniqueTerminalsGatherer gatherer = new UniqueTerminalsGatherer();
                    int numUniqueTerminalsSeq = seqRule.getGPTree().child.numNodes(gatherer);
                    int seqRuleSize = seqRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);

                    int numUniqueTerminalsRout = 0;
                    int routRuleSize = 0;
                    if (numPopulations == 2) {
                        gatherer = new UniqueTerminalsGatherer();
                        numUniqueTerminalsRout = routRule.getGPTree().child.numNodes(gatherer);
                        routRuleSize = routRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);
                    }

                    int numUniqueTerminalsOrd = 0;
                    int ordRuleSize = 0;
                    if (numPopulations == 3) {
                        gatherer = new UniqueTerminalsGatherer();
                        numUniqueTerminalsRout = routRule.getGPTree().child.numNodes(gatherer);
                        routRuleSize = routRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);

                        numUniqueTerminalsOrd = ordRule.getGPTree().child.numNodes(gatherer);
                        ordRuleSize = ordRule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL);
                    }

                    //with bad run
                /*    if (objectives.size() == 1) {
                    	 writer.write(i + "," + j + "," +
                                seqRuleSize + "," +
                                numUniqueTerminalsSeq + "," +
                                routRuleSize +"," +
                                numUniqueTerminalsRout +",0," +
                                trainFit.fitness() + "," +
                                testFit.fitness()+ ","+ 
                                result.getGenerationalTime(j)+ ","+ 
                                result.getGenerationalBadRun(j)); //
                        writer.newLine();
                    }*/

                    //with average rule size
                    if (objectives.size() == 1) {
                        if (numPopulations == 2) {
                            writer.write(i + "," + j + "," +
                                    seqRuleSize + "," +
                                    numUniqueTerminalsSeq + "," +
                                    routRuleSize + "," +
                                    numUniqueTerminalsRout + "," +
                                    0 + "," +
                                    0 + ",0," +
                                    trainFit.fitness() + "," +
                                    testFit.fitness() + "," +
                                    result.getGenerationalTime(j) + "," +
                                    //result.getGenerationalBadRun(j) + "," +
                                    result.getGenerationalAveSequencingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveRoutingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveOrderingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveRuleSizeStatStat(j)
                            );
                        } else if (numPopulations == 3) {
                            writer.write(i + "," + j + "," +
                                    seqRuleSize + "," +
                                    numUniqueTerminalsSeq + "," +
                                    routRuleSize + "," +
                                    numUniqueTerminalsRout + "," +
                                    ordRuleSize + "," +
                                    numUniqueTerminalsOrd + ",0," +
                                    trainFit.fitness() + "," +
                                    testFit.fitness() + "," +
                                    result.getGenerationalTime(j) + "," +
                                    //result.getGenerationalBadRun(j) + "," +
                                    result.getGenerationalAveSequencingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveRoutingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveOrderingRuleSizeStatStat(j) + "," +
                                    result.getGenerationalAveRuleSizeStatStat(j)
                            );
                        }
                        writer.newLine();
                    }

                    //with the number of cleared inds
                /*    if (objectives.size() == 1) {
                        writer.write(i + "," + j + "," +
                                seqRuleSize + "," +
                                numUniqueTerminalsSeq + "," +
                                routRuleSize + "," +
                                numUniqueTerminalsRout + ",0," +
                                trainFit.fitness() + "," +
                                testFit.fitness() + "," +
                                result.getGenerationalTime(j) + "," +
                                //result.getGenerationalBadRun(j) + "," +
                                result.getGenerationalAveSequencingRuleSizeStatStat(j) + "," +
                                result.getGenerationalAveRoutingRuleSizeStatStat(j) + "," +
                                result.getGenerationalAveRuleSizeStatStat(j) + "," +
                                result.getsetGenerationalSeqClearedIndsStat(j) + "," +
                                result.getsetGenerationalRouClearedIndsStat(j)
                        ); //
                        writer.newLine();
                    }*/

                    //without badrun/average rule size
                  /*  if (objectives.size() == 1) {
                   	 writer.write(i + "," + j + "," +
                               seqRuleSize + "," +
                               numUniqueTerminalsSeq + "," +
                               routRuleSize +"," +
                               numUniqueTerminalsRout +",0," +
                               trainFit.fitness() + "," +
                               testFit.fitness()+ ","+ 
                               result.getGenerationalTime(j)); //
                       writer.newLine();
                   }*/
                    else {
//                        writer.write(i + "," + j + "," +
//                                rule.getGPTree().child.numNodes(GPNode.NODESEARCH_ALL) + "," +
//                                numUniqueTerminals + ",");

                        for (int k = 0; k < objectives.size(); k++) {
                            writer.write(k + "," +
                                    trainFit.getObjective(k) + "," +
                                    testFit.getObjective(k) + ",");
                        }
                        writer.newLine();
                    }
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    /**
     * Call this main method with several parameters
     * <p>
     * /Users/dyska/Desktop/Uni/COMP489/GPJSS/grid_results/dynamic/raw/coevolution-fixed/0.85-max-flowtime/
     * simple-rule
     * 30
     * dynamic-job-shop
     * missing-0.85-4.0
     * 2
     * 1
     * max-flowtime
     */
    public static void main(String[] args) {
        int idx = 0;
        String trainPath = args[idx];
        idx++;
        RuleType ruleType = RuleType.get(args[idx]);
        idx++;
        int numRuns = Integer.parseInt(args[idx]); //30
        idx++;
        String testScenario = args[idx]; //dynamic
        idx++;
        String testSetName = args[idx]; //missing-0.85-4.0
        idx++;
        int numPopulations = Integer.parseInt(args[idx]); //2
        idx++;
        int numObjectives = Integer.parseInt(args[idx]); //1
        idx++;
        RuleTest ruleTest = new RuleTest(trainPath, ruleType, numRuns, testScenario, testSetName, numPopulations);
        for (int i = 0; i < numObjectives; i++) {
            ruleTest.addObjective(args[idx]);
            idx++;
        }

        ruleTest.writeToCSV();
    }
}
