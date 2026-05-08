package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.util.Checkpoint;
import ec.util.Parameter;
import ec.util.SortComparatorL;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import smile.clustering.HierarchicalClustering;
import smile.clustering.linkage.CompleteLinkage;
import yimei.jss.algorithm.featureConstruction.FeatureConstructionCrossoverPipeline;
import yimei.jss.algorithm.surrogateCaseLS.SilhouetteScore;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.Objective;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.niching.*;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.AbstractRuleHelper;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.basic.SPT;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.rule.operation.weighted.WSPT;
import yimei.jss.rule.workcenter.basic.WIQ;
import yimei.jss.ruleevaluation.MultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import static yimei.jss.algorithm.surrogateAccuracy.surrogateClearingMultitreeEvaluatorV1.SamePCNum;
import static yimei.jss.gp.GPRun.out_dir;

/**
 * only preserve 10% individuals
 * normalize using two reference rules
 *
 * @author luyao
 */

public class GPRuleEvolutionStateLifelongGPV0 extends GPRuleEvolutionState {

    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;

    public double thresholdSeq;
    public double thresholdRou;

    ArrayList<Integer> numSeqArchive = new ArrayList<>();
    ArrayList<Integer> numRouArchive = new ArrayList<>();

    ArrayList<int[]> seqArchiveDepth = new ArrayList<>();
    ArrayList<int[]> rouArchiveDepth = new ArrayList<>();

    public ArrayList<GPTree> seqFeatureArchive = new ArrayList<>();
    public ArrayList<GPTree> rouFeatureArchive = new ArrayList<>();      // Save selected high-level features ( will be updated in every generation)

    public ArrayList<Individual[][]> elitesEveryGen = new ArrayList<>();

    public int[] seqFeatureDepth;
    public int[] rouFeatureDepth;

    public int pcDistance;

    public static final String P_REPLICATIONS = "num-Rep";
    public int numRep;

    public static final String P_ELITES = "num-elites";
    public static int numElites;

    public double coefficientCurrentTask;

    public double seedingRatio;

    public static int generationPerTask;

    ArrayList<Double> fromLastTaskRatio = new ArrayList<>();

    //define PC
    PhenoCharacterisation[] pc = new PhenoCharacterisation[2];

    RuleOptimizationProblem problem;

    @Override
    public void setup(EvolutionState state, Parameter base) {
        Parameter p;
        //fzhang 2018.11.8 I need to do this to be able to load seed values in the AbtractRule class.
        AbstractRuleHelper.state = this;

        // Get the job seed.
        p = new Parameter("seed").push("" + 0);
        jobSeed = parameters.getLongWithDefault(p, null, 0);

        setupTerminals();

        super.setup(this, base);

        generationPerTask =  parameters.getIntWithDefault(new Parameter("generationPerTask"), null, 100);

        seedingRatio = parameters.getDoubleWithDefault(new Parameter("seedingRatio"), null, 0.1);

        pc = scoreMultiPopCoevolutionaryClearingEvaluator.getPhenoCharacterisation();

        //define problem
        problem = (RuleOptimizationProblem) evaluator.p_problem;

//        RuleOptimizationProblem problem = (RuleOptimizationProblem) state.evaluator.p_problem;
//        phenoCharacterisation[0] = SequencingPhenoCharacterisation.currentGenPhenoCharacterisation((GPRuleEvolutionState) state, problem);
//        phenoCharacterisation[1] = RoutingPhenoCharacterisation.currentGenPhenoCharacterisation((GPRuleEvolutionState) state, problem);

    }

    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        //calculate the diversity
        double[][][] indsCharListsMultiTree = surrogateClearing.clearPopulation(this, pc); //3. calculate the phenotypic characteristic
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree[0]);
        diversityValue[1] = PopulationUtils.entropy(indsCharListsMultiTree[1]);
        entropyDiversity.add(diversityValue);

        // EVALUATION
        statistics.preEvaluationStatistics(this);

        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population

        statistics.postEvaluationStatistics(this);

        //After evaluate all individuals, we record feature information about 5 elites
        // SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }

        // SHOULD WE QUIT?
        if (generation == numGenerations - 1) {

            writeDiversityToFile(entropyDiversity);
            writeFromLastTaskRatioToFile(fromLastTaskRatio);

//            writeAccuracyToFile(MSEGen,SpearmanCorrelationGen,SamePCNum,AveragePCDistance);

//            writeTerimalOccuranceToFile(seqFeatureOccurrences, 0, "seq");
//            writeTerimalOccuranceToFile(rouFeatureOccurrences, 1, "rou");
//            if (ordFeatureOccurrences.size() >= 1)
//                writeTerimalOccuranceToFile(ordFeatureOccurrences, 2, "ord");

//            writeSpearmanToFile(spearmanEstimatedFitnessCurrentGen, "estimatedCurrent");

//            writeSpearmanToFile(spearmanRealFitnessCurrentGen, "realCurrent");

//            writeSpearmanToFile(spearmanEstimatedRealFitness, "estimatedReal");

//            writeFirstSelectedCasesNumToFile(firstSelectedCasesNum);

            generation++; // in this way, the last generation value will be printed properly.  fzhang 28.3.2018
            return R_FAILURE;
        }

        // PRE-BREEDING EXCHANGING
        statistics.prePreBreedingExchangeStatistics(this);

        population = exchanger.preBreedingExchangePopulation(this);  /** Simply returns state.population. */
        statistics.postPreBreedingExchangeStatistics(this);

        String exchangerWantsToShutdown = exchanger.runComplete(this);  /** Always returns null */
        if (exchangerWantsToShutdown != null) {
            output.message(exchangerWantsToShutdown);
            /*
             * Don't really know what to return here.  The only place I could
             * find where runComplete ever returns non-null is
             * IslandExchange.  However, that can return non-null whether or
             * not the ideal individual was found (for example, if there was
             * a communication error with the server).
             *
             * Since the original version of this code didn't care, and the
             * result was initialized to R_SUCCESS before the while loop, I'm
             * just going to return R_SUCCESS here.
             */

            return R_SUCCESS;
        }

        // BREEDING
        statistics.preBreedingStatistics(this);

        if(generation == generationPerTask-1 || generation == generationPerTask*2-1  || generation == generationPerTask*3-1 ){

            //try to only save 10% individuals to the next task

            ArrayList<Individual> savedTopIndividuals = new ArrayList<>();
            PopulationUtils.sort(population);

            for (int i=0; i<population.subpops[0].individuals.length * seedingRatio; i++ ) {
                savedTopIndividuals.add(population.subpops[0].individuals[i]);
            }

            //---------------------begin to collect subtrees------------------------------
            numElites = this.parameters.getIntWithDefault(new Parameter("selectedEliteNum1"), null, 10);
            Individual elites[][] = new Individual[this.population.subpops.length][numElites];
            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                for (int i = 0; i < numElites; i++) {
                    elites[sub][i] = (Individual) this.population.subpops[sub].individuals[i].clone();
                }
            }
            elitesEveryGen.add(elites);

            thresholdSeq = this.parameters.getDoubleWithDefault(new Parameter("thresholdSeq"), null, 0.8);
            thresholdRou = this.parameters.getDoubleWithDefault(new Parameter("thresholdRou"), null, 0.8);

            featureConstruction(elitesEveryGen);

            //-------------------------------end------------------------------------------------


            //-------------------------------begin to add features to terminals--------------------

/*            GPNode[][] newTerminals = new GPNode[this.population.subpops.length][];

            for (int sub=0; sub < this.population.subpops.length; sub++) {

                int newTerminalsNumber = 0;

                if(sub == 0) {
                    newTerminalsNumber = terminals[sub].length + seqFeatureArchive.get(generation/generationPerTask).size();
                } else {
                    newTerminalsNumber = terminals[sub].length + rouFeatureArchive.get(generation/generationPerTask).size();
                }

                newTerminals[sub] = new GPNode[newTerminalsNumber];

                for (int i = 0; i < newTerminals[sub].length; i++) {
                    if(i<terminals[sub].length){
                        newTerminals[sub][i] = terminals[sub][i];
                    } else {
                        GPNode feature;
                        if(sub == 0) {
                            feature = seqFeatureArchive.get(generation/generationPerTask).get(i-terminals[sub].length).child;
                        } else {
                            feature = rouFeatureArchive.get(generation/generationPerTask).get(i-terminals[sub].length).child;
                        }
                        newTerminals[sub][i] = feature;
                    }

                }

            }

            terminals = newTerminals;*/


            //-------------------------------end--------------------------------------------------

            population.clear();

            population = initializer.initialPopulation(this, 0);

            for (int sub = 0; sub < this.population.subpops.length; sub++) {
                for (int replace = 0; replace<savedTopIndividuals.size(); replace++) {
                    population.subpops[sub].individuals[replace] = savedTopIndividuals.get(replace);
                    ((GPIndividual)population.subpops[sub].individuals[replace]).fromLastTask = true;
                }
            }

        } else {
            population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N
        }
//        population = preselection();


        // POST-BREEDING EXCHANGING
        statistics.postBreedingStatistics(this);   //position 1  here, a new pop has been generated.

        // POST-BREEDING EXCHANGING
        statistics.prePostBreedingExchangeStatistics(this);

        population = exchanger.postBreedingExchangePopulation(this);   /** Simply returns state.population. */
        statistics.postPostBreedingExchangeStatistics(this);  //position 2


        // Generate new instances if needed
        if (problem.getEvaluationModel().isRotatable()) {
            problem.rotateEvaluationModel();
        }


        // INCREMENT GENERATION AND CHECKPOINT
        generation++;
        if (checkpoint && generation % checkpointModulo == 0) {
            output.message("Checkpointing");
            statistics.preCheckpointStatistics(this);
            Checkpoint.setCheckpoint(this);
            statistics.postCheckpointStatistics(this);
        }

        return R_NOTDONE;
    }

    private double[] calculateReferenceRuleFitness(DynamicSimulation simulation) {
        // 两个参考 rule 组合：SPT+WIQ 和 WSPT+WIQ
        AbstractRule referenceSeqRule1 = new SPT(RuleType.SEQUENCING);
        AbstractRule referenceRouRule1 = new WIQ(RuleType.ROUTING);
        AbstractRule referenceSeqRule2 = new WSPT(RuleType.SEQUENCING);
        AbstractRule referenceRouRule2 = new WIQ(RuleType.ROUTING);

        double[] fitness = new double[2];

        // 统一从参数中取一次 objective
        String objectiveName = parameters.getStringWithDefault(
                new Parameter("eval.problem.eval-model.objectives.0"), null, "");
        Objective objective = Objective.get(objectiveName);

        // 计算两个参考组合的 fitness
        fitness[0] = evaluateRulePair(simulation, referenceSeqRule1, referenceRouRule1, objective);
        fitness[1] = evaluateRulePair(simulation, referenceSeqRule2, referenceRouRule2, objective);

        return fitness;
    }

    /**
     * 使用给定的排序/路由规则运行一次仿真，返回该组合在指定 objective 下的值。
     * 如果仿真没有正常结束（clockTime == Double.MAX_VALUE），抛出异常。
     */
    private double evaluateRulePair(DynamicSimulation simulation,
                                    AbstractRule sequencingRule,
                                    AbstractRule routingRule,
                                    Objective objective) {

        simulation.setSequencingRule(sequencingRule);
        simulation.setRoutingRule(routingRule);
        simulation.run();

        try {
            if (simulation.getSystemState().getClockTime() == Double.MAX_VALUE) {
                // 说明这个 rule 组合导致系统跑挂了，这里直接报错
                throw new IllegalStateException(
                        "Reference rule pair " + sequencingRule.getClass().getSimpleName() +
                                " + " + routingRule.getClass().getSimpleName() +
                                " failed to finish the simulation (clockTime == Double.MAX_VALUE).");
            }

            return simulation.objectiveValue(objective);
        } finally {
            // 无论成功还是失败，都保证重置 simulation 状态
            simulation.reset();
        }
    }



    private void setupTerminals() {
        Parameter p;

        //Need to know how many populations we're expecting here, as will need
        //one terminal set per population
        int numSubPops = parameters.getInt(new Parameter("pop.subpops"), null);
        int numTrees = parameters.getInt(new Parameter("pop.subpop.0.species.ind.numtrees"), null);


        int num = Math.max(numSubPops, numTrees);

        if (num == 1) {

            p = new Parameter(P_TERMINALS_FROM);

            terminalsFrom = new String[]{parameters.getStringWithDefault(p,
                    null, "relative")};

            p = new Parameter(P_INCLUDE_ERC);
            //includeErc seems like does not have influence.
            includeErc = new boolean[]{parameters.getBoolean(p, null, false)};
            initTerminalSet();
        } else if (num == 2) {
            terminalsFrom = new String[num];
            includeErc = new boolean[num];
            int subPopNum = 0;

            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop1TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop1TerminalSet == null) {
                //might have provided other value by mistake, we should check for this
                p = new Parameter(P_TERMINALS_FROM);
                subPop1TerminalSet = parameters.getStringWithDefault(p,
                        null, "relative");
                output.warning("No terminal set for subpopulation 1 specified - using " + subPop1TerminalSet + ".");

            }
//            terminalsFrom[subPopNum] = subPop1TerminalSet;
            terminalsFrom[subPopNum] = subPop1TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop2TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop2TerminalSet == null) {
                //use whatever we settled on for first population
                subPop2TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 2 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop2TerminalSet;
            //TODO: Add support for erc - will be false by default

            initTerminalSet(); //right
        } else {
            terminalsFrom = new String[num];
            includeErc = new boolean[num];
            int subPopNum = 0;

            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop1TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop1TerminalSet == null) {
                //might have provided other value by mistake, we should check for this
                p = new Parameter(P_TERMINALS_FROM);
                subPop1TerminalSet = parameters.getStringWithDefault(p,
                        null, "relative");
                output.warning("No terminal set for subpopulation 1 specified - using " + subPop1TerminalSet + ".");

            }
            terminalsFrom[subPopNum] = subPop1TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop2TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop2TerminalSet == null) {
                //use whatever we settled on for first population
                subPop2TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 2 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop2TerminalSet;

            subPopNum++;
            p = new Parameter(P_TERMINALS_FROM + "." + subPopNum);
            String subPop3TerminalSet = parameters.getStringWithDefault(p,
                    null, null);
            if (subPop3TerminalSet == null) {
                //use whatever we settled on for first population
                subPop3TerminalSet = subPop1TerminalSet;
                output.warning("No terminal set for subpopulation 3 specified - using terminal set for subpopulation 1.");
            }
            terminalsFrom[subPopNum] = subPop3TerminalSet;

            //TODO: Add support for erc - will be false by default

            initTerminalSet(); //right
        }
    }


    //2021.4.16 fzhang save the diversity value to csv
    public void writeDiversityToFile(ArrayList<double[]> entropyDiversity) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".diversity.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!entropyDiversity.isEmpty()) {
                for (int j = 0; j < entropyDiversity.get(0).length; j++) {
                    header.append(",diversitySubpop").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < entropyDiversity.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                for (double value : entropyDiversity.get(i)) {
                    line.append(", ").append(value);
                }
                writer.write(line.toString());
                writer.newLine();
            }

            entropyDiversity.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeFromLastTaskRatioToFile(ArrayList<Double> fromLastTaskRatio) {
        // fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".fromLastTaskRatio.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            // Dynamically create the header based on the size of entropyDiversity elements
            StringBuilder header = new StringBuilder("Gen");
            if (!fromLastTaskRatio.isEmpty()) {
                header.append(",fromLastTaskRatio");

            }
            writer.write(header.toString());
            writer.newLine();

            // Write the data
            for (int i = 0; i < fromLastTaskRatio.size(); i++) {
                StringBuilder line = new StringBuilder(i + "");
                line.append(", ").append(fromLastTaskRatio.get(i));
                writer.write(line.toString());
                writer.newLine();
            }

            fromLastTaskRatio.clear();
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    static class EliteComparator implements SortComparatorL {
        Individual[] inds;

        public EliteComparator(Individual[] inds) {
            super();
            this.inds = inds;
        }

        public boolean lt(long a, long b) {
            return inds[(int) b].fitness.betterThan(inds[(int) a].fitness);
        }

        public boolean gt(long a, long b) {
            return inds[(int) a].fitness.betterThan(inds[(int) b].fitness);
        }
    }

    //calculate the occurrance of one terminal
    public int countSubstring(String rule, String terminal) {
        int count = 0;
        int index = 0;

        while ((index = rule.indexOf(terminal, index)) != -1) {
            count++;
            index += terminal.length();
        }

        if (terminal.equals("W")) {
            return count
                    - countSubstring(rule, "WIQ")
                    - countSubstring(rule, "WKR")
                    - countSubstring(rule, "MWT")
                    - countSubstring(rule, "OWT")
                    - countSubstring(rule, "AWIS")
                    - countSubstring(rule, "MWIS")
                    - countSubstring(rule, "BWKR")
                    - countSubstring(rule, "AMWTF")
                    - countSubstring(rule, "AWIQF");
        } else if (terminal.equals("PT")) {
            return count
                    - countSubstring(rule, "NPT")
                    - countSubstring(rule, "TPT")
                    - countSubstring(rule, "APTF");
        } else if (terminal.equals("R")) {
            return count
                    - countSubstring(rule, "WKR")
                    - countSubstring(rule, "NOR")
                    - countSubstring(rule, "RDD")
                    - countSubstring(rule, "BNOR")
                    - countSubstring(rule, "BWKR");
        } else if (terminal.equals("WKR")) {
            return count
                    - countSubstring(rule, "BWKR");  // 防止 BWKR 统计两次 WKR
        } else if (terminal.equals("NOR")) {
            return count
                    - countSubstring(rule, "BNOR");  // 防止 BNOR 统计两次 NOR
        } else if (terminal.equals("MWT")) {
            return count
                    - countSubstring(rule, "AMWTF");  // 防止 BWKR 统计两次 WKR
        } else if (terminal.equals("WIQ")) {
            return count
                    - countSubstring(rule, "AWIQF");  // 防止 BNOR 统计两次 NOR
        } else if (terminal.equals("NIQ")) {
            return count
                    - countSubstring(rule, "ANIQF");  // 防止 BNOR 统计两次 NOR
        } else {
            return count;
        }
    }

    //2021.4.16 fzhang save the diversity value to csv
    public void writeTerimalOccuranceToFile(ArrayList<int[]> featureOccurance, int t, String s) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + s + "FeatureOccurance.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < terminals[t].length; i++) {
                sb.append(terminals[t][i].toStringForHumans());
                if (i < terminals[t].length - 1) {
                    sb.append(","); // 添加逗号分隔符
                }
            }

            writer.write("Gen" + "," + sb.toString());
            writer.newLine();
            for (int i = 0; i < featureOccurance.size(); i++) { //every two into one generation
                //writer.newLine();
                String[] stringArray = Arrays.stream(featureOccurance.get(i))
                        .mapToObj(String::valueOf)
                        .toArray(String[]::new);

                String content = String.join(",", stringArray);

                writer.write(i + "," + content + "\n");
            }

            featureOccurance.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // 将 int[] 转换为 List<Integer>，确保可以用作 HashMap 的键
    private static List<Integer> arrayToList(int[] arr) {
        List<Integer> list = new ArrayList<>();
        for (int num : arr) {
            list.add(num);
        }
        return list;
    }

    public static int findKneePoint(List<Double> values) {
        if (values.size() < 3) return -1; // 至少需要 3 个点

        // 找到最小值和最大值的索引
        int minIndex = 0, maxIndex = 0;
        for (int i = 0; i < values.size(); i++) {
            if (values.get(i) < values.get(minIndex)) minIndex = i;
            if (values.get(i) > values.get(maxIndex)) maxIndex = i;
        }

        // 直线参数 (P1: minIndex, P2: maxIndex)
        double x1 = minIndex, y1 = values.get(minIndex);
        double x2 = maxIndex, y2 = values.get(maxIndex);

        // 计算点到直线的最大垂直距离
        int inflectionIndex = -1;
        double maxDistance = -1;

        for (int i = 0; i < values.size(); i++) {
            if (i == minIndex || i == maxIndex) continue; // 跳过最小和最大点

            double x0 = i, y0 = values.get(i);
            double distance = pointToLineDistance(x1, y1, x2, y2, x0, y0);

            if (distance > maxDistance) {
                maxDistance = distance;
                inflectionIndex = i;
            }
        }

        return inflectionIndex;
    }

    public static double pointToLineDistance(double x1, double y1, double x2, double y2, double x0, double y0) {
        double numerator = Math.abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1);
        double denominator = Math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
        return numerator / denominator;
    }

    public void writeFirstSelectedCasesNumToFile(ArrayList<Integer> casesNum) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".firstSelectedCasesNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,casesNum");
            writer.newLine();
            for (int i = 0; i < casesNum.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + casesNum.get(i) + "\n");
            }
            casesNum.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] computeRanks(double[] doubleArray) {
        int n = doubleArray.length;
        Integer[] indices = new Integer[n];

// 初始化索引数组
        for (int i = 0; i < n; i++) {
            indices[i] = i;
        }

// 按 doubleArray 的值排序索引数组
        Arrays.sort(indices, Comparator.comparingDouble(i -> doubleArray[i]));

// 生成排名数组
        double[] rankArray = new double[n];
        for (int rank = 0; rank < n; rank++) {
            rankArray[indices[rank]] = rank + 1; // 排名从 1 开始
        }

        return rankArray;
    }

    // 计算 Spearman 相关性矩阵
    public static double[][] computeSpearmanMatrix(double[][] A) {
        int numVectors = A.length;
        double[][] spearmanMatrix = new double[numVectors][numVectors];

        SpearmansCorrelation spearman = new SpearmansCorrelation();

        for (int i = 0; i < numVectors; i++) {
            for (int j = 0; j < numVectors; j++) {
                if (i == j) {
                    spearmanMatrix[i][j] = 1.0; // 自相关性 = 1
                } else {
                    spearmanMatrix[i][j] = spearman.correlation(A[i], A[j]);
                }
            }
        }
        return spearmanMatrix;
    }

    // 计算最佳 K（使用轮廓系数）
    public static int findBestK(double[][] distanceMatrix, int minK, int maxK) {
        int bestK = minK;
        double bestScore = -1;

        for (int k = minK; k <= maxK; k++) {
            int[] labels = hierarchicalClustering(distanceMatrix, k);
            double score = SilhouetteScore.computeSilhouette(distanceMatrix, labels);
            System.out.println("K = " + k + ", 轮廓系数 = " + score);

            if (score > bestScore) {
                bestScore = score;
                bestK = k;
            }
        }
        return bestK;
    }

    // 进行层次聚类
    public static int[] hierarchicalClustering(double[][] distanceMatrix, int numClusters) {
        // 使用 Complete Linkage 进行层次聚类
        HierarchicalClustering hc = new HierarchicalClustering(new CompleteLinkage(distanceMatrix));

        // 分割成 numClusters 个簇
        return hc.partition(numClusters);
    }

    // 计算距离矩阵 (1 - 相关性)
    public static double[][] computeDistanceMatrix(double[][] similarityMatrix) {
        int n = similarityMatrix.length;
        double[][] distanceMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                distanceMatrix[i][j] = 1 - similarityMatrix[i][j]; // 相关性转距离
            }
        }
        return distanceMatrix;
    }

    // 计算均值
    private static double mean(List<Double> values) {
        return values.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
    }

    // 计算标准差
    private static double stddev(List<Double> values, double mean) {
        double variance = values.stream()
                .mapToDouble(v -> Math.pow(v - mean, 2))
                .average()
                .orElse(0.0);
        return Math.sqrt(variance);
    }

    // Z-score 归一化
    private static List<Double> zScoreNormalize(List<Double> values) {
        double mean = mean(values);
        double std = stddev(values, mean);
        double epsilon = 1e-9;  // 避免标准差为 0
        return values.stream()
                .map(v -> (v - mean) / (std + epsilon))
                .collect(Collectors.toList());
    }

    private static List<Double> minMaxNormalize(List<Double> values) {
        double min = values.stream().min(Double::compareTo).orElse(0.0);
        double max = values.stream().max(Double::compareTo).orElse(1.0);
        double epsilon = 1e-9;
        return values.stream()
                .map(v -> (v - min) / (max - min + epsilon))
                .collect(Collectors.toList());
    }

    // 对 fitness 进行 Log 变换
    private static List<Double> logTransform(List<Double> values) {
        double epsilon = 1e-9;
        return values.stream()
                .map(v -> Math.log(v + epsilon))  // 避免 log(0)
                .collect(Collectors.toList());
    }

    public void writeSpearmanToFile(ArrayList<double[]> spearman, String s) {
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + s + "Spearman.csv");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // 写入标题行（header）
            StringBuilder header = new StringBuilder("Gen");
            if (!spearman.isEmpty()) {
                for (int j = 0; j < spearman.get(0).length; j++) {
                    header.append(",").append(j);
                }
            }
            writer.write(header.toString());
            writer.newLine();

            // 写入每一行数据
            for (int i = 0; i < spearman.size(); i++) {
                StringBuilder line = new StringBuilder();
                line.append(i + 10); // 第 i 代（generation）

                double[] row = spearman.get(i);
                for (double value : row) {
                    line.append(",").append(String.format("%.4f", value)); // 保留4位小数
                }

                writer.write(line.toString());
                writer.newLine();
            }

            spearman.clear(); // 清空列表
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void writeAccuracyToFile(ArrayList<Double> mseGen, ArrayList<Double> spearmanCorrelationGen,ArrayList<Integer> samePCNum, ArrayList<Double> avePCDistance) {

        File weightFile = new File(out_dir + "/job." + jobSeed + "." + "SurrogateAccuracy.csv");
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile))) {

            // 写入标题行（header）
            StringBuilder header = new StringBuilder("Gen");
            header.append(",").append("MSE").append(",").append("SpearmanCorrelation").append(",").append("SamePCNum").append(",").append("AvePCDistance");
            writer.write(header.toString());
            writer.newLine();

            // 写入每一行数据
            for (int i = 0; i < mseGen.size(); i++) {
                StringBuilder line = new StringBuilder();
                line.append(i); // 第 i 代（generation）

                line.append(",").append(String.format("%.4f", mseGen.get(i))).append(",").append(String.format("%.4f", spearmanCorrelationGen.get(i)))
                        .append(",").append(String.format("%d", SamePCNum.get(i))).append(",").append(String.format("%.4f", avePCDistance.get(i)));


                writer.write(line.toString());
                writer.newLine();
            }
            mseGen.clear();
            spearmanCorrelationGen.clear(); // 清空列表
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void featureConstruction(ArrayList<Individual[][]> elitesEveryGen) {

        PhenoCharacterisation[] pcTask = new PhenoCharacterisation[2];

        DynamicSimulation simulation = (DynamicSimulation) ((MultipleRuleEvaluationModel)problem.getEvaluationModel()).getSchedulingSet().getSimulations().get(generation/generationPerTask);

        // then calculate phenoCharacterisation and obtain some decision points
        pcTask[0] = SequencingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation);
        pcTask[1] = RoutingPhenoCharacterisation.currentTaskPhenoCharacterisation(simulation);


        //delete the features that are not related to the new elites in the archive
//        if (this.generation > 0) {
//            deleteUnimportantFeature(pc);
//        }
//        if(this.generation%5 == 0)
//            System.out.println("1");
        //select features into the archive

        for (int sub = 0; sub < elitesEveryGen.get(generation/generationPerTask).length; sub++) {
            for (int ind = 0; ind < elitesEveryGen.get(generation/generationPerTask)[sub].length; ind++) {
                if (sub == 0) {
                    decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pcTask[sub]).decisionSituations;
                } else {
                    decisionSituationsRouting = ((RoutingPhenoCharacterisation) pcTask[sub]).decisionSituations;
                }
                AbstractRule rule = null;
                GPTree tree = ((GPIndividual) (elitesEveryGen.get(generation/generationPerTask)[sub][ind])).trees[0];
                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                double[] score = new double[nonterminals];
                double[] frequency = new double[nonterminals];
                if (nonterminals > 1) {
                    if (sub == 0) {
                        double[][] score_matrix = new double[nonterminals][];
                        for (int i = 0; i < decisionSituationsSequencing.size(); i++) {

                            int decisionSize = decisionSituationsSequencing.get(i).getQueue().size();

                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                                score_matrix[numSubtree] = new double[decisionSize];

                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.SEQUENCING, nodeToTree);

                                SequencingDecisionSituation situation = decisionSituationsSequencing.get(i);
                                List<OperationOption> queue = situation.getQueue();
                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    queue.get(candiateNum).setPriority(rule.priority(queue.get(candiateNum), situation.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }

/*                                if(i == 0){ //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for(int b=0; b<elitesEveryGen.get(this.generation)[sub].length; b++){
                                        GPTree elite = ((GPIndividual)elitesEveryGen.get(this.generation)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num=0; num < subtreeNum; num++){
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if(eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
                                }*/

                            }

                            //2019.10.24 use the rank rather than the priority value directly
                            //========================================start==========================================
                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }
                            //========================================end============================================

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }

                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsSequencing.size();
                        }

                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {

                            if (Math.abs(score[numScore]) >= thresholdSeq) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(sequencingFeature.size() < 30){
//                                    sequencingFeature.add(nodeToTree);
//                                }
                                seqFeatureArchive.add(nodeToTree);
                            }
                        }

                    } else {

                        double[][] score_matrix = new double[nonterminals][];
                        //for (int i = 0; i < 1; i++) {
                        for (int i = 0; i < decisionSituationsRouting.size(); i++) {

                            int decisionSize = decisionSituationsRouting.get(i).getQueue().size();

                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {

                                score_matrix[numSubtree] = new double[decisionSize];

                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.ROUTING, nodeToTree);

                                RoutingDecisionSituation situation = decisionSituationsRouting.get(i);
                                List<OperationOption> queue = situation.getQueue();

                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    OperationOption operationOption = queue.get(candiateNum);
                                    queue.get(candiateNum).setPriority(rule.priority(operationOption, operationOption.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }

/*                                if (i == 0) { //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for (int b = 0; b < elitesEveryGen.get(this.generation)[sub].length; b++) {
                                        GPTree elite = ((GPIndividual) elitesEveryGen.get(this.generation)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num = 0; num < subtreeNum; num++) {
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if (eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
                                }*/


                            }

                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }
                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsRouting.size();
                        }
                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {

                            if (Math.abs(score[numScore]) >= thresholdRou) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(routingFeature.size() < 30){
//                                    routingFeature.add(nodeToTree);
//                                }
                                rouFeatureArchive.add(nodeToTree);
                            }
                        }

                    }
                }
            }
        }


        //delete duplicate features
        deleteDuplicateFeatures();

//        if (this.generation < 30)
        deleteBigRedundantFeatures();
//        else {
//            deleteDuplicateFeatures();
//        }

        // record the depth of the features
        seqFeatureDepth = new int[6];
        rouFeatureDepth = new int[6];
        int sumSeqFeature = 0;

        int gen = generation/generationPerTask;
        if (!(seqFeatureArchive.size() == 0)) {
//            System.out.println("the seqFeature number in generation " + (gen + 1) + " is " + seqFeatureArchive.get(gen).size());
            for (int i = 0; i < seqFeatureArchive.size(); i++) {
//                System.out.println(seqFeatureArchive.get(gen).get(i).child.makeLispTree());
                sumSeqFeature++;
                if (seqFeatureArchive.get(i).child.depth() == 2)
                    seqFeatureDepth[0]++;
                if (seqFeatureArchive.get(i).child.depth() == 3)
                    seqFeatureDepth[1]++;
                if (seqFeatureArchive.get(i).child.depth() == 4)
                    seqFeatureDepth[2]++;
                if (seqFeatureArchive.get(i).child.depth() == 5)
                    seqFeatureDepth[3]++;
                if (seqFeatureArchive.get(i).child.depth() == 6)
                    seqFeatureDepth[4]++;
                if (seqFeatureArchive.get(i).child.depth() == 7)
                    seqFeatureDepth[5]++;
            }
        }

        seqArchiveDepth.add(seqFeatureDepth);
        numSeqArchive.add(sumSeqFeature);
        System.out.println("The number of feature in sequencing Archive is  " + sumSeqFeature);

        int sumRouFeature = 0;

        if (!(rouFeatureArchive.size() == 0)) {
//            System.out.println("the rouFeature number in generation " + (gen + 1) + " is " + rouFeatureArchive.get(gen).size());
            for (int i = 0; i < rouFeatureArchive.size(); i++) {
//                System.out.println(rouFeatureArchive.get(gen).get(i).child.makeLispTree());
                sumRouFeature++;
                if (rouFeatureArchive.get(i).child.depth() == 2)
                    rouFeatureDepth[0]++;
                if (rouFeatureArchive.get(i).child.depth() == 3)
                    rouFeatureDepth[1]++;
                if (rouFeatureArchive.get(i).child.depth() == 4)
                    rouFeatureDepth[2]++;
                if (rouFeatureArchive.get(i).child.depth() == 5)
                    rouFeatureDepth[3]++;
                if (rouFeatureArchive.get(i).child.depth() == 6)
                    rouFeatureDepth[4]++;
                if (rouFeatureArchive.get(i).child.depth() == 7)
                    rouFeatureDepth[5]++;
            }
        }

        rouArchiveDepth.add(rouFeatureDepth);
        numRouArchive.add(sumRouFeature);
        System.out.println("The number of feature in routing Archive is  " + sumRouFeature);


    }

    private void deleteBigRedundantFeatures() {

        for (int a = 0; a < seqFeatureArchive.size(); a++) {
            GPTree base = (GPTree) seqFeatureArchive.get(a).clone();
            for (int b = a + 1; b < seqFeatureArchive.size(); b++) {
                GPTree count = (GPTree) seqFeatureArchive.get(b).clone();
                if (FeatureConstructionCrossoverPipeline.LispContains(base.child, count.child)) {
                    seqFeatureArchive.remove(seqFeatureArchive.get(a));
                    a--;
                    break;
                }
            }
        }

        for (int a = 0; a < rouFeatureArchive.size(); a++) {
            GPTree base = (GPTree) rouFeatureArchive.get(a).clone();
            for (int b = a + 1; b < rouFeatureArchive.size(); b++) {
                GPTree count = (GPTree) rouFeatureArchive.get(b).clone();
                if (FeatureConstructionCrossoverPipeline.LispContains(base.child, count.child)) {
                    rouFeatureArchive.remove(rouFeatureArchive.get(a));
                    a--;
                    break;
                }
            }
        }

    }

    private void deleteDuplicateFeatures() {

        // for sequencing rule
/*        for (int gen1 = 0; gen1 < seqFeatureArchive.size(); gen1++) {
            for (int seqGenFea1 = 0; seqGenFea1 < seqFeatureArchive.get(gen1).size(); seqGenFea1++) {
                GPTree base = (GPTree) seqFeatureArchive.get(gen1).get(seqGenFea1).clone();

                for (int gen2 = 0; gen2 < seqFeatureArchive.size(); gen2++ ) {
                    for (int seqGenFea2 = 0; seqGenFea2 < seqFeatureArchive.get(gen2).size(); seqGenFea2++) {
                        GPTree count = (GPTree) seqFeatureArchive.get(gen2).get(seqGenFea2).clone();
                        if (!((gen1 == gen2) && (seqGenFea1 == seqGenFea2))) {
                            if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                                seqFeatureArchive.get(gen2).remove(seqFeatureArchive.get(gen2).get(seqGenFea2));
                                seqGenFea2--;
                            }
                        }
                    }
                }
            }
        }*/


        for (int a = 0; a < seqFeatureArchive.size(); a++) {
            GPTree base = (GPTree) seqFeatureArchive.get(a).clone();
            for (int b = a + 1; b < seqFeatureArchive.size(); b++) {
                GPTree count = (GPTree) seqFeatureArchive.get(b).clone();
                if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                    seqFeatureArchive.remove(seqFeatureArchive.get(b));
                    b--;
                }
            }
        }

        // for routing rule
/*        for (int gen1 = 0; gen1 < rouFeatureArchive.size(); gen1++) {
            for (int rouGenFea1 = 0; rouGenFea1 < rouFeatureArchive.get(gen1).size(); rouGenFea1++) {
                GPTree base = (GPTree) rouFeatureArchive.get(gen1).get(rouGenFea1).clone();

                for (int gen2 = 0; gen2 < rouFeatureArchive.size(); gen2++) {
                    for (int rouGenFea2 = 0; rouGenFea2 < rouFeatureArchive.get(gen2).size(); rouGenFea2++) {
                        GPTree count = (GPTree) rouFeatureArchive.get(gen2).get(rouGenFea2).clone();
                        if (!((gen1 == gen2) && (rouGenFea1 == rouGenFea2))) {
                            if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                                rouFeatureArchive.get(gen2).remove(rouFeatureArchive.get(gen2).get(rouGenFea2));
                                rouGenFea2--;
                            }
                        }
                    }
                }
            }
        }*/

        for (int a = 0; a < rouFeatureArchive.size(); a++) {
            GPTree base = (GPTree) rouFeatureArchive.get(a).clone();
            for (int b = a + 1; b < rouFeatureArchive.size(); b++) {
                GPTree count = (GPTree) rouFeatureArchive.get(b).clone();
                if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                    rouFeatureArchive.remove(rouFeatureArchive.get(b));
                    b--;
                }
            }
        }
    }



}


