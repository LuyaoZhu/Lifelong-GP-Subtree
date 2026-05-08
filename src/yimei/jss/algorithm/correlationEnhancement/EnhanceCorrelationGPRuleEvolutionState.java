package yimei.jss.algorithm.correlationEnhancement;

import ec.Individual;
import ec.gp.*;
import ec.util.Checkpoint;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.niching.*;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static yimei.jss.gp.GPRun.out_dir;

/**
 * The evolution state of evolving dispatching rules with GP.
 *
 * @author yimei
 */

public class EnhanceCorrelationGPRuleEvolutionState extends GPRuleEvolutionState {

    /**
     * Read the file to specify the terminals.
     */

    public Individual[][] selIndis = new Individual[2][250];   //selected individuals by knee-point

    public ArrayList<int[]> seletedIndsKneePointEveryGen = new ArrayList<>();   //record the selected individuals number after Using Knee-point

    public int[] seletedIndsKneePoint;

    public ArrayList<int[]> indsChangedFunctionEveryGen = new ArrayList<>();    //record the selected individuals number after Using Knee-point

    public int[] indsChangedFunction;

    public ArrayList<int[][]> caseNumEveryGen = new ArrayList<>();         //record the number of 6 cases occurred

    public int[][] caseNum;

    public boolean firstChange;

    public boolean changed;

    public double[][] beforeChangedAllFitnessGen = new double[2][250];  //record all individuals' fitness before changing

    public double[][][] beforeAfterChangedFitnessGen;   //subpop ,individuals, before and after fitness

    public boolean changeElites = true;


    @Override
    public int evolve() {
        if (generation > 0)
            output.message("Generation " + generation);

        //record the diversity
        PhenoCharacterisation[] pc = scoreMultiPopCoevolutionaryClearingEvaluatorV1.getPhenoCharacterisation();
        double[][][] indsCharListsMultiTree = surrogateClearing.clearPopulation(this, pc); //3. calculate the phenotypic characteristic
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree[0]);
        diversityValue[1] = PopulationUtils.entropy(indsCharListsMultiTree[1]);
        entropyDiversity.add(diversityValue);
//
//		System.out.println(diversityValue[0] + ",," + diversityValue[1]);

        //System.out.println("generation "+generation);

        // EVALUATION
        statistics.preEvaluationStatistics(this);

        if (generation == 10 || generation == 30 || generation == 50 || generation == 70 || generation == 90) {

            changeElites = false;
            evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population
            statistics.postEvaluationStatistics(this);
            if (generation >= 0) {
                for (int subIndex = 0; subIndex < population.subpops.length; subIndex++) {
                    Individual[] individuals = population.subpops[subIndex].individuals;
                    for (int indIndex = 0; indIndex < population.subpops[subIndex].individuals.length; indIndex++) {
                        //System.out.println(individuals[indIndex].fitness.fitness());
//                    if (individuals[indIndex].fitness.fitness() != Double.POSITIVE_INFINITY & individuals[indIndex].fitness.fitness() != Double.MAX_VALUE)
                        beforeChangedAllFitnessGen[subIndex][indIndex] = individuals[indIndex].fitness.fitness();
                    }
                }
            }
        }
        //---------------------------Enhance correlation-----------------------------------------------

        //step1: decide to select which individuals for enhancing correlation(use knee-point to select promising individuals)

        caseNum = new int[2][6];
        seletedIndsKneePoint = new int[2];
        indsChangedFunction = new int[2];

        //step2: for each selected individuals, calculate their subtrees' correlations

        List<SequencingDecisionSituation> decisionSituationsSequencing = null;
        List<RoutingDecisionSituation> decisionSituationsRouting = null;

        for (int i = 0; i < population.subpops.length; i++) {
            selIndis[i] = population.subpops[i].individuals.clone();
        }

        ArrayList<ArrayList<Integer>> changedIndIndex = new ArrayList<>();
        ArrayList<Integer> changedIndIndexPop0 = new ArrayList<>();
        ArrayList<Integer> changedIndIndexPop1 = new ArrayList<>();

        for (int sub = 0; sub < selIndis.length; sub++) {
            if (sub == 0) {
                decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pc[sub]).decisionSituations;
            } else {
                decisionSituationsRouting = ((RoutingPhenoCharacterisation) pc[sub]).decisionSituations;
            }
            for (int ind = 0; ind < selIndis[sub].length; ind++) {
                AbstractRule rule = null;
                GPTree tree = ((GPIndividual) (selIndis[sub][ind])).trees[0];
                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                double[] score = new double[nonterminals];
                if (nonterminals > 1) {
                    if (sub == 0) {
                        int decisionSize = decisionSituationsSequencing.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        for (int i = 0; i < decisionSituationsSequencing.size(); i++) {
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {

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

                        //begin to enhance correlation

                        firstChange = true;
                        changed = false;
                        enhanceCorrelation((GPNode) ((GPIndividual) selIndis[sub][ind]).trees[0].child.clone(), score, sub);
                        if (changed)
                            changedIndIndexPop0.add(ind);

                    } else {
                        int decisionSize = decisionSituationsRouting.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        for (int i = 0; i < decisionSituationsRouting.size(); i++) {
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
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
                        //begin to enhance correlation
                        firstChange = true;
                        changed = false;
                        enhanceCorrelation(((GPIndividual) selIndis[sub][ind]).trees[0].child, score, sub);
                        if (changed)
                            changedIndIndexPop1.add(ind);
                    }
                }
            }
        }

        changedIndIndex.add(changedIndIndexPop0);
        changedIndIndex.add(changedIndIndexPop1);

        System.out.println("------After Changing the functions------");
        evaluator.evaluatePopulation(this);  //// here, after this we evaluate the population
        statistics.postEvaluationStatistics(this);


        if (generation == 10 || generation == 30 || generation == 50 || generation == 70 || generation == 90) {
//        if (generation >= 0) {

            beforeAfterChangedFitnessGen = new double[population.subpops.length][][];

            for (int subIndex = 0; subIndex < population.subpops.length; subIndex++) {

                beforeAfterChangedFitnessGen[subIndex] = new double[changedIndIndex.get(subIndex).size()][2];

                Individual[] individuals = population.subpops[subIndex].individuals;

                for (int ind=0; ind<changedIndIndex.get(subIndex).size(); ind++){
                    beforeAfterChangedFitnessGen[subIndex][ind][0] = beforeChangedAllFitnessGen[subIndex][changedIndIndex.get(subIndex).get(ind)];
                    beforeAfterChangedFitnessGen[subIndex][ind][1] = individuals[changedIndIndex.get(subIndex).get(ind)].fitness.fitness();
                }
            }
            writeSeqChangedBeforeAfterFitnessToFile(beforeAfterChangedFitnessGen[0]);
            writeRouChangedBeforeAfterFitnessToFile(beforeAfterChangedFitnessGen[1]);
        }


//        seletedIndsKneePointEveryGen.add(seletedIndsKneePoint);
        caseNumEveryGen.add(caseNum);
        indsChangedFunction[0] = changedIndIndexPop0.size();
        indsChangedFunction[1] = changedIndIndexPop1.size();
        indsChangedFunctionEveryGen.add(indsChangedFunction);


        // SHOULD WE QUIT?
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }
        // SHOULD WE QUIT?
        if (generation == numGenerations - 1) {
            writeDiversityToFile(entropyDiversity);
            writeKneePointSelectedIndsToFile(seletedIndsKneePointEveryGen);
            writeFunctionChangedIndsToFile(indsChangedFunctionEveryGen);
            writeDifferentCasesToFile(caseNumEveryGen);
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

        population = breeder.breedPopulation(this); //!!!!!!   return newpop;  if it is NSGA-II, the population here is 2N

        // POST-BREEDING EXCHANGING
        statistics.postBreedingStatistics(this);   //position 1  here, a new pop has been generated.

        // POST-BREEDING EXCHANGING
        statistics.prePostBreedingExchangeStatistics(this);
        population = exchanger.postBreedingExchangePopulation(this);   /** Simply returns state.population. */
        statistics.postPostBreedingExchangeStatistics(this);  //position 2

        // Generate new instances if needed
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
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

    //step3: based on different cases, replace functions
    private void enhanceCorrelation(GPNode rootNode, double[] score, int subpop) {

        //If this subtree only has terminal, we ignore it temporarily (maybe we can improve it later!!!) // or converse the function??
        if (rootNode.numNodes(GPNode.NODESEARCH_NONTERMINALS) == 1)
            return;

        GPInitializer initializer = ((GPInitializer) this.initializer);
        GPTreeConstraints constraints = ((GPIndividual) this.population.subpops[0].individuals[0]).trees[0].constraints(initializer);
        GPFunctionSet set = constraints.functionset;
        GPNode[] nodes = new GPNode[set.nodes[0].length];
        for (int i = 0; i < set.nodes[0].length; i++) {
            nodes[i] = set.nodes[0][i].lightClone();
        }


//        for (int n=0; n<nodes.length; n++){
//            nodes[n].children = new GPNode[2];
//            nodes[n].parent = null;
//        }

        //need to discuss different cases

        //calculate the index of the right subtree on level1;
        int rightNonTerminalIndex = rootNode.children[0].numNodes(GPNode.NODESEARCH_NONTERMINALS) + 1; //the value is left subtree's nonterminal + 1

        //if the left node is a terminal
        if (rootNode.children[0].numNodes(GPNode.NODESEARCH_NONTERMINALS) == 0) {
            if (firstChange)
                caseNum[subpop][4]++;

            if (score[1] >= 0) {    //if the right node is positive
                if (rootNode.toStringForHumans().equals(nodes[2].toStringForHumans())) {   //if rootnode is -, we change it to +;
                    rootNode.replaceWith(nodes[1]);
                    changed = true;
                } else if (rootNode.toStringForHumans().equals(nodes[4].toStringForHumans())) { //if rootnode is /, we change it to *;
                    rootNode.replaceWith(nodes[3]);
                    changed = true;
                }
            } else if (score[1] < 0) {    //if the right node is negative
                if (rootNode.toStringForHumans().equals(nodes[1].toStringForHumans())) {   //if rootnode is +, we change it to -;
                    rootNode.replaceWith(nodes[2]);
                    changed = true;
                } else if (rootNode.toStringForHumans().equals(nodes[3].toStringForHumans())) { //if rootnode is *, we change it to /;
                    rootNode.replaceWith(nodes[4]);
                    changed = true;
                }
            }
        }
        //if the right node is a terminal
        else if (rootNode.children[1].numNodes(GPNode.NODESEARCH_NONTERMINALS) == 0) {
            if (firstChange)
                caseNum[subpop][5]++;
            //if the left node is positive, we do nothing.
            if (score[0] <= 0) {    //if the right node is negative
                firstChange = false;
                enhanceCorrelation(rootNode.children[0], Arrays.copyOfRange(score, 1, rightNonTerminalIndex), subpop);
                changed = true;
            }
        }

        //if both left and right subtrees are positive  (P&P)
        else if (score[1] >= 0 && score[rightNonTerminalIndex] >= 0) {
            if (firstChange)
                caseNum[subpop][0]++;
            if (rootNode.toStringForHumans().equals(nodes[2].toStringForHumans())) {   //if rootnode is -, we change it to +;
                rootNode.replaceWith(nodes[1]);
                changed = true;
            } else if (rootNode.toStringForHumans().equals(nodes[4].toStringForHumans())) { //if rootnode is /, we change it to *;
                rootNode.replaceWith(nodes[3]);
                changed = true;
            }
        }

        //if left subtree is positive and the right subtree is negative (P&N)
        else if (score[1] >= 0 && score[rightNonTerminalIndex] <= 0) {
            if (firstChange)
                caseNum[subpop][1]++;
            if (rootNode.toStringForHumans().equals(nodes[1].toStringForHumans())) {   //if rootnode is +, we change it to -;
                rootNode.replaceWith(nodes[2]);
                changed = true;
            } else if (rootNode.toStringForHumans().equals(nodes[3].toStringForHumans())) { //if rootnode is *, we change it to /;
                rootNode.replaceWith(nodes[4]);
                changed = true;
            }
        }

        //if left subtree is negative and the right subtree is positive (N&P)      --> (P&P)
        else if (score[1] <= 0 && score[rightNonTerminalIndex] >= 0) {
            if (firstChange)
                caseNum[subpop][2]++;
            if (rootNode.toStringForHumans().equals(nodes[2].toStringForHumans())) {   //if rootnode is -, we change it to +;
                rootNode.replaceWith(nodes[1]);
                changed = true;
            } else if (rootNode.toStringForHumans().equals(nodes[4].toStringForHumans())) { //if rootnode is /, we change it to *;
                rootNode.replaceWith(nodes[3]);
                changed = true;
            }
            firstChange = false;
            enhanceCorrelation(rootNode.children[0], Arrays.copyOfRange(score, 1, rightNonTerminalIndex), subpop);
        }

        // if neither left and right subtrees are negative (N&N)           --> (P&N)
        else if (score[1] <= 0 && score[rightNonTerminalIndex] <= 0) {
            if (firstChange)
                caseNum[subpop][3]++;
            if (rootNode.toStringForHumans().equals(nodes[1].toStringForHumans())) {   //if rootnode is +, we change it to -;
                rootNode.replaceWith(nodes[2]);
                changed = true;
            } else if (rootNode.toStringForHumans().equals(nodes[3].toStringForHumans())) { //if rootnode is *, we change it to /;
                rootNode.replaceWith(nodes[4]);
                changed = true;
            }
            firstChange = false;
            enhanceCorrelation(rootNode.children[0], Arrays.copyOfRange(score, 1, rightNonTerminalIndex), subpop);
        }

    }


    //2021.4.16 fzhang save the diversity value to csv
    public void writeDiversityToFile(ArrayList<double[]> entropyDiversity) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".diversity.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,diversitySubpop0,diversitySubpop1");
            writer.newLine();
            for (int i = 0; i < entropyDiversity.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + entropyDiversity.get(i)[0] + ", " + entropyDiversity.get(i)[1] + "\n");
            }
            entropyDiversity.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeKneePointSelectedIndsToFile(ArrayList<int[]> selectedIndsKneePointEveryGen) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".kneePointSelectedInds.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,Subpop0,Subpop1");
            writer.newLine();
            for (int i = 0; i < selectedIndsKneePointEveryGen.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + selectedIndsKneePointEveryGen.get(i)[0] + ", " + selectedIndsKneePointEveryGen.get(i)[1] + "\n");
            }
            selectedIndsKneePointEveryGen.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeFunctionChangedIndsToFile(ArrayList<int[]> indsChangedFunctionEveryGen) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".functionChangedInds.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,Subpop0,Subpop1");
            writer.newLine();
            for (int i = 0; i < indsChangedFunctionEveryGen.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + indsChangedFunctionEveryGen.get(i)[0] + ", " + indsChangedFunctionEveryGen.get(i)[1] + "\n");
            }
            indsChangedFunctionEveryGen.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeDifferentCasesToFile(ArrayList<int[][]> caseNum) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".differentCaseNum.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,SeqCase0,SeqCase1,SeqCase2,SeqCase3,SeqCase4,SeqCase5,RouCase0,RouCase1,RouCase2,RouCase3,RouCase4,RouCase5");
            writer.newLine();
            for (int i = 0; i < caseNum.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + caseNum.get(i)[0][0] + ", " + caseNum.get(i)[0][1] + ", " + caseNum.get(i)[0][2] + ", " + caseNum.get(i)[0][3] + ", " + caseNum.get(i)[0][4] + ", " + caseNum.get(i)[0][5] +
                        ", " + caseNum.get(i)[1][0] + ", " + caseNum.get(i)[1][1] + ", " + caseNum.get(i)[1][2] + ", " + caseNum.get(i)[1][3] + ", " + caseNum.get(i)[1][4] + ", " + caseNum.get(i)[1][5] + "\n");
            }
            caseNum.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeSeqChangedBeforeAfterFitnessToFile(double[][] beforeAfterChangedFitnessGen) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + generation + ".seqChangedbeforeAfterFitness.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Ind,SeqBefore,SeqAfter");
            writer.newLine();
            for (int i = 0; i < beforeAfterChangedFitnessGen.length; i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + beforeAfterChangedFitnessGen[i][0] + ", " + beforeAfterChangedFitnessGen[i][1]  + "\n");
            }
//            Arrays.fill(beforeAfterChangedFitnessGen,0);
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeRouChangedBeforeAfterFitnessToFile(double[][] beforeAfterChangedFitnessGen) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + "." + generation + ".rouChangedbeforeAfterFitness.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Ind,RouBefore,RouAfter");
            writer.newLine();
            for (int i = 0; i < beforeAfterChangedFitnessGen.length; i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + beforeAfterChangedFitnessGen[i][0] + ", " + beforeAfterChangedFitnessGen[i][1]  + "\n");
            }
//            Arrays.fill(beforeAfterChangedFitnessGen,0);
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}


