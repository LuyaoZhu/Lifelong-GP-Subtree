package yimei.jss.algorithm.featureweighted;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;
import ec.util.Parameter;
import ec.util.RandomChoice;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.scoreMultiPopCoevolutionaryClearingEvaluator;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

//fzhang 2019.9.8 choose crossover point and try to label the function node --- not finished.
public class weightedKozaNodeSelectorV6 extends KozaNodeSelector {
    //=======================================================================
    public static final String P_PRE_GENERATIONS = "pre-generations";
    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;
    static double chosenScore = 0.0;
    static Boolean parent1Func = false; //whether the parent 1 choose a function node

    //fzhang 2019.10.22 record the probability difference
    static ArrayList<Double> proDifference1 = new ArrayList<>(); //flag is true, good individual, choose the bad point
    static ArrayList<Double> proDifference2 = new ArrayList<>(); //flag is false, bad individual, choose the good point
    //=======================================================================

    public GPNode pickNode(final EvolutionState s,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree,
                           Boolean flag) {
        double rnd = s.random[thread].nextDouble(); //probability  (0,1)
        //====================================================================
        int preGenerations = s.parameters.getIntWithDefault(
                new Parameter(P_PRE_GENERATIONS), null, -1);  //50
        PhenoCharacterisation[] pc = scoreMultiPopCoevolutionaryClearingEvaluator.getPhenoCharacterisation();

        if (subpopulation == 0) {
            decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pc[subpopulation]).decisionSituations;
        } else {
            decisionSituationsRouting = ((RoutingPhenoCharacterisation) pc[subpopulation]).decisionSituations;
        }

        //====================================================================

        if (rnd > nonterminalProbability + terminalProbability + rootProbability)  // pick anyone  nonterminalProbability = 0.9  terminalProbability=0.1
        // rnd > 0.9+0.1
        //nonterminalProbability + terminalProbability + rootProbability = 1  this will not happen
        {
            parent1Func = false;
            if (nodes == -1) nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
            //including the terminals    all the possible positions
            {
                return tree.child.nodeInPosition(s.random[thread].nextInt(nodes), GPNode.NODESEARCH_ALL);
                //randomly choose a node
            }
        } else if (rnd > nonterminalProbability + terminalProbability)  // pick the root
        {
            parent1Func = false;
            return tree.child; //for example: nonterminalProbability = 0.8 terminalProbability = 0.1, of rnd = 0.9. will choose the root
        } else if (rnd > nonterminalProbability)  // pick terminals  //nonterminalProbability = 0.9
        {
            parent1Func = false;
            if (terminals == -1) terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);
            return tree.child.nodeInPosition(s.random[thread].nextInt(terminals), GPNode.NODESEARCH_TERMINALS);
            //choose the terminals
        } else  // pick nonterminals if you can
        {
            if (nonterminals == -1) nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            if (nonterminals > 0) // there are some nonterminals
            {
                double[] score = new double[nonterminals];
                //fzhang 2019.10.23 for save the index information of ndoes
                ArrayList<HashMap<Integer, Integer[]>> nodeIndex = new ArrayList<HashMap<Integer, Integer[]>>();

                if (nonterminals == 1) {
                    parent1Func = true;
                    chosenScore = 1;
                    return tree.child.nodeInPosition(0, GPNode.NODESEARCH_NONTERMINALS);
                }
                //==========================================start=======================================================
                if (s.generation < preGenerations) {
                    return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                } else {
                    AbstractRule rule = null;
                    if (subpopulation == 0) {

                        //2019.10.23 give a label to each funciton node
                        //===============================start===========================================
                        GPTree treeCloneIndex = (GPTree) tree.clone();
                        HashMap<GPNode, Integer> treeWithIndex = treeIndex(treeCloneIndex, nonterminals);
                        int[] nodeLabel = labelNodes(treeCloneIndex, treeWithIndex, nonterminals);
                        //===============================end=============================================

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

                                    //fzhang 2019.10.23 if the node has a negative label, convert the priority value to -value (the larger the value, the better)
                                    //if the node has a positive value, the smaller the better, use the original value
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority() * nodeLabel[numSubtree];
                                }
                            }
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                if (!correlation.isNaN()) {
                                    score[count] += new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                } else {
                                    score[count] += 0.5; // not sure, it is related or not, so give a middle value
                                }
                            }
                        }

                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsSequencing.size();
                        }

                    } else {
                        int decisionSize = decisionSituationsRouting.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        //for (int i = 0; i < 1; i++) {
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
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                if (!correlation.isNaN()) {
                                    score[count] += new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                } else {
                                    score[count] += 0.5;
                                }
                            }
                        }
                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsRouting.size();
                        }
                    }

                    //only for mutation 2019.9.5
                    //ignore the positive and negative correlation first---the large the score, the important the subtree
                    double[] scoreClone = score.clone();
                    int countZeroScore = 0;
                    if (flag == true) {//this is parent 1
                        parent1Func = true;
                        for (int neg = 0; neg < scoreClone.length; neg++) {
                            scoreClone[neg] = 1 / Math.abs(scoreClone[neg]); //choose the small score, no matter it is positive or negative
                        }
                    } else {
                        if (parent1Func == true) {
                            if (chosenScore > 0) {//the selected point > 0, choose the positive one to do crossover
                                for (int neg = 0; neg < scoreClone.length; neg++) {
                                    if (scoreClone[neg] <= 0) {
                                        scoreClone[neg] = 0;//that mean if chosen "+", the "-" ones do not have chance to be chosen
                                        countZeroScore++;
                                    }
                                }
                            }

                            if (chosenScore < 0) { //the selected point < 0, choose the negative one to do crossover
                                for (int neg = 0; neg < scoreClone.length; neg++) {
                                    if (scoreClone[neg] >= 0) {
                                        scoreClone[neg] = 0;
                                        countZeroScore++;
                                    } else {
                                        scoreClone[neg] = Math.abs(scoreClone[neg]);
                                    }
                                }
                            }

                            if (chosenScore == 0) { //teh selected point equals 0, it does not matter to choose which one
                                for (int neg = 0; neg < scoreClone.length; neg++) {
                                    scoreClone[neg] = Math.abs(scoreClone[neg]);
                                }
                            }
                        } else //p1 does not choose functions, when p1 choose functions, just random choose one.
                        {
                            return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                        }
                    }

                    //if(countZeroScore == scoreClone.length || parent1Func == false){ //have problem
                    if (countZeroScore == scoreClone.length) { //have problem
                        return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                    } else {
                        RandomChoice.organizeDistribution(scoreClone);
                        int index = RandomChoice.pickFromDistribution(scoreClone, s.random[0].nextDouble());
                        //2019.10.22 fzhang for analyses, calculate the probability difference: the probability that the node is chosen - 1/nonterminals (uniformly chosen probability)
                        //==================start===========================
                        if (s.generation == 1 || s.generation == 25 || s.generation == 45) {
                            if (nonterminals != 1) { //if nonterminals == 1, then the only one function will be selected. It does not make sense to record the probability difference---always 0
                                // and it will misleading the count the probability difference as 0, means the proposed algorithm have no leading for choosing the node point
                                if (flag == true) {
                                    proDifference1.add(score[index] - 1.0 / nonterminals); // should put 1.0 here, otherwise, they will return interger value, for example, 1/7=0
                                } else {
                                    proDifference2.add(score[index] - 1.0 / nonterminals);
                                }
                            }
                        }

                        if (flag == true) {
                            chosenScore = score[index];
                        }
                        return tree.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);
                    }
                    //======================================end=============================================================
                }
            } else // there ARE no nonterminals!  It must be the root node
            {
                return tree.child;
            }
        }
    }


    //fzhang 2019.10.23 get the index of each node. The goal is to label each node.
    //========================================================================================================================================================================
/*    static ArrayList<HashMap<Integer,Integer[]>> cal_child_index(GPNode node, double index, ArrayList<HashMap<Integer,Integer[]>> nodeIndex, Boolean flag) {
        if (node == null) {
            return null;
        }

        if (node.children == null || node.children.length == 0) {  //1. a node does not have child is a terminal
            ArrayList<HashMap<Integer,Integer>>[] leafScoreIdx = new ArrayList<HashMap<Integer,Integer[]>>;
            leafScoreIdx[1] = index;
            return leafScoreIdx;
        }
        else{
            GPNode left_node = node.children[0];
            GPNode right_node = node.children[1];
            double[] leftScore_rightIdx = new double[2];
            double[] rightScore_nextIdx = new double[2];

            leftScore_rightIdx = cal_child_index(left_node, index+1, nodeIndex, flag); //index can be changed in the parameters
            rightScore_nextIdx = cal_child_index(right_node, leftScore_rightIdx[1], nodeIndex, flag); //index can be changed in the parameters


            ArrayList<HashMap<Integer,Integer[]>> toReturn = new ArrayList<HashMap<Integer,Integer[]>>{node_score, rightScore_nextIdx[1]};
            return toReturn;
        }
    }*/
    //=======================================================================================================================================================================

    HashMap<GPNode, Integer> treeIndex(GPTree tree, int nonterminals) {
        HashMap<GPNode, Integer> hm = new HashMap<GPNode, Integer>();
        GPNode node = null;
        for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
            node = tree.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
            hm.put(node, numSubtree);
        }
        return hm;
    }


    int[] labelNodes(GPTree tree, HashMap<GPNode, Integer> hm, int nonterminals) {
        int[] labeledNode = new int[nonterminals];
        GPNode node = null;

        //the root node always has positive effect for the whole tree --- the small the priority value, the better.  1: positive; -1: negative
        labeledNode[0] = 1;

        for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
            node = tree.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
            int leftIndex;
            int rightIndex;

            if (hm.get(node.children[0]) != null) {
                leftIndex = hm.get(node.children[0]);
                labeledNode[leftIndex] = labeledNode[numSubtree];
            }

            if (hm.get(node.children[1]) != null) {
                rightIndex = hm.get(node.children[1]);
                if (node.toString() == "+" || node.toString() == "*" || node.toString() == "Max" || node.toString() == "Min") {
                    labeledNode[rightIndex] = labeledNode[numSubtree];
                } else {
                    labeledNode[rightIndex] = -labeledNode[numSubtree];
                }
            }
        }

        return labeledNode;
    }
}
