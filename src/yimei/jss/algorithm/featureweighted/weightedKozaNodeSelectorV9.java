package yimei.jss.algorithm.featureweighted;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;
import ec.util.Parameter;
import ec.util.RandomChoice;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
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
import java.util.Arrays;
import java.util.List;

//fzhang 2019.9.8 choose crossover point
public class weightedKozaNodeSelectorV9 extends KozaNodeSelector {
    //=======================================================================
    public static final String P_PRE_GENERATIONS = "pre-generations";
    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;
    static double chosenScore = 0.0;
    //static Boolean parent1Func = false; //whether the parent 1 choose a function node

    //fzhang 2019.10.22 record the probability difference
    static ArrayList<Double> proDifference1 = new ArrayList<>(); //flag is true, good individual, choose the bad point
    static ArrayList<Double> proDifference2 = new ArrayList<>(); //flag is false, bad individual, choose the good point
    static ArrayList<Double> coefficient1 = new ArrayList<>(); //flag is true, good individual, choose the bad point
    static ArrayList<Double> coefficient2 = new ArrayList<>(); //flag is false, bad individual, choose the good point
    static ArrayList<Double> depthRatio1 = new ArrayList<>(); //flag is true, good individual, choose the bad point
    static ArrayList<Double> depthRatio2 = new ArrayList<>();
    static int indexForParent2 = -1;
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
            //parent1Func = false;
            if (nodes == -1) nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
            //including the terminals    all the possible positions
            {
                return tree.child.nodeInPosition(s.random[thread].nextInt(nodes), GPNode.NODESEARCH_ALL);
                //randomly choose a node
            }
        } else if (rnd > nonterminalProbability + terminalProbability)  // pick the root
        {
            //parent1Func = false;
            return tree.child; //for example: nonterminalProbability = 0.8 terminalProbability = 0.1, of rnd = 0.9. will choose the root
        } else if (rnd > nonterminalProbability)  // pick terminals  //nonterminalProbability = 0.9
        {
            //parent1Func = false;
            if (terminals == -1) terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);
            return tree.child.nodeInPosition(s.random[thread].nextInt(terminals), GPNode.NODESEARCH_TERMINALS);
            //choose the terminals
        } else  // pick nonterminals if you can
        {
             if (nonterminals == -1) nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            if (nonterminals > 0) // there are some nonterminals
            {
                double[] score = new double[nonterminals];
                if(nonterminals == 1){
                   // parent1Func = true;
                    chosenScore = 1;
                    return tree.child.nodeInPosition(0, GPNode.NODESEARCH_NONTERMINALS);
                }
                //==========================================start=======================================================
                if (s.generation < preGenerations) {
                    return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                } else {
                    AbstractRule rule = null;
                    if (subpopulation == 0) {
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
                            double rank =1;
                            for(int numNodes = 0; numNodes < score_matrix.length; numNodes++){
                                double [] ranks = new double[score_matrix[numNodes].length];
                                for(int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++){
                                    for(int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++){
                                        if(score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]){
                                            rank++;
                                    }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for(int r = ranks.length-1; r >= 0; r--){
                                    for(int j = 0; j < r; j++){
                                        if (ranks[j] == ranks[r]){
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

                    } else
                        {
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

                                double rank =1;
                                for(int numNodes = 0; numNodes < score_matrix.length; numNodes++){
                                    double [] ranks = new double[score_matrix[numNodes].length];
                                    for(int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++){
                                        for(int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++){
                                            if(score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]){
                                                rank++;
                                            }
                                        }
                                        ranks[queueSize1] = rank;
                                        rank = 1;
                                    }

                                    int count = 0;
                                    for(int r = ranks.length-1; r >= 0; r--){
                                        for(int j = 0; j < r; j++){
                                            if (ranks[j] == ranks[r]){
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
                        }

                    //only for mutation 2019.9.5
                    //ignore the positive and negative correlation first---the large the score, the important the subtree
                    double[] scoreClone = score.clone();
                    double[] scoreCloneforParent2  = score.clone();
                    int countZeroScore = 0;
                    if (flag == true) {//this is parent 1
                        //parent1Func = true;
                        for (int neg = 0; neg < scoreClone.length; neg++) {
                            /*if(scoreClone[neg] == 0){
                                scoreClone[neg] = 0.01; //fix the division problem.  when the value is 0, not important, has higher change to be selected. Set it a smaller value.
                            }
                            scoreClone[neg] = 1 / Math.abs(scoreClone[neg]); //choose the small score, no matter it is positive or negative*/
                            scoreClone[neg] = 1- Math.abs(scoreClone[neg]); //in this way. We also do not need to consider the division by zero problem.
                            scoreCloneforParent2[neg] = Math.abs(scoreCloneforParent2[neg]);
                        }
                    } else {
                            //do not consider positive or negative
                            for (int neg = 0; neg < scoreClone.length; neg++) {
                                scoreClone[neg] = Math.abs(scoreClone[neg]);
                                scoreCloneforParent2[neg] = 1 - Math.abs(scoreCloneforParent2[neg]);
                            }

                            /*if (chosenScore > 0) {//the selected point > 0, choose the positive one to do crossover
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
                                    }
                                    else{
                                        scoreClone[neg] = Math.abs(scoreClone[neg]);
                                    }
                                }
                            }

                            if(chosenScore == 0){ //teh selected point equals 0, it does not matter to choose which one
                                for (int neg = 0; neg < scoreClone.length; neg++) {
                                    scoreClone[neg] = Math.abs(scoreClone[neg]);
                                }
                            }*/
                         //p1 does not choose functions, when p1 choose functions, just random choose one.
                    }

                    //if(countZeroScore == scoreClone.length || parent1Func == false){ //have problem
                    if(countZeroScore == scoreClone.length || Arrays.stream(scoreClone).distinct().count() == 1 || Arrays.stream(scoreCloneforParent2).distinct().count() == 1){ //have problem
                        return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                    }
                    else {
                        RandomChoice.organizeDistribution(scoreClone);
                        int index = RandomChoice.pickFromDistribution(scoreClone, s.random[0].nextDouble());

                        RandomChoice.organizeDistribution(scoreCloneforParent2);
                        indexForParent2 = RandomChoice.pickFromDistribution(scoreCloneforParent2, s.random[0].nextDouble());

                        //2019.10.22 fzhang for analyses, calculate the probability difference: the probability that the node is chosen - 1/nonterminals (uniformly chosen probability)
                        //==================start===========================
                        if(s.generation == 1 || s.generation == 25 || s.generation == 45){
                            if(nonterminals != 1){ //if nonterminals == 1, then the only one function will be selected. It does not make sense to record the probability difference---always 0
                                // and it will misleading the count the probability difference as 0, means the proposed algorithm have no leading for choosing the node point
                                if(flag == true){
                                    proDifference1.add(scoreClone[index] - 1.0/nonterminals); // should put 1.0 here, otherwise, they will return interger value, for example, 1/7=0
                                    coefficient1.add(score[index]);
                                }
                                else{
                                    proDifference2.add(scoreClone[index] - 1.0/nonterminals);
                                    coefficient2.add(score[index]);
                                }
                            }
                        }

                        //2019.11.11 record the depth of subtree that was chosen
                        GPTree treeClone = (GPTree) tree.clone();
                        GPNode selectedNode = treeClone.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);
                        GPNode rootNode = treeClone.child.nodeInPosition(0, GPNode.NODESEARCH_NONTERMINALS);

                        double currentDepth  = selectedNode.atDepth();
                        double depthRatio = currentDepth / rootNode.depth();
                        if(flag == true){
                            depthRatio1.add(depthRatio);
                        }
                        else{
                            depthRatio2.add(depthRatio);
                        }

                      /*  if(flag == true){
                            chosenScore = score[index];
                        }*/
                        return tree.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);
                    }
                    //======================================end=============================================================
                }
            }
            else // there ARE no nonterminals!  It must be the root node
            {
                return tree.child;
            }
        }
    }
}
