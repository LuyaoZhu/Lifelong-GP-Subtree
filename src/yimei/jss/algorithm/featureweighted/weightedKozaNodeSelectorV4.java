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

import java.util.List;

//fzhang 2019.9.8 choose mutation point

public class weightedKozaNodeSelectorV4 extends KozaNodeSelector {
    //=======================================================================
    public static final String P_PRE_GENERATIONS = "pre-generations";
    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;

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
            if (nodes == -1) nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
            //including the terminals    all the possible positions
            {
                return tree.child.nodeInPosition(s.random[thread].nextInt(nodes), GPNode.NODESEARCH_ALL);
                //randomly choose a node
            }
        } else if (rnd > nonterminalProbability + terminalProbability)  // pick the root
        {
            return tree.child; //for example: nonterminalProbability = 0.8 terminalProbability = 0.1, of rnd = 0.9. will choose the root
        } else if (rnd > nonterminalProbability)  // pick terminals  //nonterminalProbability = 0.9
        {
            if (terminals == -1) terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);
            return tree.child.nodeInPosition(s.random[thread].nextInt(terminals), GPNode.NODESEARCH_TERMINALS);
            //choose the terminals
        } else  // pick nonterminals if you can
        {
            if (nonterminals == -1) nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            if (nonterminals > 0) // there are some nonterminals
            {
                if(nonterminals == 1){
                    return tree.child.nodeInPosition(0, GPNode.NODESEARCH_NONTERMINALS);
                }
                //==========================================start=======================================================
                if (s.generation < preGenerations) {
                    return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                } else {
                    double[] score = new double[nonterminals];
                    AbstractRule rule = null;
                    //int[] validCases = new int[nonterminals];

                    if (subpopulation == 0) {
                        int decisionSize = decisionSituationsSequencing.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        //for (int i = 0; i < 1; i++) {
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
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                if(!correlation.isNaN()){
                                    //validCases[count] ++;
//                                    score[count] += Math.abs(new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]));
                                    score[count] += new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                }else
                                    {
                                     score[count] += 0.5; // not sure, it is related or not, so give a middle value
                                }
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
                                for (int count = 0; count < score_matrix.length; count++) {
                                    Double correlation = new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                    if(!correlation.isNaN()){
                                        //validCases[count] ++;
                                        score[count] += new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]);
//                                        score[count] += Math.abs(new PearsonsCorrelation().correlation(score_matrix[count], score_matrix[0]));
                                    }else
                                    {
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
                    for (int neg = 0; neg < score.length; neg++) {
                        if (flag){
                            score[neg] = 1/Math.abs(score[neg]);
                        }
                        else{
                            score[neg] = Math.abs(score[neg]);
                        }

                    }

//                    System.out.println("==============================");
//                    for(int a = 0; a < score.length; a++){
//                        System.out.println("score " + score[a]);
//                    }
                    RandomChoice.organizeDistribution(score);
                    int index = RandomChoice.pickFromDistribution(score, s.random[0].nextDouble());
                    return tree.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);

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
