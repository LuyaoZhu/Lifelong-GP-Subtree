package yimei.jss.algorithm.featureweighted;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;
import ec.util.Parameter;
import ec.util.RandomChoice;
import yimei.jss.gp.terminal.TerminalERCUniform;

import java.util.ArrayList;

public class weightedKozaNodeSelectorV2 extends KozaNodeSelector {
    //=======================================================================
    public static final String P_PRE_GENERATIONS = "pre-generations";
    ArrayList<double[]> weightsValue = null;
    double[] weightOfTerminals = null;
    //=======================================================================

    public GPNode pickNode(final EvolutionState s,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree,
                           Boolean flag)
    {
        double rnd = s.random[thread].nextDouble(); //probability  (0,1)

        //====================================================================
        int preGenerations = s.parameters.getIntWithDefault(
                new Parameter(P_PRE_GENERATIONS), null, -1);  //50
        GPNode[][] terminalSet = ((FreGPRuleEvolutionState)s).getTerminals();
        if (s.generation >= preGenerations){
            weightsValue = ((FreGPRuleEvolutionState) s).getWeights();
            if(subpopulation == 0){
                weightOfTerminals = weightsValue.get(weightsValue.size()-2); //all the weights of all the terminals
            }
            else{
                weightOfTerminals = weightsValue.get(weightsValue.size()-1);
            }

        }
        //====================================================================

        if (rnd > nonterminalProbability + terminalProbability + rootProbability)  // pick anyone  nonterminalProbability = 0.9  terminalProbability=0.1
        // rnd > 0.9+0.1
        //nonterminalProbability + terminalProbability + rootProbability = 1  this will not happen
        {
            if (nodes==-1) nodes=tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
            //including the terminals    all the possible positions
            {
                return tree.child.nodeInPosition(s.random[thread].nextInt(nodes), GPNode.NODESEARCH_ALL);
                //randomly choose a node
            }
        }
        else if (rnd > nonterminalProbability + terminalProbability)  // pick the root
        {
            return tree.child; //for example: nonterminalProbability = 0.8 terminalProbability = 0.1, of rnd = 0.9. will choose the root
        }
        else if (rnd > nonterminalProbability)  // pick terminals  //nonterminalProbability = 0.9
        {
            if (terminals==-1) terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);
            return tree.child.nodeInPosition(s.random[thread].nextInt(terminals), GPNode.NODESEARCH_TERMINALS);
            //choose the terminals
        }
        else  // pick nonterminals if you can
        {
            if (nonterminals==-1) nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            double[] score = new double[nonterminals]; //save the score (the importance of subtree)
            //the number of non-terminals
            if (nonterminals > 0) // there are some nonterminals
            {
                //==========================================start=======================================================
                if (s.generation < preGenerations){
                    return tree.child.nodeInPosition(s.random[thread].nextInt(nonterminals), GPNode.NODESEARCH_NONTERMINALS);
                }
                else{
                    double index = 0;
                    GPNode node = ind.trees[0].child;
                    cal_score(node, index, score, terminalSet, weightOfTerminals, flag);
                }

                RandomChoice.organizeDistribution(score);
                int index = RandomChoice.pickFromDistribution(score, s.random[0].nextDouble());
                return tree.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);

                //======================================end=============================================================
            }
            else // there ARE no nonterminals!  It must be the root node
            {
                return tree.child;
            }
        }
    }


//========================================================================================================================================================================
    static double[] cal_score(GPNode node, double index, double[] score, GPNode[][] terminalSet, double [] weightOfTerminals, Boolean flag) {
        if (node == null) {
            return null;
        }

        if (node.children == null || node.children.length == 0) {  //1. a node does not have child is a terminal
            //2. the length of node's child is 0 (empty array)---it is a terminal;
            String terminalName = ((TerminalERCUniform)node).getTerminal().name();
            int idxTerminal = 0;
            for(int i = 0; i < terminalSet[0].length; i++){
                if(terminalSet[0][i].toString().equals(terminalName)){
                    idxTerminal = i;
                    break;
                }
            }
            double weightOfTerminalName = weightOfTerminals[idxTerminal];

            double[] leafScoreIdx = new double[2];
            //there are some terminals whose weights are zero, then there will be some error about the probability distribution, so here give them a very small numbers
            if (flag == true){ //it is a good individual with smaller fitness---choose the part with small values (means not too much "important features" appear here)
                if(weightOfTerminalName == 0.0){
                    weightOfTerminalName = 1;
                }
                leafScoreIdx[0] = 1/weightOfTerminalName;
            }
            else{
                { //it is a bad individual with largest fitness---choose the part with larger values to give others (means a lot "important features" appear here)
                    if(weightOfTerminalName == 0.0){
                        weightOfTerminalName = 0.001;
                    }
                    leafScoreIdx[0] = weightOfTerminalName;
                }
            }
            leafScoreIdx[1] = index;
            return leafScoreIdx;
        }
        else{
            GPNode left_node = node.children[0];
            GPNode right_node = node.children[1];
            double[] leftScore_rightIdx = new double[2];
            double[] rightScore_nextIdx = new double[2];

            leftScore_rightIdx = cal_score(left_node, index+1, score, terminalSet, weightOfTerminals, flag); //index can be changed in the parameters
            rightScore_nextIdx = cal_score(right_node, leftScore_rightIdx[1], score, terminalSet, weightOfTerminals, flag); //index can be changed in the parameters

            double node_score = (leftScore_rightIdx[0] + rightScore_nextIdx[0])/2;
            score[(int)(index)] = node_score;
            double[] toReturn = new double[]{node_score, rightScore_nextIdx[1]};
            return toReturn;
        }
    }
    //=======================================================================================================================================================================
}
