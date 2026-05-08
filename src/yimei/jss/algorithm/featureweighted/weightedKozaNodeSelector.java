package yimei.jss.algorithm.featureweighted;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;
import ec.util.Parameter;
import yimei.jss.gp.terminal.TerminalERCUniform;

import java.util.ArrayList;
import java.util.stream.DoubleStream;

public class weightedKozaNodeSelector extends KozaNodeSelector {
    public static final String P_PRE_GENERATIONS = "pre-generations";

    public GPNode pickNode(final EvolutionState s,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree)
    {
        double rnd = s.random[thread].nextDouble(); //probability  (0,1)

        int preGenerations = s.parameters.getIntWithDefault(
                new Parameter(P_PRE_GENERATIONS), null, -1);  //50

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
                    GPNode node = ((GPIndividual)(ind)).trees[0].child;
                    weightsNonterminals(s, node, subpopulation, score);
//                    score[count] = scoreOfFunction;
//                    scoreOfFunction = 0.0;
                }
                int bestP = 0;
                double smallest = Double.MAX_VALUE;
                for(int p = 1; p < score.length; p++){
                    if(score[p] < smallest){
                        smallest = score[p];
                        bestP = p;
                    }
                }
               return tree.child.nodeInPosition(bestP, GPNode.NODESEARCH_NONTERMINALS);
                //======================================end=============================================================
            }
            else // there ARE no nonterminals!  It must be the root node
            {
                return tree.child;
            }
        }
    }

    static int count = 0;
    static double scoreOfFunction = 0.0;
    static int number = 0;
    static void weightsNonterminals(final EvolutionState s, GPNode node, int subpopulation, double[] score) {
        if (node == null) {
            return;
        }

        if (node.children == null || node.children.length == 0) {  //1. a node does not have child is a terminal
            //2. the length of node's child is 0 (empty array)---it is a terminal;

            GPNode[][] terminals = ((FreGPRuleEvolutionState)s).getTerminals();
            String terminalName = ((TerminalERCUniform)node).getTerminal().name();

            int idxTerminal = 0;
            for(int i = 0; i < terminals[0].length; i++){
                if(terminals[0][i].toString().equals(terminalName)){
                    idxTerminal = i;
                    break;
                }
            }

            ArrayList<double[]> weightsValue = ((FreGPRuleEvolutionState) s).getWeights();
            double[] weightOfTerminals = weightsValue.get(subpopulation); //all the weights of all the terminals
            double sum = DoubleStream.of(weightOfTerminals).sum(); //not necessary

            double weightOfTerminalName = weightOfTerminals[idxTerminal];
            scoreOfFunction += weightOfTerminalName/sum;
            number++;

            if (number == 2){
                score[count] = scoreOfFunction;
                scoreOfFunction = 0.0;
                number = 0;
            }
            return;
        }

        for(GPNode child : node.children) { //repeat to check the terminals
            if(child.children != null & child.children.length!= 0){
                count++;
            }
            weightsNonterminals(s, child, subpopulation, score);
        }
    }
}
