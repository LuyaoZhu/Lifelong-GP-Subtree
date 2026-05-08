/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package yimei.jss.algorithm.featureConstruction;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;
import yimei.jss.helper.PopulationUtils;

import java.util.ArrayList;

/*
 * Protect selected features from destroying by crossover
 * Created: 19:53 25/4/2023
 * By: Luyao Zhu
 */

public class protectedKozaNodeSelectorFrequency extends KozaNodeSelector {
//    public int useTime;
    @Override
    public GPNode pickNode(final EvolutionState s,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree) {
        //first find the subtree equal to the feature in archive
        ArrayList<GPTree> protectedSeqGPTree = new ArrayList<>();
        ArrayList<GPTree> protectedRouGPTree = new ArrayList<>();
        if (tree.child.depth() > 1) {
            nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                GPTree treeClone = (GPTree) tree.clone();
                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                GPTree baseTree = PopulationUtils.GPNodetoGPTree(node);

                //compare with the features in the archive
                if (subpopulation == 0) {
                    for (int gen = 0; gen < ((FrequencyGPRuleEvolutionState)s).seqFeatureArchive.size(); gen++) {
                        for (int f = 0; f < ((FrequencyGPRuleEvolutionState)s).seqFeatureArchive.get(gen).size(); f++) {
                            GPTree countTTree = (GPTree) ((FrequencyGPRuleEvolutionState)s).seqFeatureArchive.get(gen).get(f).clone();
                            if (baseTree.child.makeLispTree().equals(countTTree.child.makeLispTree()))
                                protectedSeqGPTree.add(baseTree);
                        }
                    }
                } else {
                    for (int gen = 0; gen < ((FrequencyGPRuleEvolutionState)s).rouFeatureArchive.size(); gen++) {
                        for (int f = 0; f < ((FrequencyGPRuleEvolutionState)s).rouFeatureArchive.get(gen).size(); f++) {
                            GPTree countTTree = (GPTree) ((FrequencyGPRuleEvolutionState)s).rouFeatureArchive.get(gen).get(f).clone();
                            if (baseTree.child.makeLispTree().equals(countTTree.child.makeLispTree()))
                                protectedRouGPTree.add(baseTree);
                        }
                    }
                }
            }
        }

        double rnd = s.random[thread].nextDouble(); //probability  (0,1)

        if (rnd > nonterminalProbability + terminalProbability + rootProbability)  // pick anyone  nonterminalProbability = 0.9  terminalProbability=0.1
        // rnd > 0.9+0.1
        //nonterminalProbability + terminalProbability + rootProbability = 1  this will not happen
        {
            if (nodes == -1)
                nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
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

            ArrayList<Integer> protectedNodeNum = new ArrayList<>();
            if (subpopulation == 0) {
                for (int a = 0; a < tree.child.numNodes(GPNode.NODESEARCH_TERMINALS); a++) {
                    GPNode baseNode = tree.child.nodeInPosition(a, GPNode.NODESEARCH_TERMINALS);
                    for (int b = 0; b < protectedSeqGPTree.size(); b++) {
                        if (protectedSeqGPTree.get(b).child.contains(baseNode)) {
                            int depth = protectedSeqGPTree.get(b).child.depth();
                            GPNode tempNode = (GPNode) baseNode.clone();
                            for (int d = 0; d < depth - 1; d++) {
                                tempNode = (GPNode) tempNode.parent;
                            }
                            if (tempNode.makeLispTree().equals(protectedSeqGPTree.get(b).child.makeLispTree())) {
                                protectedNodeNum.add(a);
                                break;
                            }
                        }
                    }
                }
            } else {
                for (int a = 0; a < tree.child.numNodes(GPNode.NODESEARCH_TERMINALS); a++) {
                    GPNode baseNode = tree.child.nodeInPosition(a, GPNode.NODESEARCH_TERMINALS);
                    for (int b = 0; b < protectedRouGPTree.size(); b++) {
                        if (protectedRouGPTree.get(b).child.contains(baseNode)) {
                            int depth = protectedRouGPTree.get(b).child.depth();
                            GPNode tempNode = (GPNode) baseNode.clone();
                            for (int d = 0; d < depth - 1; d++) {
                                tempNode = (GPNode) tempNode.parent;
                            }
                            if (tempNode.makeLispTree().equals(protectedRouGPTree.get(b).child.makeLispTree())) {
                                protectedNodeNum.add(a);
                                break;
                            }
                        }
                    }
                }
            }
            int pickedNum = s.random[thread].nextInt(terminals);
            for (int a = 0; a < protectedNodeNum.size(); a++) {
                int temp = 0;
                while (protectedNodeNum.get(a) == pickedNum) {
                    pickedNum = s.random[thread].nextInt(nonterminals);
                    temp = 1;
                }
                if (temp == 1)
                    a = -1;
//                    if(pickedNum == protectedNodeNum.get(a))
//                        break;
            }

            //record useTime
//            if (protectedNodeNum.size() > 0)
//                useTime++;

            GPNode pickedNode = tree.child.nodeInPosition(pickedNum, GPNode.NODESEARCH_TERMINALS);
            return pickedNode;
            //choose the terminals
        } else  // pick nonterminals if you can
        {
            if (nonterminals == -1) nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
            //the number of non-terminals
            if (nonterminals > 0) // there are some nonterminals
            {

                ArrayList<Integer> protectedNodeNum = new ArrayList<>();
                int begin = 0; //record the orders of the protected trees
                int end = 0;//record the orders of the protected trees
                if (subpopulation == 0) {
                    for (int a = 0; a < protectedSeqGPTree.size(); a++) {
                        GPNode baseNode = (GPNode) protectedSeqGPTree.get(a).child.clone();
                        for (int b = 0; b < tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS); b++) {
                            GPNode countNode = tree.child.nodeInPosition(b, GPNode.NODESEARCH_NONTERMINALS);
                            if (countNode.makeLispTree().equals(baseNode.makeLispTree()))
                                begin = b;
                        }
                        end = begin + baseNode.numNodes(GPNode.NODESEARCH_NONTERMINALS) - 1;
                        for (int i = begin + 1; i < end + 1; i++) {
                            protectedNodeNum.add(i);
                        }
                    }
                } else {
                    for (int a = 0; a < protectedRouGPTree.size(); a++) {
                        GPNode baseNode = (GPNode) protectedRouGPTree.get(a).child.clone();
                        for (int b = 0; b < tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS); b++) {
                            GPNode countNode = tree.child.nodeInPosition(b, GPNode.NODESEARCH_NONTERMINALS);
                            if (countNode.makeLispTree().equals(baseNode.makeLispTree()))
                                begin = b;
                        }
                        end = begin + baseNode.numNodes(GPNode.NODESEARCH_NONTERMINALS) - 1;
                        for (int i = begin + 1; i < end + 1; i++) {
                            protectedNodeNum.add(i);
                        }
                    }
                }
                int pickedNum = s.random[thread].nextInt(nonterminals);
                for (int a = 0; a < protectedNodeNum.size(); a++) {
                    int temp = 0;
                    while (protectedNodeNum.get(a) == pickedNum) {
                        pickedNum = s.random[thread].nextInt(nonterminals);
                        temp = 1;
                    }
                    if (temp == 1)
                        a = -1;
//                    if(pickedNum == protectedNodeNum.get(a))
//                        break;
                }
//                if (protectedNodeNum.size() > 0)
//                    useTime++;
                GPNode pickedNode = tree.child.nodeInPosition(pickedNum, GPNode.NODESEARCH_NONTERMINALS);
                return pickedNode;
                //choose nodes
            } else // there ARE no nonterminals!  It must be the root node
            {
                return tree.child;
            }
        }
    }


}
