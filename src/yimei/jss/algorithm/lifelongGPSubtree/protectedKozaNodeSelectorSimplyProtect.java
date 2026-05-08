/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.gp.koza.KozaNodeSelector;

import java.util.ArrayList;

/*
 * KozaNodeSelector.java
 *
 * Created: Tue Oct 12 17:21:28 1999
 * By: Sean Luke
 */

/**
 * KozaNodeSelector is a GPNodeSelector which picks nodes in trees a-la Koza I,
 * with the addition of having a probability of always picking the root.
 * The method divides the range 0.0...1.0 into four probability areas:
 *
 * <ul>
 * <li>One area specifies that the selector must pick a terminal.
 * <li>Another area specifies that the selector must pick a nonterminal (if there is one, else a terminal).
 * <li>The third area specifies that the selector pick the root node.
 * <li>The fourth area specifies that the selector pick any random node.
 * </ul>
 *
 * <p>The KozaNodeSelector chooses by probability between these four situations.
 * Then, based on the situation it has picked, it selects either a random
 * terminal, nonterminal, root, or arbitrary node from the tree and returns it.
 *
 * <p>As the selector picks a node, it builds up some statistics information
 * which makes it able to pick a little faster in subsequent passes.  Thus
 * if you want to reuse this selector on another tree, you need to call
 * reset() first.
 *
 *
 * <p><b>Parameters</b><br>
 * <table>
 * <tr><td valign=top><i>base</i>.<tt>terminals</tt><br>
 * <font size=-1>0.0 &lt;= double &lt;= 1.0,<br>
 * nonterminals + terminals + root <= 1.0</font></td>
 * <td valign=top>(the probability we must pick a terminal)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>nonterminals</tt><br>
 * <font size=-1>0.0 &lt;= double &lt;= 1.0,<br>
 * nonterminals + terminals + root <= 1.0</font></td>
 * <td valign=top>(the probability we must pick a nonterminal if possible)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>root</tt><br>
 * <font size=-1>0.0 &lt;= double &lt;= 1.0,<br>
 * nonterminals + terminals + root <= 1.0</font></td>
 * <td valign=top>(the probability we must pick the root)</td></tr>
 *
 * </table>
 *
 * <p><b>DefaultBase</b><br>
 * gp.koza.ns
 *
 * @author Sean Luke
 * @version 1.0
 */

public class protectedKozaNodeSelectorSimplyProtect extends KozaNodeSelector {

    @Override
    public GPNode pickNode(EvolutionState s, int subpopulation, int thread, GPIndividual ind, GPTree tree, Boolean flag) {
        return null;
    }

    @Override
    public GPNode pickNode(final EvolutionState s,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree) {

        double rnd = s.random[thread].nextDouble();

        int parentID = FeatureConstructionCrossoverPipelineSimplyProtect.parent1Flag ? 0 : 1;
        int treeID = FeatureConstructionCrossoverPipelineSimplyProtect.currentTreeID;

        // === CASE 1: pick ANY node ===
        if (rnd > nonterminalProbability + terminalProbability + rootProbability) {

            if (nodes == -1)
                nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL);

            return tree.child.nodeInPosition(
                    s.random[thread].nextInt(nodes),
                    GPNode.NODESEARCH_ALL
            );
        }

        // === CASE 2: pick ROOT ===
        else if (rnd > nonterminalProbability + terminalProbability) {
            return tree.child;
        }

        // === CASE 3: pick TERMINAL ===
        else if (rnd > nonterminalProbability) {

            if (terminals == -1)
                terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);

            ArrayList<Integer> protectedRanks =
                    FeatureConstructionCrossoverPipelineSimplyProtect
                            .protectedTerminalsRank[parentID][treeID];

            return pickWithProtection(
                    s, thread,
                    tree,
                    terminals,
                    protectedRanks,
                    GPNode.NODESEARCH_TERMINALS
            );
        }

        // === CASE 4: pick NONTERMINAL ===
        else {

            if (nonterminals == -1)
                nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);

            if (nonterminals == 0)
                return tree.child;

            ArrayList<Integer> protectedRanks =
                    FeatureConstructionCrossoverPipelineSimplyProtect
                            .protectedNonTerminalsRank[parentID][treeID];

            return pickWithProtection(
                    s, thread,
                    tree,
                    nonterminals,
                    protectedRanks,
                    GPNode.NODESEARCH_NONTERMINALS
            );
        }
    }

    private GPNode pickWithProtection(
            EvolutionState s,
            int thread,
            GPTree tree,
            int total,
            ArrayList<Integer> protectedRanks,
            int searchType
    ) {

        ArrayList<Integer> candidates = new ArrayList<>();

        for (int i = 0; i < total; i++) {
            if (!protectedRanks.contains(i)) {
                candidates.add(i);
            }
        }

        //record the ratio of building blocks on the trees
        

        int pickIndex;

        if (candidates.isEmpty()) {
            // fallback（避免死锁）
            pickIndex = s.random[thread].nextInt(total);
        } else {
            pickIndex = candidates.get(
                    s.random[thread].nextInt(candidates.size())
            );
        }

        return tree.child.nodeInPosition(pickIndex, searchType);
    }


}
