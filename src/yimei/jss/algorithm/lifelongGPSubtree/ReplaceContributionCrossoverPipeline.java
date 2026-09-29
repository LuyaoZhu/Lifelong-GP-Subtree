/*
  Copyright 2006 by Sean Luke
  Licensed under the Academic Free License version 3.0
  See the file "LICENSE" for more information
*/


package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPInitializer;
import ec.gp.GPNode;
import ec.gp.GPTree;
import yimei.jss.algorithm.multipletreegp.AllIndexCrossoverPipeline;

import java.util.*;

/*
 * CrossoverPipeline.java
 *
 * Created: Mon Aug 30 19:15:21 1999
 * By: Sean Luke
 */


/**
 * CrossoverPipeline is a GPBreedingPipeline which performs a strongly-typed
 * version of
 * Koza-style "Subtree Crossover".  Two individuals are selected,
 * then a single tree is chosen in each such that the two trees
 * have the same GPTreeConstraints.  Then a random node is chosen
 * in each tree such that each node's return type is type-compatible
 * with the argument type of the slot in the parent which contains
 * the other node.
 * If by swapping subtrees at these nodes the two trees will not
 * violate maximum depth constraints, then the trees perform the
 * swap, otherwise, they repeat the hunt for random nodes.
 *
 * <p>The pipeline tries at most <i>tries</i> times to a pair
 * of random nodes BOTH with valid swap constraints.  If it
 * cannot find any such pairs after <i>tries</i> times, it
 * uses the pair of its last attempt.  If either of the nodes in the pair
 * is valid, that node gets substituted with the other node.  Otherwise
 * an individual invalid node isn't changed at all (it's "reproduced").
 *
 * <p><b>Compatibility with constraints.</b>
 * Since Koza-I/II only tries 1 time, and then follows this policy, this is
 * compatible with Koza.  lil-gp either tries 1 time, or tries forever.
 * Either way, this is compatible with lil-gp.  My hacked
 * <a href="http://www.cs.umd.edu/users/seanl/gp/">lil-gp kernel</a>
 * either tries 1 time, <i>n</i> times, or forever.  This is compatible
 * as well.
 *
 * <p>This pipeline typically produces up to 2 new individuals (the two newly-
 * swapped individuals) per produce(...) call.  If the system only
 * needs a single individual, the pipeline will throw one of the
 * new individuals away.  The user can also have the pipeline always
 * throw away the second new individual instead of adding it to the population.
 * In this case, the pipeline will only typically
 * produce 1 new individual per produce(...) call.
 *
 * <p><b>Typical Number of Individuals Produced Per <tt>produce(...)</tt> call</b><br>
 * 2 * minimum typical number of individuals produced by each source, unless tossSecondParent
 * is set, in which case it's simply the minimum typical number.
 *
 * <p><b>Number of Sources</b><br>
 * 2
 *
 * <p><b>Parameters</b><br>
 * <table>
 * <tr><td valign=top><i>base</i>.<tt>tries</tt><br>
 * <font size=-1>int &gt;= 1</font></td>
 * <td valign=top>(number of times to try finding valid pairs of nodes)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>maxdepth</tt><br>
 * <font size=-1>int &gt;= 1</font></td>
 * <td valign=top>(maximum valid depth of a crossed-over subtree)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>maxsize</tt><br>
 * <font size=-1>int &gt;= 1</font></td>
 * <td valign=top>(maximum valid size, in nodes, of a crossed-over subtree)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>tree.0</tt><br>
 * <font size=-1>0 &lt; int &lt; (num trees in individuals), if exists</font></td>
 * <td valign=top>(first tree for the crossover; if parameter doesn't exist, tree is picked at random)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>tree.1</tt><br>
 * <font size=-1>0 &lt; int &lt; (num trees in individuals), if exists</font></td>
 * <td valign=top>(second tree for the crossover; if parameter doesn't exist, tree is picked at random.  This tree <b>must</b> have the same GPTreeConstraints as <tt>tree.0</tt>, if <tt>tree.0</tt> is defined.)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>ns.</tt><i>n</i><br>
 * <font size=-1>classname, inherits and != GPNodeSelector,<br>
 * or String <tt>same<tt></font></td>
 * <td valign=top>(GPNodeSelector for parent <i>n</i> (n is 0 or 1) If, for <tt>ns.1</tt> the value is <tt>same</tt>, then <tt>ns.1</tt> a copy of whatever <tt>ns.0</tt> is.  Note that the default version has no <i>n</i>)</td></tr>
 *
 * <tr><td valign=top><i>base</i>.<tt>toss</tt><br>
 * <font size=-1>bool = <tt>true</tt> or <tt>false</tt> (default)</font>/td>
 * <td valign=top>(after crossing over with the first new individual, should its second sibling individual be thrown away instead of adding it to the population?)</td></tr>
 * </table>
 *
 * <p><b>Default Base</b><br>
 * gp.koza.xover
 *
 * <p><b>Parameter bases</b><br>
 * <table>
 * <tr><td valign=top><i>base</i>.<tt>ns.</tt><i>n</i><br>
 * <td>nodeselect<i>n</i> (<i>n</i> is 0 or 1)</td></tr>
 * </table>
 *
 * @author Sean Luke
 * @version 1.0
 */

// Only simply protect the feature itself without its parent roots and only adopt it in pop1

public class ReplaceContributionCrossoverPipeline extends AllIndexCrossoverPipeline {

    int currentTreeID;

    public int produce(final int min,
                       final int max,
                       final int start,
                       final int subpopulation,
                       final Individual[] inds,
                       final EvolutionState state,
                       final int thread)

    {
        // how many individuals should we make?
        int n = typicalIndsProduced();
        if (n < min) n = min;
        if (n > max) n = max;

        // should we bother?
        if (!state.random[thread].nextBoolean(likelihood))
            return reproduce(n, start, subpopulation, inds, state, thread, true);  // DO produce children from source -- we've not done so already


        GPInitializer initializer = ((GPInitializer) state.initializer);

        for (int q = start; q < n + start; /* no increment */)  // keep on going until we're filled up
        {
            // grab two individuals from our sources
            if (sources[0] == sources[1])  // grab from the same source
                sources[0].produce(2, 2, 0, subpopulation, parents, state, thread);
            else // grab from different sources
            {
                sources[0].produce(1, 1, 0, subpopulation, parents, state, thread);
                sources[1].produce(1, 1, 1, subpopulation, parents, state, thread);
            }

            // at this point, parents[] contains our two selected individuals

            int length = parents[0].trees.length;
            if (tree1 == TREE_UNFIXED && tree2 == TREE_UNFIXED && (parents[0].trees.length == parents[1].trees.length)) {
                GPNode[][] p1 = new GPNode[length][2];     //each rule has two crossover points
                GPNode[][] p2 = new GPNode[length][2];
                for (int t = 0; t < length; t++) {

                    currentTreeID = t;

                    // prepare the nodeselectors
                    nodeselect1.reset();
                    nodeselect2.reset();


                    // pick some nodes


                    for (int x = 0; x < numTries; x++) {

                        GPNode p11;
                        GPNode p12;
                        GPNode p21;
                        GPNode p22;
                        // validity results...
                        boolean res1;
                        boolean res2;
                        // pick an unimportant node for all tasks in parent 0
                        p11 = nodeselect1.pickNode(state, subpopulation, thread, parents[0], parents[0].trees[t]);
                        // pick an important node for the current task in parent 0
                        p12 = nodeselect1.pickNode(state, subpopulation, thread, parents[0], parents[0].trees[t]);

                        // pick an unimportant node for all tasks in parent 1
                        p21 = nodeselect2.pickNode(state, subpopulation, thread, parents[1], parents[1].trees[t]);
                        // pick an important node for the current in parent 1
                        p22 = nodeselect2.pickNode(state, subpopulation, thread, parents[1], parents[1].trees[t]);

                        // check for depth and swap-compatibility limits
                        //  System.err.println(maxDepth + " " + maxSize);
                        res1 = verifyPoints(initializer, p22, p11);  // p2 can fill p1's spot -- order is important!
                        if (n - (q - start) < 2 || tossSecondParent) res2 = true;
                        else
                            res2 = verifyPoints(initializer, p12, p21);  // p1 can fill p2's spot -- order is important!

                        // did we get something that had both nodes verified?
                        // we reject if EITHER of them is invalid.  This is what lil-gp does.
                        // Koza only has numTries set to 1, so it's compatible as well.
                        if (res1 && res2) {
                            p1[t][0] = p11;
                            p1[t][1] = p12;
                            p2[t][0] = p21;
                            p2[t][1] = p22;
                            break;
                        }
                    }


                }

                // Create some new individuals based on the old ones -- since
                // GPTree doesn't deep-clone, this should be just fine.  Perhaps we
                // should change this to proto off of the main species prototype, but
                // we have to then copy so much stuff over; it's not worth it.

                GPIndividual j1 = parents[0].lightClone();
                GPIndividual j2 = null;
                if (n - (q - start) >= 2 && !tossSecondParent) j2 = parents[1].lightClone();

                // Fill in various tree information that didn't get filled in there
                j1.trees = new GPTree[parents[0].trees.length];
                if (n - (q - start) >= 2 && !tossSecondParent) j2.trees = new GPTree[parents[1].trees.length];

                // at this point, p1 or p2, or both, may be null.
                // If not, swap one in.  Else just copy the parent.

                for (int x = 0; x < j1.trees.length; x++) {
                    if (p1[x] != null)  // we've got a tree with a kicking cross position!
                    {
                        j1.trees[x] = parents[0].trees[x].lightClone();
                        j1.trees[x].owner = j1;
                        j1.trees[x].child = parents[0].trees[x].child.cloneReplacing(p2[x][1], p1[x][0]);
                        j1.trees[x].child.parent = j1.trees[x];
                        j1.trees[x].child.argposition = 0;
                        j1.evaluated = false;
                    }  // it's changed
                    else {
                        j1.trees[x] = parents[0].trees[x].lightClone();
                        j1.trees[x].owner = j1;
                        j1.trees[x].child = (GPNode) (parents[0].trees[x].child.clone());
                        j1.trees[x].child.parent = j1.trees[x];
                        j1.trees[x].child.argposition = 0;
                    }
                }

                if (n - (q - start) >= 2 && !tossSecondParent)
                    for (int x = 0; x < j2.trees.length; x++) {
                        if (p2[x] != null)  // we've got a tree with a kicking cross position!
                        {
                            j2.trees[x] = parents[1].trees[x].lightClone();
                            j2.trees[x].owner = j2;
                            j2.trees[x].child = parents[1].trees[x].child.cloneReplacing(p1[x][1], p2[x][0]);
                            j2.trees[x].child.parent = j2.trees[x];
                            j2.trees[x].child.argposition = 0;
                            j2.evaluated = false;
                        } // it's changed
                        else {
                            j2.trees[x] = parents[1].trees[x].lightClone();
                            j2.trees[x].owner = j2;
                            j2.trees[x].child = (GPNode) (parents[1].trees[x].child.clone());
                            j2.trees[x].child.parent = j2.trees[x];
                            j2.trees[x].child.argposition = 0;
                        }
                    }

                // add the individuals to the population
                inds[q] = j1;
                q++;
                if (q < n + start && !tossSecondParent) {
                    inds[q] = j2;
                    q++;
                }


            } else {
                state.output.fatal("GP AllIndexCrossover Pipeline: two individuals chosen for crossover have DIFFERENT numbers of trees! This is not supported -- you may wish to extend this method if you require this behaviour.");

            }

        }
        return n;
    }



}

