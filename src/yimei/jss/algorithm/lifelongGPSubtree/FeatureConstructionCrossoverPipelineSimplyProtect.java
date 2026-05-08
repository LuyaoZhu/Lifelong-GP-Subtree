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

public class FeatureConstructionCrossoverPipelineSimplyProtect extends AllIndexCrossoverPipeline {


    public static final int SEQ_TREE_ID = 0;
    public static final int ROU_TREE_ID = 1;
    public static final int NUM_PARENTS = 2;
    public static final int NUM_TREES = 2;

    // [parentID][treeID]
    public static ArrayList<Integer>[][] protectedNonTerminalsRank = new ArrayList[NUM_PARENTS][NUM_TREES];
    public static ArrayList<Integer>[][] protectedTerminalsRank = new ArrayList[NUM_PARENTS][NUM_TREES];

    static {
        for (int p = 0; p < NUM_PARENTS; p++) {
            for (int t = 0; t < NUM_TREES; t++) {
                protectedNonTerminalsRank[p][t] = new ArrayList<>();
                protectedTerminalsRank[p][t] = new ArrayList<>();
            }
        }
    }

    public static boolean parent1Flag;
    public static int currentTreeID = 0;


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

            //now identify the blocks
            featureIdentify(state);

            int length = parents[0].trees.length;
            if (tree1 == TREE_UNFIXED && tree2 == TREE_UNFIXED && (parents[0].trees.length == parents[1].trees.length)) {
                GPNode[] p1 = new GPNode[length];
                GPNode[] p2 = new GPNode[length];
                for (int t = 0; t < length; t++) {

                    currentTreeID = t;

                    // prepare the nodeselectors
                    nodeselect1.reset();
                    nodeselect2.reset();


                    // pick some nodes


                    for (int x = 0; x < numTries; x++) {

                        GPNode p11;
                        GPNode p21;
                        // validity results...
                        boolean res1;
                        boolean res2;
                        // pick a node in individual 1
                        parent1Flag = true;
                        p11 = nodeselect1.pickNode(state, subpopulation, thread, parents[0], parents[0].trees[t]);

                        // pick a node in individual 2
                        parent1Flag = false;
                        p21 = nodeselect2.pickNode(state, subpopulation, thread, parents[1], parents[1].trees[t]);

                        // check for depth and swap-compatibility limits
                        //  System.err.println(maxDepth + " " + maxSize);
                        res1 = verifyPoints(initializer, p21, p11);  // p2 can fill p1's spot -- order is important!
                        if (n - (q - start) < 2 || tossSecondParent) res2 = true;
                        else
                            res2 = verifyPoints(initializer, p11, p21);  // p1 can fill p2's spot -- order is important!

                        // did we get something that had both nodes verified?
                        // we reject if EITHER of them is invalid.  This is what lil-gp does.
                        // Koza only has numTries set to 1, so it's compatible as well.
                        if (res1 && res2) {
                            p1[t] = p11;
                            p2[t] = p21;
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
                        j1.trees[x].child = parents[0].trees[x].child.cloneReplacing(p2[x], p1[x]);
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
                            j2.trees[x].child = parents[1].trees[x].child.cloneReplacing(p1[x], p2[x]);
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

    private void featureIdentify(EvolutionState state) {

        clearProtectedRanks();

        GPRuleEvolutionStateLifelongGP ls =
                (GPRuleEvolutionStateLifelongGP) state;

        if (state.generation < ls.generationPerTask) {
            return;
        }

        for (int parentID = 0; parentID < parents.length; parentID++) {

            for (int treeID = 0; treeID < parents[parentID].trees.length; treeID++) {

                GPTree parentTree = parents[parentID].trees[treeID];

                List<GPTree> archiveBlocks = getArchiveBlocks(ls, treeID);

                List<MatchedBlock> matches = findMatchedBlocks(parentTree, archiveBlocks);

                recordProtectedRanks(parentID, treeID, parentTree, matches);
            }
        }
    }

    /**
     * If you only want to protect the block itself, keep this false.
     * If you also want to protect all ancestors from the block root to the full-tree root, set true.
     */
    private static final boolean PROTECT_ANCESTORS = false;

    /**
     * Store one matched building block.
     */
    private static class MatchedBlock {
        GPNode matchedRoot;
        GPTree archiveBlock;

        MatchedBlock(GPNode matchedRoot, GPTree archiveBlock) {
            this.matchedRoot = matchedRoot;
            this.archiveBlock = archiveBlock;
        }
    }

    private void clearProtectedRanks() {
        for (int p = 0; p < NUM_PARENTS; p++) {
            for (int t = 0; t < NUM_TREES; t++) {
                protectedNonTerminalsRank[p][t].clear();
                protectedTerminalsRank[p][t].clear();
            }
        }
    }

    private List<GPTree> getArchiveBlocks(GPRuleEvolutionStateLifelongGP state, int treeID) {

        List<GPTree> blocks = new ArrayList<>();

        if (state.TaskSpecificBuildingBlocks == null) {
            return blocks;
        }

        for (int taskID = 0; taskID < state.TaskSpecificBuildingBlocks.size(); taskID++) {

            ArrayList<ArrayList<GPTree>> blocksCurrentTask =
                    state.TaskSpecificBuildingBlocks.get(taskID);

            if (blocksCurrentTask == null || blocksCurrentTask.size() <= treeID) {
                continue;
            }

            ArrayList<GPTree> blocksForThisTree = blocksCurrentTask.get(treeID);

            if (blocksForThisTree != null) {
                blocks.addAll(blocksForThisTree);
            }
        }

        return blocks;
    }

    /**
     * Find all archive blocks that exactly appear in the given parent tree.
     */
    private List<MatchedBlock> findMatchedBlocks(GPTree parentTree, List<GPTree> archiveBlocks) {

        List<MatchedBlock> matches = new ArrayList<>();

        if (parentTree == null || parentTree.child == null || archiveBlocks == null || archiveBlocks.isEmpty()) {
            return matches;
        }

        int numNonTerminals = parentTree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);

        for (int pos = 0; pos < numNonTerminals; pos++) {

            GPNode candidateRoot =
                    parentTree.child.nodeInPosition(pos, GPNode.NODESEARCH_NONTERMINALS);

            for (GPTree block : archiveBlocks) {

                if (block == null || block.child == null) {
                    continue;
                }

                if (sameLispTree(candidateRoot, block.child)) {
                    matches.add(new MatchedBlock(candidateRoot, block));
                    break;
                }
            }
        }

        return matches;
    }

    /**
     * Current exact subtree matching.
     * Later you can replace this with semantic matching.
     */
    private boolean sameLispTree(GPNode a, GPNode b) {
        return a.makeLispTree().equals(b.makeLispTree());
    }

    /**
     * Record protected ranks for one parent and one tree.
     */
    private void recordProtectedRanks(
            int parentID,
            int treeID,
            GPTree fullTree,
            List<MatchedBlock> matches
    ) {
        if (matches == null || matches.isEmpty()) {
            return;
        }

        Set<Integer> protectedNonTerminalRanks = new LinkedHashSet<>();
        Set<Integer> protectedTerminalRanks = new LinkedHashSet<>();

        for (MatchedBlock match : matches) {

            Set<GPNode> protectedNodes = new HashSet<>();

            collectAllNodes(match.matchedRoot, protectedNodes);

            if (PROTECT_ANCESTORS) {
                collectAncestors(match.matchedRoot, protectedNodes);
            }

            collectNodeRanks(
                    fullTree,
                    protectedNodes,
                    GPNode.NODESEARCH_NONTERMINALS,
                    protectedNonTerminalRanks
            );

            collectNodeRanks(
                    fullTree,
                    protectedNodes,
                    GPNode.NODESEARCH_TERMINALS,
                    protectedTerminalRanks
            );
        }

        addUnique(protectedNonTerminalsRank[parentID][treeID], protectedNonTerminalRanks);
        addUnique(protectedTerminalsRank[parentID][treeID], protectedTerminalRanks);
    }

    /**
     * Collect all nodes under a matched block root.
     */
    private void collectAllNodes(GPNode root, Set<GPNode> nodes) {

        if (root == null) {
            return;
        }

        nodes.add(root);

        for (GPNode child : root.children) {
            collectAllNodes(child, nodes);
        }
    }

    /**
     * Optional: protect parent nodes from the matched block root to the full-tree root.
     */
    private void collectAncestors(GPNode node, Set<GPNode> nodes) {

        GPNode current = node;

        while (current != null && current.parent instanceof GPNode) {
            current = (GPNode) current.parent;
            nodes.add(current);
        }
    }

    /**
     * Convert protected GPNode references to ECJ node ranks.
     */
    private void collectNodeRanks(
            GPTree fullTree,
            Set<GPNode> protectedNodes,
            int searchType,
            Set<Integer> ranks
    ) {
        int numNodes = fullTree.child.numNodes(searchType);

        for (int i = 0; i < numNodes; i++) {
            GPNode node = fullTree.child.nodeInPosition(i, searchType);

            if (protectedNodes.contains(node)) {
                ranks.add(i);
            }
        }
    }

    private void addUnique(ArrayList<Integer> target, Set<Integer> source) {
        for (Integer value : source) {
            if (!target.contains(value)) {
                target.add(value);
            }
        }
    }

    /*private void featureIdentify(EvolutionState state) {
        ArrayList<ArrayList<GPTree>> protectedParent1Subtrees = new ArrayList<>();
        ArrayList<ArrayList<GPTree>> protectedParent2Subtrees = new ArrayList<>();

        ArrayList<ArrayList<Integer>> protectedParent1SubtreesRank = new ArrayList<>();
        ArrayList<ArrayList<Integer>> protectedParent2SubtreesRank = new ArrayList<>();


        if (state.generation >= ((GPRuleEvolutionStateLifelongGP) state).generationPerTask) {
            //record the feature trees and its ranks in the parents

            for (int parent = 0; parent < parents.length; parent++) {  //for two parents

                for (int treeID = 0; treeID < parents[parent].trees.length; treeID++) {

                    GPTree tree = parents[parent].trees[treeID];

                    ArrayList<GPTree> protectedSubtreesForOneRule = new ArrayList<>();
                    ArrayList<Integer> protectedSubtreesRankForOneRule = new ArrayList<>();

                    if (treeID == 0) {

                        if (tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS) > 1) {
                            int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree baseTree = PopulationUtils.GPNodetoGPTree(node);
                                //compare with the features in the archive
                                if (parent == 0) {
                                    for(int taskID=0; taskID<((GPRuleEvolutionStateLifelongGP) state).TaskSpecificBuildingBlocks.size(); taskID++) {
                                        for (int f = 0; f < ((GPRuleEvolutionStateLifelongGP) state).TaskSpecificBuildingBlocks.get(parent).size(); f++) {
                                            GPTree countTree = (GPTree) ((GPRuleEvolutionStateLifelongGP) state).TaskSpecificBuildingBlocks.get(parent).get(f).clone();
                                            if (baseTree.child.makeLispTree().equals(countTree.child.makeLispTree())) {
                                                protectedParent1SubtreesRank.add(numSubtree);
                                                protectedParent1Subtrees.add(countTree);
                                            }
                                        }
                                    }

                                } else {
                                    for (int f = 0; f < ((GPRuleEvolutionStateLifelongGP) state).seqFeatureArchive.size(); f++) {
                                        GPTree countTree = (GPTree) ((GPRuleEvolutionStateLifelongGP) state).seqFeatureArchive.get(f).clone();
                                        if (baseTree.child.makeLispTree().equals(countTree.child.makeLispTree())) {
                                            protectedParent2SubtreesRank.add(numSubtree);
                                            protectedParent2Subtrees.add(countTree);
                                        }
                                    }
                                }
                            }
                        }
                    } else // subpop == 1
                    {
                        if (tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS) > 1) {
                            int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree baseTree = PopulationUtils.GPNodetoGPTree(node);
                                //compare with the features in the archive
                                if (parent == 0) {
                                    for (int f = 0; f < ((GPRuleEvolutionStateLifelongGPV0) state).rouFeatureArchive.size(); f++) {
                                        GPTree countTree = (GPTree) ((GPRuleEvolutionStateLifelongGPV0) state).rouFeatureArchive.get(f).clone();
                                        if (baseTree.child.makeLispTree().equals(countTree.child.makeLispTree())) {
                                            protectedParent1SubtreesRank.add(numSubtree);
                                            protectedParent1Subtrees.add(countTree);
                                        }
                                    }
                                } else {
                                    for (int f = 0; f < ((GPRuleEvolutionStateLifelongGPV0) state).rouFeatureArchive.size(); f++) {
                                        GPTree countTree = (GPTree) ((GPRuleEvolutionStateLifelongGPV0) state).rouFeatureArchive.get(f).clone();
                                        if (baseTree.child.makeLispTree().equals(countTree.child.makeLispTree())) {
                                            protectedParent2SubtreesRank.add(numSubtree);
                                            protectedParent2Subtrees.add(countTree);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

//            double rnd = state.random[thread].nextDouble();
//            if (rnd >= 2) {
//                protectedParent1SubtreesRank.clear();
//                protectedParent2SubtreesRank.clear();
//                protectedParent1Subtrees.clear();
//                protectedParent2Subtrees.clear();
//            }


//            featureNumInOneTree[protectedParent1SubtreesRank.size()]++;
//            featureNumInOneTree[protectedParent2SubtreesRank.size()]++;

        //we need to record the rank of protected subtrees
        if (protectedParent1Subtrees.size() > 0) {
            int fullTreeDepth = parents[0].trees[0].child.depth();
            for (int i = 0; i < protectedParent1SubtreesRank.size(); i++) {
                int subtreeNonterminals = protectedParent1Subtrees.get(i).child.numNodes(GPNode.NODESEARCH_NONTERMINALS);

                //record the child nonTerminal nodes
                for (int num = protectedParent1SubtreesRank.get(i); num < protectedParent1SubtreesRank.get(i) + subtreeNonterminals; num++) {
                    if (!parent1ProtectedNonTerminalsRank.contains(num))
                        parent1ProtectedNonTerminalsRank.add(num);
                }

                //record the parents nodes
                    *//*int featureDepth = protectedParent1Subtrees.get(i).child.depth();
                    int fulltreeNonterminals = parents[0].trees[0].child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
*//**//*                    for (int depth = fullTreeDepth; depth > featureDepth; depth--) {
                        for (int numSubtree = 0; numSubtree < fulltreeNonterminals; numSubtree++) {
                            GPTree treeClone = (GPTree) parents[0].trees[0].clone();
                            GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                            GPTree baseTree = PopulationUtils.GPNodetoGPTree(node);
                            if (baseTree.child.depth() == depth) {
                                for (int x = 0; x < baseTree.child.children.length; x++) {
                                    if (baseTree.child.children[x].makeLispTree().equals(protectedParent1Subtrees.get(i).child.makeLispTree())) {
                                        parent1ProtectedNonTerminalsRank.add(numSubtree);
                                    }
                                }
                            }
                        }
                    }*//**//*

                    for (int depth = fullTreeDepth; depth > featureDepth; depth--) {
                        for (int numSubtree = 0; numSubtree < fulltreeNonterminals; numSubtree++) {
                            GPNode node = parents[0].trees[0].child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                            if (node.depth() == depth) {
                                if (LispContains(node, protectedParent1Subtrees.get(i).child)) {
                                    if (!parent1ProtectedNonTerminalsRank.contains(numSubtree))
                                        parent1ProtectedNonTerminalsRank.add(numSubtree);
                                }
                            }
                        }
                    }*//*

                //record the child terminal nodes
                int subtreeTerminals = protectedParent1Subtrees.get(i).child.numNodes(GPNode.NODESEARCH_TERMINALS);
                int fulltreeTerminals = parents[0].trees[0].child.numNodes(GPNode.NODESEARCH_TERMINALS);
                for (int a = 0; a < fulltreeTerminals; a++) {
                    int sameNum = 0;
                    for (int b = 0; b < subtreeTerminals; b++) {
                        GPNode baseNode = parents[0].trees[0].child.nodeInPosition(a + b, GPNode.NODESEARCH_TERMINALS);
                        GPNode countNode = protectedParent1Subtrees.get(i).child.nodeInPosition(b, GPNode.NODESEARCH_TERMINALS);
                        if (LispContains(baseNode, countNode)) {
                            sameNum++;
                        }
                    }
                    if (sameNum == subtreeTerminals) {
                        for (int num = a; num < a + subtreeTerminals; num++) {
                            if (!parent1ProtectedTerminalsRank.contains(num))
                                parent1ProtectedTerminalsRank.add(num);
                        }
                        break;
                    }
                }
            }
            //record the child terminal nodes
//                for (int a = 0; a < parents[0].trees[0].child.numNodes(GPNode.NODESEARCH_TERMINALS); a++) {
//                    if
//                    GPNode baseNode = parents[0].trees[0].child.nodeInPosition(a, GPNode.NODESEARCH_TERMINALS);
//                    for (int b = 0; b < protectedParent1Subtrees.size(); b++) {
//                        if (LispContains(protectedParent1Subtrees.get(b).child, baseNode)) {
//                            int depth = protectedParent1Subtrees.get(b).child.depth();
//                            GPNode tempNode = (GPNode) baseNode.clone();
//                            for (int d = 0; d < depth - 1; d++) {
//                                tempNode = (GPNode) tempNode.parent;
//                            }
//                            if (tempNode.makeLispTree().equals(protectedParent1Subtrees.get(b).child.makeLispTree())) {
//                                parent1ProtectedTerminalsRank.add(a);
//                                break;
//                            }
//                        }
//                    }
//                }

        }


        if (protectedParent2Subtrees.size() > 0) {

            int fullTreeDepth = parents[1].trees[0].child.depth();
            for (int i = 0; i < protectedParent2SubtreesRank.size(); i++) {
                int subtreeNonterminals = protectedParent2Subtrees.get(i).child.numNodes(GPNode.NODESEARCH_NONTERMINALS);

                //record the child nonTerminal nodes
                for (int num = protectedParent2SubtreesRank.get(i); num < protectedParent2SubtreesRank.get(i) + subtreeNonterminals; num++) {
                    if (!parent2ProtectedNonTerminalsRank.contains(num))
                        parent2ProtectedNonTerminalsRank.add(num);
                }
                //record the parents nodes
                int featureDepth = protectedParent2Subtrees.get(i).child.depth();
                int fulltreeNonterminals = parents[1].trees[0].child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
//                    for (int depth = fullTreeDepth; depth > featureDepth; depth--) {
//                        for (int numSubtree = 0; numSubtree < fulltreeNonterminals; numSubtree++) {
//                            GPTree treeClone = (GPTree) parents[1].trees[0].clone();
//                            GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
//                            GPTree baseTree = PopulationUtils.GPNodetoGPTree(node);
//                            if (baseTree.child.depth() == depth) {
//                                for (int x = 0; x < baseTree.child.children.length; x++) {
//                                    if (baseTree.child.children[x].makeLispTree().equals(protectedParent2Subtrees.get(i).child.makeLispTree())) {
//                                        parent2ProtectedNonTerminalsRank.add(numSubtree);
//                                    }
//                                }
//                            }
//                        }
//                    }

                for (int depth = fullTreeDepth; depth > featureDepth; depth--) {
                    for (int numSubtree = 0; numSubtree < fulltreeNonterminals; numSubtree++) {
                        GPNode node = parents[1].trees[0].child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                        if (node.depth() == depth) {
                            if (LispContains(node, protectedParent2Subtrees.get(i).child)) {
                                if (!parent2ProtectedNonTerminalsRank.contains(numSubtree))
                                    parent2ProtectedNonTerminalsRank.add(numSubtree);
                            }
                        }
                    }
                }


                //record the child terminal nodes
                int subtreeTerminals = protectedParent2Subtrees.get(i).child.numNodes(GPNode.NODESEARCH_TERMINALS);
                int fulltreeTerminals = parents[1].trees[0].child.numNodes(GPNode.NODESEARCH_TERMINALS);
                for (int a = 0; a < fulltreeTerminals; a++) {
                    int sameNum = 0;
                    for (int b = 0; b < subtreeTerminals; b++) {
                        GPNode baseNode = parents[1].trees[0].child.nodeInPosition(a + b, GPNode.NODESEARCH_TERMINALS);
                        GPNode countNode = protectedParent2Subtrees.get(i).child.nodeInPosition(b, GPNode.NODESEARCH_TERMINALS);
                        if (LispContains(baseNode, countNode)) {
                            sameNum++;
                        }
                    }
                    if (sameNum == subtreeTerminals) {
                        for (int num = a; num < a + subtreeTerminals; num++) {
                            if (!parent2ProtectedTerminalsRank.contains(num))
                                parent2ProtectedTerminalsRank.add(num);
                        }
                        break;
                    }
                }

            }


*//*                //record the child terminal nodes
                for (int a = 0; a < parents[1].trees[0].child.numNodes(GPNode.NODESEARCH_TERMINALS); a++) {
                    GPNode baseNode = parents[1].trees[0].child.nodeInPosition(a, GPNode.NODESEARCH_TERMINALS);
                    for (int b = 0; b < protectedParent2Subtrees.size(); b++) {
                        if (LispContains(protectedParent2Subtrees.get(b).child, baseNode)) {
                            int depth = protectedParent2Subtrees.get(b).child.depth();
                            GPNode tempNode = (GPNode) baseNode.clone();
                            for (int d = 0; d < depth - 1; d++) {
                                tempNode = (GPNode) tempNode.parent;
                            }
                            if (tempNode.makeLispTree().equals(protectedParent2Subtrees.get(b).child.makeLispTree())) {
                                parent2ProtectedTerminalsRank.add(a);
                                break;
                            }
                        }
                    }
                }*//*

        }
    }*/


    public static boolean LispContains(final GPNode node, final GPNode subnode) {
        if (subnode.makeLispTree().equals(node.makeLispTree())) return true;
        for (int x = 0; x < node.children.length; x++)
            if (LispContains(node.children[x], subnode)) return true;
        return false;
    }

}

