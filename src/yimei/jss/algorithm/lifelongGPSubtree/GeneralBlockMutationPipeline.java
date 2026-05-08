package yimei.jss.algorithm.lifelongGPSubtree;

import ec.*;
import ec.gp.*;
import ec.gp.koza.MutationPipeline;

import java.util.ArrayList;

public class GeneralBlockMutationPipeline extends MutationPipeline {

    private static final long serialVersionUID = 1L;

    @Override
    public int produce(final int min,
                       final int max,
                       final int start,
                       final int subpopulation,
                       final Individual[] inds,
                       final EvolutionState state,
                       final int thread) {

        int n = sources[0].produce(min, max, start, subpopulation, inds, state, thread);

        if (!state.random[thread].nextBoolean(likelihood)) {
            return reproduce(n, start, subpopulation, inds, state, thread, false);
        }

        GPInitializer initializer = (GPInitializer) state.initializer;
        GPRuleEvolutionStateLifelongGP ls = (GPRuleEvolutionStateLifelongGP) state;

        for (int q = start; q < n + start; q++) {

            GPIndividual i = (GPIndividual) inds[q];

            if (tree != TREE_UNFIXED && (tree < 0 || tree >= i.trees.length)) {
                state.output.fatal(
                        "GeneralBlockMutationPipeline attempted to fix tree.0 to an invalid tree index."
                );
            }

            int t;
            if (tree == TREE_UNFIXED) {
                t = i.trees.length > 1 ? state.random[thread].nextInt(i.trees.length) : 0;
            } else {
                t = tree;
            }

            boolean res = false;

            nodeselect.reset();

            GPNode p1 = null;  // selected mutation point
            GPNode p2 = null;  // selected general building block

            for (int x = 0; x < numTries; x++) {

                p1 = nodeselect.pickNode(state, subpopulation, thread, i, i.trees[t]);

                int size = GPNodeBuilder.NOSIZEGIVEN;
                if (equalSize)
                    size = p1.numNodes(GPNode.NODESEARCH_ALL);

                // 先尝试从 general blocks 取
                p2 = randomGeneralBlockNode(ls, t, state, thread);

                // ❗如果没有 general block → fallback 到随机生成
                if (p2 == null) {
                    p2 = builder.newRootedTree(state,
                            p1.parentType(initializer),
                            thread,
                            p1.parent,
                            i.trees[t].constraints(initializer).functionset,
                            p1.argposition,
                            size);
                }

                // 类型检查（很重要）
                if (!p2.swapCompatibleWith(initializer, p1)) {
                    continue;
                }

                // 深度 / size 检查
                res = verifyPoints(p2, p1);

                if (res) {
                    break;
                }

                p2.parent = p1.parent;
                p2.argposition = p1.argposition;

                if (!p2.swapCompatibleWith(initializer, p1)) {
                    continue;
                }

                res = verifyPoints(p2, p1);

                if (res) {
                    break;
                }
            }

            GPIndividual j;

            if (sources[0] instanceof BreedingPipeline) {

                j = i;

                if (res) {
                    p2.parent = p1.parent;
                    p2.argposition = p1.argposition;

                    if (p2.parent instanceof GPNode) {
                        ((GPNode) p2.parent).children[p2.argposition] = p2;
                    } else {
                        ((GPTree) p2.parent).child = p2;
                    }

                    j.evaluated = false;
                }

            } else {

                j = (GPIndividual) i.lightClone();
                j.trees = new GPTree[i.trees.length];

                for (int x = 0; x < j.trees.length; x++) {

                    if (x == t && res) {

                        j.trees[x] = (GPTree) i.trees[x].lightClone();
                        j.trees[x].owner = j;

                        p2.parent = p1.parent;
                        p2.argposition = p1.argposition;

                        j.trees[x].child =
                                i.trees[x].child.cloneReplacingNoSubclone(p2, p1);

                        j.trees[x].child.parent = j.trees[x];
                        j.trees[x].child.argposition = 0;

                        j.evaluated = false;

                    } else {

                        j.trees[x] = (GPTree) i.trees[x].lightClone();
                        j.trees[x].owner = j;
                        j.trees[x].child = (GPNode) i.trees[x].child.clone();
                        j.trees[x].child.parent = j.trees[x];
                        j.trees[x].child.argposition = 0;
                    }
                }
            }

            inds[q] = j;
        }

        return n;
    }

    private GPNode randomGeneralBlockNode(final GPRuleEvolutionStateLifelongGP state,
                                          final int treeID,
                                          final EvolutionState evoState,
                                          final int thread) {

        if (state.GeneralBuildingBlocks == null) {
            return null;
        }

        if (treeID < 0 || treeID >= state.GeneralBuildingBlocks.size()) {
            return null;
        }

        ArrayList<GPTree> blocks = state.GeneralBuildingBlocks.get(treeID);

        if (blocks == null || blocks.isEmpty()) {
            return null;
        }

        GPTree selectedBlock = blocks.get(
                evoState.random[thread].nextInt(blocks.size())
        );

        if (selectedBlock == null || selectedBlock.child == null) {
            return null;
        }

        return (GPNode) selectedBlock.child.clone();
    }

    @Override
    public int produceFrac(int min,
                           int max,
                           int start,
                           int subpopulation,
                           Individual[] inds,
                           EvolutionState state,
                           int thread) {
        return 0;
    }
}