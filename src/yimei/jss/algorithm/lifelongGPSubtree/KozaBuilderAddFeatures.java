package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.gp.*;
import ec.gp.koza.HalfBuilder;

public class KozaBuilderAddFeatures extends HalfBuilder {

    public GPNode newRootedTree(final EvolutionState state,
                                final GPType type,
                                final int thread,
                                final GPNodeParent parent,
                                final GPFunctionSet set,
                                final int argposition,
                                final int requestedSize)
    {
        if (state.random[thread].nextDouble() < pickGrowProbability)
            return growNode(state,0,state.random[thread].nextInt(maxDepth-minDepth+1) + minDepth,type,thread,parent,argposition,set);
        else
            return fullNode(state,0,state.random[thread].nextInt(maxDepth-minDepth+1) + minDepth,type,thread,parent,argposition,set);
    }


    protected GPNode growNode(final EvolutionState state,
                              final int current,
                              final int max,
                              final GPType type,
                              final int thread,
                              final GPNodeParent parent,
                              final int argposition,
                              final GPFunctionSet set) {
        boolean triedTerminals = false;

        int t = type.type;
        GPNode[] terminals = set.terminals[t];
        GPNode[] nodes = set.nodes[t];

        if (nodes.length == 0)
            errorAboutNoNodeWithType(type, state);

        if ((current + 1 >= max) &&
                (triedTerminals = true) &&
                terminals.length != 0)
        {
            GPNode proto = terminals[state.random[thread].nextInt(terminals.length)];

            // 关键改动：深拷贝整个 feature 子树
            GPNode n = (GPNode)(proto.clone());

            n.argposition = (byte) argposition;
            n.parent = parent;
            n.resetNode(state, thread, set);

            return n;
        }
        else {
            if (triedTerminals)
                warnAboutNoTerminalWithType(type, false, state);

            GPNode proto = nodes[state.random[thread].nextInt(nodes.length)];

            // 如果你怀疑 nodes 里也可能混入已有子树，也可以统一改 clone()
            GPNode n = (GPNode)(proto.lightClone());
            n.resetNode(state, thread, set);
            n.argposition = (byte) argposition;
            n.parent = parent;

            GPType[] childtypes = n.constraints(((GPInitializer) state.initializer)).childtypes;
            for (int x = 0; x < childtypes.length; x++)
                n.children[x] = growNode(state, current + 1, max, childtypes[x], thread, n, x, set);

            return n;
        }
    }


    protected GPNode fullNode(final EvolutionState state,
                              final int current,
                              final int max,
                              final GPType type,
                              final int thread,
                              final GPNodeParent parent,
                              final int argposition,
                              final GPFunctionSet set) {
        boolean triedTerminals = false;

        int t = type.type;
        GPNode[] terminals = set.terminals[t];
        GPNode[] nonterminals = set.nonterminals[t];
        GPNode[] nodes = set.nodes[t];

        if (nodes.length == 0)
            errorAboutNoNodeWithType(type, state);

        if ((current + 1 >= max ||
                warnAboutNonterminal(nonterminals.length == 0, type, false, state)) &&
                (triedTerminals = true) &&
                terminals.length != 0)
        {
            GPNode proto = terminals[state.random[thread].nextInt(terminals.length)];

            // 关键改动：用深拷贝，不用 lightClone()
            GPNode n = (GPNode)(proto.clone());

            // clone() 会把子节点 parent/argposition 设好，
            // 但根节点自己的 parent/argposition 还需要这里补上
            n.argposition = (byte) argposition;
            n.parent = parent;

            // 如有 ERC，可只对根 reset；是否递归 reset 取决于你的设计
            n.resetNode(state, thread, set);

            return n;
        }
        else {
            if (triedTerminals)
                warnAboutNoTerminalWithType(type, false, state);

            GPNode[] nodesToPick = set.nonterminals[type.type];
            if (nodesToPick == null || nodesToPick.length == 0)
                nodesToPick = set.terminals[type.type];

            GPNode n = (GPNode)(nodesToPick[state.random[thread].nextInt(nodesToPick.length)].lightClone());
            n.resetNode(state, thread);
            n.argposition = (byte) argposition;
            n.parent = parent;

            GPType[] childtypes = n.constraints(((GPInitializer) state.initializer)).childtypes;
            for (int x = 0; x < childtypes.length; x++)
                n.children[x] = fullNode(state, current + 1, max, childtypes[x], thread, n, x, set);

            return n;
        }
    }


}
