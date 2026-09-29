package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.gp.*;
import ec.gp.koza.KozaNodeSelector;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.CalcPriorityProblem;
import yimei.jss.gp.data.DoubleData;
import yimei.jss.gp.terminal.DoubleERC;
import yimei.jss.gp.terminal.TerminalERC;
import yimei.jss.gp.terminal.TerminalERCUniform;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;
import yimei.jss.simulation.state.SystemState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TaskwiseImportanceNodeSelector extends KozaNodeSelector {

    public enum ImportanceMode {
        CURRENT_TASK,
        HISTORICAL_TASKS,
        ALL_TASKS
    }

    private java.util.ArrayList<PhenoCharacterisation[]> decisionSituationsEachTask;

    private boolean pickHighImportance = true;

    private ImportanceMode mode = ImportanceMode.ALL_TASKS;

    private int currentTaskIndex = -1;

    private int treeIndex = 0;

    private double lastSelectedImportance = 0.0;

    private double lastSelectedBehaviourChange = 0.0;

    private TerminalERCUniform terminalPrototype = null;

    public static final java.util.HashMap<String, double[]> importanceCache = new java.util.HashMap<>();

    private String currentGroupLabel = null; // "HG" / "CS" / "CG"，由外部（pipeline）在调用 pickNode 前设置

    public void setCurrentGroupLabel(String currentGroupLabel) {
        this.currentGroupLabel = currentGroupLabel;
    }

    private static int cachedGeneration = -1;

    private static final double TOP_IMPORTANCE_RATIO = 0.2; // 前 20%，改这个常量就能调比例


    public void setDecisionSituationsEachTask(
            java.util.ArrayList<PhenoCharacterisation[]> decisionSituationsEachTask) {
        this.decisionSituationsEachTask = decisionSituationsEachTask;
    }

    public void setPickHighImportance(boolean pickHighImportance) {
        this.pickHighImportance = pickHighImportance;
    }

    public double getLastSelectedBehaviourChange() {

        return lastSelectedBehaviourChange;
    }

    public void setMode(ImportanceMode mode) {
        this.mode = mode;
    }

    public double getLastSelectedImportance() {

        return lastSelectedImportance;
    }

    public void setCurrentTaskIndex(int currentTaskIndex) {
        this.currentTaskIndex = currentTaskIndex;
    }

    public void setTreeIndex(int treeIndex) {
        this.treeIndex = treeIndex;
    }

    private String buildCacheKey(GPIndividual ind) {

        String tag;

        if (mode == ImportanceMode.HISTORICAL_TASKS) {

            tag = "HG_" + currentTaskIndex;

        } else if (mode == ImportanceMode.CURRENT_TASK) {

            tag = "CS_" + currentTaskIndex;

        } else {

            tag = "CG_" + cachedGeneration;
        }

        return System.identityHashCode(ind)
                + "_"
                + treeIndex
                + "_"
                + tag;
    }

    private void refreshGeneration(EvolutionState state) {

        if (cachedGeneration != state.generation) {

            cachedGeneration = state.generation;

            importanceCache.entrySet().removeIf(
                    e -> e.getKey().contains("_CG_"));
        }
    }


    @Override
    public GPNode pickNode(final EvolutionState state,
                           final int subpopulation,
                           final int thread,
                           final GPIndividual ind,
                           final GPTree tree) {

        double rnd = state.random[thread].nextDouble(); //probability  (0,1)

        if (rnd > 1)  // pick anyone  nonterminalProbability = 0.9  terminalProbability=0.1
        // rnd > 0.9+0.1
        //nonterminalProbability + terminalProbability + rootProbability = 1  this will not happen
        {
            if (nodes == -1) nodes = tree.child.numNodes(GPNode.NODESEARCH_ALL); //nodes: the number of node in the tree
            //including the terminals    all the possible positions
            {
                return tree.child.nodeInPosition(state.random[thread].nextInt(nodes), GPNode.NODESEARCH_ALL);
                //randomly choose a node
            }
        } else if (rnd > 1)  // pick the root
        {
            return tree.child; //for example: nonterminalProbability = 0.8 terminalProbability = 0.1, of rnd = 0.9. will choose the root
        } else if (rnd > 0.9)  // pick terminals  //nonterminalProbability = 0.9
        {
            if (terminals == -1) terminals = tree.child.numNodes(GPNode.NODESEARCH_TERMINALS);
            return tree.child.nodeInPosition(state.random[thread].nextInt(terminals), GPNode.NODESEARCH_TERMINALS);
            //choose the terminals
        } else  // pick nonterminals if you can
        {
            refreshGeneration(state);

            String cacheKey = buildCacheKey(ind);

            double[] cachedScores = importanceCache.get(cacheKey);

            if (tree == null || tree.child == null) {
                return null;
            }

            if (cachedScores == null) {


                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);

                if (nonterminals <= 0) {
                    return tree.child;
                }

                cachedScores = new double[nonterminals];
                cachedScores[0] = 1;
                List<double[]> importanceVectorsThisTree = new ArrayList<>();
                // 遍历所有 nonterminal subtree
                for (int nodeIndex = 0; nodeIndex < nonterminals; nodeIndex++) {

                    GPTree treeClone = (GPTree) tree.clone();
                    GPNode node;
                    GPNode originalNode;

                    node = treeClone.child.nodeInPosition(nodeIndex, GPNode.NODESEARCH_NONTERMINALS);

                    originalNode = tree.child.nodeInPosition(nodeIndex, GPNode.NODESEARCH_NONTERMINALS);

                    if (node == null || node.parent instanceof GPTree) {
                        continue;
                    }

                    if (terminalPrototype == null && node.children != null) {

                        for (GPNode child : node.children) {
                            if (child instanceof TerminalERCUniform) {
                                terminalPrototype = (TerminalERCUniform) child;
                                break;
                            }
                        }
                    }

                    double[] importanceVector =
                            getReplaceImportanceVector(
                                    node,
                                    originalNode,
                                    tree,
                                    treeIndex);

                    cachedScores[nodeIndex] = aggregateImportance(importanceVector);

                    importanceVectorsThisTree.add(importanceVector); // 新写法：先收集

                }

                importanceCache.put(cacheKey, cachedScores);
                recordSubtreeImportance(state, importanceVectorsThisTree); // 整棵树算完后统一调用一次，内部排序取前 20%
            }

            return sampleByImportance(
                    state,
                    thread,
                    tree,
                    cachedScores);

        }
    }

    private double[] getReplaceImportanceVector(
            GPNode node,
            GPNode originalNode,
            GPTree fullTree,
            int treeIndex) {

        if (decisionSituationsEachTask == null
                || decisionSituationsEachTask.isEmpty()) {
            return new double[]{0.0};
        }

        int taskNum =
                decisionSituationsEachTask.size();

        double[] importanceVector =
                new double[taskNum];

        for (int task = 0; task < taskNum; task++) {

//            switch (mode) {
//
//                case ALL_TASKS:
//                    break;
//
//                case HISTORICAL_TASKS:
//                    if (task == taskNum - 1) continue;
//                    break;
//
//                case CURRENT_TASK:
//                    if (task != taskNum - 1) continue;
//                    break;
//            }

            PhenoCharacterisation[] pc =
                    decisionSituationsEachTask.get(task);

            importanceVector[task] =
                    calculateTaskImportance(
                            node,
                            originalNode,
                            fullTree,
                            pc,
                            treeIndex);
        }

        return importanceVector;
    }

    private void recordSubtreeImportance(EvolutionState state, List<double[]> importanceVectors) {

        if (currentGroupLabel == null) {
            return;
        }

        if (!(state instanceof GPRuleEvolutionStateLifelongGPReplace)) {
            return;
        }

        if (importanceVectors == null || importanceVectors.isEmpty()) {
            return;
        }

        int rankIndex;

        switch (currentGroupLabel) {
            case "HG": rankIndex = 0; break; // previous
            case "CS": rankIndex = 1; break; // current
            case "CG": rankIndex = 2; break; // all
            default: return;
        }

        // 每个节点的 importanceVector -> {previous, current, all}
        List<double[]> triples = new ArrayList<>(importanceVectors.size());

        for (double[] importanceVector : importanceVectors) {

            int taskNum = importanceVector.length;

            if (taskNum == 0) {
                continue;
            }

            double sumAll = 0.0;
            for (double v : importanceVector) {
                sumAll += v;
            }

            double current = importanceVector[taskNum - 1];
            double all = sumAll / taskNum;
            double previous = (taskNum == 1)
                    ? importanceVector[0]
                    : (sumAll - current) / (taskNum - 1);

            triples.add(new double[]{ previous, current, all });
        }

        if (triples.isEmpty()) {
            return;
        }

        // 按该 group 对应的口径从高到低排序
        triples.sort((a, b) -> Double.compare(b[rankIndex], a[rankIndex]));

        int keepCount = (int) Math.ceil(triples.size() * TOP_IMPORTANCE_RATIO);
        keepCount = Math.max(1, Math.min(keepCount, triples.size()));

        GPRuleEvolutionStateLifelongGPReplace s = (GPRuleEvolutionStateLifelongGPReplace) state;

        for (int i = 0; i < keepCount; i++) {
            double[] t = triples.get(i);
            s.addSubtreeImportanceSample(currentGroupLabel, t[0], t[1], t[2]);
        }
    }
    private double aggregateImportance(double[] vector) {

        if (vector == null || vector.length == 0) {
            return 0.0;
        }

        if (mode == ImportanceMode.CURRENT_TASK) {
            int idx = currentTaskIndex;

            if (idx < 0 || idx >= vector.length) {
                idx = vector.length - 1;
            }

            return vector[idx];
        }

        if (mode == ImportanceMode.HISTORICAL_TASKS) {
            int end = currentTaskIndex;

            if (end <= 0) {
                return vector[0];
            }

            double sum = 0.0;
            int count = 0;

            for (int i = 0; i < end && i < vector.length; i++) {
                sum += vector[i];
                count++;
            }

            return count == 0 ? 0.0 : sum / count;
        }

        double sum = 0.0;

        for (double v : vector) {
            sum += v;
        }

        return sum / vector.length;
    }

    private double calculateTaskImportance(
            GPNode node,
            GPNode originalNode,
            GPTree fullTree,
            PhenoCharacterisation[] pc,
            int treeIndex) {

        try {

            List<?> decisionSituations;

            RuleType ruleType;

            if (treeIndex == 0) {

                decisionSituations =
                        ((SequencingPhenoCharacterisation) pc[0])
                                .decisionSituations;

                ruleType =
                        RuleType.SEQUENCING;

            } else {

                decisionSituations =
                        ((RoutingPhenoCharacterisation) pc[1])
                                .decisionSituations;

                ruleType =
                        RuleType.ROUTING;
            }

            if (decisionSituations == null
                    || decisionSituations.isEmpty()) {
                return 0.0;
            }

//             subtree mean
            double subtreeMean =
                    calculateSubtreeMean(
                            node,
                            decisionSituations,
                            ruleType);

//            double subtreeMean = 1.0;

            // replace subtree
            GPTree modifiedTree =
                    replaceNodeWithMeanConstant(
                            fullTree,
                            node,
                            originalNode,
                            subtreeMean);

            double weightedSimilarity = 0.0;

            for (Object situation : decisionSituations) {

                double[] originalRanking = evaluateRanking(fullTree, situation, ruleType);
                double[] modifiedRanking = evaluateRanking(modifiedTree, situation, ruleType);

                if (originalRanking.length <= 1
                        || modifiedRanking.length <= 1) {
                    continue;
                }

                double corr =
                        new SpearmansCorrelation()
                                .correlation(
                                        originalRanking,
                                        modifiedRanking);

                if (Double.isNaN(corr)) {
                    corr = 0.0;
                }

                weightedSimilarity += corr;
            }

            double similarity =
                    weightedSimilarity / decisionSituations.size();

            double importance =
                    1.0 - similarity;

            lastSelectedBehaviourChange =
                    importance;

            if (Double.isNaN(importance)
                    || Double.isInfinite(importance)) {
                return 0.0;
            }

            return Math.max(0.0, importance);

        } catch (Exception e) {

            return 0.0;
        }
    }

    private GPNode sampleByImportance(
            EvolutionState state,
            int thread,
            GPTree tree,
            double[] scores) {

        double[] probs = new double[scores.length];

        double sum = 0.0;

        double maxScore = Arrays.stream(scores).max().getAsDouble();

        for (int i = 1; i < scores.length; i++) {

            double score = scores[i];

            double prob;

            if (pickHighImportance) {

                prob = Math.max(score, 1e-10);

            } else {

                prob = maxScore - score + 1e-10;
            }

            probs[i] = prob;

            sum += prob;
        }

        double r = state.random[thread].nextDouble() * sum;

        double cumulative = 0.0;

        for (int i = 0; i < probs.length; i++) {

            cumulative += probs[i];

            if (r <= cumulative) {

                lastSelectedImportance = scores[i];

                return tree.child.nodeInPosition(i, GPNode.NODESEARCH_NONTERMINALS);
            }
        }

        lastSelectedImportance = scores[scores.length - 1];

        return tree.child.nodeInPosition(scores.length - 1, GPNode.NODESEARCH_NONTERMINALS);
    }


    private double calculateSubtreeMean(GPNode node,
                                        List<?> decisionSituations,
                                        RuleType ruleType) {

        GPTree subtreeTree =
                PopulationUtils.GPNodetoGPTree((GPNode) node.clone());

        AbstractRule rule =
                new GPRule(ruleType, subtreeTree);

        double sum = 0.0;
        int count = 0;

        for (Object obj : decisionSituations) {

            if (ruleType == RuleType.SEQUENCING) {
                SequencingDecisionSituation situation =
                        (SequencingDecisionSituation) obj;

                List<OperationOption> queue =
                        situation.getQueue();

                for (OperationOption op : queue) {
                    sum += rule.priority(
                            op,
                            situation.getWorkCenter(),
                            situation.getSystemState());

                    count++;
                }
            } else {
                RoutingDecisionSituation situation =
                        (RoutingDecisionSituation) obj;

                List<OperationOption> queue =
                        situation.getQueue();

                for (OperationOption op : queue) {
                    sum += rule.priority(
                            op,
                            op.getWorkCenter(),
                            situation.getSystemState());

                    count++;
                }
            }
        }

        return sum / Math.max(count, 1);
    }

    private GPTree replaceNodeWithMeanConstant(
            GPTree fullTree,
            GPNode targetNode,
            GPNode originalNode,
            double meanValue) {

        GPTree modifiedTree =
                (GPTree) fullTree.clone();

        GPNode constantNode =
                createConstantNode(meanValue);

        modifiedTree.child =
                fullTree.child.cloneReplacing(
                        constantNode,
                        originalNode);

        modifiedTree.child.parent =
                modifiedTree;

        modifiedTree.child.argposition =
                0;

        return modifiedTree;
    }


    private double[] evaluateRanking(GPTree tree,
                                     Object situation,
                                     RuleType ruleType) {

        AbstractRule rule =
                new GPRule(ruleType, tree);

        List<OperationOption> queue;
        WorkCenter workCenter;
        SystemState systemState;

        if (ruleType == RuleType.SEQUENCING) {
            SequencingDecisionSituation s =
                    (SequencingDecisionSituation) situation;

            queue = s.getQueue();
            workCenter = s.getWorkCenter();
            systemState = s.getSystemState();

        } else {
            RoutingDecisionSituation s =
                    (RoutingDecisionSituation) situation;

            queue = s.getQueue();

            if (queue == null || queue.isEmpty()) {
                return new double[0];
            }

            workCenter = queue.get(0).getWorkCenter();
            systemState = s.getSystemState();
        }

        double[] priorities =
                new double[queue.size()];

        for (int i = 0; i < queue.size(); i++) {
            OperationOption op = queue.get(i);

            if (ruleType == RuleType.SEQUENCING) {
                priorities[i] =
                        rule.priority(op, workCenter, systemState);
            } else {
                priorities[i] =
                        rule.priority(op, op.getWorkCenter(), systemState);
            }
        }

        return prioritiesToRanks(priorities);
    }

    private double[] prioritiesToRanks(double[] priorities) {
        double[] ranks =
                new double[priorities.length];

        for (int i = 0; i < priorities.length; i++) {
            double rank = 1.0;

            for (int j = 0; j < priorities.length; j++) {
                if (priorities[j] < priorities[i]) {
                    rank++;
                }
            }

            ranks[i] = rank;
        }

        return ranks;
    }

    private GPNode createConstantNode(double value) {

        try {

            // create NEW wrapper
            TerminalERCUniform newNode =
                    new TerminalERCUniform();

            newNode.children =
                    new GPNode[0];

            newNode.parent =
                    null;

            newNode.argposition =
                    0;

            // clone internal ERC
            DoubleERC erc =
                    new DoubleERC();

            erc.children =
                    new GPNode[0];

            erc.value =
                    value;

            // set protected field terminal
            java.lang.reflect.Field field =
                    TerminalERC.class
                            .getDeclaredField("terminal");

            field.setAccessible(true);

            field.set(newNode, erc);

            return newNode;

        } catch (Exception e) {

            throw new RuntimeException(
                    "Failed to create constant node.",
                    e);
        }
    }



}