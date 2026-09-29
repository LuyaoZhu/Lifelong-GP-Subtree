package yimei.jss.algorithm.lifelongGPSubtree;

import ec.EvolutionState;
import ec.Individual;
import ec.gp.*;
import yimei.jss.algorithm.multipletreegp.AllIndexCrossoverPipeline;
import yimei.jss.niching.PhenoCharacterisation;

import java.util.ArrayList;

import static yimei.jss.rule.AbstractRuleHelper.state;

public class ReplaceBasedAllCrossoverPipeline extends AllIndexCrossoverPipeline {

    protected TaskwiseImportanceNodeSelector replaceSelectorP1;
    protected TaskwiseImportanceNodeSelector replaceSelectorP2;

    protected TaskwiseImportanceNodeSelector donorSelectorP1;
    protected TaskwiseImportanceNodeSelector donorSelectorP2;

    protected ArrayList<PhenoCharacterisation[]> decisionSituationsEachTask;

    protected ArrayList<GPIndividual> historicalGeneralists;
    protected ArrayList<GPIndividual> currentSpecialists;
    protected ArrayList<GPIndividual> currentGeneralists;

    protected int currentTaskIndex = -1;
    protected int tournamentSize = 7;
    protected int lastArm = -1;

    protected TaskwiseImportanceNodeSelector.ImportanceMode replaceMode =
            TaskwiseImportanceNodeSelector.ImportanceMode.ALL_TASKS;

    protected TaskwiseImportanceNodeSelector.ImportanceMode donorMode =
            TaskwiseImportanceNodeSelector.ImportanceMode.ALL_TASKS;

    public void setDecisionSituationsEachTask(
            ArrayList<PhenoCharacterisation[]> decisionSituationsEachTask) {
        this.decisionSituationsEachTask = decisionSituationsEachTask;
    }

    public void setHistoricalGeneralists(ArrayList<GPIndividual> historicalGeneralists) {
        this.historicalGeneralists = historicalGeneralists;
    }

    public void setCurrentSpecialists(ArrayList<GPIndividual> currentSpecialists) {
        this.currentSpecialists = currentSpecialists;
    }

    public void setCurrentGeneralists(ArrayList<GPIndividual> currentGeneralists) {
        this.currentGeneralists = currentGeneralists;
    }

    public void setCurrentTaskIndex(int currentTaskIndex) {
        this.currentTaskIndex = currentTaskIndex;
    }

    public void setTournamentSize(int tournamentSize) {
        this.tournamentSize = Math.max(1, tournamentSize);
    }

    public void setReplaceMode(TaskwiseImportanceNodeSelector.ImportanceMode replaceMode) {
        this.replaceMode = replaceMode;
    }

    public void setDonorMode(TaskwiseImportanceNodeSelector.ImportanceMode donorMode) {
        this.donorMode = donorMode;
    }

    public int getLastArm() {
        return lastArm;
    }


    @Override
    public int produce(final int min,
                       final int max,
                       final int start,
                       final int subpopulation,
                       final Individual[] inds,
                       final EvolutionState state,
                       final int thread) {

        int n = typicalIndsProduced();

        if (n < min) {
            n = min;
        }

        if (n > max) {
            n = max;
        }

        if (!state.random[thread].nextBoolean(likelihood)) {
            return reproduce(n, start, subpopulation, inds, state, thread, true);
        }

        GPInitializer initializer = (GPInitializer) state.initializer;

        boolean useAdaptiveSelection =
                shouldUseAdaptiveParentSelection(state);

        if (useAdaptiveSelection) {
            if (decisionSituationsEachTask == null || decisionSituationsEachTask.isEmpty()) {
                state.output.fatal(
                        "decisionSituationsEachTask is null or empty in ReplaceBasedAllCrossoverPipeline.");
            }

            initialiseSelectors();
        }

        for (int q = start; q < n + start; ) {

            GPIndividual parent1;
            GPIndividual parent2;

            int arm = -1;

            if (useAdaptiveSelection) {

                arm = AdaptiveParentSelection.selectArm(state, thread);
                lastArm = arm;

                replaceSelectorP1.setMode(modeOfParentFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT1));
                replaceSelectorP1.setCurrentGroupLabel(parentGroupFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT1));

                replaceSelectorP2.setMode(modeOfParentFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT2));
                replaceSelectorP2.setCurrentGroupLabel(parentGroupFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT2));

                donorSelectorP1.setMode(modeOfParentFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT1));
                donorSelectorP1.setCurrentGroupLabel(parentGroupFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT1));

                donorSelectorP2.setMode(modeOfParentFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT2));
                donorSelectorP2.setCurrentGroupLabel(parentGroupFromArm(arm, AdaptiveParentSelection.ParentRole.PARENT2));

                parent1 =
                        (GPIndividual) AdaptiveParentSelection.selectParentFromArm(
                                arm,
                                AdaptiveParentSelection.ParentRole.PARENT1,
                                historicalGeneralists,
                                currentSpecialists,
                                currentGeneralists,
                                state,
                                thread,
                                tournamentSize);

                parent2 =
                        (GPIndividual) AdaptiveParentSelection.selectParentFromArm(
                                arm,
                                AdaptiveParentSelection.ParentRole.PARENT2,
                                historicalGeneralists,
                                currentSpecialists,
                                currentGeneralists,
                                state,
                                thread,
                                tournamentSize);

                if (parent1 == null || parent2 == null) {
                    parent1 = null;
                    parent2 = null;
                    arm = -1;

                    replaceSelectorP1.setCurrentGroupLabel(null);
                    replaceSelectorP2.setCurrentGroupLabel(null);
                    donorSelectorP1.setCurrentGroupLabel(null);
                    donorSelectorP2.setCurrentGroupLabel(null);
                }

            } else {
                parent1 = null;
                parent2 = null;
            }

            // =====================================================
            // Standard parent selection fallback / phase-1 behaviour
            // =====================================================

            if (parent1 == null || parent2 == null) {

                if (sources[0] == sources[1]) {
                    sources[0].produce(
                            2,
                            2,
                            0,
                            subpopulation,
                            parents,
                            state,
                            thread);
                } else {
                    sources[0].produce(
                            1,
                            1,
                            0,
                            subpopulation,
                            parents,
                            state,
                            thread);

                    sources[1].produce(
                            1,
                            1,
                            1,
                            subpopulation,
                            parents,
                            state,
                            thread);
                }

                parent1 = (GPIndividual) parents[0];
                parent2 = (GPIndividual) parents[1];
                arm = -1;
            }

            if (tree1 != TREE_UNFIXED
                    || tree2 != TREE_UNFIXED
                    || parent1.trees.length != parent2.trees.length) {
                state.output.fatal(
                        "GP AllIndexCrossover Pipeline: two individuals chosen for crossover have DIFFERENT numbers of trees, or fixed-tree crossover is not supported.");
            }

            int length = parent1.trees.length;

            GPNode[] replace1 = new GPNode[length];
            GPNode[] donor1 = new GPNode[length];

            GPNode[] replace2 = new GPNode[length];
            GPNode[] donor2 = new GPNode[length];

            // =====================================================
            // Select crossover points
            //
            // Phase 1 / task 0:
            //   standard AllIndexCrossoverPipeline node selection
            //
            // Phase 2:
            //   importance-guided replacement/donor node selection
            // =====================================================

            for (int t = 0; t < length; t++) {

                nodeselect1.reset();
                nodeselect2.reset();

                if (useAdaptiveSelection) {
                    replaceSelectorP1.setTreeIndex(t);
                    donorSelectorP1.setTreeIndex(t);
                    replaceSelectorP2.setTreeIndex(t);
                    donorSelectorP2.setTreeIndex(t);

                    replaceSelectorP1.reset();
                    donorSelectorP1.reset();
                    replaceSelectorP2.reset();
                    donorSelectorP2.reset();

                }

                for (int x = 0; x < numTries; x++) {

                    if (useAdaptiveSelection) {

                        // =====================================
                        // child1
                        //
                        // parent2 donates to parent1
                        // =====================================

                        replace1[t] =
                                replaceSelectorP1.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent1,
                                        parent1.trees[t]);

                        double replaceImportance1 =
                                replaceSelectorP1
                                        .getLastSelectedImportance();

/*                        donor1[t] =
                                nodeselect2.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent2,
                                        parent2.trees[t]);

                        double donorImportance1 = 0;*/

                        donor1[t] =
                                donorSelectorP2.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent2,
                                        parent2.trees[t]);

                        double donorImportance1 =
                                donorSelectorP2
                                        .getLastSelectedImportance();

                        // =====================================
                        // child2
                        //
                        // parent1 donates to parent2
                        // =====================================

                        replace2[t] =
                                replaceSelectorP2.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent2,
                                        parent2.trees[t]);

                        donor2[t] =
                                donorSelectorP1.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent1,
                                        parent1.trees[t]);

/*                        donor2[t] =
                                nodeselect1.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent1,
                                        parent1.trees[t]);*/

                        // =====================================
                        // semantic crossover statistics
                        // =====================================

                        GPRuleEvolutionStateLifelongGPReplace s =
                                (GPRuleEvolutionStateLifelongGPReplace) state;

                        // -------------------------------------
                        // Tree 1 / Tree 2 donor
                        // -------------------------------------

                        s.importanceSum[t * 2] +=
                                donorImportance1;

                        s.importanceCount[t * 2]++;

                        s.depthSum[t * 2] +=
                                donor1[t].depth();

                        s.depthCount[t * 2]++;

                        // -------------------------------------
                        // Tree 1 / Tree 2 replace
                        // -------------------------------------

                        s.importanceSum[t * 2 + 1] +=
                                replaceImportance1;

                        s.importanceCount[t * 2 + 1]++;

                        s.depthSum[t * 2 + 1] +=
                                replace1[t].depth();

                        s.depthCount[t * 2 + 1]++;

                    } else {

                        // =====================================
                        // standard crossover
                        // =====================================

                        replace1[t] =
                                nodeselect1.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent1,
                                        parent1.trees[t]);

                        donor1[t] =
                                nodeselect2.pickNode(
                                        state,
                                        subpopulation,
                                        thread,
                                        parent2,
                                        parent2.trees[t]);

                        replace2[t] = donor1[t];

                        donor2[t] = replace1[t];

                        GPRuleEvolutionStateLifelongGPReplace s =
                                (GPRuleEvolutionStateLifelongGPReplace) state;

                        s.depthSum[t * 2] +=
                                Math.max( donor1[t].depth(), replace1[t].depth());

                        s.depthCount[t * 2]++;

                        s.depthSum[t * 2 + 1] +=
                                Math.min( donor1[t].depth(), replace1[t].depth());

                        s.depthCount[t * 2 + 1]++;

                    }

                    boolean res1 =
                            verifyPoints(
                                    initializer,
                                    donor1[t],
                                    replace1[t]);

                    boolean res2 = true;

                    if (n - (q - start) >= 2
                            && !tossSecondParent) {

                        res2 =
                                verifyPoints(
                                        initializer,
                                        donor2[t],
                                        replace2[t]);
                    }

                    if (res1 && res2) {
                        break;
                    }

                    if (x == numTries - 1) {

                        replace1[t] = null;
                        donor1[t] = null;

                        replace2[t] = null;
                        donor2[t] = null;
                    }
                }
            }

            GPIndividual child1 =
                    createChild(
                            parent1,
                            parent2,
                            replace1,
                            donor1,
                            length,
                            arm);

            inds[q] = child1;
            q++;

            if (q < n + start && !tossSecondParent) {
                GPIndividual child2 =
                        createChild(
                                parent2,
                                parent1,
                                replace2,
                                donor2,
                                length,
                                arm);

                inds[q] = child2;
                q++;
            }
        }

        return n;
    }

    protected boolean shouldUseAdaptiveParentSelection(EvolutionState state) {

        if (!(state instanceof GPRuleEvolutionStateLifelongGPReplace)) {
            return false;
        }

        GPRuleEvolutionStateLifelongGPReplace s =
                (GPRuleEvolutionStateLifelongGPReplace) state;

        // First task or current-task-only phase:
        // use traditional parent selection and traditional node selection.
        if (s.onlyCurrentTaskPhase) {
            return false;
        }

        if (s.decisionSituationsEachTask == null
                || s.decisionSituationsEachTask.size() <= 1) {
            return false;
        }

        if (historicalGeneralists == null || historicalGeneralists.isEmpty()) {
            return false;
        }

        if (currentSpecialists == null || currentSpecialists.isEmpty()) {
            return false;
        }

        return true;
    }

    private void initialiseSelectors() {

        if (replaceSelectorP1 == null) {
            replaceSelectorP1 = new TaskwiseImportanceNodeSelector();
            replaceSelectorP1.setPickHighImportance(false);
        }
        if (replaceSelectorP2 == null) {
            replaceSelectorP2 = new TaskwiseImportanceNodeSelector();
            replaceSelectorP2.setPickHighImportance(false);
        }
        if (donorSelectorP1 == null) {
            donorSelectorP1 = new TaskwiseImportanceNodeSelector();
            donorSelectorP1.setPickHighImportance(true);
        }
        if (donorSelectorP2 == null) {
            donorSelectorP2 = new TaskwiseImportanceNodeSelector();
            donorSelectorP2.setPickHighImportance(true);
        }
        TaskwiseImportanceNodeSelector[] all = {
                replaceSelectorP1, replaceSelectorP2,
                donorSelectorP1, donorSelectorP2
        };

        for (TaskwiseImportanceNodeSelector s : all) {
            s.setDecisionSituationsEachTask(decisionSituationsEachTask);
            s.setCurrentTaskIndex(currentTaskIndex);
        }
    }

    private GPIndividual createChild(GPIndividual recipientParent,
                                     GPIndividual donorParent,
                                     GPNode[] replaceNodes,
                                     GPNode[] donorNodes,
                                     int length,
                                     int arm) {

        GPIndividual child =
                (GPIndividual) recipientParent.lightClone();
        recipientParent.clone();

        child.trees = new GPTree[length];

        for (int t = 0; t < length; t++) {
            child.trees[t] =
                    (GPTree) recipientParent.trees[t].lightClone();

            child.trees[t].owner = child;

            if (replaceNodes[t] != null && donorNodes[t] != null) {
                child.trees[t].child =
                        recipientParent.trees[t]
                                .child
                                .cloneReplacing(
                                        donorNodes[t],
                                        replaceNodes[t]);

                child.evaluated = false;
            } else {
                child.trees[t].child =
                        (GPNode) recipientParent.trees[t]
                                .child
                                .clone();
            }

            child.trees[t].child.parent = child.trees[t];
            child.trees[t].child.argposition = 0;
        }

        child.context = new GPIndividual[2];
        child.context[0] = recipientParent;
        child.context[1] = donorParent;

        child.producedByArm = arm;

        if (arm >= 0 && arm < AdaptiveParentSelection.generatedByArm.length) {
            AdaptiveParentSelection.generatedByArm[arm] += 1;
        }

        child.creationType = CreationType.CROSSOVER;

        child.birthGeneration =
                state.generation;

        return child;
    }

    private TaskwiseImportanceNodeSelector.ImportanceMode modeOfParentFromArm(
            int arm,
            AdaptiveParentSelection.ParentRole role) {

        // 这里按你 AdaptiveParentSelection 里的 arm 定义改
        // 示例：
        // arm 0: HG x CS
        // arm 1: HG x CG
        // arm 2: CS x CG
        // arm 3: CG x CG

        if (arm == 0) {
            return role == AdaptiveParentSelection.ParentRole.PARENT1
                    ? TaskwiseImportanceNodeSelector.ImportanceMode.HISTORICAL_TASKS
                    : TaskwiseImportanceNodeSelector.ImportanceMode.CURRENT_TASK;
        }

        if (arm == 1) {
            return role == AdaptiveParentSelection.ParentRole.PARENT1
                    ? TaskwiseImportanceNodeSelector.ImportanceMode.HISTORICAL_TASKS
                    : TaskwiseImportanceNodeSelector.ImportanceMode.ALL_TASKS;
        }

        if (arm == 2) {
            return role == AdaptiveParentSelection.ParentRole.PARENT1
                    ? TaskwiseImportanceNodeSelector.ImportanceMode.CURRENT_TASK
                    : TaskwiseImportanceNodeSelector.ImportanceMode.ALL_TASKS;
        }

        return TaskwiseImportanceNodeSelector.ImportanceMode.ALL_TASKS;
    }

    private String parentGroupFromArm(
            int arm,
            AdaptiveParentSelection.ParentRole role) {

        if (arm == 0) {

            return role ==
                    AdaptiveParentSelection.ParentRole.PARENT1
                    ? "HG"
                    : "CS";
        }

        if (arm == 1) {

            return role ==
                    AdaptiveParentSelection.ParentRole.PARENT1
                    ? "HG"
                    : "CG";
        }

        if (arm == 2) {

            return role ==
                    AdaptiveParentSelection.ParentRole.PARENT1
                    ? "CS"
                    : "CG";
        }

        return "CG";
    }


}