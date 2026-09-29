package yimei.jss.algorithm.lifelongGPSubtree;

import ec.*;
import ec.simple.SimpleBreeder;
import ec.gp.GPIndividual;

public class LifelongSimpleBreeder
        extends SimpleBreeder {

    @Override
    protected void breedPopChunk(
            Population newpop,
            EvolutionState state,
            int[] numinds,
            int[] from,
            int threadnum) {

        for (int subpop = 0;
             subpop < newpop.subpops.length;
             subpop++) {

            if (!shouldBreedSubpop(
                    state,
                    subpop,
                    threadnum)) {

                for (int ind = from[subpop];
                     ind < numinds[subpop] - from[subpop];
                     ind++) {

                    newpop.subpops[subpop]
                            .individuals[ind] =
                            state.population
                                    .subpops[subpop]
                                    .individuals[ind];
                }

            } else {

                BreedingPipeline bp;

                if (clonePipelineAndPopulation) {

                    bp =
                            (BreedingPipeline)
                                    newpop.subpops[subpop]
                                            .species
                                            .pipe_prototype
                                            .clone();

                } else {

                    bp =
                            (BreedingPipeline)
                                    newpop.subpops[subpop]
                                            .species
                                            .pipe_prototype;
                }

                // =====================================
                // configure lifelong pipeline
                // =====================================

                configureLifelongPipeline(
                        bp,
                        state);

                // =====================================
                // normal ECJ logic
                // =====================================

                if (!bp.produces(
                        state,
                        newpop,
                        subpop,
                        threadnum)) {

                    state.output.fatal(
                            "Breeding pipeline does not produce expected species.");
                }

                bp.prepareToProduce(
                        state,
                        subpop,
                        threadnum);

                int x = from[subpop];

                int upperbound =
                        from[subpop]
                                + numinds[subpop];

                while (x < upperbound) {

                    x += bp.produce(
                            1,
                            upperbound - x,
                            x,
                            subpop,
                            newpop.subpops[subpop]
                                    .individuals,
                            state,
                            threadnum);
                }

                if (x > upperbound) {

                    state.output.fatal(
                            "Breeding pipeline overwrote another pipeline.");
                }

                bp.finishProducing(
                        state,
                        subpop,
                        threadnum);
            }
        }
    }

    // =====================================================
    // recursively configure lifelong pipelines
    // =====================================================

    protected void configureLifelongPipeline(
            BreedingSource source,
            EvolutionState state) {

        if (source == null) {
            return;
        }

        // =================================================
        // our lifelong crossover
        // =================================================

        if (source instanceof
                ReplaceBasedAllCrossoverPipeline) {

            ReplaceBasedAllCrossoverPipeline p =
                    (ReplaceBasedAllCrossoverPipeline)
                            source;

            GPRuleEvolutionStateLifelongGPReplace s =
                    (GPRuleEvolutionStateLifelongGPReplace)
                            state;

            // =============================================
            // decision situations
            // =============================================

            p.setDecisionSituationsEachTask(
                    s.decisionSituationsEachTask);

            // =============================================
            // current task
            // =============================================

            p.setCurrentTaskIndex(
                    s.decisionSituationsEachTask.size() - 1);

            // =============================================
            // archives
            // =============================================

            p.setHistoricalGeneralists(s.historicalGeneralists);

            p.setCurrentSpecialists(s.currentSpecialists);

            p.setCurrentGeneralists(s.currentGeneralists);

            // =============================================
            // importance mode
            // =============================================

            p.setReplaceMode(
                    TaskwiseImportanceNodeSelector
                            .ImportanceMode.ALL_TASKS);

            p.setDonorMode(
                    TaskwiseImportanceNodeSelector
                            .ImportanceMode.ALL_TASKS);
        }

        // =================================================
        // recursively configure child pipelines
        // =================================================

        if (source instanceof BreedingPipeline) {

            BreedingPipeline bp =
                    (BreedingPipeline) source;

            if (bp.sources != null) {

                for (int i = 0;
                     i < bp.sources.length;
                     i++) {

                    configureLifelongPipeline(
                            bp.sources[i],
                            state);
                }
            }
        }
    }

    // =====================================================
    // helper
    // =====================================================

    protected java.util.ArrayList<GPIndividual>
    castIndividuals(
            java.util.ArrayList<? extends Individual>
                    inds) {

        java.util.ArrayList<GPIndividual>
                result =
                new java.util.ArrayList<>();

        if (inds == null) {
            return result;
        }

        for (Individual ind : inds) {

            if (ind instanceof GPIndividual) {

                result.add(
                        (GPIndividual) ind);
            }
        }

        return result;
    }
}