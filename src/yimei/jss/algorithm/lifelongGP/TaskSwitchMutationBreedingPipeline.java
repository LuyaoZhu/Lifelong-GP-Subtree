package yimei.jss.algorithm.lifelongGP;

import ec.BreedingSource;
import ec.EvolutionState;
import ec.Individual;
import ec.SelectionMethod;
import ec.breed.MultiBreedingPipeline;

public class TaskSwitchMutationBreedingPipeline extends MultiBreedingPipeline {

    @Override
    public int produce(final int min,
                       final int max,
                       final int start,
                       final int subpopulation,
                       final Individual[] inds,
                       final EvolutionState state,
                       final int thread)

    {
        BreedingSource s;

        // ===== 全突变代 =====
        // 假设 sources[1] 是 mutation pipeline
        if(state.generation == ((GPRuleEvolutionStateLifelongGPV10N1)state).generationPerTask-1)
        {
            s = sources[1];
        }
        // ===== 正常代 =====
        else
        {
            s = sources[BreedingSource.pickRandom(
                    sources,
                    state.random[thread].nextDouble())];
        }

        int total;

        if (generateMax)
        {
            if (maxGeneratable == 0)
                maxGeneratable = maxChildProduction();

            int n = maxGeneratable;

            if (n < min) n = min;
            if (n > max) n = max;

            total = s.produce(
                    n,
                    n,
                    start,
                    subpopulation,
                    inds,
                    state,
                    thread);
        }
        else
        {
            total = s.produce(
                    min,
                    max,
                    start,
                    subpopulation,
                    inds,
                    state,
                    thread);
        }

        // clone if necessary
        if (s instanceof SelectionMethod)
            for(int q = start; q < total + start; q++)
                inds[q] = (Individual)(inds[q].clone());

        return total;
    }

}
