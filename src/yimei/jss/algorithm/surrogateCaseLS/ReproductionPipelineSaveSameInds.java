package yimei.jss.algorithm.surrogateCaseLS;

import ec.EvolutionState;
import ec.Individual;
import ec.SelectionMethod;
import ec.breed.ReproductionPipeline;

public class ReproductionPipelineSaveSameInds extends ReproductionPipeline {


    public int produce(
            final int min,
            final int max,
            final int start,
            final int subpopulation,
            final Individual[] inds,
            final EvolutionState state,
            final int thread)
    {
        // grab individuals from our source and stick 'em right into inds.
        // we'll modify them from there
        int n = sources[0].produce(min,max,start,subpopulation,inds,state,thread);

        if (mustClone || sources[0] instanceof SelectionMethod)
            for(int q=start; q < n+start; q++){
                inds[q] = (Individual)(inds[q].clone());
                ((GPRuleEvolutionStateSavedSurrogateCasesV6)state).individualBeforeFitness.add(inds[q].fitness.fitness());
                ((GPRuleEvolutionStateSavedSurrogateCasesV6)state).sameIndividualsIndex.add(q);
            }
        return n;
    }

}
