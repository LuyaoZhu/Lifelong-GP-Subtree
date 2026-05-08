package yimei.jss.algorithm.surrogateCaseLS;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;

import java.util.List;

public class PCEvaluator extends SimpleEvaluator {

    protected void evalPopChunk(EvolutionState state, int[] numinds, int[] from,
                                int threadnum, SimpleProblemForm p)
    {
        ((ec.Problem)p).prepareToEvaluate(state,threadnum);

        Subpopulation[] subpops = state.population.subpops;
        int len = subpops.length;

        for(int pop=0;pop<len;pop++)
        {
            // start evaluatin'!
            int fp = from[pop];
            int upperbound = fp+numinds[pop];
            Individual[] inds = subpops[pop].individuals;

            if(state.generation >= ((GPRuleEvolutionStateSavedSurrogateCasesV6)state).switchGen) {
                //real evaluate one ind with the same PC
                for (List<Integer> list : ((GPRuleEvolutionStateSavedSurrogateCasesV6) state).pcToIndexMap.values()) {
                    p.evaluate(state, inds[list.get(0)], pop, threadnum);

                    //then assign the fitness to other inds with the same PC
                    double[] fitnesses = new double[1];
                    fitnesses[0] = inds[list.get(0)].fitness.fitness();

                    if (list.size() > 1) {
                        for (int ind : list) {
                            ((MultiObjectiveFitness) inds[ind].fitness).setObjectives(state, fitnesses);
                            inds[ind].evaluated = true;
                        }
                    }
                }
            } else {
                    for (int x=fp;x<upperbound;x++)
                        p.evaluate(state,inds[x], pop, threadnum);
            }
        }

        ((ec.Problem)p).finishEvaluating(state,threadnum);
    }

}
