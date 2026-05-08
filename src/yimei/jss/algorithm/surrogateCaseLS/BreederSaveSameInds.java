package yimei.jss.algorithm.surrogateCaseLS;

import ec.EvolutionState;
import ec.Individual;
import ec.Population;
import ec.simple.SimpleBreeder;
import ec.util.Parameter;
import ec.util.QuickSort;

public class BreederSaveSameInds extends SimpleBreeder {

    protected void loadElites(EvolutionState state, Population newpop)
    {
        // are our elites small enough?
        for(int x=0;x<state.population.subpops.length;x++)
        {
            if (numElites(state, x)>state.population.subpops[x].individuals.length)
                state.output.error("The number of elites for subpopulation " + x + " exceeds the actual size of the subpopulation",
                        new Parameter(EvolutionState.P_BREEDER).push(P_ELITE).push(""+x));
            if (numElites(state, x)==state.population.subpops[x].individuals.length)
                state.output.warning("The number of elites for subpopulation " + x + " is the actual size of the subpopulation",
                        new Parameter(EvolutionState.P_BREEDER).push(P_ELITE).push(""+x));
        }
        state.output.exitIfErrors();

        //=================find out the index of the best indivudual=======================
        // we assume that we're only grabbing a small number (say <10%), so
        // it's not being done multithreaded
        int[] bestIndex = new int[state.population.subpops.length];
        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            int best = 0;
            Individual[] oldinds = state.population.subpops[sub].individuals;
            for (int x = 1; x < oldinds.length; x++) {
                if (oldinds[x].fitness.betterThan(oldinds[best].fitness)) {
                    best = x;
                }
            }
            bestIndex[sub] = best;
            //the index of best individuals in each subpopulation are saved in bestIndex[sub], like bestIndex[0] = 0, bestIndex[1]=145
        }

        for(int sub=0;sub<state.population.subpops.length;sub++) //0  1
        {
            //System.out.println(!shouldBreedSubpop(state, sub, 0));  //nothing is printed out.
            //skip here
            if (!shouldBreedSubpop(state, sub, 0))  // don't load the elites for this one, we're not doing breeding of it  true
            {
                continue;
            }

            //System.out.println("numElites(state, sub)  "+numElites(state, sub) );  //always eauql two, seems like no related to eval.num-elites and breed.elite.0
            // if the number of elites is 1, then we handle this by just finding the best one.
            if (numElites(state, sub) == 1) {
                Individual[] oldinds = state.population.subpops[sub].individuals; // sub = 0,   1
                Individual[] inds = newpop.subpops[sub].individuals; // null
                if (state.population.subpops.length > 1) {
                    int otherSubPop = (sub+1)%2;  // otherSubPop = 1,  0
                    Individual[] oldindsOtherSubpop = state.population.subpops[otherSubPop].individuals;
                    //want to also insert context of best individual   save the best
                    //Auxiliary variable, used by coevolutionary processes, to store the individuals
                    //involved in producing this given Fitness value.
                    Individual otherCollab = (Individual) oldindsOtherSubpop[bestIndex[otherSubPop]].fitness.getContext()[sub].clone();
                    inds[inds.length-2] = otherCollab;
                }
                Individual elite = (Individual)(oldinds[bestIndex[sub]].clone());
                inds[inds.length-1] = elite;
            }
            else if (numElites(state, sub)>0)  // we'll need to sort
            {
                //define int[] orderPop, length = 512 and its elements are from 0 to 511
                int[] orderedPop = new int[state.population.subpops[sub].individuals.length];
                for(int x=0;x<state.population.subpops[sub].individuals.length;x++)
                    orderedPop[x] = x;
                //orderPop[0]= 0, orderPop[1]= 1, orderPop[2]= 2....orderPop[511]= 511

                // sort the best so far where "<" means "not as fit as"
                QuickSort.qsort(orderedPop, new EliteComparator(state.population.subpops[sub].individuals));
                // load the top N individuals

                Individual[] inds = newpop.subpops[sub].individuals; // has not value
                Individual[] oldinds = state.population.subpops[sub].individuals; //has values
                for(int x=inds.length-numElites(state, sub);x<inds.length;x++) {//start from 510, because numElites(state,sub)
                    inds[x] = (Individual) (oldinds[orderedPop[x]].clone());
                    ((GPRuleEvolutionStateSavedSurrogateCasesV6)state).individualBeforeFitness.add(inds[x].fitness.fitness());
                    ((GPRuleEvolutionStateSavedSurrogateCasesV6)state).sameIndividualsIndex.add(x);
                }
            }
        }

        // optionally force reevaluation
        unmarkElitesEvaluated(state, newpop);
    }


}
