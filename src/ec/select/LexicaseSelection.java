package ec.select;

import ec.EvolutionState;
import ec.Individual;
import ec.SelectionMethod;
import ec.util.MersenneTwisterFast;
import ec.util.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

/**
 *
 * @author Eric 'Siggy' Scott
 */
public class LexicaseSelection extends SelectionMethod
    {
    private static final long serialVersionUID = 1;
    
    public static final String P_LEXICASESELECT = "lexicaseselect";

    @Override
    public Parameter defaultBase()
        {
        return SelectDefaults.base().push(P_LEXICASESELECT);
        }

    @Override
    public int produce(final int subpopulation, final EvolutionState state, final int thread)
        {
        assert(state != null);
        assert(subpopulation >= 0);
        assert(subpopulation < state.population.subpops.length);
        assert(state.population.subpops[subpopulation] != null);
        assert(state.population.subpops[subpopulation].individuals.length > 0);
        
        final ArrayList<Individual> pop = new ArrayList<Individual> (Arrays.asList(state.population.subpops[subpopulation].individuals));
        
        // Initialize the candidates to the entire population
        final ArrayList<Integer> candidates = new ArrayList<Integer>();
        for (int i = 0; i < pop.size(); i++)
            candidates.add(i); //fzhang: give each element's value as its index
        
        // Shuffle test cases
        assert(pop.get(candidates.get(0)).fitness.trials != null);
        final int numCases = pop.get(0).fitness.trials.size(); //numCases means how many fitness/objective values
        if (numCases == 0)
            state.output.fatal(String.format("Attempted to use %s on an individual with an empty list of trials.", this.getClass().getSimpleName()));
        final int[] caseOrder = new int[numCases];
        for (int i = 0; i < numCases; i++)
            caseOrder[i] = i;
        shuffle(state, caseOrder); //fzhang: change the order of elements in caseOrder randomly
        
        for (int i = 0; i < caseOrder.length; i++)
            {
            final int currentCase = caseOrder[i];
            
            // Find the best value of the current test case
            //Fitness best = (Fitness) pop.get(candidates.get(0)).fitness.trials.get(currentCase);
            Individual best = pop.get(candidates.get(0));
            double bestFitness = (double) pop.get(candidates.get(0)).fitness.trials.get(currentCase);
            for (int j = 1; j < candidates.size(); j++)
                {
                assert(pop.get(candidates.get(j)).fitness.trials != null);
                assert(pop.get(candidates.get(j)).fitness.trials.size() == numCases);
                //final Fitness caseFitness = (Fitness) pop.get(candidates.get(j)).fitness.trials.get(currentCase);
                final Individual caseFitness = pop.get(candidates.get(j));
                //if (caseFitness.fitness.betterThan(best.fitness))
                if (caseFitness.fitness.fitness() < best.fitness.fitness())
                    //best = caseFitness;
                    bestFitness = caseFitness.fitness.fitness();
                }
            
            // Reduce candidates to the subset that performs best on the current test case
            final Iterator<Integer> it = candidates.iterator();
            while (it.hasNext())
                {
                //final Fitness caseFitness = (Fitness) pop.get(it.next()).fitness.trials.get(currentCase);
                final Individual caseFitness = pop.get(it.next());
                //if (caseFitness.compareTo(best) > 0) // if strictly worse than best
                if (caseFitness.fitness.fitness() > bestFitness)
                    it.remove();
                }
            
            // If only one individual is left, return it
            if (candidates.size() == 1)
                {
                return candidates.get(0);
                }
            
            // If this was the last test case, return a random candidate
            if (i == caseOrder.length - 1)
                {
                return candidates.get(state.random[0].nextInt(candidates.size()));
                }
            }
        throw new IllegalStateException();
        }
    

    private void shuffle(final EvolutionState state, final int[] a)
        {
        final MersenneTwisterFast mtf = state.random[0];
        for(int x = a.length - 1; x >= 1; x--)
            {
            int rand = mtf.nextInt(x+1);
            int obj = a[x];
            a[x] = a[rand];
            a[rand] = obj;
            }
        }
    
    }
