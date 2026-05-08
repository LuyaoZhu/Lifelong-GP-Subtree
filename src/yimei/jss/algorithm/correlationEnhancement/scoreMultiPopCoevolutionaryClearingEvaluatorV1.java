package yimei.jss.algorithm.correlationEnhancement;

import ec.EvolutionState;
import ec.coevolve.GroupedProblemForm;
import ec.coevolve.MultiPopCoevolutionaryEvaluator;
import ec.util.Parameter;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;

/**
 * Created by fzhang on 2019.9.7.
 * get the decision situation from phenoCharacterisation
 * <p>
 * !!!!!!!!!This class only change one thing. If there are two times evaluation in the same generation, we'll not change the elites.
 */
//a class can not extend from more than one class
public class scoreMultiPopCoevolutionaryClearingEvaluatorV1 extends MultiPopCoevolutionaryEvaluator {

    //fzhang 2018.10.9 to get the pre-generation value
    public static final String P_PRE_GENERATIONS = "pre-generations";

    protected boolean clear = true;

    protected static PhenoCharacterisation[] phenoCharacterisation;

    public static PhenoCharacterisation[] getPhenoCharacterisation() {
        return phenoCharacterisation;
    }

    public PhenoCharacterisation getPhenoCharacterisation(int index) {
        return phenoCharacterisation[index];
    }

    public void setup(final EvolutionState state, final Parameter base) {
        super.setup(state, base);

        String filePath = state.parameters.getString(new Parameter("filePath"), null);
        //It's a little tricky to know whether we have 1 or 2 populations here, so we will assume
        //2 for the purpose of the phenoCharacterisation, and ignore the second object if only
        //1 is used
        phenoCharacterisation = new PhenoCharacterisation[2];
        if (filePath == null) {
            //dynamic simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation();
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation();
        } else {
            //static simulation
            phenoCharacterisation[0] =
                    SequencingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
            phenoCharacterisation[1] =
                    RoutingPhenoCharacterisation.defaultPhenoCharacterisation(filePath);
        }
    }

    public void setClear(boolean clear) {
        this.clear = clear;
    }


    public void evaluatePopulation(final EvolutionState state)
    {
        // determine who needs to be evaluated
        boolean[] preAssessFitness = new boolean[state.population.subpops.length];
        boolean[] postAssessFitness = new boolean[state.population.subpops.length];
        for(int i = 0; i < state.population.subpops.length; i++)
        {
            postAssessFitness[i] = shouldEvaluateSubpop(state, i, 0);
            //System.out.println(shouldEvaluateSubpop(state, i, 0));  //true  true
            preAssessFitness[i] = postAssessFitness[i] || (state.generation == 0);  // always prepare (set up trials) on generation 0
        }

        // do evaluation
        beforeCoevolutionaryEvaluation( state, state.population, (GroupedProblemForm)p_problem );

        ((GroupedProblemForm)p_problem).preprocessPopulation(state,state.population, preAssessFitness, false);
        performCoevolutionaryEvaluation( state, state.population, (GroupedProblemForm)p_problem );

        //change the elites for cooperative coevolution
        ((GroupedProblemForm)p_problem).postprocessPopulation(state, state.population, postAssessFitness, false);


        if (((EnhanceCorrelationGPRuleEvolutionStateTest) state).changeElites)
            afterCoevolutionaryEvaluation(state, state.population, (GroupedProblemForm) p_problem);//change the elites
        ((EnhanceCorrelationGPRuleEvolutionStateTest) state).changeElites = true;
    }


}
