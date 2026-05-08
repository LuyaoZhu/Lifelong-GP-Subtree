package yimei.jss.niching;

import ec.EvolutionState;
import ec.coevolve.MultiPopCoevolutionaryEvaluator;
import ec.util.Parameter;

/**
 * Created by fzhang on 2019.9.7.
 * get the decision situation from phenoCharacterisation
 */
//a class can not extend from more than one class
public class scoreMultiPopCoevolutionaryClearingEvaluator extends MultiPopCoevolutionaryEvaluator {

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

    @Override
    public void evaluatePopulation(final EvolutionState state) {
        super.evaluatePopulation(state); //fzhang  all the evolution process is the same. for one individual, get a fitness value
    }

    public void setClear(boolean clear) {
        this.clear = clear;
    }
}
