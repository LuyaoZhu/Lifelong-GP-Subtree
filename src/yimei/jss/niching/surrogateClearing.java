package yimei.jss.niching;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import ec.gp.GPIndividual;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;

/**
2019.8.24 calculate the (phenotypy characteristic) PC of ech individuals in each population
 not used in weighted feature idea
 */
public class surrogateClearing {
    static double[][][] indsCharLists = null; //save the PC information

    public static double[][][] clearPopulation(final EvolutionState state,
                                       PhenoCharacterisation[] pc) {

        RuleType[] ruleTypes = new RuleType[pc.length];
        if(state.population.subpops.length == 2) {
            ruleTypes = new RuleType[]{RuleType.SEQUENCING, RuleType.ROUTING}; //ruleType is an array}
        }else if(state.population.subpops.length == 3) {
            ruleTypes = new RuleType[]{RuleType.SEQUENCING, RuleType.ROUTING, RuleType.ORDERING}; //ruleType is an array}
        }
        indsCharLists = new double[state.population.subpops.length][(int)(state.population.subpops[0].individuals.length)][];

        for (int subpopNum = 0; subpopNum < state.population.subpops.length; subpopNum++) {//each population

            Subpopulation subpop = state.population.subpops[subpopNum];
            RuleType ruleType = ruleTypes[subpopNum];  //ruleType is a rule type---ruleType[0] = SEQUENCING  ruleType[1] = ROUTING
            PhenoCharacterisation phenoCharacterisation = pc[subpopNum];//fzhang 2018.10.02  define two phenotype characteristic---phenoCharacterisation
            Individual[] inds = subpop.individuals; //should according to the fitness value
            phenoCharacterisation.setReferenceRule(new GPRule(ruleType,((GPIndividual)inds[0]).trees[0]));

            //each individuals
            for(int ind = 0; ind < (int)(inds.length); ind ++){
                int[] charList = phenoCharacterisation.characterise(  //.characterise: calculate the distance
                        new GPRule(ruleType,((GPIndividual)inds[ind]).trees[0]));
                indsCharLists[subpopNum][ind] = new double[charList.length];

                //each PC information---convert int[] to int[][]
                for(int numFeature = 0; numFeature < charList.length; numFeature ++){
                    indsCharLists[subpopNum][ind][numFeature] = charList[numFeature];
                }
            }
        }
        return indsCharLists;
    }
}
