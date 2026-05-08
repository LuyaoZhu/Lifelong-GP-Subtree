package yimei.jss.niching;

import ec.EvolutionState;
import ec.Individual;
import ec.Subpopulation;
import yimei.jss.gp.GPRuleEvolutionState;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static yimei.jss.gp.GPRun.out_dir;

/**
 * The clearing method for niching.
 *
 * Created by fzhang on 2019.9.12.
 */
public class multitreeClearingFirst {

    static ArrayList<Integer> clearedIndsArray = new ArrayList<>();
    protected static long jobSeed;
	// delete the poor individuals      radius: control the range of each niche   capacity: determines the number of individuals in each niche.
    public static void clearPopulation(final EvolutionState state,
                                       double radius, int capacity,
                                       PhenoCharacterisation[] pc) {

        double[][] phenotypicOfIntermedidatePop = phenotypicForSurrogate.phenotypicPopulation(state, pc, false); //fzhang 2019.9.26  add true here to make consistent


        for (int subpopNum = 0; subpopNum < state.population.subpops.length; subpopNum++) {
                int clearedInds = 0;
                Subpopulation subpop = state.population.subpops[subpopNum];
                Individual[] sortedPop = subpop.individuals;

            //fzhang 2018.10.2 calculate the distance of each individual and the reference rule
                List<double[]> sortedPopCharLists = new ArrayList<>();
                for (int i = 0;  i < phenotypicOfIntermedidatePop.length; i++) {
                    //where the examined rule set the chosen operation by reference rule  for example: 3. means examined rule set the chosen operation by reference rule as the third one
                    double[] charList = phenotypicOfIntermedidatePop[i];
                    sortedPopCharLists.add(charList);
                }

                // clear this subpopulation
                for (int i = 0; i < sortedPop.length; i++) {
                    // skip the cleared individuals
                    if (((Clearable)sortedPop[i].fitness).isCleared()) {
                        continue;
                    }

                    int numWinners = 1;
                    for (int j = i+1; j < sortedPop.length; j++) {
                        // skip the cleared individuals
                        if (((Clearable)sortedPop[j].fitness).isCleared()) {
                            continue;
                        }

                        // calculate the distance between individuals i and j
                        double distance = PhenoCharacterisation.distance(
                                sortedPopCharLists.get(i), sortedPopCharLists.get(j));
                        if (distance > radius) {
                            // Individual j is not in the niche
                            continue; //if distance, means two individuals are the same, clear (below) the individual, get out of current loop
                        }

                        if (numWinners < capacity) { //when set capacity to 1, the code will never go to here
                            numWinners ++;
                        }
                        else {
                            // Clear the fitness of individual j
                            ((Clearable)sortedPop[j].fitness).clear();
                            //fzhang 2019.9.11
                            sortedPop[j].evaluated = true;
                            clearedInds++;
                        }
                    }
                }
                clearedIndsArray.add(clearedInds);
        }

        if(state.generation == state.numGenerations-2){ //do not need to breed in the last generation
            jobSeed = ((GPRuleEvolutionState)state).getJobSeed();
            writeToFile(jobSeed, state.numGenerations);
        }
    }

    public static void clearPopulationRanks(final EvolutionState state,
                                       double radius, int capacity,
                                       PhenoCharacterisation[] pc) {

        int[][][] phenotypicOfIntermedidatePop = phenotypicForSurrogate.muchBetterPhenotypicPopulationRanks(state, pc); //fzhang 2019.9.26  add true here to make consistent

        for (int subpopNum = 0; subpopNum < state.population.subpops.length; subpopNum++) {
            int clearedInds = 0;
            Subpopulation subpop = state.population.subpops[subpopNum];
            Individual[] sortedPop = subpop.individuals;

            //fzhang 2018.10.2 calculate the distance of each individual and the reference rule
            List<int[][]> sortedPopCharLists = new ArrayList<>();
            for (int i = 0;  i < phenotypicOfIntermedidatePop.length; i++) {
                //where the examined rule set the chosen operation by reference rule  for example: 3. means examined rule set the chosen operation by reference rule as the third one
                int[][] charList = phenotypicOfIntermedidatePop[i];
                sortedPopCharLists.add(charList);
            }

            // clear this subpopulation
            for (int i = 0; i < sortedPop.length; i++) {
                // skip the cleared individuals
                if (((Clearable)sortedPop[i].fitness).isCleared()) {
                    continue;
                }

                int numWinners = 1;
                for (int j = i+1; j < sortedPop.length; j++) {
                    // skip the cleared individuals
                    if (((Clearable)sortedPop[j].fitness).isCleared()) {
                        continue;
                    }

                    // calculate the distance between individuals i and j
                    double distance = PhenoCharacterisation.spearmanDistance(
                            sortedPopCharLists.get(i), sortedPopCharLists.get(j));
                    if (distance > radius) {
                        // Individual j is not in the niche
                        continue; //if distance, means two individuals are the same, clear (below) the individual, get out of current loop
                    }

                    if (numWinners < capacity) { //when set capacity to 1, the code will never go to here
                        numWinners ++;
                    }
                    else {
                        // Clear the fitness of individual j
                        ((Clearable)sortedPop[j].fitness).clear();
                        //fzhang 2019.9.11
                        sortedPop[j].evaluated = true;
                        clearedInds++;
                    }
                }
            }
            clearedIndsArray.add(clearedInds);
        }

        if(state.generation == state.numGenerations-2){ //do not need to breed in the last generation
            jobSeed = ((GPRuleEvolutionState)state).getJobSeed();
            writeToFile(jobSeed, state.numGenerations);
        }
    }

    public static void writeToFile(long jobSeed, int numGenerations) {
        //fzhang 2019.5.21 save the weight values
        File weightFile = new File(out_dir+"/job." + jobSeed + ".clearedInds.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,ClearedInds");
            writer.newLine();
            for (int i = 0; i < clearedIndsArray.size(); i += 1) { //every two into one generation
                //writer.newLine();
                writer.write(i + "," + clearedIndsArray.get(i) + "\n");
            }
            writer.write(numGenerations-1 + "," + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    //===================================end=================================
}
