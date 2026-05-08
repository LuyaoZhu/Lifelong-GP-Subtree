package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Fitness;
import ec.multiobjective.MultiObjectiveFitness;
import yimei.jss.jobshop.WorkCenter;
import yimei.jss.rule.AbstractRule;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.simulation.Simulation;

import java.util.ArrayList;
import java.util.List;

public class MultipleTreeMultipleRuleEvaluationModelSubCasesV5 extends MultipleTreeMultipleRuleEvaluationModel {

    public int caseNum;

    @Override
    public void evaluate(List<Fitness> currentFitnesses,
                         List<AbstractRule> rules,
                         EvolutionState state) {
        //expecting 2 rules here - one routing rule and one sequencing rule
/*        if (rules.size() != 2) {
            System.out.println("Rule evaluation failed!");
            System.out.println("Expecting 2 rules, only 1 found.");
            return;
        }*/

        //System.out.println(rules.size()); //2 repeat
        countInd++;

        AbstractRule sequencingRule = rules.get(0); // for each arraylist in list, they have two elements, the first one is sequencing rule and the second one is routing rule
        AbstractRule routingRule = rules.get(1);

        //System.out.println(sequencingRule);  //"GPRule"  repeat
        //System.out.println(routingRule);   //"GPRule"  repeat

        //System.out.println(objectives.size()); //1  repeat
        //code taken from Abstract Rule
        double[] fitnesses = new double[objectives.size()];

        ArrayList<Double> fitnessesList = new ArrayList<>();

        List<Simulation> simulations = schedulingSet.getSimulations();

        int col = 0;

        if(state.generation >= ((GPRuleEvolutionStateDownSampling)state).switchGen) {
            caseNum = (((GPRuleEvolutionStateDownSamplingV5) state).caseIndex + 1);
        } else {
            caseNum = ((GPRuleEvolutionStateDownSamplingV5) state).caseInOneInstance;
        }


        //System.out.println(simulations.size()); // 1 repeat
        //System.out.println(schedulingSet.getReplications().get(0)); //1 repeat

        for (int j = 0; j < simulations.size(); j++) {
            Simulation simulation = simulations.get(j);
//            ((DynamicSimulation)simulation).reseed(list.get(state.generation).getKey());

                simulation.numBatchesRecorded = caseNum * ((GPRuleEvolutionStateDownSampling)state).batchInOneCase;

                simulation.setSequencingRule(sequencingRule); //indicate different individuals
                simulation.setRoutingRule(routingRule);

            if (rules.size() == 3) {
                AbstractRule orderingRule = rules.get(2);
                simulation.setOrderingRule(orderingRule);
            }
            //System.out.println(simulation);
            simulation.run();

            for (int i = 0; i < objectives.size(); i++) {
                //2018.10.23  cancel normalized process
//                double normObjValue = simulation.objectiveValue(objectives.get(i))  // this line: the value of makespan
//                        / schedulingSet.getObjectiveLowerBound(i, col);

                double ObjValue = simulation.objectiveValue(objectives.get(i)); // this line: the value of makespan
                double[] caseObjective = simulation.objectiveValueMultiCase(objectives.get(i),caseNum);

                //in essence, here is useless. because if w.numOpsInQueue() > 100, the simulation has been canceled in run(). here is a double check
                for (WorkCenter w : simulation.getSystemState().getWorkCenters()) {
                    if (w.numOpsInQueue() > 100) {
                        //this was a bad run

                        //fzhang cancel normalized process
//                      normObjValue = Double.MAX_VALUE;
                        if(objectives.get(0).getName().endsWith("profit"))
                            ObjValue = -Double.MAX_VALUE;
                        else
                            ObjValue = Double.MAX_VALUE;

                        //System.out.println(systemState.getJobsInSystem().size());
                        //System.out.println(systemState.getJobsCompleted().size());

                        //normObjValue = normObjValue*(systemState.getJobsInSystem().size()/systemState.getJobsCompleted().size());
                        countBadrun++;
                        break;
                    }
                }

                //2018.10.23  cancel normalized process
//                fitnesses[i] += normObjValue;  //the value of fitness is the normalization of the objective value
                fitnesses[i] += ObjValue;
                for (int caseIndex =0; caseIndex < caseObjective.length; caseIndex++) {
                    fitnessesList.add(caseObjective[caseIndex]);  //add the case fitness to this list
                }
            }
            col++;

            simulation.reset();
        }


            for (int i = 0; i < fitnesses.length; i++) {
                fitnesses[i] /= col;
            }

        //System.out.println(currentFitnesses.size()); //1
        for (Fitness fitness : currentFitnesses) {
            MultiObjectiveFitness f = (MultiObjectiveFitness) fitness;
            f.setObjectives(state, fitnesses);
        }

        currentFitnesses.get(0).trials = fitnessesList;


        //modified by fzhang 23.5.2018  save bad run information for one population
        //if(countInd % 1024 == 0) {
   /*     if(countInd % state.population.subpops[0].individuals.length == 0) {
            genNumBadRun.add(countBadrun);
            countBadrun = 0;
           }

        //if(countInd == 1024*51)
        if(countInd == state.population.subpops[0].individuals.length*state.numGenerations)
            WriteCountBadrun(state,null);*/
    }
}
