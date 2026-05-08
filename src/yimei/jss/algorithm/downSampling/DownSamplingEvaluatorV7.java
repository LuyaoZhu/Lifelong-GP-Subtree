package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Fitness;
import ec.Individual;
import ec.Subpopulation;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import ec.simple.SimpleEvaluator;
import ec.simple.SimpleProblemForm;
import yimei.jss.jobshop.Objective;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleevaluation.AbstractEvaluationModel;
import yimei.jss.ruleevaluation.MultipleTreeMultipleRuleEvaluationModel;
import yimei.jss.ruleoptimisation.MultipleTreeRuleOptimizationProblem;
import yimei.jss.simulation.DynamicSimulation;
import yimei.jss.simulation.Simulation;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class DownSamplingEvaluatorV7 extends SimpleEvaluator {

    protected void evalPopChunk(EvolutionState state, int[] numinds, int[] from,
                                int threadnum, SimpleProblemForm p)
    {
        ((ec.Problem)p).prepareToEvaluate(state,threadnum);

        Subpopulation[] subpops = state.population.subpops;
        int len = subpops.length;

        for(int pop=0;pop<len;pop++) {
            // start evaluatin'!
            int fp = from[pop];
            int upperbound = fp + numinds[pop];
            Individual[] inds = subpops[pop].individuals;

            if (state.generation < ((GPRuleEvolutionStateDownSampling) state).switchGen) {
                for (int x=fp;x<upperbound;x++)
                    p.evaluate(state,inds[x], pop, threadnum);
            } else {
                List<Simulation> commonSimulation = new ArrayList<>();
                for (Map.Entry<List<Integer>, List<Simulation>> entry : ((GPRuleEvolutionStateDownSamplingV1) state).simulationsForIndividualsInOneCluster.entrySet()) {
                    if (entry.getValue().size() > 0) {
                        commonSimulation = entry.getValue();
                        break;
                    }
                }
                for (Map.Entry<List<Integer>, List<Simulation>> entry : ((GPRuleEvolutionStateDownSamplingV1) state).simulationsForIndividualsInOneCluster.entrySet()) {
                    for (Integer ind : entry.getKey()) {
                        if (entry.getValue().size() > 0) {
                            evaluate(state, inds[ind], pop, threadnum, entry.getValue());
                        } else {
                            evaluate(state, inds[ind], pop, threadnum, commonSimulation);
                        }
                    }
                }
            }
        }

        ((ec.Problem)p).finishEvaluating(state,threadnum);
    }

    public void evaluate(EvolutionState state, Individual indi, int subpopulation, int threadnum, List<Simulation> simulations) {

        //GPRule rule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) indi).trees[0]);

        //modified by fzhang 23.5.2018  read two rules from one individual
        GPRule sequencingRule = new GPRule(RuleType.SEQUENCING, ((GPIndividual) indi).trees[0]);
        GPRule routingRule = new GPRule(RuleType.ROUTING, ((GPIndividual) indi).trees[1]);

        List rules = new ArrayList();
        List fitnesses = new ArrayList();

        //rules.add(rule);
        //modified by fzhang  to save two rules for evaluating from one individual
        rules.add(sequencingRule);
        rules.add(routingRule);

        if(((GPIndividual) indi).trees.length == 3){
            GPRule orderingRule = new GPRule(RuleType.ORDERING, ((GPIndividual) indi).trees[2]);
            rules.add(orderingRule);
        }

        fitnesses.add(indi.fitness);

        AbstractEvaluationModel evaluationModel = ((MultipleTreeRuleOptimizationProblem)state.evaluator.p_problem).getEvaluationModel();

        evaluate(fitnesses, rules, state, evaluationModel, simulations);

        indi.evaluated = true;
    }

    public void evaluate(List<Fitness> currentFitnesses,
                         List<AbstractRule> rules,
                         EvolutionState state,
                         AbstractEvaluationModel evaluationModel,
                         List<Simulation> finalSelectedSimulations) {
        //expecting 2 rules here - one routing rule and one sequencing rule
/*        if (rules.size() != 2) {
            System.out.println("Rule evaluation failed!");
            System.out.println("Expecting 2 rules, only 1 found.");
            return;
        }*/

        //System.out.println(rules.size()); //2 repeat

        List<Objective> objectives = evaluationModel.getObjectives();

        AbstractRule sequencingRule = rules.get(0); // for each arraylist in list, they have two elements, the first one is sequencing rule and the second one is routing rule
        AbstractRule routingRule = rules.get(1);

        //System.out.println(sequencingRule);  //"GPRule"  repeat
        //System.out.println(routingRule);   //"GPRule"  repeat

        //System.out.println(objectives.size()); //1  repeat
        //code taken from Abstract Rule
        double[] fitnesses = new double[objectives.size()];

        ArrayList<Double> fitnessesList = new ArrayList<>();

        List<Simulation> simulations = new ArrayList<>();

        if(state.generation < ((GPRuleEvolutionStateDownSampling)state).switchGen){
            simulations = ((MultipleTreeMultipleRuleEvaluationModel)evaluationModel).getSchedulingSet().getSimulations();
        } else {
            for (int i = 0; i < finalSelectedSimulations.size(); i++) {
                simulations.add(((DynamicSimulation) finalSelectedSimulations.get(i)).deepCloneAllSame());
            }
        }

        int col = 0;

        //System.out.println(simulations.size()); // 1 repeat
        //System.out.println(schedulingSet.getReplications().get(0)); //1 repeat

        for (int j = 0; j < simulations.size(); j++) {
            Simulation simulation = simulations.get(j);
//            ((DynamicSimulation)simulation).reseed(list.get(state.generation).getKey());

            //========================change here======================================
            if(state.generation >= ((GPRuleEvolutionStateDownSampling)state).switchGen){
                simulation.setBothRuleWithoutReset(sequencingRule,routingRule);
            } else {
                simulation.setSequencingRule(sequencingRule); //indicate different individuals
                simulation.setRoutingRule(routingRule);
            }

            if (rules.size() == 3) {
                AbstractRule orderingRule = rules.get(2);
                simulation.setOrderingRule(orderingRule);
            }
            //System.out.println(simulation);
            simulation.run();

/*            if ((!(simulation.getClockTime() == Double.MAX_VALUE)) && simulations.size() > 1) {
                for (Batch batch : simulation.getSystemState().getBatchesCompleted()) {
                    if (batch.getId() < 0.2 * simulation.getBatchesRecorded() + simulation.getWarmupBatches()) {
                        for (Job job : batch.getJobs()) {
                            simulation.getSystemState().getJobsCompleted().remove(job);
                        }
                    }
                }
            }*/

//            if (!(simulation.getClockTime() == Double.MAX_VALUE)) {
//                double utilisation = 0;
//                for (WorkCenter machine : simulation.getSystemState().getWorkCenters()) {
//                    utilisation += machine.getBusyTime() / simulation.getClockTime();
//                }
//                System.out.println("utilisation: " + utilisation / 10);
//                System.out.println("The number of completed jobs" + " is " + simulation.getSystemState().getJobsCompleted().size());
//            }

            for (int i = 0; i < objectives.size(); i++) {
                //2018.10.23  cancel normalized process
//                double normObjValue = simulation.objectiveValue(objectives.get(i))  // this line: the value of makespan
//                        / schedulingSet.getObjectiveLowerBound(i, col);

                double ObjValue = simulation.objectiveValue(objectives.get(i)); // this line: the value of makespan



                //in essence, here is useless. because if w.numOpsInQueue() > 100, the simulation has been canceled in run(). here is a double check
                    if (simulation.getSystemState().getClockTime() ==  Double.MAX_VALUE) {
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
                    }


                //2018.10.23  cancel normalized process
//                fitnesses[i] += normObjValue;  //the value of fitness is the normalization of the objective value
                fitnesses[i] += ObjValue;
                fitnessesList.add(ObjValue);  //add the case fitness to this list

            }
            col++;

            simulation.reset();
        }


            for (int i = 0; i < fitnesses.length; i++) {
                fitnesses[i] /= col;
                if (fitnesses[i] >= Double.POSITIVE_INFINITY || fitnesses[i] <= Double.NEGATIVE_INFINITY){
                    fitnesses[i] = Double.MAX_VALUE;
                }
            }




        //System.out.println(currentFitnesses.size()); //1
        for (Fitness fitness : currentFitnesses) {
            MultiObjectiveFitness f = (MultiObjectiveFitness) fitness;
            f.setObjectives(state, fitnesses);
        }

        currentFitnesses.get(0).trials = fitnessesList;

    }

}
