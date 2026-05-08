package yimei.jss.algorithm.downSampling;

import ec.EvolutionState;
import ec.Fitness;
import ec.Individual;
import ec.gp.GPIndividual;
import ec.multiobjective.MultiObjectiveFitness;
import yimei.jss.jobshop.Objective;
import yimei.jss.niching.PhenoCharacterisation;
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

public class MultipleTreeRuleOptimizationProblemSubSimulations extends MultipleTreeRuleOptimizationProblem {

	public void evaluate(EvolutionState state, Individual indi, int subpopulation, int threadnum) {

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

		evaluate(fitnesses, rules, state, evaluationModel, ((GPIndividual) indi).PC);

		indi.evaluated = true;
	}

	public void evaluate(List<Fitness> currentFitnesses,
						 List<AbstractRule> rules,
						 EvolutionState state,
						 AbstractEvaluationModel evaluationModel,
						 int[] PC) {


		AbstractRule sequencingRule = rules.get(0); // for each arraylist in list, they have two elements, the first one is sequencing rule and the second one is routing rule
		AbstractRule routingRule = rules.get(1);

		List<Objective> objectives = evaluationModel.getObjectives();

		double[] fitnesses = new double[objectives.size()];

		ArrayList<Double> fitnessesList = new ArrayList<>();

		List<Simulation> simulations = new ArrayList<>();

		if(state.generation < ((GPRuleEvolutionStateDownSampling)state).switchGen){
			simulations = ((MultipleTreeMultipleRuleEvaluationModel)evaluationModel).getSchedulingSet().getSimulations();
		} else {

			double minDistance = Double.MAX_VALUE;
			List<Simulation> finalSelectedSimulations = new ArrayList<>();
			for (Map.Entry<int[], List<Simulation>> entry : ((GPRuleEvolutionStateDownSamplingV7)state).PCSimulations.entrySet()) {
				int[] indsCharList = entry.getKey();
				double distance = PhenoCharacterisation.distance(PC, indsCharList);
				if (distance == 0) {
					finalSelectedSimulations = entry.getValue();
					break;
				} else if (distance < minDistance) {
					minDistance = distance;
					finalSelectedSimulations = entry.getValue();
				}
			}

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

			if(state.generation < ((GPRuleEvolutionStateDownSampling)state).switchGen){
				simulation.setSequencingRule(sequencingRule); //indicate different individuals
				simulation.setRoutingRule(routingRule);
			} else {
				simulation.setBothRuleWithoutReset(sequencingRule,routingRule);

			}

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

				double ObjValue;
				if(simulation.getSystemState().getClockTime() == Double.MAX_VALUE) { //means this is a bad run
					if(objectives.get(0).getName().endsWith("profit")) {
						ObjValue = -Double.MAX_VALUE;
					}
					else {
						ObjValue = Double.MAX_VALUE;
					}
				} else {
					ObjValue = simulation.objectiveValue(objectives.get(i)); // this line: the value of makespan
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
