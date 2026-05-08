package yimei.jss.algorithm.featureweighted;

import ec.EvolutionState;
import ec.gp.GPNode;
import ec.util.Checkpoint;
import ec.util.Parameter;
import yimei.jss.feature.FeatureIgnorable;
import yimei.jss.feature.FeatureUtil;
import yimei.jss.feature.ignore.Ignorer;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.gp.TerminalsChangable;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.OptionalDouble;

/**
 * Created by fzhang on 22.5.2019.
 * set the weights of featurs according to the frequency of features (based on top-k individuals)
 */
public class FreGPRuleEvolutionState extends GPRuleEvolutionState implements TerminalsChangable, FeatureIgnorable {

	public static final String P_IGNORER = "ignorer";
	public static final String P_PRE_GENERATIONS = "pre-generations";
	public static final String P_POP_ADAPT_FRAC_ELITES = "pop-adapt-frac-elites";
	public static final String P_POP_ADAPT_FRAC_ADAPTED = "pop-adapt-frac-adapted";
	public static final String P_DO_ADAPT = "feature-selection-adapt-population";
	public static final String P_NUM_TOP_INDS = "num-topinds";

	private Ignorer ignorer;
	private int preGenerations;
	private double fracElites;
	private double fracAdapted;
	private boolean doAdapt;
	private int topInds;

	ArrayList<Double> aveDepthRatio = new ArrayList<>();
	ArrayList<Integer> failSwapTimes = new ArrayList<>();

	@Override
	public Ignorer getIgnorer() {
		return ignorer;
	}

	@Override
	public void setIgnorer(Ignorer ignorer) {
		this.ignorer = ignorer;
	}

	@Override
	public void setup(EvolutionState state, Parameter base) {
		super.setup(state, base);

		ignorer = (Ignorer)(state.parameters.getInstanceForParameter(
				new Parameter(P_IGNORER), null, Ignorer.class)); //ignorer = yimei.jss.feature.ignore.SimpleIgnorer
		preGenerations = state.parameters.getIntWithDefault(
				new Parameter(P_PRE_GENERATIONS), null, -1);  //50
		fracElites = state.parameters.getDoubleWithDefault(
				new Parameter(P_POP_ADAPT_FRAC_ELITES), null, 0.0); //0.0
		fracAdapted = state.parameters.getDoubleWithDefault(
				new Parameter(P_POP_ADAPT_FRAC_ADAPTED), null, 1.0); //0.1
		doAdapt = state.parameters.getBoolean(new Parameter(P_DO_ADAPT),
				null, true);
		topInds = state.parameters.getInt(new Parameter(P_NUM_TOP_INDS), null);
	}

	//fzhang 2019.5.19 pick terminals based on weighting power
	double[][] weights =  null;
	ArrayList<double[]>  saveWeights = new ArrayList<>();
	@Override
	public GPNode pickTerminalRandom(int subPopNum) {

	/*	int index = -1; //random[0].nextInt(terminals[subPopNum].lenth)

		//1. if no weights is set for features, then we will choose features uniformly.
		if (weights == null || weights[subPopNum] == null) {
			//output.warning("weights are null");
			index = random[0].nextInt(terminals[subPopNum].length);
		}
		else //2. otherwise, choose features based on their weighting power
		{
			//System.out.println(subPopNum);
			index = RandomChoice.pickFromDistribution(weights[subPopNum], random[0].nextDouble());
		}*/
		//System.out.println("The index of chosen features: " + index);

		//randomly choose the subtree
		int index = random[0].nextInt(terminals[subPopNum].length);

		return terminals[subPopNum][index];
	}


	@Override
	public int evolve() {
		//if (generation > 0)
		output.message("Generation " + generation);

		// EVALUATION
		statistics.preEvaluationStatistics(this);
		evaluator.evaluatePopulation(this); //Feature selection, firstly, evaluate population as usual; then clear population



		//removeInfInds();

		//=======================//fzhang 2019.8.9 in order to calculate weights based on all previous generations=======================================
		/*if (generation >= preGenerations)
		{
			//fzhang 2019.8.9 in order to calculate weights based on all previous generations
			if(generation == preGenerations){
				weights = new double[this.population.subpops.length][terminals[0].length];
				copyWeights = new double[this.population.subpops.length][terminals[0].length];
			}

			ArrayList<HashMap<String, Integer>> stats = PopulationUtils.Frequency(population, topInds); //stats contains two values, one is terminal name
			//and the other is its frequency
			//stats.toString();

//			weights = new double[this.population.subpops.length][];
			for(int subpop = 0; subpop < this.population.subpops.length; subpop++)
			{
//				weights[subpop] = new double[terminals[subpop].length];
				weights[subpop] = copyWeights[subpop].clone();
				for(int i = 0; i < terminals[0].length;i++)
				{
					String name = (terminals[0][i]).name();//the terminals in each population is same.
					for (int indIdx = subpop*topInds; indIdx < subpop*topInds+topInds; indIdx ++) {
						if(stats.get(indIdx).containsKey(name))
						{
							weights[subpop][i] += stats.get(indIdx).get(name);
						}
					}

				}
				//save the weights values in each generation
				//saveWeights.add(weights[subpop]) this is a java style, the weights will be changed later
				copyWeights[subpop] = weights[subpop].clone();
				saveWeights.add(weights[subpop].clone()); //need to use clone to copy the array
				RandomChoice.organizeDistribution(weights[subpop]); // organizeDistribution will change the input, that what we want, we need this to choose terminals
			} // for(int subpop = 0; ...

		}*/

		//fzhang 2019.5.19 check the frequency of terminals in each generation and set them as weighting power---better than taking all the previous generations into account (above code)
	/*	if (generation >= preGenerations)
		{
			ArrayList<HashMap<String, Integer>> stats = PopulationUtils.Frequency(population, topInds); //stats contains two values, one is terminal name
			//and the other is its frequency
			//stats.toString();
			weights = new double[this.population.subpops.length][];
			for(int subpop = 0; subpop < this.population.subpops.length; subpop++)
			{
				weights[subpop] = new double[terminals[subpop].length];
				for(int i = 0; i < terminals[0].length;i++)
				{
					String name = (terminals[0][i]).name();//the terminals in each population is same.
					for (int indIdx = subpop*topInds; indIdx < subpop*topInds+topInds; indIdx ++) {
						if(stats.get(indIdx).containsKey(name))
						{
							weights[subpop][i] += stats.get(indIdx).get(name);
						}
					}
				}
				//save the weights values in each generation
				//saveWeights.add(weights[subpop]) this is a java style, the weights will be changed later
				saveWeights.add(weights[subpop].clone()); //need to use clone to copy the array
				RandomChoice.organizeDistribution(weights[subpop]); // organizeDistribution will change the input --- weights[subpop], this is what we want. It will be used for pick the terminals
			}
		}*/

		//clearing, niching  MultiPopCoevolutionaryEvaluator
		statistics.postEvaluationStatistics(this);

		// SHOULD WE QUIT?
		if (evaluator.runComplete(this) && quitOnRunComplete)
		{
			output.message("Found Ideal Individual");
			return R_SUCCESS;
		}

		// SHOULD WE QUIT?
		if (generation == numGenerations-1)
		{
            //writetoFile();
			writeDepthRatioToFile();
			failSwapTimes();
			generation++;
			return R_FAILURE;
		}


		// PRE-BREEDING EXCHANGING
		statistics.prePreBreedingExchangeStatistics(this);
		population = exchanger.preBreedingExchangePopulation(this);
		statistics.postPreBreedingExchangeStatistics(this);

		String exchangerWantsToShutdown = exchanger.runComplete(this);
		if (exchangerWantsToShutdown!=null)
		{
			output.message(exchangerWantsToShutdown);
			/*
			 * Don't really know what to return here.  The only place I could
			 * find where runComplete ever returns non-null is
			 * IslandExchange.  However, that can return non-null whether or
			 * not the ideal individual was found (for example, if there was
			 * a communication error with the server).
			 *
			 * Since the original version of this code didn't care, and the
			 * result was initialized to R_SUCCESS before the while loop, I'm
			 * just going to return R_SUCCESS here.
			 */

			return R_SUCCESS;
		}

		// BREEDING
		statistics.preBreedingStatistics(this);

		population = breeder.breedPopulation(this);

		// POST-BREEDING EXCHANGING
		statistics.postBreedingStatistics(this);

		// POST-BREEDING EXCHANGING
		statistics.prePostBreedingExchangeStatistics(this);
		population = exchanger.postBreedingExchangePopulation(this);
		statistics.postPostBreedingExchangeStatistics(this);

		// Generate new instances if needed
		RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
		if (problem.getEvaluationModel().isRotatable()) {
			problem.rotateEvaluationModel();
		}

		//fzhang 2019.10.22  write the probability differenence information into csv
		//===========start=======
		if(this.generation == 1 || this.generation == 25 || this.generation == 45){
			writeProDifferenceTrue();
			writeProDifferenceFalse();
			weightedKozaNodeSelectorV9.proDifference1 = new ArrayList<>();
			weightedKozaNodeSelectorV9.proDifference2 = new ArrayList<>();

			writeCoefficientTrue();
			writeCoefficientFalse();
			weightedKozaNodeSelectorV9.coefficient1 = new ArrayList<>();
			weightedKozaNodeSelectorV9.coefficient2 = new ArrayList<>();
		}
		//==========end========

		//fzhang 2019.11.11 save the depthRatio (average value in each generation)
		OptionalDouble aveDepthRatioTrue = weightedKozaNodeSelectorV9.depthRatio1
				.stream()
				.mapToDouble(Double::doubleValue)
				.average();

		OptionalDouble aveDepthRatioFalse = weightedKozaNodeSelectorV9.depthRatio2
				.stream()
				.mapToDouble(Double::doubleValue)
				.average();

		aveDepthRatio.add(aveDepthRatioTrue.orElse(-1));
		aveDepthRatio.add(aveDepthRatioFalse.orElse(-1));
		weightedKozaNodeSelectorV9.depthRatio1 = new ArrayList<>();
		weightedKozaNodeSelectorV9.depthRatio2 = new ArrayList<>();

		//fzhang failSwapTimes
		failSwapTimes.add(weightedCrossoverPipelineV3.toolarge);
		weightedCrossoverPipelineV3.toolarge = 0;

		// INCREMENT GENERATION AND CHECKPOINT
		generation++;
		if (checkpoint && generation%checkpointModulo == 0)
		{
			output.message("Checkpointing");
			statistics.preCheckpointStatistics(this);
			Checkpoint.setCheckpoint(this);
			statistics.postCheckpointStatistics(this);
		}

		return R_NOTDONE;
	}

/*	private void removeInfInds() {
		// TODO Auto-generated method stub

	}*/

	@Override
	public void adaptPopulation(int subPopNum) {
		FeatureUtil.adaptPopulationThreeParts(this, fracElites, fracAdapted, subPopNum);
	}

	public void writetoFile() {
		//fzhang 2019.5.21 save the weight values
		File weightFile = new File("job." + jobSeed + ".weight.csv"); // jobSeed = 0
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
			writer.write("Gen, sNIQ, sWIQ, sMWT, sPT, sNPT, sOWT, sWKR, sNOR, sW, sTIS, "
					     + "rNIQ, rWIQ, rMWT, rPT, rNPT, rOWT, rWKR, rNOR, rW, rTIS");
			writer.newLine();
			double[] zero = new double[terminals[0].length];
//			zero

			for (int i = 0; i < preGenerations; i++) {
				//writer.newLine();
				writer.write(i + ", " + Arrays.toString(zero).replaceAll("\\[|\\]", ""));
				writer.write(", " + Arrays.toString(zero).replaceAll("\\[|\\]", "")  + "\n");
				//Returns a string representation of the contents of the specified array.
			}
			for (int i = 0; i < saveWeights.size(); i += 2) { //every two into one generation
				//writer.newLine();
				writer.write(i/2+preGenerations + ", " + Arrays.toString(saveWeights.get(i)).replaceAll("\\[|\\]", ""));
				writer.write(", " + Arrays.toString(saveWeights.get(i+1)).replaceAll("\\[|\\]", "")  + "\n");
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeProDifferenceTrue() {
		//fzhang 2019.5.21 save the weight values
		File proDiff1 = new File("job." + jobSeed + "." + generation + ".proDifferenceTrue.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(proDiff1));
			write.write("Gen, proDifferenceTrue");
			write.newLine();

			ArrayList<Double> temp = weightedKozaNodeSelectorV9.proDifference1;
			for (int i = 0; i < temp.size(); i ++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + temp.get(i) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeProDifferenceFalse() {
		//fzhang 2019.5.21 save the weight values
		File proDiff2 = new File("job." + jobSeed + "." + generation + ".proDifferenceFalse.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(proDiff2));
			write.write("Gen, proDifferenceFalse");
			write.newLine();
			ArrayList<Double> temp = weightedKozaNodeSelectorV9.proDifference2;
			for (int i = 0; i < temp.size(); i ++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + temp.get(i) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeCoefficientTrue() {
		//fzhang 2019.5.21 save the weight values
		File coefficient1 = new File("job." + jobSeed + "." + generation + ".coefficientTrue.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(coefficient1));
			write.write("Gen, coefficientTrue");
			write.newLine();

			ArrayList<Double> temp = weightedKozaNodeSelectorV9.coefficient1;
			for (int i = 0; i < temp.size(); i ++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + temp.get(i) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeCoefficientFalse() {
		//fzhang 2019.5.21 save the weight values
		File coefficient2 = new File("job." + jobSeed + "." + generation + ".coefficientFalse.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(coefficient2));
			write.write("Gen, coefficientFalse");
			write.newLine();
			ArrayList<Double> temp = weightedKozaNodeSelectorV9.coefficient2;
			for (int i = 0; i < temp.size(); i ++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + temp.get(i) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	//fzhang 2019.11.11 write depth ratio to file
	public void writeDepthRatioToFile() {
		File depthRatioFile = new File("job." + jobSeed + ".depthRatio.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(depthRatioFile));
			write.write("Gen, aveDepthRatioTrue, aveDepthRatioFalse");
			write.newLine();
			for (int i = 0; i < aveDepthRatio.size()/2; i++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + aveDepthRatio.get(2*i) + ", " + aveDepthRatio.get(2*i + 1) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	//fzhang 2019.11.11 write depth ratio to file
	public void failSwapTimes() {
		File failSwapTimesFile = new File("job." + jobSeed + ".failSwapTimes.csv"); // jobSeed = 0
		try {
			BufferedWriter write = new BufferedWriter(new FileWriter(failSwapTimesFile));
			write.write("Gen, failSwapTimes");
			write.newLine();
			for (int i = 0; i < failSwapTimes.size(); i++) { //every two into one generation
				//writer.newLine();
				write.write(i + ", " + failSwapTimes.get(i) + "\n");
			}
			write.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public ArrayList<double[]> getWeights(){
		return saveWeights;
	}
}
