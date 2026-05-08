package yimei.jss.algorithm.simplification;

import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.util.Checkpoint;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.niching.PhenoCharacterisation;
import yimei.jss.niching.RoutingPhenoCharacterisation;
import yimei.jss.niching.SequencingPhenoCharacterisation;
import yimei.jss.niching.scoreMultiPopCoevolutionaryClearingEvaluator;
import yimei.jss.rule.AbstractRule;
import yimei.jss.rule.RuleType;
import yimei.jss.rule.operation.evolved.GPRule;
import yimei.jss.ruleoptimisation.RuleOptimizationProblem;
import yimei.jss.simulation.RoutingDecisionSituation;
import yimei.jss.simulation.SequencingDecisionSituation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static yimei.jss.gp.GPRun.out_dir;


/**
 * Created by luyao on 22.8.2023.
 * simplify GP Rule according to the correlation between subtrees and the fullTree
 * 1.replace the fullTree with the most relevant subtree;
 * 2.replace the unimportant subtrees with important 2-layer subtrees (relevant the elites)
 */
public class SimplifyGPRuleEvolutionState extends GPRuleEvolutionState {

    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;

//    public double thresholdSeq;
//
//    public double thresholdRou;

//    ArrayList<Integer> numSeqArchive = new ArrayList<>();
//    ArrayList<Integer> numRouArchive = new ArrayList<>();
//
//    ArrayList<int[]> seqArchiveDepth = new ArrayList<>();
//    ArrayList<int[]> rouArchiveDepth = new ArrayList<>();

    //    public  ArrayList<ArrayList<GPTree>> seqFeatureArchive = new ArrayList<>();
//    public ArrayList<ArrayList<GPTree>> rouFeatureArchive = new ArrayList<>();      // Save selected high-level features ( will be updated in every generation)
//
//    public ArrayList<Individual[][]> elitesEveryGen = new ArrayList<>();
//
    public int numElites;
//
//    public int[] seqFeatureDepth;
//    public int[] rouFeatureDepth;


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
        if (evaluator.runComplete(this) && quitOnRunComplete) {
            output.message("Found Ideal Individual");
            return R_SUCCESS;
        }

        // SHOULD WE QUIT?
        if (generation == numGenerations - 1) {
            //writetoFile();
            generation++;
            return R_FAILURE;
        }

        // ----begin-----feature construction based on the importance of the subtree in elites
//        PopulationUtils.sort(this.population);
//        numElites = this.parameters.getIntWithDefault(new Parameter("selectedEliteNum"), null, 1);
//        Individual elites[][] = new Individual[this.population.subpops.length][numElites];
//        for (int sub = 0; sub < this.population.subpops.length; sub++) {
//            for (int i = 0; i < numElites; i++) {
//                elites[sub][i] = (Individual) this.population.subpops[sub].individuals[i].clone();
//            }
//        }
//        elitesEveryGen.add(elites);
//        thresholdSeq = this.parameters.getDoubleWithDefault(new Parameter("thresholdSeq"), null, 0.8);
//        thresholdRou = this.parameters.getDoubleWithDefault(new Parameter("thresholdRou"), null, 0.8);


        //-----end------


        // PRE-BREEDING EXCHANGING
        statistics.prePreBreedingExchangeStatistics(this);
        population = exchanger.preBreedingExchangePopulation(this);
        statistics.postPreBreedingExchangeStatistics(this);

        String exchangerWantsToShutdown = exchanger.runComplete(this);
        if (exchangerWantsToShutdown != null) {
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

//        protectedKozaNodeSelector.useTime = 0;

        population = breeder.breedPopulation(this);

//        System.out.println(protectedKozaNodeSelector.useTime);

        // POST-BREEDING EXCHANGING
        statistics.postBreedingStatistics(this);

        // POST-BREEDING EXCHANGING
        statistics.prePostBreedingExchangeStatistics(this);
        population = exchanger.postBreedingExchangePopulation(this);
        statistics.postPostBreedingExchangeStatistics(this);

        //After breed new population, simplify these new individuals

        double rnd = this.random[0].nextDouble();
        if (rnd > 0.0) {
            simplifyGPRule(this);
        }


        // Generate new instances if needed
        RuleOptimizationProblem problem = (RuleOptimizationProblem) evaluator.p_problem;
        if (problem.getEvaluationModel().isRotatable()) {
            problem.rotateEvaluationModel();
        }


        // INCREMENT GENERATION AND CHECKPOINT
        generation++;
        if (checkpoint && generation % checkpointModulo == 0) {
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

    public void simplifyGPRule(GPRuleEvolutionState state) {
        PhenoCharacterisation[] pc = scoreMultiPopCoevolutionaryClearingEvaluator.getPhenoCharacterisation();

//        ArrayList<GPTree> sequencingFeature = new ArrayList<>();
//        ArrayList<GPTree> routingFeature = new ArrayList<>();
        for (int sub = 0; sub < state.population.subpops.length; sub++) {
            for (int ind = 0; ind < state.population.subpops[sub].individuals.length; ind++) {
                if (sub == 0) {
                    decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pc[sub]).decisionSituations;
                } else {
                    decisionSituationsRouting = ((RoutingPhenoCharacterisation) pc[sub]).decisionSituations;
                }
                AbstractRule rule = null;
                GPTree tree = ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0];
                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                double[] score = new double[nonterminals];
//                double[] frequency = new double[nonterminals];
                if (nonterminals > 1) {
                    if (sub == 0) {
                        int decisionSize = decisionSituationsSequencing.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        for (int i = 0; i < decisionSituationsSequencing.size(); i++) {
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {

                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.SEQUENCING, nodeToTree);

                                SequencingDecisionSituation situation = decisionSituationsSequencing.get(i);
                                List<OperationOption> queue = situation.getQueue();
                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    queue.get(candiateNum).setPriority(rule.priority(queue.get(candiateNum), situation.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }

                                // record the frequency of subtrees
//                                if(i == 0){ //只需记录一次
//                                    int occurrence = 0; //记录次特征在elites中出现的次数
//                                    for(int b=0; b<elitesEveryGen.get(state.generation)[sub].length; b++){
//                                        GPTree elite = ((GPIndividual)elitesEveryGen.get(state.generation)[sub][b]).trees[0];
//                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
//                                        for (int num=0; num < subtreeNum; num++){
//                                            GPTree eliteClone = (GPTree) elite.clone();
//                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
//                                            if(eliteNode.makeLispTree().equals(node.makeLispTree()))
//                                                occurrence++;
//                                        }
//                                    }
//                                    frequency[numSubtree] = occurrence;
//                                }

                            }

                            //2019.10.24 use the rank rather than the priority value directly
                            //========================================start==========================================
                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }
                            //========================================end============================================

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }

                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsSequencing.size();
                        }

                        //记录下最高的相关性
                        double maxCorrelation = 0;
                        for (int numScore = 1; numScore < score.length; numScore++) {
                            if (Math.abs(score[numScore]) > maxCorrelation) {
                                maxCorrelation = Math.abs(score[numScore]);
                            }
                        }
                        //也许相关性最高的有多个子树
                        ArrayList<Integer> selectedIndexes = new ArrayList<>();
                        for (int numScore = 1; numScore < score.length; numScore++) {
                            if (Math.abs(score[numScore]) == maxCorrelation)
                                selectedIndexes.add(numScore);
                        }
                        //从中选出最小节点的子树,确定其索引
                        int index = 0;
                        int minNumNodes = Integer.MAX_VALUE;
                        for (int a = 0; a < selectedIndexes.size(); a++) {
                            GPTree treeClone = (GPTree) tree.clone();
                            GPNode node = treeClone.child.nodeInPosition(selectedIndexes.get(a), GPNode.NODESEARCH_NONTERMINALS);
                            GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                            int sumNodes = nodeToTree.child.numNodes(GPNode.NODESEARCH_TERMINALS) + nodeToTree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                            if (sumNodes < minNumNodes) {
                                minNumNodes = sumNodes;
                                index = selectedIndexes.get(a);
                            }
                        }
                        System.out.println(minNumNodes);

                        //start replacing
                        GPTree treeClone = (GPTree) tree.clone();
                        GPNode node = treeClone.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);
                        GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0] = (GPTree) nodeToTree.clone();
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].constraints = 1;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].printTerminalsAsVariablesInC = true;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].printTwoArgumentNonterminalsAsOperatorsInC = true;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].owner = (GPIndividual) state.population.subpops[sub].individuals[ind];

                    } else {
                        int decisionSize = decisionSituationsRouting.get(0).getQueue().size();
                        double[][] score_matrix = new double[nonterminals][decisionSize];
                        //for (int i = 0; i < 1; i++) {
                        for (int i = 0; i < decisionSituationsRouting.size(); i++) {
                            for (int numSubtree = 0; numSubtree < nonterminals; numSubtree++) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numSubtree, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                                rule = new GPRule(RuleType.ROUTING, nodeToTree);

                                RoutingDecisionSituation situation = decisionSituationsRouting.get(i);
                                List<OperationOption> queue = situation.getQueue();

                                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                                    OperationOption operationOption = queue.get(candiateNum);
                                    queue.get(candiateNum).setPriority(rule.priority(operationOption, operationOption.getWorkCenter(), situation.getSystemState()));
                                    score_matrix[numSubtree][candiateNum] = queue.get(candiateNum).getPriority();
                                }
                            }

                            double rank = 1;
                            for (int numNodes = 0; numNodes < score_matrix.length; numNodes++) {
                                double[] ranks = new double[score_matrix[numNodes].length];
                                for (int queueSize1 = 0; queueSize1 < score_matrix[numNodes].length; queueSize1++) {
                                    for (int queueSize2 = 0; queueSize2 < score_matrix[numNodes].length; queueSize2++) {
                                        if (score_matrix[numNodes][queueSize2] < score_matrix[numNodes][queueSize1]) {
                                            rank++;
                                        }
                                    }
                                    ranks[queueSize1] = rank;
                                    rank = 1;
                                }

                                int count = 0;
                                for (int r = ranks.length - 1; r >= 0; r--) {
                                    for (int j = 0; j < r; j++) {
                                        if (ranks[j] == ranks[r]) {
                                            count++;
                                        }
                                    }
                                    ranks[r] += count;
                                    count = 0;
                                }

                                score_matrix[numNodes] = ranks;
                            }

                            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
                            for (int count = 0; count < score_matrix.length; count++) {
                                Double correlation = new SpearmansCorrelation().correlation(score_matrix[count], score_matrix[0]);
                                score[count] += correlation;
                            }
                        }
                        for (int numScore = 0; numScore < score.length; numScore++) {
                            score[numScore] = score[numScore] / decisionSituationsRouting.size();
                        }

                        //记录下最高的相关性
                        double maxCorrelation = 0;
                        for (int numScore = 1; numScore < score.length; numScore++) {
                            if (Math.abs(score[numScore]) > maxCorrelation) {
                                maxCorrelation = Math.abs(score[numScore]);
                            }
                        }
                        //也许相关性最高的有多个子树
                        ArrayList<Integer> selectedIndexes = new ArrayList<>();
                        for (int numScore = 1; numScore < score.length; numScore++) {
                            if (Math.abs(score[numScore]) == maxCorrelation)
                                selectedIndexes.add(numScore);
                        }
                        //从中选出最小节点的子树,确定其索引
                        int index = 0;
                        int minNumNodes = Integer.MAX_VALUE;
                        for (int a = 0; a < selectedIndexes.size(); a++) {
                            GPTree treeClone = (GPTree) tree.clone();
                            GPNode node = treeClone.child.nodeInPosition(selectedIndexes.get(a), GPNode.NODESEARCH_NONTERMINALS);
                            GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                            int sumNodes = nodeToTree.child.numNodes(GPNode.NODESEARCH_TERMINALS) + nodeToTree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                            if (sumNodes < minNumNodes) {
                                minNumNodes = sumNodes;
                                index = selectedIndexes.get(a);
                            }
                        }

                        //start replacing
                        GPTree treeClone = (GPTree) tree.clone();
                        GPNode node = treeClone.child.nodeInPosition(index, GPNode.NODESEARCH_NONTERMINALS);
                        GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0] = (GPTree) nodeToTree.clone();
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].constraints = 0;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].printTerminalsAsVariablesInC = true;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].printTwoArgumentNonterminalsAsOperatorsInC = true;
                        ((GPIndividual) (state.population.subpops[sub].individuals[ind])).trees[0].owner = (GPIndividual) state.population.subpops[sub].individuals[ind];

                    }
                }
            }
        }


        // record the depth of the features
//        seqFeatureDepth = new int[6];
//        rouFeatureDepth = new int[6];
//        int sumSeqFeature = 0;
//        for (int gen = 0; gen < seqFeatureArchive.size(); gen++) {
//            if (!(seqFeatureArchive.get(gen).size() == 0)) {
//                System.out.println("the seqFeature number in generation " + (gen + 1) + " is " + seqFeatureArchive.get(gen).size());
//                for (int i = 0; i < seqFeatureArchive.get(gen).size(); i++) {
//                    sumSeqFeature++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 2)
//                        seqFeatureDepth[0]++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 3)
//                        seqFeatureDepth[1]++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 4)
//                        seqFeatureDepth[2]++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 5)
//                        seqFeatureDepth[3]++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 6)
//                        seqFeatureDepth[4]++;
//                    if (seqFeatureArchive.get(gen).get(i).child.depth() == 7)
//                        seqFeatureDepth[5]++;
//                }
//            }
//        }
//        seqArchiveDepth.add(seqFeatureDepth);
//        numSeqArchive.add(sumSeqFeature);
//        System.out.println("The number of feature in sequencing Archive is  " + sumSeqFeature);
//
//        int sumRouFeature = 0;
//        for (int gen = 0; gen < rouFeatureArchive.size(); gen++) {
//            if (!(rouFeatureArchive.get(gen).size() == 0)) {
//                System.out.println("the rouFeature number in generation " + (gen + 1) + " is " + rouFeatureArchive.get(gen).size());
//                for (int i = 0; i < rouFeatureArchive.get(gen).size(); i++) {
//                    sumRouFeature++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 2)
//                        rouFeatureDepth[0]++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 3)
//                        rouFeatureDepth[1]++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 4)
//                        rouFeatureDepth[2]++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 5)
//                        rouFeatureDepth[3]++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 6)
//                        rouFeatureDepth[4]++;
//                    if (rouFeatureArchive.get(gen).get(i).child.depth() == 7)
//                        rouFeatureDepth[5]++;
//                }
//            }
//        }
//        rouArchiveDepth.add(rouFeatureDepth);
//        numRouArchive.add(sumRouFeature);
//        System.out.println("The number of feature in routing Archive is  " + sumRouFeature);


    }

    public void writeArchiveInformationToFile(ArrayList<Integer> numSeqArchive, ArrayList<Integer> numRouArchive, ArrayList<int[]> seqArchiveDepth, ArrayList<int[]> rouArchiveDepth) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File(out_dir + "/job." + jobSeed + ".archiveInformation.csv");
//        File weightFile = new File( "/job." + jobSeed + ".archiveInformation.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,numSeqArchive,numRouArchive,seqDepth2,seqDepth3,seqDepth4,seqDepth5,seqDepth6,seqDepth7,rouDepth2,rouDepth3,rouDepth4,rouDepth5,rouDepth6,rouDepth7");
            writer.newLine();
            for (int i = 0; i < numSeqArchive.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + numSeqArchive.get(i) + ", " + numRouArchive.get(i) + ", "
                        + seqArchiveDepth.get(i)[0] + ", " + seqArchiveDepth.get(i)[1] + ", " + seqArchiveDepth.get(i)[2] + ", "
                        + seqArchiveDepth.get(i)[3] + ", " + seqArchiveDepth.get(i)[4] + ", " + seqArchiveDepth.get(i)[5] + ", "
                        + rouArchiveDepth.get(i)[0] + ", " + rouArchiveDepth.get(i)[1] + ", " + rouArchiveDepth.get(i)[2] + ", "
                        + rouArchiveDepth.get(i)[3] + ", " + rouArchiveDepth.get(i)[4] + ", " + rouArchiveDepth.get(i)[5] + "\n");
            }
            numSeqArchive.clear();
            numRouArchive.clear();
            seqArchiveDepth.clear();
            rouArchiveDepth.clear();
//			writer.write(numGenerations -1 + ", " + 0 + "\n");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


//    public boolean isEqual(GPTree t1, GPTree t2) {
//        if(t1.child == t2.child)
//        return doAdapt;
//    }
}
