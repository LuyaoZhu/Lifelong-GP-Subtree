package yimei.jss.algorithm.featureConstruction;

import ec.Individual;
import ec.gp.GPIndividual;
import ec.gp.GPNode;
import ec.gp.GPTree;
import ec.util.Checkpoint;
import ec.util.Parameter;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import yimei.jss.gp.GPRuleEvolutionState;
import yimei.jss.helper.PopulationUtils;
import yimei.jss.jobshop.OperationOption;
import yimei.jss.niching.*;
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
 * Created by fzhang on 22.5.2019.
 * select features from elites at the current generation
 */
public class FreGPRuleEvolutionStateV1 extends GPRuleEvolutionState{

    List<SequencingDecisionSituation> decisionSituationsSequencing = null;
    List<RoutingDecisionSituation> decisionSituationsRouting = null;

    public double thresholdSeq;
    public double thresholdRou;

    ArrayList<Integer> numSeqArchive = new ArrayList<>();
    ArrayList<Integer> numRouArchive = new ArrayList<>();

    ArrayList<int[]> seqArchiveDepth = new ArrayList<>();
    ArrayList<int[]> rouArchiveDepth = new ArrayList<>();

    public  ArrayList<ArrayList<GPTree>> seqFeatureArchive = new ArrayList<>();
    public ArrayList<ArrayList<GPTree>> rouFeatureArchive = new ArrayList<>();      // Save selected high-level features ( will be updated in every generation)

    public ArrayList<Individual[][]> elitesEveryGen = new ArrayList<>();

    public int numElites;

    public int[] seqFeatureDepth;
    public int[] rouFeatureDepth;

    ArrayList<double[]> entropyDiversity = new ArrayList<>();



    @Override
    public int evolve() {
        //if (generation > 0)
        output.message("Generation " + generation);

        //record the diversity
        PhenoCharacterisation[] pc = scoreMultiPopCoevolutionaryClearingEvaluator.getPhenoCharacterisation();
        double[][][] indsCharListsMultiTree = surrogateClearing.clearPopulation(this,pc); //3. calculate the phenotypic characteristic
        double[] diversityValue = new double[this.population.subpops.length];
        diversityValue[0] = PopulationUtils.entropy(indsCharListsMultiTree[0]);
        diversityValue[1] = PopulationUtils.entropy(indsCharListsMultiTree[1]);
        entropyDiversity.add(diversityValue);

        System.out.println(diversityValue[0] + ",," + diversityValue[1]);

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
            writeArchiveInformationToFile(numSeqArchive,numRouArchive,seqArchiveDepth,rouArchiveDepth);
            return R_FAILURE;
        }

        elitesEveryGen.clear();
        rouFeatureArchive.clear();
        seqFeatureArchive.clear();
        // ----begin-----feature construction based on the importance of the subtree in elites
        PopulationUtils.sort(this.population);
        numElites = this.parameters.getIntWithDefault(new Parameter("selectedEliteNum"), null, 1);
        Individual elites[][] = new Individual[this.population.subpops.length][numElites];
        for (int sub = 0; sub < this.population.subpops.length; sub++) {
            for (int i = 0; i < numElites; i++) {
                elites[sub][i] = (Individual) this.population.subpops[sub].individuals[i].clone();
            }
        }
        elitesEveryGen.add(elites);
        thresholdSeq = this.parameters.getDoubleWithDefault(new Parameter("thresholdSeq"), null, 0.8);
        thresholdRou = this.parameters.getDoubleWithDefault(new Parameter("thresholdRou"), null, 0.8);

        featureConstruction(this, elitesEveryGen);

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

        System.out.println((double)protectedKozaNodeSelectorV1.useTime/protectedKozaNodeSelectorV1.sum);

        protectedKozaNodeSelectorV1.sum = 0;
        protectedKozaNodeSelectorV1.useTime = 0;

//        System.out.println(protectedKozaNodeSelector.useTime);

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

    public void featureConstruction(GPRuleEvolutionState state, ArrayList<Individual[][]> elitesEveryGen) {
        PhenoCharacterisation[] pc = scoreMultiPopCoevolutionaryClearingEvaluator.getPhenoCharacterisation();
        //delete the features that are not related to the new elites in the archive
//        if (this.generation > 0) {
//            deleteUnimportantFeature(pc);
//        }
        if(this.generation%5 == 0)
            System.out.println("1");
        //select features into the archive
        ArrayList<GPTree> sequencingFeature = new ArrayList<>();
        ArrayList<GPTree> routingFeature = new ArrayList<>();
        for (int sub = 0; sub < elitesEveryGen.get(0).length; sub++) {
            for (int ind = 0; ind < elitesEveryGen.get(0)[sub].length; ind++) {
                if (sub == 0) {
                    decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pc[sub]).decisionSituations;
                } else {
                    decisionSituationsRouting = ((RoutingPhenoCharacterisation) pc[sub]).decisionSituations;
                }
                AbstractRule rule = null;
                GPTree tree = ((GPIndividual) (elitesEveryGen.get(0)[sub][ind])).trees[0];
                int nonterminals = tree.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                double[] score = new double[nonterminals];
                double[] frequency = new double[nonterminals];
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

                                if(i == 0){ //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for(int b=0; b<elitesEveryGen.get(0)[sub].length; b++){
                                        GPTree elite = ((GPIndividual)elitesEveryGen.get(0)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num=0; num < subtreeNum; num++){
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if(eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
                                }

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

                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {
                            GPTree treeClone = (GPTree) tree.clone();
                            GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
//                            System.out.println("correlation is " + Math.abs(score[numScore]) + ", depth ratio is " + ((double)node.depth()/treeClone.child.depth()));

                            if (Math.abs(score[numScore]) > thresholdSeq) {
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(nodeToTree.child.depth() < 5){
//                                    sequencingFeature.add(nodeToTree);
//                                }
                                sequencingFeature.add(nodeToTree);
                            }
                        }

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

                                if(i == 0){ //只需记录一次
                                    int occurrence = 0; //记录次特征在elites中出现的次数
                                    for(int b=0; b<elitesEveryGen.get(0)[sub].length; b++){
                                        GPTree elite = ((GPIndividual)elitesEveryGen.get(0)[sub][b]).trees[0];
                                        int subtreeNum = elite.child.numNodes(GPNode.NODESEARCH_NONTERMINALS);
                                        for (int num=0; num < subtreeNum; num++){
                                            GPTree eliteClone = (GPTree) elite.clone();
                                            GPNode eliteNode = eliteClone.child.nodeInPosition(num, GPNode.NODESEARCH_NONTERMINALS);
                                            if(eliteNode.makeLispTree().equals(node.makeLispTree()))
                                                occurrence++;
                                        }
                                    }
                                    frequency[numSubtree] = occurrence;
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
                        //save the correlation > threshold to the feature archive

                        for (int numScore = 1; numScore < score.length; numScore++) {
                            if (Math.abs(score[numScore]) > thresholdRou) {
                                GPTree treeClone = (GPTree) tree.clone();
                                GPNode node = treeClone.child.nodeInPosition(numScore, GPNode.NODESEARCH_NONTERMINALS);
                                GPTree nodeToTree = PopulationUtils.GPNodetoGPTree(node);
//                                if(nodeToTree.child.depth() < 5){
//                                    routingFeature.add(nodeToTree);
//                                }
                                routingFeature.add(nodeToTree);
                            }
                        }

                    }
                }
            }
        }
        seqFeatureArchive.add(sequencingFeature);
        rouFeatureArchive.add(routingFeature);


        //delete duplicate features
        deleteDuplicateFeatures();


        // record the depth of the features
        seqFeatureDepth = new int[6];
        rouFeatureDepth = new int[6];
        int sumSeqFeature = 0;
        for (int gen = 0; gen < seqFeatureArchive.size(); gen++) {
            if(!(seqFeatureArchive.get(gen).size() == 0)){
                System.out.println("the seqFeature number in generation "+(gen+1)+" is "+seqFeatureArchive.get(gen).size());
                for (int i=0; i<seqFeatureArchive.get(gen).size(); i++){
                    sumSeqFeature++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 2)
                        seqFeatureDepth[0]++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 3)
                        seqFeatureDepth[1]++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 4)
                        seqFeatureDepth[2]++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 5)
                        seqFeatureDepth[3]++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 6)
                        seqFeatureDepth[4]++;
                    if(seqFeatureArchive.get(gen).get(i).child.depth() == 7)
                        seqFeatureDepth[5]++;
                }
            }
        }
        seqArchiveDepth.add(seqFeatureDepth);
        numSeqArchive.add(sumSeqFeature);
        System.out.println("The number of feature in sequencing Archive is  " + sumSeqFeature);

        int sumRouFeature = 0;
        for (int gen = 0; gen < rouFeatureArchive.size(); gen++) {
            if(!(rouFeatureArchive.get(gen).size() == 0)){
                System.out.println("the rouFeature number in generation "+(gen+1)+" is "+rouFeatureArchive.get(gen).size());
                for (int i=0; i<rouFeatureArchive.get(gen).size(); i++){
                    sumRouFeature++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 2)
                        rouFeatureDepth[0]++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 3)
                        rouFeatureDepth[1]++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 4)
                        rouFeatureDepth[2]++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 5)
                        rouFeatureDepth[3]++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 6)
                        rouFeatureDepth[4]++;
                    if(rouFeatureArchive.get(gen).get(i).child.depth() == 7)
                        rouFeatureDepth[5]++;
                }
            }
        }
        rouArchiveDepth.add(rouFeatureDepth);
        numRouArchive.add(sumRouFeature);
        System.out.println("The number of feature in routing Archive is  " + sumRouFeature);


    }

    private void deleteDuplicateFeatures() {

        // for sequencing rule
        for (int gen1 = 0; gen1 < seqFeatureArchive.size(); gen1++) {
            for (int seqGenFea1 = 0; seqGenFea1 < seqFeatureArchive.get(gen1).size(); seqGenFea1++) {
                GPTree base = (GPTree) seqFeatureArchive.get(gen1).get(seqGenFea1).clone();

                for (int gen2 = 0; gen2 < seqFeatureArchive.size(); gen2++ ) {
                    for (int seqGenFea2 = 0; seqGenFea2 < seqFeatureArchive.get(gen2).size(); seqGenFea2++) {
                        GPTree count = (GPTree) seqFeatureArchive.get(gen2).get(seqGenFea2).clone();
                        if (!((gen1 == gen2) && (seqGenFea1 == seqGenFea2))) {
                            if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                                seqFeatureArchive.get(gen2).remove(seqFeatureArchive.get(gen2).get(seqGenFea2));
                                seqGenFea2--;
                            }
                        }
                    }
                }
            }
        }

        // for routing rule
        for (int gen1 = 0; gen1 < rouFeatureArchive.size(); gen1++) {
            for (int rouGenFea1 = 0; rouGenFea1 < rouFeatureArchive.get(gen1).size(); rouGenFea1++) {
                GPTree base = (GPTree) rouFeatureArchive.get(gen1).get(rouGenFea1).clone();

                for (int gen2 = 0; gen2 < rouFeatureArchive.size(); gen2++) {
                    for (int rouGenFea2 = 0; rouGenFea2 < rouFeatureArchive.get(gen2).size(); rouGenFea2++) {
                        GPTree count = (GPTree) rouFeatureArchive.get(gen2).get(rouGenFea2).clone();
                        if (!((gen1 == gen2) && (rouGenFea1 == rouGenFea2))) {
                            if (base.child.makeLispTree().equals(count.child.makeLispTree())) {
                                rouFeatureArchive.get(gen2).remove(rouFeatureArchive.get(gen2).get(rouGenFea2));
                                rouGenFea2--;
                            }
                        }
                    }
                }
            }
        }
    }

    private void deleteUnimportantFeature(PhenoCharacterisation[] pc) {
        AbstractRule rule = null;

        // sequencing rule
        int seqFeatureNum = 0;
        for (int gen = 0; gen < seqFeatureArchive.size(); gen++) {
            for (int seqGenFea = 0; seqGenFea < seqFeatureArchive.get(gen).size(); seqGenFea++) {
                seqFeatureNum++;
            }
        }

        double[][] scoreSeq = new double[seqFeatureNum][numElites];
        GPTree allSeqFeature[] = new GPTree[seqFeatureNum];

        decisionSituationsSequencing = ((SequencingPhenoCharacterisation) pc[0]).decisionSituations;
        int decisionSize = decisionSituationsSequencing.get(0).getQueue().size();
        double[][] score_matrixSeqFeature = new double[seqFeatureNum][decisionSize];
        double[][] score_matrixSeqElites = new double[numElites][decisionSize];
        for (int i = 0; i < decisionSituationsSequencing.size(); i++) {
            int num = 0;  //record the num st feature in the archive
            for (int gen = 0; gen < seqFeatureArchive.size(); gen++) {
                for (int seqGenFea = 0; seqGenFea < seqFeatureArchive.get(gen).size(); seqGenFea++) {

                    rule = new GPRule(RuleType.SEQUENCING, seqFeatureArchive.get(gen).get(seqGenFea));
                    allSeqFeature[num] = (GPTree) seqFeatureArchive.get(gen).get(seqGenFea).clone();
                    SequencingDecisionSituation situation = decisionSituationsSequencing.get(i);
                    List<OperationOption> queue = situation.getQueue();
                    for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                        queue.get(candiateNum).setPriority(rule.priority(queue.get(candiateNum), situation.getWorkCenter(), situation.getSystemState()));
                        score_matrixSeqFeature[num][candiateNum] = queue.get(candiateNum).getPriority();
                    }
                    num++;
                }
            }
            double rank = 1;
            for (int numNodes = 0; numNodes < score_matrixSeqFeature.length; numNodes++) {
                double[] ranks = new double[score_matrixSeqFeature[numNodes].length];
                for (int queueSize1 = 0; queueSize1 < score_matrixSeqFeature[numNodes].length; queueSize1++) {
                    for (int queueSize2 = 0; queueSize2 < score_matrixSeqFeature[numNodes].length; queueSize2++) {
                        if (score_matrixSeqFeature[numNodes][queueSize2] < score_matrixSeqFeature[numNodes][queueSize1]) {
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

                score_matrixSeqFeature[numNodes] = ranks;
            }

            //calculate matrix of elites
            for (int elite = 0; elite < elitesEveryGen.get(this.generation)[0].length; elite++) {
                GPTree eliteTree = (GPTree) ((GPIndividual) elitesEveryGen.get(this.generation)[0][elite]).trees[0].clone();
                rule = new GPRule(RuleType.SEQUENCING, eliteTree);
                SequencingDecisionSituation situation = decisionSituationsSequencing.get(i);
                List<OperationOption> queue = situation.getQueue();
                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                    queue.get(candiateNum).setPriority(rule.priority(queue.get(candiateNum), situation.getWorkCenter(), situation.getSystemState()));
                    score_matrixSeqElites[elite][candiateNum] = queue.get(candiateNum).getPriority();
                }
            }
            rank = 1;
            for (int numNodes = 0; numNodes < score_matrixSeqElites.length; numNodes++) {
                double[] ranks = new double[score_matrixSeqElites[numNodes].length];
                for (int queueSize1 = 0; queueSize1 < score_matrixSeqElites[numNodes].length; queueSize1++) {
                    for (int queueSize2 = 0; queueSize2 < score_matrixSeqElites[numNodes].length; queueSize2++) {
                        if (score_matrixSeqElites[numNodes][queueSize2] < score_matrixSeqElites[numNodes][queueSize1]) {
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

                score_matrixSeqElites[numNodes] = ranks;
            }


            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
            for (int count = 0; count < score_matrixSeqFeature.length; count++) {
                for (int elite = 0; elite < score_matrixSeqElites.length; elite++) {
                    Double correlation = new SpearmansCorrelation().correlation(score_matrixSeqFeature[count], score_matrixSeqElites[elite]);
                    scoreSeq[count][elite] += correlation;
                }
            }
        }
        double[] scoreSeqAve = new double[score_matrixSeqFeature.length];
        for (int count = 0; count < score_matrixSeqFeature.length; count++) {
            for (int elite = 0; elite < score_matrixSeqElites.length; elite++) {
                scoreSeq[count][elite] = scoreSeq[count][elite] / decisionSituationsSequencing.size();
                scoreSeqAve[count] += scoreSeq[count][elite];
            }
            scoreSeqAve[count] = scoreSeqAve[count] / numElites;
            if (Math.abs(scoreSeqAve[count]) < thresholdSeq) {
                for (int gen = 0; gen < seqFeatureArchive.size(); gen++) {
                    for (int seqGenFea = 0; seqGenFea < seqFeatureArchive.get(gen).size(); seqGenFea++) {
                        if ((seqFeatureArchive.get(gen).get(seqGenFea).child.makeLispTree()).equals(allSeqFeature[count].child.makeLispTree())) {
                            seqFeatureArchive.get(gen).remove(seqFeatureArchive.get(gen).get(seqGenFea));
                            seqGenFea--;
                        }
                    }
                }
            }
        }


        //routing rule
        int rouFeatureNum = 0;
        for (int gen = 0; gen < rouFeatureArchive.size(); gen++) {
            for (int rouGenFea = 0; rouGenFea < rouFeatureArchive.get(gen).size(); rouGenFea++) {
                rouFeatureNum++;
            }
        }

        double[][] scoreRou = new double[rouFeatureNum][numElites];
        GPTree allRouFeature[] = new GPTree[rouFeatureNum];

        decisionSituationsRouting = ((RoutingPhenoCharacterisation) pc[1]).decisionSituations;
        decisionSize = decisionSituationsRouting.get(0).getQueue().size();
        double[][] score_matrixRouFeature = new double[rouFeatureNum][decisionSize];
        double[][] score_matrixRouElites = new double[numElites][decisionSize];
        for (int i = 0; i < decisionSituationsRouting.size(); i++) {
            int num = 0;  //record the num st feature in the archive
            for (int gen = 0; gen < rouFeatureArchive.size(); gen++) {
                for (int rouGenFea = 0; rouGenFea < rouFeatureArchive.get(gen).size(); rouGenFea++) {

                    rule = new GPRule(RuleType.ROUTING, rouFeatureArchive.get(gen).get(rouGenFea));
                    allRouFeature[num] = (GPTree) rouFeatureArchive.get(gen).get(rouGenFea).clone();
                    RoutingDecisionSituation situation = decisionSituationsRouting.get(i);
                    List<OperationOption> queue = situation.getQueue();

                    for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                        OperationOption operationOption = queue.get(candiateNum);
                        queue.get(candiateNum).setPriority(rule.priority(operationOption, operationOption.getWorkCenter(), situation.getSystemState()));
                        score_matrixRouFeature[num][candiateNum] = queue.get(candiateNum).getPriority();
                    }
                    num++;
                }
            }
            double rank = 1;
            for (int numNodes = 0; numNodes < score_matrixRouFeature.length; numNodes++) {
                double[] ranks = new double[score_matrixRouFeature[numNodes].length];
                for (int queueSize1 = 0; queueSize1 < score_matrixRouFeature[numNodes].length; queueSize1++) {
                    for (int queueSize2 = 0; queueSize2 < score_matrixRouFeature[numNodes].length; queueSize2++) {
                        if (score_matrixRouFeature[numNodes][queueSize2] < score_matrixRouFeature[numNodes][queueSize1]) {
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

                score_matrixRouFeature[numNodes] = ranks;
            }

            //calculate matrix of elites
            for (int elite = 0; elite < elitesEveryGen.get(this.generation)[1].length; elite++) {
                GPTree eliteTree = (GPTree) ((GPIndividual) elitesEveryGen.get(this.generation)[1][elite]).trees[0].clone();
                rule = new GPRule(RuleType.ROUTING, eliteTree);
                RoutingDecisionSituation situation = decisionSituationsRouting.get(i);
                List<OperationOption> queue = situation.getQueue();
                for (int candiateNum = 0; candiateNum < queue.size(); candiateNum++) {
                    OperationOption operationOption = queue.get(candiateNum);
                    queue.get(candiateNum).setPriority(rule.priority(operationOption, operationOption.getWorkCenter(), situation.getSystemState()));
                    score_matrixRouElites[elite][candiateNum] = queue.get(candiateNum).getPriority();
                }
            }
            rank = 1;
            for (int numNodes = 0; numNodes < score_matrixRouElites.length; numNodes++) {
                double[] ranks = new double[score_matrixRouElites[numNodes].length];
                for (int queueSize1 = 0; queueSize1 < score_matrixRouElites[numNodes].length; queueSize1++) {
                    for (int queueSize2 = 0; queueSize2 < score_matrixRouElites[numNodes].length; queueSize2++) {
                        if (score_matrixRouElites[numNodes][queueSize2] < score_matrixRouElites[numNodes][queueSize1]) {
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

                score_matrixRouElites[numNodes] = ranks;
            }


            //but for ranks based, there is no NaN values, all the ranks are as 1,3,4,5,2
            for (int count = 0; count < score_matrixRouFeature.length; count++) {
                for (int elite = 0; elite < score_matrixRouElites.length; elite++) {
                    Double correlation = new SpearmansCorrelation().correlation(score_matrixRouFeature[count], score_matrixRouElites[elite]);
                    scoreRou[count][elite] += correlation;
                }
            }
        }
        double[] scoreRouAve = new double[score_matrixRouFeature.length];
        for (int count = 0; count < score_matrixRouFeature.length; count++) {
            for (int elite = 0; elite < score_matrixRouElites.length; elite++) {
                scoreRou[count][elite] = scoreRou[count][elite] / decisionSituationsRouting.size();
                scoreRouAve[count] += scoreRou[count][elite];
            }
            scoreRouAve[count] = scoreRouAve[count] / numElites;
            if (Math.abs(scoreRouAve[count]) < thresholdRou) {
                for (int gen = 0; gen < rouFeatureArchive.size(); gen++) {
                    for (int rouGenFea = 0; rouGenFea < rouFeatureArchive.get(gen).size(); rouGenFea++) {
                        if ((rouFeatureArchive.get(gen).get(rouGenFea).child.makeLispTree()).equals(allRouFeature[count].child.makeLispTree())) {
                            rouFeatureArchive.get(gen).remove(rouFeatureArchive.get(gen).get(rouGenFea));
                            rouGenFea--;
                        }
                    }
                }
            }
        }


    }

    public void writeArchiveInformationToFile(ArrayList<Integer> numSeqArchive, ArrayList<Integer> numRouArchive, ArrayList<int[]> seqArchiveDepth, ArrayList<int[]> rouArchiveDepth) {
        //fzhang 2019.5.21 save the number of cleared individuals
        File weightFile = new File( out_dir + "/job." + jobSeed +  ".archiveInformation.csv");
//        File weightFile = new File( "/job." + jobSeed + ".archiveInformation.csv"); // jobSeed = 0
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(weightFile));
            writer.write("Gen,numSeqArchive,numRouArchive,seqDepth2,seqDepth3,seqDepth4,seqDepth5,seqDepth6,seqDepth7,rouDepth2,rouDepth3,rouDepth4,rouDepth5,rouDepth6,rouDepth7");
            writer.newLine();
            for (int i = 0; i < numSeqArchive.size(); i++) { //every two into one generation
                //writer.newLine();
                writer.write(i + ", " + numSeqArchive.get(i) + ", " + numRouArchive.get(i) + ", "
                        + seqArchiveDepth.get(i)[0] + ", " + seqArchiveDepth.get(i)[1] + ", "+ seqArchiveDepth.get(i)[2] + ", "
                        + seqArchiveDepth.get(i)[3] + ", " + seqArchiveDepth.get(i)[4] + ", "+ seqArchiveDepth.get(i)[5] + ", "
                        + rouArchiveDepth.get(i)[0] + ", " + rouArchiveDepth.get(i)[1] + ", "+ rouArchiveDepth.get(i)[2] + ", "
                        + rouArchiveDepth.get(i)[3] + ", " + rouArchiveDepth.get(i)[4] + ", "+ rouArchiveDepth.get(i)[5] +"\n");
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
