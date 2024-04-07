package main;
import java.util.Random;

import utility.L1distance;

import anonymization.Mechanisms;
import anonymization.SyntheticDataMechanism;

import datahandler.DataGeneration;


public class CheckMechanism {
	
	public enum Mechanism {
		Laplace, Geometric, SyntheticData
	}
	
	public enum DistributionType {
		Uniform, Bimodal, Geometric
	}
	
	public static void main(String[] args) {
		DistributionType distType = DistributionType.Uniform;
		Mechanism mechanism = Mechanism.SyntheticData;
		
		//int[] numBlocksArray = {5, 10, 50, 100, 250};
		//testVaryingNumBlocks(numBlocksArray, distType, mechanism);
		
		//double[] epsilonArray = {0.1, 0.2, 0.5, 1, 2, 5, 10};
		//testVaryingEpsilon(epsilonArray, distType, mechanism);
		int[] vals = {3, 5, 7, 10};
		for (int val: vals) {
			int sum = sum(val);
			System.out.println(val + ": " + sum);
		}
		
	}
	
	public static int sum(int level) {
		int sum = 0;
		for (int i=0; i<=level; i++) 
			sum += (1 << i);
		for (int i=level+1; i<=20; i++)
			sum += (1 << i)/2;
		return sum;
	}
	public static void testVaryingEpsilon(double[] epsilonArray, DistributionType distType, Mechanism mechanism) {
		int seed = 123456;
		Random rand = new Random(seed);
		int numRuns = 20;
		// Laplace 
		double sensitivity = 2;
		// synthetic data
		double delta = 1E-5;
		double alpha = 500000;
		// input data
		int totalNumPeople = 500;
		int numBlocks = 5;
		int blockSize = totalNumPeople/numBlocks;
		double geoProb = 0.4;
		
		int[] inputArray = new int[numBlocks];
		if (distType == DistributionType.Uniform) {	
			for (int i=0 ; i<numBlocks; i++)
				inputArray[i] = blockSize;
		}
		else if (distType == DistributionType.Geometric) {
			inputArray = DataGeneration.genGeometricArray(numBlocks, totalNumPeople, geoProb);
		}
		
		for (double epsilon:epsilonArray) {
			//int blockSize = 50;
			//totalNumPeople = blockSize*numBlocks;
			double laplaceStdDev = sensitivity/epsilon;
				
			double[] outputArray = new double[numBlocks];
		
			boolean isPDP = false;
			double[] priorArray = new double[numBlocks];
			if (mechanism == Mechanism.SyntheticData) {
				for (int i=0; i<numBlocks; i++) {
					priorArray[i] = alpha;
				}
				isPDP = SyntheticDataMechanism.checkSyntheticData(inputArray, priorArray, epsilon, delta);
				//System.out.println(numBlocks + ", " + isPDP);
			}
			
			double totalDist = 0;
			
			for (int k=0; k<numRuns; k++) {
				if (mechanism == Mechanism.Laplace) {
					for (int i=0; i<numBlocks; i++) {
						outputArray[i] = Mechanisms.getLaplace(inputArray[i], laplaceStdDev, rand);
					}
					
					double l1Dist = L1distance.getL1Dist(inputArray, outputArray);
					//System.out.println("Distance=" + l1Dist);
					totalDist += l1Dist;
				}
				else if (mechanism == Mechanism.SyntheticData) {
					if (isPDP) {
						int[] output = SyntheticDataMechanism.synthesize(inputArray, priorArray, totalNumPeople, seed*k);
						double l1Dist  = L1distance.getL1Dist(inputArray, output);
						totalDist += l1Dist;
					}
				}
			}
			
			if (mechanism == Mechanism.Laplace)
				System.out.println("epsilon=" + epsilon + ", dist=" + totalDist/numRuns);
			else if (mechanism == Mechanism.SyntheticData)
				System.out.println("epsilon=" + epsilon + ", dist=" + totalDist/numRuns + ", " + isPDP);
		}
	}
	
	// fix total number of people, vary number of origin blocks (or block size)
	// use 1 destination only
	public static void testVaryingNumBlocks(int[] numBlocksArray, DistributionType distType, Mechanism mechanism) {
		int seed = 123456;
		Random rand = new Random(seed);
		int numRuns = 20;
		// Laplace 
		double epsilon = 1;
		double sensitivity = 2;
		double laplaceStdDev = sensitivity/epsilon;
		// synthetic data
		double delta = 1E-5;
		double alpha = 50;
		// input data
		int totalNumPeople = 500;
		double geoProb = 0.4;
		
		for (int numBlocks: numBlocksArray) {
			int blockSize = totalNumPeople/numBlocks;
			//int blockSize = 50;
			//totalNumPeople = blockSize*numBlocks;
			
			alpha = Math.max(10, 4*totalNumPeople/numBlocks);
			
			int[] inputArray = new int[numBlocks];
			double[] outputArray = new double[numBlocks];
		
			if (distType == DistributionType.Uniform) {	
				for (int i=0 ; i<numBlocks; i++)
					inputArray[i] = blockSize;
			}
			else if (distType == DistributionType.Geometric) {
				inputArray = DataGeneration.genGeometricArray(numBlocks, totalNumPeople, geoProb);
			}
			
			boolean isPDP = false;
			double[] priorArray = new double[numBlocks];
			if (mechanism == Mechanism.SyntheticData) {
				for (int i=0; i<numBlocks; i++) {
					priorArray[i] = alpha;
				}
				isPDP = SyntheticDataMechanism.checkSyntheticData(inputArray, priorArray, epsilon, delta);
				//System.out.println(numBlocks + ", " + isPDP);
			}
			
			double totalDist = 0;
			
			for (int k=0; k<numRuns; k++) {
				if (mechanism == Mechanism.Laplace) {
					for (int i=0; i<numBlocks; i++) {
						outputArray[i] = Mechanisms.getLaplace(inputArray[i], laplaceStdDev, rand);
					}
					
					double l1Dist = L1distance.getL1Dist(inputArray, outputArray);
					//System.out.println("Distance=" + l1Dist);
					totalDist += l1Dist;
				}
				else if (mechanism == Mechanism.SyntheticData) {
					if (isPDP) {
						int[] output = SyntheticDataMechanism.synthesize(inputArray, priorArray, totalNumPeople, seed*k);
						double l1Dist  = L1distance.getL1Dist(inputArray, output);
						totalDist += l1Dist;
					}
				}
			}
			if (mechanism == Mechanism.Laplace)
				System.out.println("numBlocks=" + numBlocks + ", dist=" + totalDist/numRuns);
			else if (mechanism == Mechanism.SyntheticData)
				System.out.println("numBlocks=" + numBlocks + ", dist=" + totalDist/numRuns + ", " + isPDP);
		}
	}
	
	// generate k-by-k contingency table, add noise according to Laplace or Geometric 
	// mechanism, compute L1-distance.
	public static void testTableData_Laplace() {
		int numBlocks = 5;
		int minPerPair = 10;
		int maxPerPair = 100;
		
		int seed = 123456;
		Random rand = new Random(seed);
		int numRuns = 20;
		Mechanism mechanism = Mechanism.Laplace;
		// for Laplace mechanism
		double epsilon = 0.5;
		double sensitivity = 2;
		double laplaceStdDev = sensitivity/epsilon;
		// for geometric mechanism
		double alpha = Math.exp(-epsilon);
		
//		double[] epsilonArr = {0.1, 0.25, 0.5, 0.75, 1};
		double[] epsilonArr = {1};
		for (int l=0; l<epsilonArr.length; l++) {
			epsilon = epsilonArr[l];
		
//		int[] maxValues = {10, 50, 100, 200, 500, 1000};
//		for (int l=0; l<maxValues.length; l++) {
//			minPerBlock = maxPerBlock = maxValues[l];
			
			if (mechanism == Mechanism.Geometric) {
				alpha = Math.exp(-epsilon);
				//System.out.println("epsilon=" + epsilon + ", alpha=" + alpha);
			}
			else if (mechanism == Mechanism.Laplace) {
				laplaceStdDev = sensitivity/epsilon;
			}
			
			double totalDist = 0;
			int numPeople = 0;
			for (int k=0; k<numRuns; k++) {
				
				int[][] inputTable = DataGeneration.genUniformTable(numBlocks, minPerPair, maxPerPair, rand);
				
				double[][] outputTable = new double[numBlocks][numBlocks];
				
				for (int i=0; i<numBlocks; i++) {
					for (int j=0; j<numBlocks; j++) {
						numPeople += inputTable[i][j];
						if (mechanism == Mechanism.Laplace)
							//outputBlocks[i][j] = (int) Math.round(Mechanisms.getLaplace(blocks[i][j], laplaceStdDev, rand));
							outputTable[i][j] = Mechanisms.getLaplace(inputTable[i][j], laplaceStdDev, rand);
						else if (mechanism == Mechanism.Geometric) {
							outputTable[i][j] = Mechanisms.getGeometric((int) inputTable[i][j], alpha, rand);
						}
						else {
							System.out.println("Unknown mechanism");
							System.exit(1);
						}
					}
				}
				/*System.out.println("Input:");
				for (int i=0; i<numBlocks; i++) {
					for (int j=0; j<numBlocks; j++) {
						System.out.print(blocks[i][j] + "\t");
					}
					System.out.print("\n");
				}	
				System.out.println("Output:");
				for (int i=0; i<numBlocks; i++) {
					for (int j=0; j<numBlocks; j++) {
						System.out.print(outputBlocks[i][j] + "\t");
					}
					System.out.print("\n");
				}*/
				
				double l1Dist = L1distance.getL1Dist(inputTable, outputTable);
				//System.out.println("Distance=" + l1Dist);
				totalDist += l1Dist;
			}
			System.out.println("epsilon=" + epsilon + ", dist=" + totalDist/numRuns + ", averagePairSize=" + numPeople/numRuns/(numBlocks*numBlocks));
		}
	}
}
