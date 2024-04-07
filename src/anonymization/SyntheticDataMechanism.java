package anonymization;

import java.util.Random;

import utility.L1distance;
import cc.mallet.util.Maths;
import distributions.Dirichlet;

// Privacy: Theory meets practice on the map (ICDE 2008)
public class SyntheticDataMechanism {
	public static final int SEED = 123456;
	// Algorithm (1): 
	// for all destination blocks d do:
	// 1. let n(d) be the histogram of origin blocks
	// 2. choose prior sample \alpha(d) with |\alpha(d)| = O(n(d))
	// 3. choose output sample size m = O(n(d))
	// 4. sample m points from multinomial distribution w/ prior D(n(d)_1+\alpha_1,...,n(d)_k+\alpha_k))
	
	/**
	 * originArr: the histogram of origin blocks corresponding to a destination
	 * priorArr: the prior distribution for this destination 
	 * numSamples: num of synthetic data points, m
	 */
	public static int[] synthesize(int[] originArr, double[] priorArr, int numSamples, int seed) {
		if (originArr.length != priorArr.length) {
			System.err.println("Origin vector and prior vector must have the same size.");
			System.exit(1);
		}
		double magnitude = 0;
		for (int i=0; i<priorArr.length; i++) {
			magnitude += (priorArr[i] + originArr[i]);
		}
		double[] params = new double[priorArr.length];
		for (int i=0; i<priorArr.length; i++) {
			params[i] = (double) (priorArr[i] + originArr[i])/magnitude;
			//System.out.println("Param " + i + ": " + params[i]);
		}
		
		Dirichlet priorDist = new Dirichlet(magnitude, params);
		priorDist.initRandom(seed);
		
		int[] syntheticData = priorDist.drawObservation(numSamples);
		return syntheticData;
	}
	
	public static boolean[] checkSyntheticData(int[][] inputTable, double[][] priorTable, double epsilon, double delta) {
		// how to satisfy the \epsilon requirement
		int numDestBlocks = inputTable.length;
		boolean[] isPDP = new boolean[numDestBlocks];
		
		for (int dest=0; dest<numDestBlocks; dest++) {
			isPDP[dest] = checkSyntheticData(inputTable[dest], priorTable[dest], epsilon, delta);
		}
		return isPDP;
	}
	
	
	public static boolean checkSyntheticData(int[] inputArray, double[] priorArray, double epsilon, double delta) {
		// how to satisfy the \epsilon requirement	
		int numOrigBlocks = inputArray.length;
		
		// numCommuter = m = n
		int numCommuters = 0;
		for (int orig=0; orig<numOrigBlocks; orig++) {
			numCommuters += inputArray[orig];
		}
		int numSyntheticPoints = numCommuters;
		
		// min(alpha(d)_i)
		double minAlpha = priorArray[0];
		double totalAlpha = priorArray[0];
		for (int orig=1; orig<numOrigBlocks; orig++) {
			minAlpha = Math.min(minAlpha, priorArray[orig]);
			totalAlpha += priorArray[orig];
		}
		
		double reference0Sample = Utility.reference0Sample(numCommuters, numSyntheticPoints, minAlpha, totalAlpha-minAlpha, Math.exp(epsilon)-1);
		//System.out.println("==== Dest " + dest);
		//System.out.println("Reference 0 Sample: " + reference0Sample);
		if (reference0Sample <= delta*(Math.exp(epsilon)-2)/(2*numOrigBlocks*Math.exp(epsilon))) {
			System.out.println("Mininum delta: " + reference0Sample/((Math.exp(epsilon)-2)/(2*numOrigBlocks*Math.exp(epsilon))));
			//System.out.println("Satisfies (" + epsilon + ", " + delta + ") p.d.p.");
			return true;
		}
		else {
			//System.out.println("NOT Satisfy (" + epsilon + ", " + delta + ") p.d.p.");
			return false;
		}
	}
	
	
	public static boolean checkSyntheticData(int[] inputArray, double alpha, double epsilon, double delta) {
		// how to satisfy the \epsilon requirement	
		int numOrigBlocks = inputArray.length;
		
		// numCommuter = m = n
		int numCommuters = 0;
		for (int orig=0; orig<numOrigBlocks; orig++) {
			numCommuters += inputArray[orig];
		}
		int numSyntheticPoints = numCommuters;
		
		// min(alpha(d)_i)
		double minAlpha = alpha;
		double totalAlpha = alpha*inputArray.length;
		
		
		double reference0Sample = Utility.reference0Sample(numCommuters, numSyntheticPoints, minAlpha, totalAlpha-minAlpha, Math.exp(epsilon)-1);
		//System.out.println("==== Dest " + dest);
		//System.out.println("Reference 0 Sample: " + reference0Sample);
		if (reference0Sample <= delta*(Math.exp(epsilon)-2)/(2*numOrigBlocks*Math.exp(epsilon))) {
			//System.out.println("Mininum delta: " + reference0Sample/((Math.exp(epsilon)-2)/(2*numOrigBlocks*Math.exp(epsilon))));
			//System.out.println("Satisfies (" + epsilon + ", " + delta + ") p.d.p.");
			return true;
		}
		else {
			//System.out.println("NOT Satisfy (" + epsilon + ", " + delta + ") p.d.p.");
			return false;
		}
	}
	
	public static double searchForAlpha(int[] inputArray, double initialGuess, double epsilon, double delta) {	
		if (checkSyntheticData(inputArray, initialGuess, epsilon, delta))
			return initialGuess;
		
		double lowerBound = initialGuess;
		double upperBound = 2*initialGuess;
		
		while (!checkSyntheticData(inputArray, upperBound, epsilon, delta)) {
			upperBound *= 2;
			if (upperBound > 10000) {
				System.out.println("cannot find alpha");
				return upperBound;
			}
		}
		
		double guess;	
		while (upperBound - lowerBound > initialGuess) {
			guess = (upperBound + lowerBound)/2;
	
			if (checkSyntheticData(inputArray, guess, epsilon, delta))
				upperBound = guess;
			else 
				lowerBound = guess;
		}
		return upperBound;
	}
	public static void main(String[] args) {		
		int[] numBlocksArray = {5, 10, 50, 100, 250};
		for (int numBlocks : numBlocksArray) {
			System.out.println("NUM BLOCKS: " + numBlocks);
			testVaryingAlpha(numBlocks, 500/numBlocks);
		}
		
	}
	
	
	// assuming distribution per destination is uniform.
	public static void testVaryingAlpha(int numBlocks, int numPerBlock) {
		int[] alphas = {10, 20, 50, 100, 200, 500};
		
		int[][] originTable = new int[1][numBlocks];
		for (int i=0; i<numBlocks; i++)
			originTable[0][i] = numPerBlock;
		
		double epsilon = 1;
		double delta = 1E-5; // 1E-5;
		double[][] priorTable = new double[originTable.length][originTable[0].length];
		for (int i=0; i<alphas.length; i++) {
			int alpha = alphas[i];
			for (int j=0; j<originTable.length; j++)
				for (int k=0; k<originTable[0].length; k++)
					priorTable[j][k] = alpha;
			
			boolean[] isPDP = checkSyntheticData(originTable, priorTable, epsilon, delta);
			
			
			Random rand = new Random(SEED);
			int numRuns = 20;
			
			for (int j=0; j<isPDP.length; j++) {
				if (isPDP[j]) {		
					double totalDistance = 0;
					for (int k=0; k<numRuns; k++) {
						int seed = rand.nextInt();
						int[] syntheticData = synthesize(originTable[j], priorTable[j], numBlocks*numPerBlock, seed);
						double distance = L1distance.getL1Dist(originTable[j], syntheticData);
						//System.out.println("Run " + j + ", Distance=" + distance);
						totalDistance += distance;
					}
				
					System.out.println("Alpha: " + alpha + ", Average distance: " + totalDistance/numRuns/isPDP.length);
				}
			}
		}
	}
	
	private static class Utility {
		// Gamma function is defined as \int_0^{+infty} x^{t-1} e^{-x} dx
		static double M;
		static double ALPHA1; 
		static double C;
		public static double reference0Sample(int n, int m, double alpha1, double alpha2, double c) {
			M = m; ALPHA1 = alpha1; C= c;
			double max = 0;
			for (int x=0; x<=n; x++) {
				double f_x = f(x);
//				System.out.println(Maths.logGamma(m+1));
//				System.out.println(Maths.logGamma(n+alpha1+alpha2));
//				System.out.println(Maths.logGamma(f_x+1));
//				System.out.println(Maths.logGamma(m-f_x+1));
//				System.out.println(Maths.logGamma(x+alpha1));
//				System.out.println(Maths.logGamma(n-x+alpha2));
				double nominator = Maths.logGamma(m+1) + Maths.logGamma(n+alpha1+alpha2) -
					Maths.logGamma(f_x+1) - Maths.logGamma(m-f_x+1) - Maths.logGamma(x+alpha1) - Maths.logGamma(n-x+alpha2);
				double denominator = Maths.logGamma(m+n+alpha1+alpha2) -
					Maths.logGamma(x+f_x+alpha1) - Maths.logGamma(n-x+m-f_x+alpha2);
				max = Math.max(max, Math.exp(nominator - denominator));
			}
			return max;
		}
		
		private static double f(double x) {
			return Math.min(M, C*(ALPHA1 + Math.max(x-1, 0)));
		}
	}
}
