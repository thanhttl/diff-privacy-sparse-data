package datastructures;


import java.util.Random;

import distributions.Geometric;


/**
 * for publishing scarse data, e.g., with many zero entries in the frequency matrix.
 * 1. add Laplacian noise for all non-zero entries.
 * 2. choose a subset of zero entries (#entries = B(n, p)) and upgrade to 
 * a value >= threshold \theta. 
 * n is the number of zeros in the original data.
 * p = Pr[ x>= \theta |x_0 = 0  ] = \frac{1-\alpha}{1+\alpha^2} \alpha^\theta
 * For \epsilon-df, \alpha = e^{-\epsilon}   
 * @author ttran, 7/2010
 *
 */
public class Filtering_Absolute extends Sampling {
	double theta;
	
	public Filtering_Absolute(double alpha, double theta, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		this.alpha = alpha;
		this.theta = theta;
		noiseRand = new Random(noiseSeed);
		sampleRand = new Random(sampleSeed);
		binomialRand = new Random(zeroPosSeed);
	}
	
//	public boolean isNonZeroSampled(double input) {
//		return true;
//	}
	
	public double sampleNonZeroEntry(double input) {
		double output = input + Geometric.getSample(noiseRand);
		if (Math.abs(output) < theta)
			output = 0;
		return output;
	}
	
	// sampling directly from the distribution P(x) = (1-\alpha)*\alpha^{b-\theta}
	public double sampleZeroEntry(double prob) {
		double sample;
		if (prob > 0.5)
		    sample = Math.ceil(Math.log(2*(1-prob))/Math.log(alpha) + theta - 1);
		else
			sample = Math.ceil(-theta - Math.log(2*prob)/Math.log(alpha));
		return sample;
	}
	
	public int getNumZerosSelected(int numZeros) {
		double thetaProb = 2*Math.pow(alpha, theta)/(1+alpha); //p
		//System.out.println("Theta Prob: " + thetaProb);
		// approximate the binomial dist w/ a Gaussian dist
		double variance = numZeros*thetaProb*(1-thetaProb);
		int numUpgrades = (int) Math.round(noiseRand.nextGaussian()*Math.sqrt(variance) + thetaProb*numZeros);
		while (numUpgrades < thetaProb*numZeros/5 || numUpgrades > numZeros)
			numUpgrades = (int) Math.round(noiseRand.nextGaussian()*Math.sqrt(variance) + thetaProb*numZeros);
		return numUpgrades;
	}
	
	public int getNumZerosSelected(long numZeros) {
		double thetaProb = 2*Math.pow(alpha, theta)/(1+alpha); //p
		
		double mean = thetaProb*numZeros;
		double stddev = Math.sqrt(numZeros*thetaProb*(1-thetaProb));
		int numUpgrades = (int) Math.round(noiseRand.nextGaussian()*stddev + mean);
		//while (numUpgrades < thetaProb*numZeros/5 || numUpgrades > numZeros)
		//	numUpgrades = (int) Math.round(noiseRand.nextGaussian()*Math.sqrt(variance) + thetaProb*numZeros);
		while (numUpgrades < mean-3*stddev || numUpgrades > mean+3*stddev)
			numUpgrades = (int) Math.round(noiseRand.nextGaussian()*stddev + mean);
		return numUpgrades;
	}
	
	
	
	public static void main(String[] args) {
		int INPUT_SIZE = 1000000;
		double[] input = new double[INPUT_SIZE];
		int i=0;
		Random rand = new Random(92929);
		for (; i< INPUT_SIZE*0.8; i++) {
			input[i] = 0;
		}
		for (; i<INPUT_SIZE; i++) {
			input[i] = 1 + rand.nextInt(100);
		}
		
		double epsilon = 1;
		double alpha = Math.exp(-epsilon);
		
		double[] thetas = {0, 5, 10, 20, 50};
		
		int seed = 123456;
		for (double theta: thetas ) {
			int numRepeats = 20;
			// compute the errors
			double L1Dist = 0;
			double L2Dist = 0;
			double totalCount = 0;
			
			for (int j=0; j<numRepeats; j++) {
				Filtering_Signed sampling = new Filtering_Signed(alpha, theta, seed, 98654, 3421);
				double[] output = new double[input.length];
				int sampleSize = sampling.anonymize_basic(input, output);
				
				for (int k=0; k<INPUT_SIZE; k++) {
					double estCount = output[k];
					//double trueCount = origData.get(i).count;
					double trueCount = input[k];
					totalCount += trueCount;
					L1Dist += Math.abs(estCount - trueCount);
					L2Dist += Math.pow(estCount - trueCount, 2);
				}
			}
			
			System.out.println("theta=" + theta + "; Size=" + INPUT_SIZE + "; TotalCount=" + totalCount/numRepeats + 
					"; L1Dist=" + L1Dist/numRepeats + "; L2Dist=" + Math.sqrt(L2Dist)/numRepeats);
		}
	}
	
}
