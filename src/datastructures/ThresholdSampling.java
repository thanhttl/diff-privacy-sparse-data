package datastructures;

import java.util.HashMap;
import java.util.Random;

import distributions.Geometric;

public class ThresholdSampling extends Sampling {
	
	// two parameters: alpha= exp(-e), 
	// D determines the sample size |S|, D= \sum M(i,j) / |S|
	final double theta;
	final double constantC;
	final double probGeqTheta;
	final double probLeqTheta;
	// store values between (-D, D) 
	HashMap<Integer, Double> cachedCDF;
	final boolean adjustingWeights;
	
	public ThresholdSampling(double alpha, double theta, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		this(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed, true);
	}
	
	public ThresholdSampling(double alpha, double theta, int noiseSeed, int sampleSeed, int zeroPosSeed, boolean adjusting) {
		this.alpha = alpha;
		this.theta = theta;
		noiseRand = new Random(noiseSeed);
		sampleRand = new Random(sampleSeed);
		binomialRand = new Random(zeroPosSeed);
		constantC = 1.0 / (2*alpha*(1-Math.pow(alpha, theta)));
		probGeqTheta = constantC*theta* Math.pow(alpha, theta)*(1- alpha);
		probLeqTheta = 0.5 + constantC*(alpha-(theta+1)*Math.pow(alpha,theta+1)+theta*Math.pow(alpha,theta+2));
		cachedCDF = new HashMap<Integer,Double>();
		adjustingWeights = adjusting;		
	}

	
	
	// allow negative values in the output
	public double sampleNonZeroEntry(double input) {
		int noise = Geometric.getSample(noiseRand);
		double output = input + noise;
		if (sampleRand.nextDouble() <= Math.abs(output)/theta) {
			if (output>0 && output<theta)
				output = (int) theta;
			else if (output<0 && output>-theta)
				output = (int) -theta;
			return output;
		}
		else
			return 0;
	}
		
	// sampling according to CDF, need a binary search over the CDF
	public double sampleZeroEntry(double prob) {
		double sample;
		if (prob <= probGeqTheta) // case 1
			sample = (int) Math.ceil(-Math.log(prob/(1-alpha)/(constantC*theta))/Math.log(alpha));
		else if (prob >= probLeqTheta) // case 4
			sample = (int) Math.ceil(Math.log((1-prob)/(constantC*theta*(1-alpha)))/Math.log(alpha) - 1);
		else {
			int lower, upper;
			if (prob >= 0.5) { // case 3, search between 0 and D
				lower = 1; upper = (int) theta;
			}
			else {
				lower = (int) -theta; upper = 1;
			}
			while (upper - lower > 1) {
				int midPoint = (lower+upper)/2;
				double probMidPoint;
				if (!cachedCDF.containsKey(midPoint)) {
					probMidPoint = cdf(midPoint);
					cachedCDF.put(midPoint, probMidPoint);
				}
				else {
					probMidPoint = cachedCDF.get(midPoint);
				}
				if (probMidPoint == prob)
					return midPoint;
				else if (probMidPoint < prob) 
					lower = midPoint+1;
				else
					upper = midPoint-1;
			}
			sample = upper;
		}
		
		// adjustingWeights for threshold, non adjusting weights for calling from priority.
		if (adjustingWeights && Math.abs(sample) < theta) {
			if (sample < 0) sample = (int) -theta;
			else sample = (int) theta;
		}
		return sample;
	}
	

	public int getNumZerosSelected(int numZeros) {
		double selectProb = 2*alpha*(1-Math.pow(alpha,theta))/ (theta*(1-alpha*alpha)); //p
		
		// approximate the binomial dist w/ a Gaussian dist
		double variance = numZeros*selectProb*(1-selectProb);
		double stddev = Math.sqrt(variance);
		double mean = selectProb*numZeros;
		int numUpgrades = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		while (numUpgrades < mean - 3*stddev || numUpgrades > mean + 3*stddev) {
			numUpgrades = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		}	
		return numUpgrades;
	}
	
	public int getNumZerosSelected(long numZeros) {
		double selectProb = 2*alpha*(1-Math.pow(alpha,theta))/ (theta*(1-alpha*alpha)); //p
		
		// approximate the binomial dist w/ a Gaussian dist
		double variance = numZeros*selectProb*(1-selectProb);
		double stddev = Math.sqrt(variance);
		double mean = selectProb*numZeros;
		int numUpgrades = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		while (numUpgrades < mean - 3*stddev || numUpgrades > mean + 3*stddev) {
			numUpgrades = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		}	
		return numUpgrades;
	}
	
	// x \in [-theta, theta]
	private double cdf(double x) {
		if (x < -theta || x > theta) {
			System.out.println("Wrong case");
			System.exit(1);
		}
		if (x <= 0) {
			return constantC*(-x*Math.pow(alpha, -x) + (x+1)*Math.pow(alpha, -x+1) - Math.pow(alpha, theta+1));
		}
		else {
			return 0.5 + constantC*alpha*(1-(x+1)*Math.pow(alpha,x)+x*Math.pow(alpha,x+1));
		}
			
	}
	
	public static void main(String[] args) {
		double alpha = Math.exp(-0.1);
		int D = 1000;
		int noiseSeed = 123456; int sampleSeed = 45323;
		ThresholdSampling wrs = new ThresholdSampling(alpha, D, noiseSeed, sampleSeed, 2345);
		double[] probs = {0.1, 0.3, 0.5, 0.7, 0.9};
		for (double prob: probs) {
			double x = wrs.sampleZeroEntry(prob);
			System.out.println(x);
		}
	}

}
