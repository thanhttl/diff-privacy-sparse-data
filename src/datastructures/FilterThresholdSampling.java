package datastructures;

import java.util.HashMap;
import java.util.Random;

import distributions.Geometric;

public class FilterThresholdSampling extends Sampling {
	
	// two parameters: alpha= exp(-e), 
	// D determines the sample size |S|, D= \sum M(i,j) / |S|
	final int theta; // theta >=1
	final double D;
	final double constantC;
	final double probGeqD;
	final double probLeqD;
	// store values between (-D, D) 
	HashMap<Integer, Double> cachedCDF;
	
	public FilterThresholdSampling(double alpha, int theta, double D, int noiseSeed, int sampleSeed) {
		this.alpha = alpha;
		this.theta = theta;
		this.D = D;
		noiseRand = new Random(noiseSeed); sampleRand = new Random(sampleSeed);
		binomialRand = new Random((int) (3*noiseSeed + 2*sampleSeed));
		constantC = 1.0 / (2*(theta*Math.pow(alpha,theta) - (theta-1)*Math.pow(alpha, theta+1) - Math.pow(alpha, D+1)));
		probGeqD = constantC*D*(1-alpha)*Math.pow(alpha, D);
		probLeqD = 1 - constantC*D*(1-alpha)* Math.pow(alpha, D+1);
		cachedCDF = new HashMap<Integer,Double>();
		
	}
	
	public double sampleNonZeroEntry(double input) {
		System.err.println("This method is empty. Cannot call it now");
		System.exit(1);
		return -1;
	}
	
	// sampling according to CDF, need a binary search over the CDF
	public double sampleZeroEntry(double prob) {
		double sample;
		if (prob <= probGeqD) // case 1
			sample = (int) Math.ceil(-Math.log(prob/(1-alpha)/(constantC*D))/Math.log(alpha));
		else if (prob >= probLeqD) // case 4
			sample = (int) Math.ceil(Math.log((1-prob)/(1-alpha)/(constantC*D))/Math.log(alpha) -1);
		else {
			int lower, upper;
			if (prob > 0.5) { // case 3, search between theta and D
				lower = theta; upper = (int) D+1;
			}
			else {
				lower = (int) -D-1; upper = -theta;
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
		
		return sample;
	}
	

	public int getNumZerosSelected(int numZeros) {
		double selectProb = 2*(theta*Math.pow(alpha, theta)- (theta-1)*Math.pow(alpha,theta+1)-Math.pow(alpha,D+1))/ (D*(1-alpha*alpha)); //p
		// check the conditions if this binomial dist can be approximated by a Gaussian
		// approximate the binomial dist w/ a Gaussian dist
		double mean = selectProb*numZeros;
		double variance = numZeros*selectProb*(1-selectProb);
		double stddev = Math.sqrt(variance);
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*Math.sqrt(variance) + selectProb*numZeros);
		while (numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev);
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		return numZerosSelected;
	}
	
	public int getNumZerosSelected(long numZeros) {
		double selectProb = 2*(theta*Math.pow(alpha, theta)- (theta-1)*Math.pow(alpha,theta+1)-Math.pow(alpha,D+1))/ (D*(1-alpha*alpha)); //p
		// check the conditions if this binomial dist can be approximated by a Gaussian
		// approximate the binomial dist w/ a Gaussian dist
		double mean = selectProb*numZeros;
		double variance = numZeros*selectProb*(1-selectProb);
		double stddev = Math.sqrt(variance);
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*Math.sqrt(variance) + selectProb*numZeros);
		while (numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev);
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		return numZerosSelected;
	}
	
	// x \in [-D, D]
	private double cdf(double x) {
		if (x < -D || x > D) {
			System.out.println("Wrong case");
			System.exit(1);
		}
		if (x <= -theta) {
			return constantC*(-x*Math.pow(alpha, -x) + (x+1)*Math.pow(alpha, -x+1) - Math.pow(alpha, D+1));
		}
		else {
			return 0.5 + constantC*(theta*Math.pow(alpha,theta) - (theta-1)*Math.pow(alpha, theta+1) -(x+1)*Math.pow(alpha,x+1)+ x*Math.pow(alpha,x+2));
		}
			
	}
}
