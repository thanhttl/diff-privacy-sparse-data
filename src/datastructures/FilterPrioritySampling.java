package datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import common.QuickSelect;

import datastructures.UtilityObject.NonZeroEntry;
import datastructures.UtilityObject.NonZeroEntry_Long;
import datastructures.UtilityObject.ZeroIndexInterval;
import datastructures.PrioritySampling.PriorityComparator;
import distributions.Geometric;

public class FilterPrioritySampling extends PrioritySampling {
	int theta;
//	int K;
//	public double omega;
//	double probGeqOmega;
//	int noiseSeed; int sampleSeed; int zeroPosSeed;
//	public static boolean useQuickSelect = true;
	
//	public FilterPrioritySampling(double alpha, int K, int noiseSeed, int sampleSeed, int zeroPosSeed) {
//		this.alpha = alpha;
//		this.K = K;
//		this.noiseSeed = noiseSeed; this.sampleSeed = sampleSeed; this.zeroPosSeed = zeroPosSeed;
//		noiseRand = new Random(noiseSeed);
//		sampleRand = new Random(sampleSeed);
//		binomialRand = new Random(zeroPosSeed);
//		omega = 0;	
//	}
	
	public FilterPrioritySampling(double alpha, int theta, int K, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		super(alpha, K, noiseSeed, sampleSeed, zeroPosSeed);
		this.theta = theta;	
	}
	
	public int anonymize(NonZeroEntry[] nonZeros, int[] zeroIndices, double[] output) {
		return anonymize(nonZeros, zeroIndices, zeroIndices.length, output);
	}
	
	public int anonymize(NonZeroEntry[] nonZeros, int[] zeroIndices, int numZeros, double[] output) {
	//	int numPassThreshold = 0;
		NonZeroEntry[] passThresholdArray = new NonZeroEntry[nonZeros.length];
		int numPassed = 0;
		for (NonZeroEntry entry: nonZeros) {
			entry.value += Geometric.getSample(noiseRand);
			if (Math.abs(entry.value) >= theta) {
				entry.rand = sampleRand.nextDouble();
				passThresholdArray[numPassed++] = entry;
			}
		}
	
		Comparator comparator = new PriorityComparator(); 
		NonZeroEntry[] passThresholdNonZeros = null; 
		
		passThresholdNonZeros = null; 
		
 		if (!useQuickSelect) {
 			passThresholdNonZeros = new NonZeroEntry[numPassed];
 			for (int k=0; k<numPassed; k++)
 				passThresholdNonZeros[k++] = passThresholdArray[k];
 	
 			Arrays.sort(passThresholdNonZeros, comparator);	
 		}
 			
 		double probGeqOmega = K*(1.0)/numZeros;
 		
 		int numZerosSelected;
 
 		if (probGeqOmega > 1- cdf(theta-1)  ) {
//			System.out.println("Expected num zeros upgraded is smaller than num required.");
			probGeqOmega = (1-cdf(theta-1)); 
			double variance = numZeros*probGeqOmega*(1-probGeqOmega);
			numZerosSelected = (int) Math.round(noiseRand.nextGaussian()*Math.sqrt(variance) + probGeqOmega*numZeros);
			if (numZerosSelected > numZeros)
				numZerosSelected = numZeros;
			omega= theta;
		}
 		else {
 			omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
 			//System.out.println("Prob>=omega: " + probGeqOmega + ", omega=" + omega);
 			numZerosSelected = getNumZerosSelected(numZeros); // this need to be >= K
 			while (numZerosSelected <= K  || numZerosSelected > 2*K)
 				numZerosSelected = getNumZerosSelected(numZeros);
 		}
		
		FilterThresholdSampling wrs = new FilterThresholdSampling(alpha, theta, (int) omega, noiseSeed, sampleSeed);
		
		NonZeroEntry[] zeroArray;
		
		HashSet<Integer> sampledIndices =null;
		if (rejectionSampling)
			sampledIndices =  new HashSet<Integer>(numZerosSelected);
		
		if (!useQuickSelect)
			zeroArray = new NonZeroEntry[numZerosSelected];
		else 
			zeroArray = new NonZeroEntry[numZerosSelected + numPassed];
		
		if (!rejectionSampling) {
			int numShuffles = zeroIndices.length/4;
			for (int i=0; i<numShuffles; i++) {
				int position = 4*noiseRand.nextInt(numShuffles);
				int tempIndex = zeroIndices[position];
				zeroIndices[position] = zeroIndices[i];
				zeroIndices[i] = tempIndex;
				
			}
		}
		
		
		for (int i=0; i<numZerosSelected; i++) {
			int index; 
			if (!rejectionSampling)
				index = zeroIndices[i];
			else {
				index = zeroIndices[noiseRand.nextInt(numZeros)];
				while (sampledIndices.contains(index))
					index = zeroIndices[noiseRand.nextInt(numZeros)];
				sampledIndices.add(index);
			}
			double sample;
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double prob = noiseRand.nextDouble();
			sample = wrs.sampleZeroEntry(prob);
			double randNum;
			if (Math.abs(sample) >= omega)
				randNum = sampleRand.nextDouble();
			else
				randNum = sampleRand.nextDouble()*Math.abs(sample)/omega;
			
			zeroArray[i] = ObjectPool.allocateNonZeroEntry();
			zeroArray[i].setPosition(index);
			zeroArray[i].setValue(sample);
			zeroArray[i].setRand(randNum);			
		}
	
		if (!useQuickSelect) {
			Arrays.sort(zeroArray, comparator);
		
			int sampleSize = Math.min(K, numZerosSelected+passThresholdNonZeros.length);
			
			// merge the two arrays based on their priorities
			NonZeroEntry[] merged = new NonZeroEntry[K];
			
			int currNonZeroIndex = passThresholdNonZeros.length-1;
			int currZeroIndex = zeroArray.length-1;
			
			
			for (int i=0; i< Math.min(K,sampleSize); i++) {
				if (currNonZeroIndex < 0)
					merged[i] = zeroArray[currZeroIndex--];
				else if (currZeroIndex < 0)
					merged[i] = passThresholdNonZeros[currNonZeroIndex--];
				else {
					if (comparator.compare(passThresholdNonZeros[currNonZeroIndex], zeroArray[currZeroIndex])  > 0) {
						merged[i] = passThresholdNonZeros[currNonZeroIndex--];
					}
					else {
						merged[i] = zeroArray[currZeroIndex--];
					}
				}		
			}
			
			if (sampleSize == K) {
				// find the (K+1)-th priority
				double priority;
				if (currNonZeroIndex < 0)
					priority = zeroArray[currZeroIndex].getPriority();
				else if (currZeroIndex < 0)
					priority = passThresholdNonZeros[currNonZeroIndex].getPriority();
				else if (zeroArray[currZeroIndex].getPriority() > passThresholdNonZeros[currNonZeroIndex].getPriority())
					priority = zeroArray[currZeroIndex].getPriority();
				else
					priority = passThresholdNonZeros[currNonZeroIndex].getPriority();
			
				
				for (int i=0; i<K; i++) {
					NonZeroEntry entry = merged[i];
					// not adjusting weights
					
					if (Math.abs(entry.value) > priority) 
						output[entry.getPosition()] = entry.value;
					else {
						if (entry.value > 0)
							output[entry.getPosition()] = priority;
						else
							output[entry.getPosition()] = -priority;
					}	
				}
			}
			else {
				for (int i=0; i<sampleSize; i++) {
					NonZeroEntry entry = merged[i];
					output[entry.getPosition()] = entry.value;
				}
			}
			return sampleSize;
		}
		else { // Quick select
			for (int i=numZerosSelected; i<zeroArray.length; i++)
				zeroArray[i] = passThresholdArray[i-numZerosSelected];
			
			int numEntries = zeroArray.length;
			if (K < numEntries) {
				QuickSelect.quickselect(zeroArray, numEntries-1-K);
				double kPlus1Priority = zeroArray[numEntries-1-K].getPriority();
				//System.out.println("K+1 Priority=" + kPlus1Priority);
				for (int i=numEntries-K; i<numEntries; i++) {
					NonZeroEntry entry = zeroArray[i];
						
					if (Math.abs(entry.value) > kPlus1Priority) 
						output[entry.getPosition()] = entry.value;
					else {
						if (entry.value > 0)
							output[entry.getPosition()] = kPlus1Priority;
						else
							output[entry.getPosition()] = -kPlus1Priority;
					}		
				}
				return K;
			}
			else {
				for (int i=0; i<numEntries; i++) {
					NonZeroEntry entry = zeroArray[i];
					output[entry.getPosition()] = entry.value;
				}
			}
			return numEntries;
		
		}
	}
	
	public int anonymize(ArrayList<NonZeroEntry> nonZeros, ArrayList<ZeroIndexInterval> zeroIndices, int inputSize, ArrayList<NonZeroEntry> output) {
		int numNonZeros = nonZeros.size();
		int numZeros = inputSize - numNonZeros;
		
		ArrayList<NonZeroEntry> selectedArray = new ArrayList<NonZeroEntry>();
		HashSet<Integer> nonZeroIndices = new HashSet<Integer>();
		
		for (int i=0; i< numNonZeros; i++) {
			NonZeroEntry entry = nonZeros.get(i);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (Math.abs(sample) >= theta) {
				NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				outEntry.setRand(sampleRand.nextDouble());
				selectedArray.add(outEntry);
			}
			nonZeroIndices.add(entry.getPosition());
		}
		
//		double probGeqOmega = K*(1.0)/numZeros;
// 		
// 		int numZerosSelected = getNumZerosSelected(numZeros);
// 
		HashSet<Integer> sampledIndices = sampleZeroIndices(nonZeroIndices, zeroIndices, inputSize, noiseRand);
		
		FilterThresholdSampling wrs = new FilterThresholdSampling(alpha, theta, (int) omega, noiseSeed, sampleSeed);
		
		Iterator<Integer> sampledIndexIterator = sampledIndices.iterator();
		while (sampledIndexIterator.hasNext()) {
			int index = sampledIndexIterator.next();			
			
			double prob = noiseRand.nextDouble();
			double sample = wrs.sampleZeroEntry(prob);
			double randNum;
			if (Math.abs(sample) >= omega)
				randNum = sampleRand.nextDouble();
			else
				randNum = sampleRand.nextDouble()*Math.abs(sample)/omega;
			
			NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			outEntry.setRand(randNum);
			selectedArray.add(outEntry);
		}		
		
		getKItems(selectedArray, output);
		
		return output.size();
	}
	
	public int anonymize_basic(int[] input, double[] output) {
		ArrayList<NonZeroEntry> tempArray = new ArrayList<NonZeroEntry>();
		for (int i=0; i<input.length; i++) {
			double val = sampleNonZeroEntry(input[i]);
			if (Math.abs(val) >= theta) {
				double rand = sampleRand.nextDouble();
				NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
				entry.setPosition(i);
				entry.setValue(val);
				entry.setRand(rand);
				tempArray.add(entry);
			}
		}
		if (tempArray.size() < K) {
			for (NonZeroEntry entry: tempArray) {
				output[entry.getPosition()] = entry.value;
			}
			return tempArray.size();
		}
		ArrayList<NonZeroEntry> outArray = new ArrayList<NonZeroEntry>();
		getKItems(tempArray, outArray);
		for (NonZeroEntry entry: outArray) {
			output[entry.getPosition()] = entry.value;
		}
		ObjectPool.returnNonZeroEntries(tempArray);
		ObjectPool.returnNonZeroEntries(outArray);
		tempArray.clear();
		outArray.clear();
		return K;
	}
	
	
	public int getNumZerosSelected(int numZeros) {
		double probGeqOmega = K*(1.0)/numZeros;
		if (probGeqOmega > 1- cdf(theta-1)  ) {
//			System.out.println("Expected num zeros upgraded is smaller than num required.");
//			probGeqOmega = (1-cdf(theta-1));
			omega= theta;
			probGeqOmega = 1- cdf(theta-1);
		}
		else {
			omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		}
		
		double mean = probGeqOmega*numZeros;
		double stddev = Math.sqrt(numZeros*probGeqOmega*(1-probGeqOmega));
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
//		while (numZerosSelected < K || numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev)
		while (numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev)
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + probGeqOmega*numZeros);
		
		return numZerosSelected;
	
	}
	
	public int getNumZerosSelected(long numZeros) {
		double probGeqOmega = K*(1.0)/numZeros;
		if (probGeqOmega > 1- cdf(theta-1)  ) {
			omega = theta;
			probGeqOmega = 1- cdf(theta-1);
		}
		else {
			omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		}
		
		double mean = probGeqOmega*numZeros;
		double stddev = Math.sqrt(numZeros*probGeqOmega*(1-probGeqOmega));
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);

		while (numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev)
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + probGeqOmega*numZeros);
		
		return numZerosSelected;
	
	}
		
	// allow negative values in the output
	public double sampleNonZeroEntry(double input) {
		int noise = Geometric.getSample(noiseRand);
		double output = input + noise;
		return output;
	}
	
	// change this..
	public double sampleZeroEntry(double prob) {
		System.out.println("should call from FilterThresholdSampling");
		System.exit(1);
		return 0;
	}
	
	public int anonymize(ArrayList<NonZeroEntry_Long> nonZeros, long inputSize, ArrayList<NonZeroEntry_Long> output) {
		int numNonZeros = nonZeros.size();
		long numZeros = inputSize - numNonZeros;
		
		ArrayList<NonZeroEntry_Long> selectedArray = new ArrayList<NonZeroEntry_Long>();
		HashSet<Long> nonZeroIndices = new HashSet<Long>();
		
		for (int i=0; i< numNonZeros; i++) {
			NonZeroEntry_Long entry = nonZeros.get(i);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (Math.abs(sample) >= theta) {
				NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				outEntry.setRand(sampleRand.nextDouble());
				selectedArray.add(outEntry);
			}
			nonZeroIndices.add(entry.getPosition());
		}
		
//		double probGeqOmega = K*(1.0)/numZeros;
// 		
// 		int numZerosSelected = getNumZerosSelected(numZeros);
// 
		HashSet<Long> sampledIndices = sampleZeroIndices(nonZeroIndices, inputSize, noiseRand);
		
		FilterThresholdSampling wrs = new FilterThresholdSampling(alpha, theta, (int) omega, noiseSeed, sampleSeed);
		
		Iterator<Long> sampledIndexIterator = sampledIndices.iterator();
		while (sampledIndexIterator.hasNext()) {
			long index = sampledIndexIterator.next();			
			
			double prob = noiseRand.nextDouble();
			double sample = wrs.sampleZeroEntry(prob);
			double randNum;
			if (Math.abs(sample) >= omega)
				randNum = sampleRand.nextDouble();
			else
				randNum = sampleRand.nextDouble()*Math.abs(sample)/omega;
			
			NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			outEntry.setRand(randNum);
			selectedArray.add(outEntry);
		}		
		
		getKItems_Long(selectedArray, output);
		
		return output.size();
	}
	
	// solve for y, given Pr(Y<=y) = prob
	private double searchForOmega(double prob) {
		double unit = 0.1; // assume y is discretized.
		double sample;
		double lower = omega;
		double upper = 2*omega;
		if (upper ==0) {
			upper = 100;
			omega = 50;
		}
		double upperCdf = cdf(upper);
		while (upperCdf < prob) {
			upper += omega;
			upperCdf = cdf(upper);
		}
			
		while (upper - lower > unit) {
			double midPoint = (lower+upper)/2;
			double probMidPoint = cdf(midPoint);
			
			if (probMidPoint == prob)
				return midPoint;
			else if (probMidPoint < prob) 
				lower = midPoint+unit;
			else
				upper = midPoint;
		}
		sample = upper;
		
		// pick r from [0, 1];
		return sample;
	}

	// Pr(Y <= y), y>=0
	protected double cdf(double y) {
		double probGeqYPlus1 = 2*(theta*Math.pow(alpha, theta)- (theta-1)*Math.pow(alpha,theta+1)-Math.pow(alpha,y+2))/ ((y+1)*(1-alpha*alpha)); //p
		//double probGeqYPlus1 = 2*alpha*(1- Math.pow(alpha, y+1))/((y+1)*(1-alpha*alpha));
		return 1 - probGeqYPlus1;
			
	}
	
//	public void getKItems(ArrayList<NonZeroEntry> selectedArray, ArrayList<NonZeroEntry> output) {
//		int numEntries = selectedArray.size();
//		NonZeroEntry[] tempArray = new NonZeroEntry[selectedArray.size()];
//		for (int i=0; i<tempArray.length; i++)
//			tempArray[i] = selectedArray.get(i);
//		
//		QuickSelect.quickselect(tempArray, numEntries-K);
//		double kPlus1Priority = tempArray[numEntries-K-1].getPriority();
////		
////		Comparator comparator = new PriorityComparator();
////		Collections.sort(selectedArray, comparator);
////		double kPlus1Priority = selectedArray.get(numEntries-K-1).getPriority();
//		
//		for (int i=0; i<numEntries; i++) {
////		for (int i=numEntries-K; i<numEntries; i++) {
//			NonZeroEntry entry = selectedArray.get(i);
//
//			if (entry.getPriority() <= kPlus1Priority)
//				continue;
//			
//				NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
//				outEntry.position = entry.getPosition();
//				
//				if (Math.abs(entry.value) >= kPlus1Priority) 
//					outEntry.setValue(entry.value);
//				else {
//					if (entry.value > 0)
//						outEntry.setValue( kPlus1Priority);
//					else
//						outEntry.setValue( -kPlus1Priority);
//				}
//				output.add(outEntry);
////			}
//		}
//	}
//	
//	public void getKItems_Long(ArrayList<NonZeroEntry_Long> selectedArray, ArrayList<NonZeroEntry_Long> output) {
//		int numEntries = selectedArray.size();
//		NonZeroEntry_Long[] tempArray = new NonZeroEntry_Long[selectedArray.size()];
//		for (int i=0; i<tempArray.length; i++)
//			tempArray[i] = selectedArray.get(i);
//		
//		QuickSelect.quickselect(tempArray, numEntries-K);
//		double kPlus1Priority = tempArray[numEntries-K-1].getPriority();
////		
////		Comparator comparator = new PriorityComparator();
////		Collections.sort(selectedArray, comparator);
////		double kPlus1Priority = selectedArray.get(numEntries-K-1).getPriority();
//		
//		for (int i=0; i<numEntries; i++) {
////		for (int i=numEntries-K; i<numEntries; i++) {
//			NonZeroEntry_Long entry = selectedArray.get(i);
//
//			if (entry.getPriority() <= kPlus1Priority)
//				continue;
//			
//				NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
//				outEntry.position = entry.getPosition();
//				
//				if (Math.abs(entry.value) >= kPlus1Priority) 
//					outEntry.setValue(entry.value);
//				else {
//					if (entry.value > 0)
//						outEntry.setValue( kPlus1Priority);
//					else
//						outEntry.setValue( -kPlus1Priority);
//				}
//				output.add(outEntry);
////			}
//		}
//	}
}
