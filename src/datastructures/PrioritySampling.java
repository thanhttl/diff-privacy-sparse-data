package datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Random;

import common.QuickSelect;

import datastructures.UtilityObject.NonZeroEntry;
import datastructures.UtilityObject.NonZeroEntry_Long;
import datastructures.UtilityObject.ZeroIndexInterval;
import distributions.Geometric;

public class PrioritySampling extends Sampling {
	int K;
	
	public double omega;
	double probGeqOmega;
	int noiseSeed; int sampleSeed; int zeroPosSeed;
	public static boolean useQuickSelect = true;
	//boolean adjusting = false;
	
	public PrioritySampling(double alpha, int K, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		this.alpha = alpha;
		this.K = K;
		this.noiseSeed = noiseSeed; this.sampleSeed = sampleSeed; this.zeroPosSeed = zeroPosSeed;
		noiseRand = new Random(noiseSeed);
		sampleRand = new Random(sampleSeed);
		binomialRand = new Random(zeroPosSeed);
		omega = 0;	
	}
	
/*	
	public int anonymize(int[] input, int[] output) {
		ArrayList<Integer> zeroIndices = new ArrayList<Integer>();
		ArrayList<NonZeroEntry> nonZeroEntries = new ArrayList<NonZeroEntry>();
		
		int numZeros = 0;
		
		for (int i=0; i<input.length; i++) {
			if (input[i] == 0) {
				numZeros ++;
				output[i] = 0;
				zeroIndices.add(i);
			}
			else {
				nonZeroEntries.add(new NonZeroEntry(i, input[i], 0));
			}
		}
	
		
		for (NonZeroEntry entry: nonZeroEntries) {
			entry.value += Geometric.getSample(alpha, noiseRand);
			entry.rand = sampleRand.nextDouble();
		}
		
		NonZeroEntry[] nonZeroArray = new NonZeroEntry[nonZeroEntries.size()];
		for (int i=0; i<nonZeroArray.length; i++)
			nonZeroArray[i] = nonZeroEntries.get(i);
		
		//System.out.println("Processed non-zero");
		Comparator comparator = new PriorityComparator();
 		Arrays.sort(nonZeroArray, comparator);	
 		
 		probGeqOmega = K*(1.0)/numZeros;
		omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		if (omega == 0)
			omega = 1;

		ThresholdSampling wrs = new ThresholdSampling(alpha, (int) omega, noiseSeed, sampleSeed, false);
		
		int numZerosSelected = getNumZerosSelected(numZeros); // this need to be > K
		while ((nonZeroArray.length + numZerosSelected) <= K || numZerosSelected > 2*K)
			numZerosSelected = getNumZerosSelected(numZeros);
		
		// System.out.println("Num zeros selected: " + numZerosSelected);
 		// int numZerosSelected = K+1;	
		// randomly sample from the list of zero indices
		//Collections.shuffle(zeroIndices, new Random(12498));
		Collections.shuffle(zeroIndices, new Random(noiseSeed+sampleSeed/3));
		
		NonZeroEntry[] zeroArray = new NonZeroEntry[numZerosSelected];
		
		for (int i=0; i<numZerosSelected; i++) {	
			int index = zeroIndices.get(i);
			
			int sample;
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double prob = noiseRand.nextDouble();
			sample = wrs.sampleZeroEntry(prob);
			double randNum;
			if (Math.abs(sample) >= omega)
				randNum = sampleRand.nextDouble();
			else
				randNum = sampleRand.nextDouble()*Math.abs(sample)/omega;
			
//			if (!noiseRand.nextBoolean())
//				sample = -sample;
			zeroArray[i] = new NonZeroEntry(index, sample, randNum);			
		}
		
		Arrays.sort(zeroArray, comparator);
		// merge the two arrays based on their priorities
		NonZeroEntry[] merged = new NonZeroEntry[K];
		
		int currNonZeroIndex = nonZeroArray.length-1;
		int currZeroIndex = zeroArray.length-1;
		
		for (int i=0; i<K; i++) {
			if (currNonZeroIndex < 0)
				merged[i] = zeroArray[currZeroIndex--];
			else if (currZeroIndex < 0)
				merged[i] = nonZeroArray[currNonZeroIndex--];
			else {
				if (comparator.compare(nonZeroArray[currNonZeroIndex], zeroArray[currZeroIndex])  > 0) {
					merged[i] = nonZeroArray[currNonZeroIndex--];
				}
				else {
					merged[i] = zeroArray[currZeroIndex--];
				}
			}
				
		}
		
		// find the (K+1)-th priority
		double priority;
		if (currNonZeroIndex < 0)
			priority = zeroArray[currZeroIndex].getPriority();
		else if (currZeroIndex < 0)
			priority = nonZeroArray[currNonZeroIndex].getPriority();
		else if (zeroArray[currZeroIndex].getPriority() > nonZeroArray[currNonZeroIndex].getPriority())
			priority = zeroArray[currZeroIndex].getPriority();
		else
			priority = nonZeroArray[currNonZeroIndex].getPriority();
	
		
		for (int i=0; i<K; i++) {
			NonZeroEntry entry = merged[i];
		
			// adjusting weights
			if (Math.abs(entry.value) > priority) 
				output[entry.getPosition()] = entry.value;
			else {
				if (entry.value > 0)
					output[entry.getPosition()] = (int) priority;
				else
					output[entry.getPosition()] = (int) -priority;
			}	
		}
		
		return K;
	}
*/	
	public int anonymize(ArrayList<NonZeroEntry> nonZeros, ArrayList<ZeroIndexInterval> zeroIndices, int inputSize, ArrayList<NonZeroEntry> output) {
		int numZeros = inputSize - nonZeros.size();
		probGeqOmega = K*(1.0)/numZeros;
		omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		if (omega == 0)
			omega = 1;
		
		
		ArrayList<NonZeroEntry> selectedArray = new ArrayList<NonZeroEntry>();
			
		HashSet<Integer> nonZeroIndices = new HashSet<Integer>();
		for (NonZeroEntry entry: nonZeros) {
			double newVal =  sampleNonZeroEntry(entry.value);
			if (newVal != 0) {
				NonZeroEntry noisyEntry = ObjectPool.allocateNonZeroEntry();
				noisyEntry.setValue(newVal);
				noisyEntry.setPosition(entry.getPosition());
				noisyEntry.setRand(sampleRand.nextDouble());
				selectedArray.add(noisyEntry);
			}
			nonZeroIndices.add(entry.getPosition());
		}
		
		int numZerosSelected = getNumZerosSelected(numZeros); // this need to be > K
				
		ThresholdSampling wrs = new ThresholdSampling(alpha, (int) omega, noiseSeed, sampleSeed, zeroPosSeed, false);
		
		HashSet<Integer> sampledIndices =  new HashSet<Integer>((int) (1.2*numZerosSelected));
		                            
		for (int i=0; i<numZerosSelected; i++) {	
			int index; 
			index = binomialRand.nextInt(inputSize);
			while (sampledIndices.contains(index) || nonZeroIndices.contains(index)) {
				index = binomialRand.nextInt(inputSize);
			}
			sampledIndices.add(index);
			
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
			
			NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
			entry.setPosition(index);
			entry.setValue(sample);
			entry.setRand(randNum);
			selectedArray.add(entry);
		}
		
		getKItems(selectedArray, output);
		ObjectPool.returnNonZeroEntries(selectedArray);
		selectedArray.clear();
		return output.size();	
	}
	
	public void getKItems(ArrayList<NonZeroEntry> selectedArray, ArrayList<NonZeroEntry> output) {
		int numEntries = selectedArray.size();
		NonZeroEntry[] tempArray = new NonZeroEntry[selectedArray.size()];
		for (int i=0; i<tempArray.length; i++)
			tempArray[i] = selectedArray.get(i);
		
		QuickSelect.quickselect(tempArray, numEntries-K);
		double kPlus1Priority = tempArray[numEntries-K-1].getPriority();
//		
//		Comparator comparator = new PriorityComparator();
//		Collections.sort(selectedArray, comparator);
//		double kPlus1Priority = selectedArray.get(numEntries-K-1).getPriority();
		
		for (int i=0; i<numEntries; i++) {
//		for (int i=numEntries-K; i<numEntries; i++) {
			NonZeroEntry entry = selectedArray.get(i);

			if (entry.getPriority() <= kPlus1Priority)
				continue;
			
				NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
				outEntry.position = entry.getPosition();
				
				if (Math.abs(entry.value) >= kPlus1Priority) 
					outEntry.setValue(entry.value);
				else {
					if (entry.value > 0)
						outEntry.setValue( kPlus1Priority);
					else
						outEntry.setValue( -kPlus1Priority);
				}
				output.add(outEntry);
//			}
		}
	}
	
	public int anonymize_basic(int[] input, double[] output) {
		ArrayList<NonZeroEntry> entries = new ArrayList<NonZeroEntry>();
		for (int i=0; i<input.length; i++) {
			NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
			entry.setPosition(i);
			entry.setValue(input[i] + Geometric.getSample(noiseRand));			
			entry.rand = sampleRand.nextDouble();
			entries.add(entry);
		}

		Comparator comparator = new PriorityComparator();
 		Collections.sort(entries, comparator);	
 		
 		// find the (K+1)-th priority
 		double kPlus1Priority = entries.get(entries.size()-K-1).getPriority();
 		
		// System.out.println("currZeroIndex= " + currZeroIndex);
		
		for (int i=0; i<K; i++) {
			NonZeroEntry entry = entries.get(entries.size()-i-1);
			
		//	 adjusting weights
			if (Math.abs(entry.value) >= kPlus1Priority) 
				output[entry.getPosition()] = entry.value;
			else {
				if (entry.value > 0)
					output[entry.getPosition()] = kPlus1Priority;
				else
					output[entry.getPosition()] = - kPlus1Priority;
 			}
		}
		
		return K;
	}
	
	public int getNumZerosSelected(int numZeros) {
		probGeqOmega = K*(1.0)/numZeros;
		omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		if (omega == 0)
			omega = 1;
		double variance = numZeros*probGeqOmega*(1-probGeqOmega);
		double stddev = Math.sqrt(variance);
		double mean = probGeqOmega*numZeros;
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		while (numZerosSelected < K || numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev)
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		
		return numZerosSelected;
	}
	
	public int getNumZerosSelected(long numZeros) {
		probGeqOmega = K*(1.0)/numZeros;
		omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		if (omega == 0)
			omega = 1;
		double variance = numZeros*probGeqOmega*(1-probGeqOmega);
		double stddev = Math.sqrt(variance);
		double mean = probGeqOmega*numZeros;
		int numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		while (numZerosSelected < mean-3*stddev || numZerosSelected > mean+3*stddev)
			numZerosSelected = (int) Math.round(binomialRand.nextGaussian()*stddev + mean);
		
		return numZerosSelected;
	}
	
	
//	public boolean isNonZeroSampled(double input) {		
//		return true;
//	}
	
	// allow negative values in the output
	public double sampleNonZeroEntry(double input) {
		int noise = Geometric.getSample(noiseRand);
		double output = input + noise;
		return output;
	}
	
	// change this..
	public double sampleZeroEntry(double prob) {
		System.out.println("should call from ThresholdSampling");
		System.exit(1);
		return 0;
	}
	
	
	public int anonymize(ArrayList<NonZeroEntry_Long> nonZeros, long inputSize, ArrayList<NonZeroEntry_Long> output) {
		long numZeros = inputSize - nonZeros.size();
		probGeqOmega = K*(1.0)/numZeros;
		omega = Math.round(searchForOmega(1-probGeqOmega)); // this is prob(Y<= omega-1);
		if (omega == 0)
			omega = 1;
		
		
		ArrayList<NonZeroEntry_Long> selectedArray = new ArrayList<NonZeroEntry_Long>();	
		HashSet<Long> nonZeroIndices = new HashSet<Long>();
		
		for (NonZeroEntry_Long entry: nonZeros) {
			double newVal =  sampleNonZeroEntry(entry.value);
			if (newVal != 0) {
				NonZeroEntry_Long noisyEntry = ObjectPool.allocateNonZeroEntry_Long();
				noisyEntry.setValue(newVal);
				noisyEntry.setPosition(entry.getPosition());
				noisyEntry.setRand(sampleRand.nextDouble());
				selectedArray.add(noisyEntry);
			}
			nonZeroIndices.add(entry.getPosition());
		}
		
		int numZerosSelected = getNumZerosSelected(numZeros); // this need to be > K
				
		ThresholdSampling wrs = new ThresholdSampling(alpha, (int) omega, noiseSeed, sampleSeed, zeroPosSeed, false);
		
		HashSet<Long> sampledIndices =  new HashSet<Long>((int) (1.2*numZerosSelected));
		                            
		for (int i=0; i<numZerosSelected; i++) {	
			long index = Math.abs(binomialRand.nextLong() % inputSize);
			while (sampledIndices.contains(index) || nonZeroIndices.contains(index)) {
				index = Math.abs(binomialRand.nextLong() % inputSize);
			}
			sampledIndices.add(index);
			
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
			
			NonZeroEntry_Long entry = ObjectPool.allocateNonZeroEntry_Long();
			entry.setPosition(index);
			entry.setValue(sample);
			entry.setRand(randNum);
			selectedArray.add(entry);
		}
		
		getKItems_Long(selectedArray, output);
		ObjectPool.returnNonZeroEntries_Long(selectedArray);
		selectedArray.clear();
		return output.size();	
	}
	
	
	public void getKItems_Long(ArrayList<NonZeroEntry_Long> selectedArray, ArrayList<NonZeroEntry_Long> output) {
		int numEntries = selectedArray.size();
		NonZeroEntry_Long[] tempArray = new NonZeroEntry_Long[selectedArray.size()];
		for (int i=0; i<tempArray.length; i++)
			tempArray[i] = selectedArray.get(i);
		
		QuickSelect.quickselect(tempArray, numEntries-K);
		double kPlus1Priority = tempArray[numEntries-K-1].getPriority();
//		
//		Comparator comparator = new PriorityComparator();
//		Collections.sort(selectedArray, comparator);
//		double kPlus1Priority = selectedArray.get(numEntries-K-1).getPriority();
		
		for (int i=0; i<numEntries; i++) {
//		for (int i=numEntries-K; i<numEntries; i++) {
			NonZeroEntry_Long entry = selectedArray.get(i);

			if (entry.getPriority() <= kPlus1Priority)
				continue;
			
				NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
				outEntry.position = entry.getPosition();
				
				if (Math.abs(entry.value) >= kPlus1Priority) 
					outEntry.setValue(entry.value);
				else {
					if (entry.value > 0)
						outEntry.setValue( kPlus1Priority);
					else
						outEntry.setValue( -kPlus1Priority);
				}
				output.add(outEntry);
//			}
		}
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
				upper = midPoint-unit;
		}
		sample = upper;
		
		// pick r from [0, 1];
		
		return sample;
	}

	// Pr(Y <= y), y>=0
	private double cdf(double y) {
		//double probGeqYPlus1 = 2*Math.pow(alpha, y+1)/(1+alpha) + 
		//1*alpha*(1+y*Math.pow(alpha, y+1)- omega*Math.pow(alpha, y))/((y+1)*(1-alpha*alpha));
		double probGeqYPlus1 = 2*alpha*(1- Math.pow(alpha, y+1))/((y+1)*(1-alpha*alpha));
		return 1 - probGeqYPlus1;
			
	}
	 
	
	
	static public class PriorityComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			if (o1 instanceof NonZeroEntry && o2 instanceof NonZeroEntry) {
				NonZeroEntry e1 = (NonZeroEntry) o1;
				NonZeroEntry e2 = (NonZeroEntry) o2;
				if (e1.getPriority() >= e2.getPriority())
					return 1;
				else 
					return -1;
			}
			return 0;
		}
	}
	
	public static void main(String[] args) {
		int inputSize = 1000;
		int sampleSize = 100;
		int numNonZeros = 850;
		
		double mean = 100; double stddev = 10;
		int numRepeats = 10;
		double absError = 0;
		
		for (int k=0; k<numRepeats; k++) {
			Random rand = new Random(3462 + 53*k);
			int[] inArray = new int[inputSize];
			for (int i=0; i<numNonZeros; i++) {
				inArray[i] = (int) Math.round(rand.nextGaussian()*stddev + mean);
			}
			double[] outArray = new double[inputSize];
			PrioritySampling sampling = new PrioritySampling(0.99, sampleSize, 2345, 2432, 3999);
			sampling.anonymize_basic(inArray, outArray);
		
			for (int i=0; i<inputSize; i++) {
				absError += Math.abs(inArray[i] - outArray[i]);
			}
		}
		System.out.println("Abs error=" + absError);
	}
}
