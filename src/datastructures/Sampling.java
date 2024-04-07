package datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import datastructures.UtilityObject.*;
import distributions.Geometric;

public abstract class Sampling {
	double alpha;
	Random noiseRand, sampleRand; 
	Random binomialRand; 
	//boolean useLinkedList = false; // for zero entries;
	public static boolean rejectionSampling = true;
	
	public abstract double sampleZeroEntry(double prob);
	public abstract int getNumZerosSelected(int numZeroes);
	public abstract double sampleNonZeroEntry(double input);
	//public abstract boolean isNonZeroSampled(double input);
	
	public abstract int getNumZerosSelected(long numZeroes);
	
	public int anonymize_basic(double[] input, double[] output) {
		int totalNumSampled = 0;
		for (int i=0; i<input.length; i++) {
			output[i] = sampleNonZeroEntry(input[i]);
			if (output[i] != 0)
				totalNumSampled++;
		}
		return totalNumSampled;
	}
	
	public int anonymize(ArrayList<NonZeroEntry> nonZeros, int inputSize, ArrayList<NonZeroEntry> output) {
		return anonymize(nonZeros, null, inputSize, output);
	}
	
	public int anonymize(ArrayList<NonZeroEntry> nonZeros, ArrayList<ZeroIndexInterval> zeroIndices, int inputSize, ArrayList<NonZeroEntry> output) {
		int numNonZeros = nonZeros.size();
		
		//int totalNumSampled = 0;
		
		HashSet<Integer> nonZeroIndices = new HashSet<Integer>();
		
		for (int i=0; i< numNonZeros; i++) {
			NonZeroEntry entry = nonZeros.get(i);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (sample != 0) {
				NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				output.add(outEntry);
		//		totalNumSampled ++;
			}
			nonZeroIndices.add(entry.getPosition());
		}
		
		//System.out.println("Num zeros selected: " + numZerosSelected);	
		// ArrayList<Integer> upgradedList = new ArrayList<Integer>();
		// shuffle the zero indices
		HashSet<Integer> sampledIndices = sampleZeroIndices(nonZeroIndices, zeroIndices, inputSize, binomialRand);
		
		//totalNumSampled += sampledIndices.size();
		
		Iterator<Integer> sampledIndexIterator = sampledIndices.iterator();
		while (sampledIndexIterator.hasNext()) {
			int index = sampledIndexIterator.next();
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double y = sampleRand.nextDouble();
			double sample = sampleZeroEntry(y);
			NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			output.add(outEntry);			
		}		
		return output.size();
	}
	
	public HashSet<Integer> sampleZeroIndices(HashSet<Integer> nonZeroIndices, ArrayList<ZeroIndexInterval> zeroIndices, int inputSize, Random rand) {
		HashSet<Integer> sampledZeroIndices = new HashSet<Integer>();
		int numZeros = inputSize - nonZeroIndices.size();
		
		int numZerosSelected = getNumZerosSelected(numZeros);
		
		for (int i=0; i<numZerosSelected; i++) {
			// randomly sample from the list of zero indices by rejection sampling
			int index;
			if (zeroIndices == null) {
				index = rand.nextInt(inputSize);
				while (nonZeroIndices.contains(index) || sampledZeroIndices.contains(index)) {
					index = rand.nextInt(inputSize);
				}
			}
			else {
				do {
					double prob = rand.nextDouble();
					// find the interval that this prob is in
					int lowerIndex = 0, upperIndex = zeroIndices.size()-1;
					int middleIndex = 0; ///zeroIndices.size()/2;
					boolean findInterval = false;
					
					while (!findInterval) {
						if (prob <= zeroIndices.get(lowerIndex).cumProb) {
							findInterval = true; middleIndex = lowerIndex;
						}
						else if (prob > zeroIndices.get(upperIndex-1).cumProb) {
							findInterval = true; middleIndex =  upperIndex;
						}
						else {
							middleIndex = (lowerIndex + upperIndex)/2; 
							if (prob > zeroIndices.get(middleIndex-1).cumProb && prob <= zeroIndices.get(middleIndex).cumProb)
								findInterval = true;
							else if (zeroIndices.get(middleIndex-1).cumProb >= prob) {
								upperIndex = middleIndex-1;
							}
							else {
								lowerIndex = middleIndex+1; 
							}
						}
					}
					double prevProb;
					if (middleIndex == 0)
						prevProb = 0;
					else 
						prevProb = zeroIndices.get(middleIndex-1).cumProb;
					
					ZeroIndexInterval interval = zeroIndices.get(middleIndex);
					
					double endProb = interval.cumProb;
					double percentage = (prob-prevProb)/(endProb-prevProb);
					index = interval.startIndex + (int) Math.round(percentage*(interval.endIndex-interval.startIndex));
				} while (sampledZeroIndices.contains(index)); 
			} 
//				boolean checking = true; int vary =1;
//				while (checking) {
//					//if (sampleRand.nextBoolean()) {
//						if (index-vary>=0 && !nonZeroIndices.contains(index-vary) && !sampledZeroIndices.contains(index-vary)) {
//							index = index-vary; break;
//						}
//						if (index+vary<inputSize && !nonZeroIndices.contains(index+vary) && !sampledZeroIndices.contains(index+vary)) {
//							index = index+vary; break;
//						}
//					}
//					else {
//						if (index+vary<inputSize && !nonZeroIndices.contains(index+vary) && !sampledZeroIndices.contains(index+vary)) {
//							index = index+vary; break;
//						}
//						if (index-vary>=0 && !nonZeroIndices.contains(index-vary) && !sampledZeroIndices.contains(index-vary)) {
//							index = index-vary; break;
//						}
//					}
//					vary += 1;
//				}
//			}
			sampledZeroIndices.add(index);
		}
		return sampledZeroIndices;	
	}
	
	// for dyadic ranges, when the map stores the non-zero indices already
	public int anonymize(HashMap<Integer,NonZeroEntry> nonZeros, int inputSize, HashMap<Integer,NonZeroEntry> output) {
		//for (NonZeroEntry entry : nonZeros.values()) {
		Iterator<Integer> entryIter = nonZeros.keySet().iterator();
		while (entryIter.hasNext()) {
			Integer index = entryIter.next();
			NonZeroEntry entry = nonZeros.get(index);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (sample != 0) {
				NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				output.put(index, outEntry);
			}
		}
		
		//System.out.println("Num zeros selected: " + numZerosSelected);	
		// ArrayList<Integer> upgradedList = new ArrayList<Integer>();
		// shuffle the zero indices
		
		HashSet<Integer> sampledIndices = sampleZeroIndices(nonZeros, inputSize, binomialRand);
		
		//totalNumSampled += sampledIndices.size();
		
		Iterator<Integer> sampledIndexIterator = sampledIndices.iterator();
		while (sampledIndexIterator.hasNext()) {
			Integer index = sampledIndexIterator.next();
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double y = sampleRand.nextDouble();
			double sample = sampleZeroEntry(y);
			NonZeroEntry outEntry = ObjectPool.allocateNonZeroEntry();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			output.put(index, outEntry);			
		}		
		return output.size();
	}
	
	public HashSet<Integer> sampleZeroIndices(HashMap<Integer,NonZeroEntry> nonZeros, int inputSize, Random rand) {	
		int numZeros = inputSize - nonZeros.size();
		
		int numZerosSelected = getNumZerosSelected(numZeros);
//System.out.println("Num zeros selected: " + numZerosSelected);		
		HashSet<Integer> sampledZeroIndices = new HashSet<Integer>(numZerosSelected);
		for (int i=0; i<numZerosSelected; i++) {
			// randomly sample from the list of zero indices by rejection sampling
			int index;
			
			index = rand.nextInt(inputSize);
			while (nonZeros.containsKey(index) || sampledZeroIndices.contains(index)) {
					index = rand.nextInt(inputSize);
			}
			sampledZeroIndices.add(index);
		}
		return sampledZeroIndices;
	}
	
//	public static void main(String[] args) {
//		double[] input = {60, 20, 0, 80, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//		double alpha = Math.exp(-0.1/2);
//		int D = 20;
//		ThresholdSampling wrs = new ThresholdSampling(alpha, D, 12345, 64324);
//		double[] output = new double[input.length];
//		int sampleSize = wrs.anonymize(input, output);
//		System.out.println("Weighted random sampling");
//		System.out.println("sample size=" + sampleSize);
//		for (double val: output) {
//			System.out.print(val + ", ");
//		}
//		System.out.println();
//		
//		double theta = 5;
//		Filtering_Signed ts = new Filtering_Signed(alpha, theta, 12345, 59709);
//		output = new double[input.length];
//		sampleSize = ts.anonymize(input, output);
//		System.out.println("Threshold sampling");
//		System.out.println("sample size=" + sampleSize);
//		for (double val: output) {
//			System.out.print(val + ", ");
//		}
//		System.out.println();
//	}
	// ************ FOR LONG INPUT SIZE
	public int anonymize(ArrayList<NonZeroEntry_Long> nonZeros, long inputSize, ArrayList<NonZeroEntry_Long> output) {
		int numNonZeros = nonZeros.size();
		
		//int totalNumSampled = 0;	
		HashSet<Long> nonZeroIndices = new HashSet<Long>();
		
		for (int i=0; i< numNonZeros; i++) {
			NonZeroEntry_Long entry = nonZeros.get(i);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (sample != 0) {
				NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				output.add(outEntry);
		//		totalNumSampled ++;
			}
			nonZeroIndices.add(entry.getPosition());
		}

//System.out.println("get zero indices");
		HashSet<Long> sampledIndices = sampleZeroIndices(nonZeroIndices, inputSize, binomialRand);
System.out.println("# zero indices = " + sampledIndices.size());		
		//totalNumSampled += sampledIndices.size();
		
		Iterator<Long> sampledIndexIterator = sampledIndices.iterator();
		
		while (sampledIndexIterator.hasNext()) {
			long index = sampledIndexIterator.next();
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double y = sampleRand.nextDouble();
			double sample = sampleZeroEntry(y);
			NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			output.add(outEntry);			
		}		
		return output.size();
	}
	
	public HashSet<Long> sampleZeroIndices(HashSet<Long> nonZeroIndices, long inputSize, Random rand) {
		long numZeros = inputSize - nonZeroIndices.size();
		
		int numZerosSelected = getNumZerosSelected(numZeros);
//		System.out.println("num zeros selected:" + numZerosSelected);		
		HashSet<Long> sampledZeroIndices = new HashSet<Long>((int) (1.05*numZerosSelected));

		for (int i=0; i<numZerosSelected; i++) {
			// randomly sample from the list of zero indices by rejection sampling
			long index = Math.abs(rand.nextLong() % inputSize);
			while (nonZeroIndices.contains(index) || sampledZeroIndices.contains(index)) {
				index = Math.abs(rand.nextLong() % inputSize);
			}
			sampledZeroIndices.add(index);
		}
		return sampledZeroIndices;	
	}
	
	// for dyadic ranges, when the map stores the non-zero indices already
	public int anonymize(HashMap<Long,NonZeroEntry_Long> nonZeros, long inputSize, HashMap<Long,NonZeroEntry_Long> output) {
		//for (NonZeroEntry_Long entry : nonZeros.values()) {
		Iterator<Long> entryIter = nonZeros.keySet().iterator();
		while (entryIter.hasNext()) {
			Long index = entryIter.next();
			NonZeroEntry_Long entry = nonZeros.get(index);
			double sample =  sampleNonZeroEntry(entry.getValue());
			if (sample != 0) {
				NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
				outEntry.setPosition(entry.getPosition());
				outEntry.setValue(sample);
				output.put(index, outEntry);
			}
		}
		
		//System.out.println("Num zeros selected: " + numZerosSelected);	
		// ArrayList<Integer> upgradedList = new ArrayList<Integer>();
		// shuffle the zero indices
		System.out.println("sample zero entries");
		HashSet<Long> sampledIndices = sampleZeroIndices(nonZeros, inputSize, binomialRand);
		System.out.println("num zero entries: " + sampledIndices.size());
		//totalNumSampled += sampledIndices.size();
		
		Iterator<Long> sampledIndexIterator = sampledIndices.iterator();
		while (sampledIndexIterator.hasNext()) {
			Long index = sampledIndexIterator.next();
			
			// sampling directly from the distribution P(M''(i,j)=v | (i,j) selected)
			// using the CDF
			double y = sampleRand.nextDouble();
			double sample = sampleZeroEntry(y);
			NonZeroEntry_Long outEntry = ObjectPool.allocateNonZeroEntry_Long();
			outEntry.setPosition(index);
			outEntry.setValue(sample);
			output.put(index, outEntry);			
		}		
		return output.size();
	}
	
	public HashSet<Long> sampleZeroIndices(HashMap<Long,NonZeroEntry_Long> nonZeros, long inputSize, Random rand) {	
		long numZeros = inputSize - nonZeros.size();
		
		int numZerosSelected = getNumZerosSelected(numZeros);
System.out.println("Num zeros selected: " + numZerosSelected);		
		HashSet<Long> sampledZeroIndices = new HashSet<Long>(numZerosSelected);
		for (int i=0; i<numZerosSelected; i++) {
			// randomly sample from the list of zero indices by rejection sampling
			long index = Math.abs(rand.nextLong() % inputSize);
			while (nonZeros.containsKey(index) || sampledZeroIndices.contains(index)) {
					index = Math.abs(rand.nextLong() % inputSize);
			}
			sampledZeroIndices.add(index);
		}
		return sampledZeroIndices;
	}
	
}
