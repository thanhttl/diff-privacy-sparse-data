package datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

import main.Run;
import main.Run.Mechanism;

import datastructures.UtilityObject.NonZeroEntry;

import distributions.Geometric;

public class DyadicRange_Sparse_1d implements DyadicRange_Sparse {
	int sensitivity = 1;
	HashMap<Integer,NonZeroEntry> dyadicRanges;
	HashMap<Integer,NonZeroEntry> anonymizedRanges;
	int inputSize; // with zeros added to have 2^k entries
	int orgInputSize; // orig input size
	int numLevels;
	ArrayList<NonZeroEntry> currEntries = new ArrayList<NonZeroEntry>();
	ArrayList<NonZeroEntry> nextEntries = new ArrayList<NonZeroEntry>();
	
	public DyadicRange_Sparse_1d(int inputSize) { 
		orgInputSize = inputSize;
		int depth = (int) (Math.log(inputSize)/Math.log(2)); 		
		if (inputSize > (1 << depth)) {
			depth += 1; // need to add zero entries to the end of the array
			this.inputSize = (1 << depth);
		}
		// top node is counted as a level
		numLevels = depth+1;
		
		dyadicRanges = new HashMap<Integer,NonZeroEntry>();
		anonymizedRanges = new HashMap<Integer,NonZeroEntry>();
		
//		
//		for (NonZeroEntry entry: inputArr) {
//			NonZeroEntry drEntry = ObjectPool.allocateNonZeroEntry();
//			int thisPosition = entry.getPosition() + this.inputSize;
//			drEntry.setPosition(thisPosition);
//			drEntry.setValue(entry.value);
//			dyadicRanges.put(thisPosition, drEntry);
//			currEntries.add(drEntry);
//		}	
//		
//		
//		for (int i=1; i<depth; i++) {
////			System.out.println(currEntries.size());
//			for (NonZeroEntry entry: currEntries) {
//				int thisPosition = entry.getPosition()/2;
//				if (dyadicRanges.containsKey(thisPosition)) {
//					dyadicRanges.get(thisPosition).value += entry.value;
//				}
//				else {
//					NonZeroEntry parent = ObjectPool.allocateNonZeroEntry();
//					parent.setPosition(thisPosition);
//					parent.setValue(entry.value);
//					dyadicRanges.put(thisPosition, parent);
//					nextEntries.add(parent);
//				}
//			}
//			ArrayList<NonZeroEntry> temp = currEntries;
//			currEntries = nextEntries;
//			nextEntries = temp;
//			nextEntries.clear();
//		}	
	}
	
	public void update(ArrayList<NonZeroEntry> inputArr) { 
		
		for (NonZeroEntry entry: inputArr) {
			NonZeroEntry drEntry = ObjectPool.allocateNonZeroEntry();
			int thisPosition = entry.getPosition() + this.inputSize;
			drEntry.setPosition(thisPosition);
			drEntry.setValue(entry.value);
			dyadicRanges.put(thisPosition, drEntry);
			currEntries.add(drEntry);
		}	
		
		
		for (int i=1; i<numLevels-1; i++) {
//			System.out.println(currEntries.size());
			for (NonZeroEntry entry: currEntries) {
				int thisPosition = entry.getPosition()/2;
				if (dyadicRanges.containsKey(thisPosition)) {
					dyadicRanges.get(thisPosition).value += entry.value;
				}
				else {
					NonZeroEntry parent = ObjectPool.allocateNonZeroEntry();
					parent.setPosition(thisPosition);
					parent.setValue(entry.value);
					dyadicRanges.put(thisPosition, parent);
					nextEntries.add(parent);
				}
			}
			ArrayList<NonZeroEntry> temp = currEntries;
			currEntries = nextEntries;
			nextEntries = temp;
			nextEntries.clear();
		}	
	}
	
	public int anonymize(Mechanism mechanism, int theta, int K, double epsilon, int noiseSeed, int sampleSeed, int zeroPosSeed) {
//		System.out.println(dyadicRanges.size());
		double alpha = Math.exp(-epsilon/(sensitivity*numLevels));
		
		int sampleSize = 0;
		ArrayList<NonZeroEntry> inArray = null;
		ArrayList<NonZeroEntry> outArray = null;
		if (mechanism != Mechanism.Geometric) {
			currEntries.clear();
			inArray = currEntries; //new ArrayList<NonZeroEntry>();
			for (NonZeroEntry entry: dyadicRanges.values()) {
				if (entry.value > 0)
					inArray.add(entry);
			}
			nextEntries.clear();
			outArray = nextEntries; //new ArrayList<NonZeroEntry>();
		}
		
		if (mechanism == Mechanism.Filter1) {			
			Filtering_Signed ts = new Filtering_Signed(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ts.anonymize(inArray, inputSize, outArray);
		}
		else if (mechanism == Mechanism.Filter2) {	
			Filtering_Absolute ts = new Filtering_Absolute(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ts.anonymize(inArray, inputSize, outArray);
		}	
		else if (mechanism == Mechanism.Threshold) {
			int absValueSum = 0;
			for (int i=0; i<inArray.size(); i++)
				absValueSum += Math.abs(inArray.get(i).getValue());
			
			absValueSum += 2*alpha/(1-alpha*alpha)*(inputSize - inArray.size());
			int D = (int) (absValueSum*(1.0)/K);
//			System.out.println("D=" + D);
			Sampling wrs;
			
			wrs = new ThresholdSampling(alpha, D, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = wrs.anonymize(inArray, inputSize, outArray);
		}				
		else if (mechanism == Mechanism.Priority) {
			PrioritySampling ps = new PrioritySampling(alpha, K, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ps.anonymize(inArray, inputSize, outArray);	
		}
		else if (mechanism == Mechanism.Filter_Priority) {
			FilterPrioritySampling ps = new FilterPrioritySampling(alpha, theta, K, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ps.anonymize(inArray, inputSize, outArray);
		}
		else if (mechanism == Mechanism.Geometric) {
			Random rand = new Random(noiseSeed);
			Geometric.setAlpha(alpha);
			for (int i=0; i<inputSize; i++) {
				NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
				entry.position = i;
				int noise = Geometric.getSample(rand);
				if (dyadicRanges.containsKey(i)) 
					entry.value = dyadicRanges.get(i).value + noise;
				else
					entry.value = noise;
				anonymizedRanges.put(i, entry);
			}
			return inputSize;
		}
		for (NonZeroEntry entry: outArray)
			anonymizedRanges.put(entry.getPosition(), entry);
		
		return sampleSize;
	}
	
	public int getOrgInputSize() {
		return orgInputSize;
	}
	
	public void reset() {
		for (NonZeroEntry entry: dyadicRanges.values()) {
			ObjectPool.returnNonZeroEntry(entry);
		}
		for (NonZeroEntry entry: anonymizedRanges.values()) {
			ObjectPool.returnNonZeroEntry(entry);
		}
		dyadicRanges.clear();
		anonymizedRanges.clear();
		
		currEntries.clear();
		nextEntries.clear();
	}
	
	public double rangeQuery(int startIndex, int endIndex, boolean isAnonymized) {
		return range(startIndex + inputSize, endIndex + inputSize, isAnonymized);
	}
	
	// both startIndex and endIndex are inclusive in the range
	private double range(int startIndex, int endIndex, boolean isAnonymized) {
		HashMap<Integer, NonZeroEntry> dr;
		if (!isAnonymized)
			dr = dyadicRanges;
		else 
			dr = anonymizedRanges;
		
		if (startIndex > endIndex)
			return 0;
		if (startIndex == endIndex) {
			if (dr.containsKey(startIndex))
				return dr.get(startIndex).getValue();
			else 
				return 0;
		}
		
		if (startIndex % 2 == 0 && endIndex % 2 ==1) {
			return range(startIndex/2, endIndex/2, isAnonymized); 
		}
		else if (startIndex % 2 == 1 && endIndex % 2 ==1) {
			double val = 0;
			if (dr.containsKey(startIndex))
				val = dr.get(startIndex).value;
			return val + range((startIndex+1)/2, endIndex/2, isAnonymized);
		}
		else if (startIndex % 2 == 0 && endIndex % 2 ==0) {
			double val = 0;
			if (dr.containsKey(endIndex))
				val = dr.get(endIndex).value;
			return range(startIndex/2, (endIndex-1)/2, isAnonymized) + val;  
		}
		else {
			double val1 = 0;
			if (dr.containsKey(startIndex))
				val1 = dr.get(startIndex).value;
			double val2 = 0;
			if (dr.containsKey(endIndex))
				val2 = dr.get(endIndex).value;
			return range((startIndex+1)/2, (endIndex-1)/2, isAnonymized) + val1 + val2;
		}
	}
	
	public static void main(String[] args) {
		double[] temp = {3, 8, 6, 0, 4 , 0, 10};
		ArrayList<NonZeroEntry> list = new ArrayList<NonZeroEntry>();
		for (int i=0; i<temp.length; i++) {
			if (temp[i] != 0) {
				NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
				entry.position = i;
				entry.value = temp[i];
				list.add(entry);
			}
		}
		DyadicRange_Sparse_1d dr = new DyadicRange_Sparse_1d(temp.length);
		dr.update(list);
		dr.anonymize(Run.Mechanism.Filter1, 10, 20, 0.1, 243632, 98433, 56296);
		System.out.println(dr.rangeQuery(1, 7, false));
		System.out.println(dr.rangeQuery(1, 7, true));
	}
}
