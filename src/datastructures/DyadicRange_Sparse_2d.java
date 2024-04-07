package datastructures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

import main.Run;
import main.Run.Mechanism;

import datastructures.UtilityObject.NonZeroEntry;
import datastructures.UtilityObject.NonZeroEntry_2D;

import distributions.Geometric;

public class DyadicRange_Sparse_2d implements DyadicRange_Sparse {
	int sensitivity = 1;
	HashMap<Integer,NonZeroEntry> dyadicRanges;
	HashMap<Integer,NonZeroEntry> anonymizedRanges;
	int xSize, ySize; // with zeros added to have 2^k entries
	int orgXSize, orgYSize; // orig input size
	int logx, logy;
	ArrayList<NonZeroEntry> currEntries = new ArrayList<NonZeroEntry>();
	ArrayList<NonZeroEntry> nextEntries = new ArrayList<NonZeroEntry>();
	int MULTIPLE;
	
	public DyadicRange_Sparse_2d(int xSize, int ySize) { 
		orgXSize = xSize; orgYSize = ySize;
		this.xSize = xSize; this.ySize = ySize;
		logx = (int) (Math.log(xSize)/Math.log(2)); 		
		if (xSize > (1 << logx)) {
			logx += 1; // need to add zero entries to the end of the array
			this.xSize = (1 << logx);
		}
		
		logy = (int) (Math.log(ySize)/Math.log(2)); 		
		if (ySize > (1 << logy)) {
			logy += 1; // need to add zero entries to the end of the array
			this.ySize = (1 << logy);
		}
		
		MULTIPLE = 2*ySize;
		// top node is counted as a level
		//numLevels = depth+1;
		
		dyadicRanges = new HashMap<Integer,NonZeroEntry>();
		anonymizedRanges = new HashMap<Integer,NonZeroEntry>();
	}
	
	public void update(ArrayList<NonZeroEntry> inputArr) { 
//		ArrayList<NonZeroEntry_2D> tempList = new ArrayList<NonZeroEntry_2D>();
//		for (NonZeroEntry entry: inputArr) {
//			NonZeroEntry_2D drEntry = ObjectPool.allocateNonZeroEntry_2D();
//			int x = entry.getPosition()/orgYSize + this.xSize;
//			int y = entry.getPosition() % orgYSize + this.ySize;
//			drEntry.setXY(x, y);
//			drEntry.setValue(entry.value);
//			//int index = x*MULTIPLE+y;
//			//dyadicRanges.put(index, drEntry);
//			tempList.add(drEntry);
//		}	
		
		//for (NonZeroEntry_2D entry : tempList) {
		for (NonZeroEntry entry : inputArr) {
			int x = entry.getPosition()/orgYSize + this.xSize;
			int y = entry.getPosition() % orgYSize + this.ySize;
			while (x>=1) {
				int tempY = y;
				while (tempY>=1) {
					int index = x*MULTIPLE + tempY;
					if (dyadicRanges.containsKey(index))
						dyadicRanges.get(index).value += entry.value;
					else {
						NonZeroEntry parent = ObjectPool.allocateNonZeroEntry();
						parent.setPosition(index);
						parent.setValue(entry.value);
						dyadicRanges.put(index, parent);
					}
					tempY = tempY/2;
				}
				x = x/2;
			}
		}
		
	}
	
	public void update2D(ArrayList<NonZeroEntry_2D> inputArr) { 
		ArrayList<NonZeroEntry_2D> tempList = new ArrayList<NonZeroEntry_2D>();
		for (NonZeroEntry_2D entry: inputArr) {
			NonZeroEntry_2D drEntry = ObjectPool.allocateNonZeroEntry_2D();
			int x = entry.getX() + this.xSize;
			int y = entry.getY() + this.ySize;
			drEntry.setXY(x, y);
			drEntry.setValue(entry.value);
			//int index = x*MULTIPLE+y;
			//dyadicRanges.put(index, drEntry);
			tempList.add(drEntry);
		}	
		
		for (NonZeroEntry_2D entry : tempList) {
			int x = entry.x;
			int y = entry.y;
			while (x>=1) {
				int tempY = y;
				while (tempY>=1) {
					int index = x*MULTIPLE + tempY;
					if (dyadicRanges.containsKey(index))
						dyadicRanges.get(index).value += entry.value;
					else {
						NonZeroEntry parent = ObjectPool.allocateNonZeroEntry();
						parent.setPosition(index);
						parent.setValue(entry.value);
						dyadicRanges.put(index, parent);
					}
					tempY = tempY/2;
				}
				x = x/2;
			}
		}
		
	}
	
	public int anonymize(Mechanism mechanism, int theta, int K, double epsilon, int noiseSeed, int sampleSeed, int zeroPosSeed) {
//		System.out.println(dyadicRanges.size());
		double alpha = Math.exp(-epsilon/(sensitivity*((logx+1)*(logy+1)-1)));
		
		int sampleSize = 0;
		ArrayList<NonZeroEntry> inArray = null;
		ArrayList<NonZeroEntry> outArray = null;
		if (mechanism != Mechanism.Geometric && mechanism != Mechanism.Filter2) {
			currEntries.clear();
			inArray = currEntries; //new ArrayList<NonZeroEntry>();
			for (NonZeroEntry entry: dyadicRanges.values()) {
				inArray.add(entry);
			}
			
		}
		
		nextEntries.clear();
		outArray = nextEntries; //new ArrayList<NonZeroEntry>();
		
		int inputSize = xSize*ySize;
		System.out.println("n=" + inArray.size() + ", m="+inputSize);
		
		if (mechanism == Mechanism.Filter1) {			
			Filtering_Signed ts = new Filtering_Signed(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ts.anonymize(inArray, inputSize, outArray);
		}
		else if (mechanism == Mechanism.Filter2) {	
			Filtering_Absolute ts = new Filtering_Absolute(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
			
			sampleSize = ts.anonymize(dyadicRanges, inputSize, anonymizedRanges);
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
			for (int i=0; i<xSize; i++) {
				for (int j=0; j<ySize; j++) {
					NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
					int index = i*MULTIPLE + j;
					entry.setPosition(index);
					int noise = Geometric.getSample(rand);
					if (dyadicRanges.containsKey(index)) 
						entry.value = dyadicRanges.get(index).value + noise;
					else
						entry.value = noise;
					anonymizedRanges.put(index, entry);
				}
			}
			return inputSize;
		}
		
			
		return sampleSize;
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
	
	public double rangeQuery(int startX, int endX, int startY, int endY, boolean isAnonymized) {
		ArrayList<Integer> xIndices = new ArrayList<Integer>();
		ArrayList<Integer> yIndices = new ArrayList<Integer>();
		range(startX, endX, xIndices);
		range(startY, endY, yIndices);
		double sum = 0;
		for (int x : xIndices) {
			for (int y: yIndices) {
				int index = x*MULTIPLE + y;
				if (isAnonymized) {
					if (anonymizedRanges.containsKey(index)) sum += anonymizedRanges.get(index).value;
				}
				else {
					if (dyadicRanges.containsKey(index)) sum += dyadicRanges.get(index).value;
				}
			}
		}
		return sum;
	}
	
	// both startIndex and endIndex are inclusive in the range
	private void range(int startIndex, int endIndex, ArrayList<Integer> indices) {
		HashMap<Integer, NonZeroEntry> dr = null;
//		if (!isAnonymized)
//			dr = dyadicRanges;
//		else 
//			dr = anonymizedRanges;
		
		if (startIndex > endIndex) 
			return;
		if (startIndex == endIndex) {
			indices.add(startIndex);
			return;
		}
		
		if (startIndex % 2 == 0 && endIndex % 2 ==1) {
			range(startIndex/2, endIndex/2, indices); 
		}
		else if (startIndex % 2 == 1 && endIndex % 2 ==1) {
			indices.add(startIndex);
			range((startIndex+1)/2, endIndex/2, indices);
		}
		else if (startIndex % 2 == 0 && endIndex % 2 ==0) {
			indices.add(endIndex);
			range(startIndex/2, (endIndex-1)/2, indices);  
		}
		else {
			indices.add(startIndex);
			indices.add(endIndex);
			range((startIndex+1)/2, (endIndex-1)/2, indices);
		}
		return;
	}
	

	public static void main(String[] args) {
		DyadicRange_Sparse_2d dr = new DyadicRange_Sparse_2d(8, 4);
		
		ArrayList<NonZeroEntry_2D> list = new ArrayList<NonZeroEntry_2D>();
		list.add(new NonZeroEntry_2D(0, 1, 5, 0));
		list.add(new NonZeroEntry_2D(2, 2, 2, 0));
		list.add(new NonZeroEntry_2D(4, 0, 3, 0));
		list.add(new NonZeroEntry_2D(4, 3, 4, 0));
		list.add(new NonZeroEntry_2D(7, 1, 1, 0));
		
		dr.update2D(list);
		
		for (NonZeroEntry entry : dr.dyadicRanges.values()) {
			System.out.println(entry);
		}
		
		ArrayList<Integer> indices = new ArrayList<Integer>();
		dr.range(10,14,indices);
		System.out.println(indices.toString());
		System.out.println(dr.rangeQuery(10,14,4,5,false));
		System.out.println(dr.rangeQuery(12,12,4,7,false));
	}
}
