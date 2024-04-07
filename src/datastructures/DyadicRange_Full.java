package datastructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import main.Run;


import datastructures.UtilityObject.NonZeroEntry;
import distributions.Geometric;

public class DyadicRange_Full {
	// input: 1-d array with 2^m entries 
	// dyadic intervals are [j*2^a, (j+1)*2^a - 1]
	// a \in [0, m]; j \in [0, 2^{m-a} -1 ] (this includes the original data)

	BinaryTree intervalTree;
	double[] intervalArray;
	int numLevels; // not including the tree root
	Random rand;
	public int sampleSize; // for the case of weighted random sampling
    // temp variables
	ArrayList<BinaryTree> treeArray = new ArrayList<BinaryTree>(); 
	ArrayList<NonZeroEntry> nonZeroArray = new ArrayList<NonZeroEntry>(); // only for measuring time experiments
	int[] tempZeroIndices;
	double[] tempArray;
	int numZeros;
	double nonZeroSum;
	public static boolean consistencyCheck = true;
	
	public DyadicRange_Full() {
	}
	
	public DyadicRange_Full(double[] inputArr, int seed) {
		this(inputArr, 1, 0, seed);
	}
	
	public DyadicRange_Full(double[] inputArr, double epsilon, double sensitivity, int seed) {
		buildGeometric(inputArr, epsilon, sensitivity, seed);
		
	}
	
	public void buildGeometric(double[] inputArr, double epsilon, double sensitivity, int seed) {
		int length = inputArr.length;
		int depth = (int) (Math.log(length)/Math.log(2));
		
		if (length > (1 << depth)) {
			depth += 1; // need to add zero entries to the end of the array
		}
		numLevels = depth;
		
		rand = new Random(seed);
		
		int newLength = 1 << depth;
		if (intervalArray == null || intervalArray.length != 2*newLength-1)
			intervalArray = new double[newLength*2 - 1]; 
		formDyadicRanges(inputArr, depth, intervalArray);
				
		//double scale = numLevels*sensitivity/epsilon;
		//double alpha = Math.exp(-epsilon/(sensitivity*depth));
		for (int i=0; i< intervalArray.length; i++) {
			intervalArray[i] = intervalArray[i] + Geometric.getSample(rand);
		}
		
		if (consistencyCheck)
			checkConsistency(intervalArray);
		
		formTree();
	}
	
	// param is sample size
	public DyadicRange_Full(double[] inputArr, double epsilon, double sensitivity, Run.Mechanism mechanism, int param, int theta, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		this(inputArr, epsilon, sensitivity, mechanism, param, theta, false, noiseSeed, sampleSeed, zeroPosSeed);
	}
	

	public DyadicRange_Full(double[] inputArr, double epsilon, double sensitivity, Run.Mechanism mechanism, int param, int theta, boolean measuringTime, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		sample(inputArr, epsilon, sensitivity, mechanism, param, theta, measuringTime, noiseSeed, sampleSeed, zeroPosSeed);
	}
	
	public void sample(double[] inputArr, double epsilon, double sensitivity, Run.Mechanism mechanism, int param, int theta, boolean measuringTime, int noiseSeed, int sampleSeed, int zeroPosSeed) {
		int length = inputArr.length;
		int depth = (int) (Math.log(length)/Math.log(2));
		
		if (length > (1 << depth)) {
			depth += 1; // need to add zero entries to the end of the array
		}
		numLevels = depth;
		
		int newLength = 1 << depth;
		if (intervalArray == null || intervalArray.length != 2*newLength-1)
			intervalArray = new double[newLength*2 - 1];
		
		formDyadicRanges(inputArr, depth, intervalArray);
		
		double alpha = Math.exp(-epsilon/(sensitivity*depth));	
		double nonZeroMean = (double) nonZeroSum/(intervalArray.length - numZeros);
//		double zeroNonZeroRatio = 1.0*numZeros/(intervalArray.length - numZeros);
		//System.out.println("NonZeroMean=" + nonZeroMean + ", ZeroNonZeroRatio=" + zeroNonZeroRatio);
		
		Sampling sampling = null;
		if (mechanism == Run.Mechanism.Threshold) {
			int D = (int) ((nonZeroMean*(intervalArray.length-numZeros) + 2*alpha/(1-alpha*alpha)*numZeros)/param);
			//System.out.println("D=" + D);
			sampling = new ThresholdSampling(alpha, D, noiseSeed, sampleSeed, zeroPosSeed);
		}
		else if (mechanism == Run.Mechanism.Priority){
			int K = (int) param;
			sampling = new PrioritySampling(alpha, K, noiseSeed, sampleSeed, zeroPosSeed);
		}
		else if (mechanism == Run.Mechanism.Filter_Priority){
			int K = (int) param;
			//System.out.println("K=" + K + ", theta=" + theta);
			if (! consistencyCheck)
				sampling = new FilterPrioritySampling(alpha, theta, K, noiseSeed, sampleSeed, zeroPosSeed);
			else
				sampling = new Filtering_Absolute(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
		}
		else if (mechanism == Run.Mechanism.Filter1) {
			int threshold = (int) Run.computeFilteringTheta(param, mechanism, alpha, numZeros, intervalArray.length-numZeros, nonZeroMean);
			//System.out.println("theta=" + threshold);
			sampling = new Filtering_Signed(alpha, threshold, noiseSeed, sampleSeed, zeroPosSeed);
		}
		else if (mechanism == Run.Mechanism.Filter2) {
			//int threshold = (int) RunSampling.computeFilteringTheta(param, mechanism, alpha, numZeros, intervalArray.length-numZeros, nonZeroMean);
			int threshold = theta;
			sampling = new Filtering_Absolute(alpha, threshold, noiseSeed, sampleSeed, zeroPosSeed);
		}
		
		if (tempArray == null || tempArray.length != intervalArray.length)
			tempArray = new double[intervalArray.length];
				
		if (!measuringTime)
			sampleSize = sampling.anonymize_basic(intervalArray, tempArray);
		else {
			ArrayList<NonZeroEntry> outputArray = new ArrayList<NonZeroEntry>();
			sampleSize = sampling.anonymize(nonZeroArray, intervalArray.length, outputArray); //tempZeroIndices, numZeros, tempArray);
			for (int i=0; i<tempArray.length; i++)
				tempArray[i] = 0;
			for (NonZeroEntry entry: outputArray) {
				tempArray[entry.getPosition()] = entry.value;
			}
		}
		
		if (mechanism == Run.Mechanism.Filter_Priority && consistencyCheck) {
//			System.out.println("Sample size before:" + sampleSize);
			sampleSize = checkConsistency(tempArray);
//			System.out.println("Sample size after:" + sampleSize);
			if (sampleSize > param) {
				PrioritySampling sampling2 = new PrioritySampling(0, (int) param, noiseSeed + 234, sampleSeed + 975, zeroPosSeed);
				for (int i=0; i< intervalArray.length; i++)
					intervalArray[i] = 0;
				sampleSize = sampling2.anonymize_basic(tempArray, intervalArray);
			}
			else {
				double[] temp = intervalArray;
				intervalArray = tempArray;
				tempArray = temp;
			}
		}
		else {
			double[] temp = intervalArray;
			intervalArray = tempArray;
			tempArray = temp;
		}
		
		formTree();	
		if (mechanism == Run.Mechanism.Filter2 && consistencyCheck)
			sampleSize = checkConsistency(intervalTree);

	}
	
	private int checkConsistency(double[] array) {
		//int numLevels = (int) (Math.log(array.length+1)/Math.log(2));
		//int currIndex = array.length-1;
		int numRetained = 0;
		for (int i=0; i<(array.length-1)/2; i++) {
		//for (; currIndex> array.length/2; currIndex--) {
			int currIndex = array.length-1-i;
//			System.out.println(currIndex + ", " + (array.length-2-2*i) + ", " + (array.length-3-2*i));
			if (array[currIndex] == 0) {
				array[array.length-2-2*i] = 0;
				array[array.length-3-2*i] = 0;
			}
			else {
				numRetained ++;
			}
		}
		for (int i=(array.length-1)/2; i<array.length; i++) {
			int currIndex = array.length-1-i;
			if (array[currIndex] != 0)
				numRetained ++;
		}
		return numRetained;
	}
	
	private int checkConsistency(BinaryTree tree) {
		if (tree == null)
			return 0;
		if (tree.value == 0) {
			reset(tree.left);
			reset(tree.right);
			return 0;
		}
		return 1 + checkConsistency(tree.left) + checkConsistency(tree.right);
	}
	
	private void reset(BinaryTree tree) {
		if (tree == null)
			return;
		tree.value = 0;
		reset(tree.left);
		reset(tree.right);
		return;
	}
	
	private void formTree() {
		int newLength = 1 << numLevels;
		
		for (int i=0; i<intervalArray.length; i++) {
			BinaryTree tree = ObjectPool.allocalteTree(); //new BinaryTree(intervalArray[i], null, null);
			tree.value = intervalArray[i];
			treeArray.add(tree);
		}
		int currNumItems = newLength;
		//ArrayList<BinaryTree> newTreeArray = new ArrayList<BinaryTree>();
		
		int currIndex = 0;
		while (currNumItems > 0) {
			for (int i=0; i < currNumItems/2; i++) {
				BinaryTree left = treeArray.get(currIndex+2*i);
				BinaryTree right = treeArray.get(currIndex+2*i+1);
				BinaryTree parent = treeArray.get(currIndex+currNumItems+i);
				parent.setLeft(left);
				parent.setRight(right);
				
			}
			currIndex += currNumItems;
			currNumItems /= 2;
		}
		
		intervalTree = treeArray.get(treeArray.size()-1);
		// assign unique ids for each tree node for hashing, level-order traversal
//		int id = 0;
//		intervalTree.id = (id++);
//		intervalTree.left.id = (id++);
//		intervalTree.right.id = (id++);
//		
//		LinkedList<BinaryTree> queue = new LinkedList<BinaryTree>();
//		queue.add(intervalTree.left);
//		queue.add(intervalTree.right);
//		while (!queue.isEmpty()) {
//			BinaryTree first = queue.remove(0);
//			
//			if (first.left != null) {
//				queue.add(first.left);
//				first.left.id = (id++);
//			}
//			if (first.right != null) {
//				queue.add(first.right);
//				first.right.id = (id++);
//			}
//		}	
	}
	
	public void reset() {
		numLevels = 0;
		sampleSize = 0;
		if (intervalArray != null) {
			for (int i=0; i<intervalArray.length; i++) {
				intervalArray[i] = 0;
			}
		}
		if (tempArray != null) {
			for (int i=0; i<tempArray.length; i++) {
				tempArray[i] = 0;
			}
		}
		ObjectPool.returnTree(intervalTree);
		treeArray.clear();
		for (NonZeroEntry entry: nonZeroArray)
			ObjectPool.returnNonZeroEntry(entry);
		nonZeroArray.clear();
		
	}
	
	// sum of entries in range [startIndex, endIndex], indices start from 0.
	public double rangeQuery(int startIndex, int endIndex) {
		if (startIndex == endIndex) {
			return intervalArray[startIndex];
		}
		// startBits = bit-string(startIndex), endBits = bit-string(endIndex)
		// for start Index, if bit = 0, sum (left, and right) , right is immediate child
		// find the least significant 1 of startBits
		String startBitString = Integer.toBinaryString(startIndex);
		String endBitString = Integer.toBinaryString(endIndex);
		
		for (int i=0; i<=1; i++) {
			String st;
			if (i==0) st = startBitString;
			else st = endBitString;
			
			int numAdded = numLevels - st.length();
			StringBuffer leadingZeroes = new StringBuffer(numAdded);
			for (int j=0; j<numAdded; j++) {
				leadingZeroes.append('0');
			}
			leadingZeroes.append(st);
			if (i==0)
				startBitString = leadingZeroes.toString();
			else
				endBitString = leadingZeroes.toString();
		}
		
		BinaryTree subTree = intervalTree;
		if (startBitString.length() == endBitString.length()) {
			int index = 0;
			while (startBitString.charAt(index) == endBitString.charAt(index)) {
//				System.out.println(startBitString.charAt(index));
				if (startBitString.charAt(index) == '0')
					subTree = subTree.left;
				else
					subTree = subTree.right;
				index ++;
			}
			startBitString = startBitString.substring(index+1);
			endBitString = endBitString.substring(index+1);
		}
		
		
		double leftVal = subTree.left.getSum(startBitString, -1);
		double rightVal = subTree.right.getSum(endBitString, 1);
		return leftVal+rightVal;
		
	}
	
//	public void addLaplacian(double epsilon, double sensitivity) {
//		double scale = numLevels*sensitivity/epsilon;
//		intervalTree.left.addLaplacian(scale, rand);
//		intervalTree.right.addLaplacian(scale, rand);
//	}
	
	public double getCount() {
		return intervalTree.value;
	}
	
	public int getNumLevels() {
		return numLevels;
	}
	
	
	public void show() {
		for (int i=0; i<intervalArray.length-1; i++) {
			System.out.print(intervalArray[i] + ", ");
		}
		System.out.println(intervalArray[intervalArray.length-1]);
		System.out.println("TREE:");
		System.out.println(intervalTree.toString());
	}
	
	private void formDyadicRanges(double[] inputArray, int depth, double[] intervalArray) {
		
		int currInterval = 0;
		numZeros = 0;
		nonZeroSum = 0;
		
		if (tempZeroIndices == null || tempZeroIndices.length < intervalArray.length);
			tempZeroIndices = new int[intervalArray.length];
		
		for (int i=currInterval; i<currInterval+inputArray.length; i++)
			intervalArray[i] = inputArray[i-currInterval];
		
		for (int a=depth-1; a>=0; a--) {
			int numIntervals = (1 << a) ;
			int currInputIndex = currInterval;
			currInterval = currInterval + 2*numIntervals;
			
			for (int j=0; j < numIntervals; j++) {
				intervalArray[currInterval+j] = intervalArray[currInputIndex] + intervalArray[currInputIndex+1];
				currInputIndex = currInputIndex+2;
				
				if (intervalArray[currInterval+j] == 0) {
					tempZeroIndices[numZeros++] = currInterval+j;
					numZeros++;
				}
				else {
					NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
					entry.setPosition(currInterval+j);
					entry.setValue((int)intervalArray[currInterval+j]);
					nonZeroArray.add(entry);
					nonZeroSum += entry.value;
				}
			}
		}		
		
	}

	public static void main(String[] args) {
//		double[] inputData = {2, 3, 8, 1, 7, 6};
//		DyadicRange di = new DyadicRange(inputData, 123);
//		di.show();
//		double rangeSum = di.rangeQuery(1, 4);
//		System.out.println("sum:" + rangeSum);
		double[] temp = {3, 5, 6, 0, 4 , 0, 3, 6, 0, 4, 6, 3, 3, 0, 4};
		DyadicRange_Full di = new DyadicRange_Full();
		di.checkConsistency(temp);
		for (int i=0; i< temp.length; i++) {
			System.out.print(temp[i] + ", ");
		}
			
	}
	
	public static class BinaryTree {
		double value;
		BinaryTree left;
		BinaryTree right;
		int id = -1;
		
		public BinaryTree(double value, BinaryTree left, BinaryTree right) {
			this.value = value;
			this.left = left;
			this.right = right;
		}
		
		public void setLeft(BinaryTree left) {
			this.left = left;
		}
		
		public void setRight(BinaryTree right) {
			this.right = right;
		}
		
		public BinaryTree getLeft() {
			return left;
		}
		
		public BinaryTree getRight() {
			return right;
		}
		
		public String toString() {
			if (left == null && right == null)
				return "ID:" + id + ", Val:" + value;
			
			return "(" + left.toString() + ", ID:" + id + ", Val:" + value + ", " + right.toString() + ")";
		}

		
		// side = -1: left, 1: right
		// the size of bitString = # of levels of tree
		public double getSum(String bitString, int side) {
			if (bitString.length() == 0)
				return value;
			if (side < 0) {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '1') 
					return right.getSum(subString, side);
				else if (subString.indexOf('1') == -1)
					return value;
				else 
					return left.getSum(bitString.substring(1), side) + right.value; 
			}
			else {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '0') 
					return left.getSum(subString, side);
				else if (subString.indexOf('0') == -1)
					return value;
				else 
					return right.getSum(bitString.substring(1), side) + left.value; 
			}
			
		}
		
		// side = -1: left, 1: right
		// the size of bitString = # of levels of tree
		public double getSum_Hashed(String bitString, int side, int level, HashMap<Integer, AMS> hashedVals) {
			boolean useMedian = true;
			double thisVal;
			if (!hashedVals.containsKey(level))
				thisVal = value;
			else 
				thisVal = hashedVals.get(level).count(id, useMedian);
			
			if (bitString.length() == 0) 
				return thisVal;
			
			double leftVal, rightVal;
			if (!hashedVals.containsKey(level+1)) {
				leftVal = left.value;
				rightVal = right.value;
			}
			else {
				leftVal = hashedVals.get(level+1).count(left.id, useMedian);
				rightVal = hashedVals.get(level+1).count(right.id, useMedian);
			}
				
			if (side < 0) {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '1') 
					return right.getSum_Hashed(subString, side, level+1, hashedVals);
				else if (subString.indexOf('1') == -1)
					return thisVal;
				else 
					return left.getSum_Hashed(subString, side, level+1, hashedVals) + rightVal; 
			}
			else {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '0') 
					return left.getSum_Hashed(subString, side, level+1, hashedVals);
				else if (subString.indexOf('0') == -1)
					return value;
				else 
					return right.getSum_Hashed(subString, side, level+1, hashedVals) + leftVal; 
			}
			
		}
		
		// side = -1: left, 1: right
		// the size of bitString = # of levels of tree
		public double getSum_Hashed(String bitString, int side, int level, int startLevel, AMS hashedVals) {
			boolean useMedian = true;
			double thisVal;
			if (level < startLevel)
				thisVal = value;
			else 
				thisVal = hashedVals.count(id, useMedian);
			
			if (bitString.length() == 0) 
				return thisVal;
			
			double leftVal, rightVal;
			if (level+1 < startLevel) {
				leftVal = left.value;
				rightVal = right.value;
			}
			else {
				leftVal = hashedVals.count(left.id, useMedian);
				rightVal = hashedVals.count(right.id, useMedian);
			}
				
			if (side < 0) {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '1') 
					return right.getSum_Hashed(subString, side, level+1, startLevel, hashedVals);
				else if (subString.indexOf('1') == -1)
					return thisVal;
				else 
					return left.getSum_Hashed(subString, side, level+1, startLevel, hashedVals) + rightVal; 
			}
			else {
				String subString = bitString.substring(1);
				if (bitString.charAt(0) == '0') 
					return left.getSum_Hashed(subString, side, level+1, startLevel, hashedVals);
				else if (subString.indexOf('0') == -1)
					return value;
				else 
					return right.getSum_Hashed(subString, side, level+1, startLevel, hashedVals) + leftVal; 
			}
			
		}
		
//		public void addLaplacian(double scale, Random rand) {
//			value += Math.round(Laplace.getSample(scale, rand));
//			if (left != null && right != null) {
//				left.addLaplacian(scale, rand);
//				right.addLaplacian(scale, rand);
//			}
//		}
	}
}
