package datastructures;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import datastructures.DyadicRange_Full.BinaryTree;


public class HashedDyadicRange {
	// dont hash for all levels <= this
	final static int FIXED_UPTO_LEVEL = 3;
	
	// hashing parameters
	boolean isOneHashTable;
	// array of <level-id, hashed array>, the AMS has depth=1 for now
	HashMap<Integer, AMS> hashedLevels;
	AMS hashTable;
	
	int maxLevel;
	Random rand;
	int[][] test;
	// TODO: just keep a reference to the original DyadicInterval here 
	// to access levels close to the root.
	DyadicRange_Full dyadicIntervals;
	
	public HashedDyadicRange(DyadicRange_Full di, int seed, boolean isOneHashTable) {
		this.isOneHashTable = isOneHashTable;
		rand = new Random(seed);
		int numLevels = di.numLevels;
		
		if (!isOneHashTable) {
			hashedLevels = new HashMap<Integer, AMS>();
			test = new int[4][numLevels];
			for (int i=0; i<4; i++) {
				for (int j=0; j<numLevels; j++) {
					test[i][j] = rand.nextInt();
					if (test[i][j] < 0)
						test[i][j] = -test[i][j];
				}
			}
		}
		else {	
			int numBuckets = 0;
			for (int i=FIXED_UPTO_LEVEL+1; i<=20; i++)
				numBuckets += (1 << i)/2;
			
			System.out.println(numBuckets);

			hashTable = new AMS(numBuckets, 1, seed);
			test = new int[4][1];
			for (int i=0; i<4; i++) {
				test[i][0] = rand.nextInt();
				if (test[i][0] < 0)
					test[i][0] = -test[i][0];
			}
		}
		dyadicIntervals = di;
		maxLevel = dyadicIntervals.numLevels;
		
		
		// traverse the tree, di, and hash based on the node id: first find the
		// level, then hash the id and insert into that level.
		BinaryTree tree = di.intervalTree;
		
		LinkedList<BinaryTree> queue = new LinkedList<BinaryTree>();
		queue.add(tree.left);
		queue.add(tree.right);
		
		
		if (isOneHashTable) {
			int numKeepUnchanged = (1 << (FIXED_UPTO_LEVEL+1)) - 1;
			while (!queue.isEmpty()) {
				BinaryTree first = queue.remove(0);
				
				if (first.left != null) {
					queue.add(first.left);
				}
				if (first.right != null) {
					queue.add(first.right);
				}
				int id = first.id;
				if (id < numKeepUnchanged) {
					continue;
				}
				else {
					hashTable.update((long) id, (int) first.value);
				}
			}
		}
		else {
			int currLevel = FIXED_UPTO_LEVEL + 1;
			int numKeepUnchanged = (1 << currLevel) - 1;
			int numCurrLevelEntries = 0;
			// TODO: take the width of sketch as an input parameter.
			// currently just half the size of the total num entries per level.
			AMS ams = new AMS(1 << (currLevel-1), 1, seed);
			while (!queue.isEmpty()) {
				BinaryTree first = queue.remove(0);
				
				if (first.left != null) {
					queue.add(first.left);
				}
				if (first.right != null) {
					queue.add(first.right);
				}
				int id = first.id;
				if (id < numKeepUnchanged) {
					continue;
				}
				else {
					numCurrLevelEntries++;
					ams.update((long) id, (int) first.value);
					if (numCurrLevelEntries == (1 << currLevel)) {
						hashedLevels.put(currLevel, ams);
						currLevel++;
						numCurrLevelEntries = 0;
						ams = new AMS(1 << (currLevel-1), 1, seed);
					}
				}
			}
		}
	}
	
	public void addLaplacian(double scale) {
		if (isOneHashTable) {
			hashTable.addLaplacian(scale);
		}
		else {
			Iterator<Integer> keyIter = hashedLevels.keySet().iterator();
			while (keyIter.hasNext()) {
				Integer key = keyIter.next();
				AMS ams = hashedLevels.get(key);
				ams.addLaplacian(scale);
			}
		}
	}
	
	public int count(int level, int id) {
		boolean useMedian = true;
		AMS ams = hashedLevels.get(level);
		int count = ams.count(id, useMedian);
		return count;
	}
	
	// sum of entries in range [startIndex, endIndex], indices start from 0.
	public double rangeQuery(int startIndex, int endIndex) {
		boolean useMedian = true;
		if (startIndex == endIndex) {
			if (isOneHashTable) {
				return hashTable.count(startIndex, useMedian);
			}
			else {
				AMS leafAms = hashedLevels.get(maxLevel);
				return leafAms.count(startIndex, useMedian);
			}
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
			
			int numAdded = maxLevel - st.length();
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
		
		BinaryTree subTree = dyadicIntervals.intervalTree;
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
		
		if (isOneHashTable) {
			double leftVal = subTree.left.getSum_Hashed(startBitString, -1, 1, FIXED_UPTO_LEVEL+1, hashTable);
			double rightVal = subTree.right.getSum_Hashed(endBitString, 1, 1, FIXED_UPTO_LEVEL+1, hashTable);
			return leftVal+rightVal;
		}
		else {
			double leftVal = subTree.left.getSum_Hashed(startBitString, -1, 1, hashedLevels);
			double rightVal = subTree.right.getSum_Hashed(endBitString, 1, 1, hashedLevels);
			//		System.out.println("Left:" + leftVal + ", Right:" + rightVal);
			return leftVal+rightVal;
		}
	}
}
