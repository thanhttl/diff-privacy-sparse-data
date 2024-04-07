package main;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Random;
import java.util.StringTokenizer;

import datahandler.DataGeneration;
import datastructures.AMS;
import datastructures.DyadicRange_Full;
import datastructures.HaarWavelet;
import datastructures.HashedDyadicRange;
import distributions.Laplace;


public class RunSketches {
	//static int MAX_NUM = 1 << 25;
	
	public static class PairData {
		//long dest;
		//long source;
		int count;
		
		public PairData(long d, long s, int c) {
			//dest = d; source = s; 
			count = c;
		}
	}
	
	public enum Mechanism {Laplace, Hash_Lap, Dyadic_Lap, Dyadic_Hash_Lap, Wavelet_Lap };
	
	public static void main(String[] args) {
		String fileName = "data/onthemap/nj/od/nj_od_main_ja_2008_1.csv";
		//runOnTheMapData_Point(fileName, Mechanism.Hash_Lap);
		double[] epsilonVals = {0.1, 0.5, 1};

		for (double epsilon: epsilonVals) {
//			runOnTheMapData_Point(fileName, Mechanism.Hash_Lap, epsilon);
			runOnTheMapData_Range(fileName, Mechanism.Dyadic_Lap, epsilon, 1000000);
//			runSynthetic_Point(Mechanism.Hash_Lap, epsilon);
		}

	}
	
	
	public static void runSynthetic_Point(Mechanism mech, double epsilon) {
		
		int seed = 123456;
		double zeroPercent = 0.9;
		//HashMap<Integer, PairData> origData = new HashMap<Integer, PairData>();
		double[] counts = new double[1000000];
		double mean = 20;
		double stdDev = 1;
		Random rand = new Random(1234);
		
		for (int i=0; i<counts.length*(1-zeroPercent); i++) {
			double sample = Math.round(rand.nextGaussian()*stdDev + mean);
			if (sample < 0)
				counts[i] = 1;
			else
				counts[i] = sample;
		}
		
		int width = 10000;
		int sensitivity = 2;
		double[] output = new double[counts.length];
		pointQueries(counts, output, mech, epsilon, sensitivity, seed, width);
	}


	public static void runOnTheMapData_Point(String fileName, Mechanism mech, double epsilon, int numEntries) {

		int seed = 123456;
		
		double[] counts = new double[numEntries];
		
		// read the file
		BufferedReader reader=null;
		try {
			reader = new BufferedReader(new FileReader(fileName));
			reader.readLine(); // first line contains attribute names
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
				
		for (int id=0; id<numEntries; id++) {
			int count = 0;
			try {
				String line = reader.readLine();
				if (line == null)
					break;
				StringTokenizer st = new StringTokenizer(line, ", ");
				
				long dest = Long.valueOf(st.nextToken());
				long source = Long.valueOf(st.nextToken());
				count = Integer.valueOf(st.nextToken());
				//origData.put(id, new PairData(dest, source, count));
				counts[id] = count;
			}
			catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			
		}
		
		int width = 10000;
		int sensitivity = 2;
		double[] output = new double[counts.length];
		pointQueries(counts, output, mech, epsilon, sensitivity, seed, width);
	}
	
	
	public static void runOnTheMapData_Range(String fileName, Mechanism mech, double epsilon, int numEntries) {
		//double epsilon = 0.1;
		int seed = 123456;
		
		//HashMap<Integer, PairData> origData = new HashMap<Integer, PairData>();
		double[] counts = new double[numEntries];
		
		// read the file
		BufferedReader reader=null;
		try {
			reader = new BufferedReader(new FileReader(fileName));
			reader.readLine(); // first line contains attribute names
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
				
		for (int id=0; id<numEntries; id++) {
			int count = 0;
			try {
				String line = reader.readLine();
				if (line == null)
					break;
				StringTokenizer st = new StringTokenizer(line, ", ");
				
				long dest = Long.valueOf(st.nextToken());
				long source = Long.valueOf(st.nextToken());
				count = Integer.valueOf(st.nextToken());
				//origData.put(id, new PairData(dest, source, count));
				counts[id] = count;
			}
			catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}		
		}
		
		if (mech == Mechanism.Laplace) {
			double scale = 2/epsilon;
			Random rand = new Random(seed);
			double L1Dist = 0;
			double L2Dist = 0;
			
			double queryNoise = 0;
			for (int i=0; i<numEntries; i++) {
				queryNoise += Laplace.getSample(scale, rand);
				L1Dist += Math.abs(queryNoise);
				L2Dist += Math.pow(queryNoise, 2);
			}
			System.out.println("epsilon=" + epsilon + "; L1Dist=" + L1Dist + "; L2Dist=" + Math.sqrt(L2Dist));
			return;
		}
//		else if (mech == Mechanism.Dyadic_Lap) {
//			DyadicRange dyadicRanges = new DyadicRange(counts, epsilon, 2, seed);
//			//dyadicRanges.addLaplacian(epsilon, 2);
//			
//			double L1Dist = 0;
//			double L2Dist = 0;
//			
//			double trueSum = 0;
//			for (int i=0; i<numEntries; i++) {
//				//System.out.println(i);
//				double estSum = dyadicRanges.rangeQuery(0, i);
//				trueSum += counts[i];
//				L1Dist += Math.abs(estSum - trueSum);
//				L2Dist += Math.pow(estSum - trueSum, 2);
//			}
//			
//			System.out.println("epsilon=" + epsilon + "; TotalCount=" + dyadicRanges.getCount() + 
//					"; L1Dist=" + L1Dist + "; L2Dist=" + Math.sqrt(L2Dist));
//		}
		else if (mech == Mechanism.Wavelet_Lap) {
			HaarWavelet wavelet = new HaarWavelet(counts.clone(), seed);
			
			int depth = (int) (Math.log(counts.length)/Math.log(2));
			if (counts.length > (1 << depth)) {
				depth += 1; // need to add zero entries to the end of the array
			}
			
			wavelet.addLaplacian(epsilon, false);
			
			double L1Dist = 0;
			double L2Dist = 0;
			
			double trueSum = 0;
			double estSum = 0;
			double[] data = wavelet.getData();
			
			for (int i=0; i< numEntries; i++) {
				//System.out.println(i);
				estSum += data[i];
				trueSum += counts[i];
				L1Dist += Math.abs(estSum - trueSum);
				L2Dist += Math.pow(estSum - trueSum, 2);
			}
			
			System.out.println("epsilon=" + epsilon + "; TotalCount=" + trueSum + 
					"; L1Dist=" + L1Dist + "; L2Dist=" + Math.sqrt(L2Dist));
		}
//		else if (mech == Mechanism.Dyadic_Hash_Lap) {
//			boolean useOneHashTable = true;
//			DyadicRange dyadicRanges = new DyadicRange(counts, seed);
//			HashedDyadicRange hashedRanges = new HashedDyadicRange(dyadicRanges, seed, useOneHashTable);
//			double scale = 2*dyadicRanges.getNumLevels()/epsilon;
//			hashedRanges.addLaplacian(scale);
//			
//			double L1Dist = 0;
//			double L2Dist = 0;
//			double trueSum = 0;
//			for (int i=0; i< numEntries; i++) {
//				// System.out.println(i);
//				double estSum = hashedRanges.rangeQuery(0, i);
//				trueSum += counts[i];
//				L1Dist += Math.abs(estSum - trueSum);
//				L2Dist += Math.pow(estSum - trueSum, 2);
//			}
//			
//			System.out.println("TotalCount=" + dyadicRanges.getCount() + 
//					"; L1Dist=" + L1Dist + "; L2Dist=" + Math.sqrt(L2Dist));
//			
//		}
	}
	
	public static void runSketchesUniformData() {
		int numBlocks = 100; 
		int minPerBlock = 50; 
		int maxPerBlock = 50;
		int seed = 1234;
		Random rand = new Random(999);
		
		int numBuckets = 1000;
		int depth = 10;
		
		double epsilon = 1;
		
		int[][] contTable = DataGeneration.genUniformTable(numBlocks, minPerBlock, maxPerBlock, rand);
		
		
		int[] depths = {2, 4, 6, 8, 10, 15, 20};
		for (int d: depths) {
			depth = d;
			
			AMS sketch = new AMS(numBuckets, depth, seed);	
			// add entries in contTable to sketch
			for (int i=0; i<numBlocks; i++) {
				for (int j=0; j<numBlocks; j++) {
					// a unique entry id, i*numBlocks + j
					int item = i*numBlocks + j;
					sketch.update(item, contTable[i][j]);
				}
			}
			
			// add Laplacian noise to the sketch
			double scale = 2*depth/epsilon; 
			sketch.addLaplacian(scale);
			//System.out.println("Sketch after adding noise");
			//sketch.show();
			
			// compute the distance
			int L1Dist = 0;
			int L2Dist = 0;
			boolean useMedian = true;
			for (int i=0; i<numBlocks; i++) {
				for (int j=0; j<numBlocks; j++) {
					int item = i*numBlocks + j;
					int estimateCount = sketch.count(item, useMedian);
					L1Dist += Math.abs(estimateCount - contTable[i][j]);
					L2Dist += Math.pow(estimateCount - contTable[i][j],2);
				}
			}
		
			System.out.println("Depth=" + d + "; Total count=" + sketch.getTotalCount() + "; L1Dist=" + L1Dist +
				"; L2Dist=" + Math.sqrt(L2Dist));
		}
			
	}
	
	public static void pointQueries(double[] counts, double[] output, Mechanism mech, double epsilon, int sensitivity, int seed, int width) {
		
		if (output.length != counts.length) {
			System.out.println("RunSketches: output has different size from input");
			System.exit(1);
		}
		
		if (mech == Mechanism.Laplace) {
			double scale = sensitivity/epsilon;
			Random rand = new Random(seed);
			double L1Dist = 0;
			double L2Dist = 0;
			double mean = 0; double squared = 0;
			
			for (int i=0; i< counts.length; i++) {
				// noise can be positive or negative
				double noise = Laplace.getSample(scale, rand);
				output[i] = counts[i] + noise;
				L1Dist += Math.abs(noise);
				L2Dist += Math.pow(noise, 2);
				mean += Math.abs(noise);
				squared += Math.pow(noise, 2);
			}
			
//			mean = mean/counts.length;	
//			double variance = squared/counts.length - mean*mean;
//			System.out.println("epsilon=" + epsilon + "; L1=" + L1Dist + "; L2=" + Math.sqrt(L2Dist) + "; Mean=" + mean + "; Var=" + variance);
			return;
		}
		
		// mech == Mechanism.Sketch_Lap
		int[] depths = {1, 3, 5, 10};
		for (int d: depths) {
			int depth = d;                                                                                            
		
		AMS sketch = new AMS(width, depth, seed);
		for (int i=0; i< counts.length; i++) {
			int count = (int) counts[i];
			sketch.update(i, count);
		}
		
		// add Laplacian noise
		double scale = 2*depth/epsilon;
		// double scale = 2/epsilon;
		sketch.addLaplacian(scale);
		
		// compute the errors
		double L1Dist = 0;
		double L2Dist = 0;
		boolean useMedian = false;
		
		for (int i=0; i<counts.length; i++) {
			double estCount = sketch.count(i, useMedian);
			output[i] = estCount;
			//double trueCount = origData.get(i).count;
			double trueCount = counts[i];
			L1Dist += Math.abs(estCount - trueCount);
			L2Dist += Math.pow(estCount - trueCount, 2);
		}
		
		System.out.println("epsilon=" + epsilon + "; Depth=" + depth + "; Total count=" + sketch.getTotalCount() + 
				"; L1Dist=" + L1Dist + "; L2Dist=" + Math.sqrt(L2Dist));
		}
		
	}
}
