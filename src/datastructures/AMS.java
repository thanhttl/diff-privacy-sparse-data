package datastructures;

import java.util.Arrays;
import java.util.Random;

import distributions.Laplace;

public class AMS {
	Random rand;
	int[][] counts;
	int count;
	int depth;
	int buckets;
	int[][] test;
	
	public AMS(int buckets, int depth, int seed) {
		rand = new Random(seed);
		this.depth = depth;
		this.buckets = buckets;
		
		counts = new int[depth][buckets];
		for (int i=0; i<depth; i++) {
			for (int j=0; j<buckets; j++) {
				counts[i][j] = 0;
			}
		}
		count = 0;
		// need to initialize the hash functions
		test = new int[4][depth];
		for (int i=0; i<4; i++) {
			for (int j=0; j<depth; j++) {
				test[i][j] = rand.nextInt();
				if (test[i][j] < 0)
					test[i][j] = -test[i][j];
			}
		}
	}
	
	  
	public void update(long item, int diff) {
		// update the sketch
		// hash to one bucket in each row
		// then multiply by {+1, -1} chosen at random
		count += diff;
		for (int i=0; i<depth; i++)  {
			int hash = (int) PRNG.hash31(test[0][i], test[1][i], item);
			hash = hash % buckets;
			//long mult = PRNG.fourwise(test[2][j], test[3][j], test[4][j], test[5][j], item);
			int mult = (int) PRNG.hash31(test[2][i], test[3][i], item);
			// System.out.println(mult & 1);
			if ((mult & 1) == 1)
				counts[i][hash] += diff;
			else
				counts[i][hash] -= diff;
			
//			int hash = Math.abs(rand.nextInt());
//			hash = hash % buckets;
//			boolean bool = rand.nextBoolean();
//			if (bool)
//				counts[j][hash] += diff;
//			else
//				counts[j][hash] -= diff;
		}
	}
	
	public int count(int item, boolean useMedian) {
		int[] estimates = new int[depth];
		for (int i=0; i<depth; i++) {
			int hash = (int) PRNG.hash31(test[0][i], test[1][i], item);
			hash = hash % buckets;
			//int mult = (int) PRNG.fourwise(test[2][i], test[3][i], test[4][i], test[5][i], item);
			int mult = (int) PRNG.hash31(test[2][i], test[3][i], item);
			if ((mult & 1) == 1)
				estimates[i] = counts[i][hash];
			else
				estimates[i] = -counts[i][hash];
			
//			int hash = Math.abs(rand.nextInt());
//			hash = hash % buckets;
//			boolean bool = rand.nextBoolean();
//			if (bool)
//				estimates[i] = counts[i][hash];
//			else
//				estimates[i] = -counts[i][hash];
		}
		if (useMedian) {
			if (depth == 1)
				return estimates[0];
			else if (depth == 2)
				return (estimates[0]+estimates[1])/2;
			else {
				Arrays.sort(estimates);
				return(estimates[depth/2]);	
			}
		}
		else { // return the mean
			double sum = 0;
			for (int i=0; i<depth; i++)
				sum += estimates[i];
			return (int) (sum/depth);
		}
	}
	
	
	public void addLaplacian(double scale) {
		for (int i=0; i<depth; i++) {
			for (int j=0; j<buckets; j++) {
				counts[i][j] += Math.round(Laplace.getSample(scale, rand));
			}
		}
	}
	
	
	public void show() {
		for (int i=0; i<depth; i++) {
			for (int j=0; j<buckets; j++) {
				System.out.print(counts[i][j] + "\t");
			}
			System.out.println();
		}
	}
	
	public int getTotalCount() {
		return count;
	}
	
	public static void main(String[] args) {
		
		AMS ams = new AMS(100, 7, 123456);
		Random rand = new Random(9999);
		
		int numItems = 1000;
		int[] inputs = new int[numItems];
		
		for (int i=0; i<numItems; i++) {
			int diff = rand.nextInt(100);
			inputs[i] = diff;
			
			ams.update(i, diff);
		}
		ams.show();
		
		int L1Dist = 0;
		int L2Dist = 0;
		for (int i=0; i<numItems; i++) {
			int count = ams.count(i, true);
			L1Dist += Math.abs(count - inputs[i]);
			L2Dist += Math.pow(count - inputs[i], 2);
		}
		
		System.out.println("Total count=" + ams.getTotalCount());
		System.out.println("L1 Dist=" + L1Dist);
		System.out.println("L2 Dist=" + Math.sqrt(L2Dist));
		
	}
}
