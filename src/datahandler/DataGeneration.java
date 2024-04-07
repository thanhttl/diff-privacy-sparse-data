package datahandler;

import java.util.Random;

public class DataGeneration {
	/**
	 * return the pairs #Commuters[destination][source]
	 * number of people per cell is uniformly drawn from a range [min, max]
	 */
	public static int[][] genUniformTable(int numBlocks, int minPerBlock, int maxPerBlock, Random rand) {
		int[][] synData = new int[numBlocks][numBlocks];
		int rangePerBlock = maxPerBlock - minPerBlock;
		for (int i=0; i<numBlocks; i++)
			for (int j=0; j<numBlocks; j++)
				synData[i][j] = minPerBlock + (int) (rand.nextDouble()*rangePerBlock);
		return synData;
	}
	
	public static int[] genGeometricArray(int numBlocks, int totalPeople, double prob) {
		int[] data = new int[numBlocks];
		double totalProb = 1 - Math.pow(1-prob, numBlocks);
		int sumPeople = 0;
		for (int i=1; i<numBlocks; i++) {
			data[i-1] = (int) Math.round(Math.pow(1-prob, i-1)*prob/totalProb * totalPeople);
			sumPeople += data[i-1];
		}
		data[numBlocks-1] = totalPeople - sumPeople;
		return data;
	}
	
	public static void main(String[] args) {
		int[] geometricData = genGeometricArray(5, 500, 0.4);
		for (int n: geometricData) {
			System.out.print(n + ", ");
		}
	}
}
