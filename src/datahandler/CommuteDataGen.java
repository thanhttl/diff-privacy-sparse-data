package datahandler;

import java.util.Random;

public class CommuteDataGen {
	/**
	 * return the pairs #Commuters[destination][source]
	 */
	public static int[][] genData(int numBlocks, int minPerBlock, int maxPerBlock, Random rand) {
		int[][] synData = new int[numBlocks][numBlocks];
		int rangePerBlock = maxPerBlock - minPerBlock;
		for (int i=0; i<numBlocks; i++)
			for (int j=0; j<numBlocks; j++)
				synData[i][j] = minPerBlock + (int) (rand.nextDouble()*rangePerBlock);
		return synData;
	}
}
