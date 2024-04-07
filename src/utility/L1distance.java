package utility;

public class L1distance {
	
	// L1(a, b) = \sum_i | a_i - b_i |
	public static double getL1Dist(int[][] t1, double[][] t2) {
		double dist = 0;
		int totalNum = 0;
		for (int i=0; i<t1.length; i++) {
			for (int j=0; j<t1[0].length; j++) {
				dist += Math.abs(t1[i][j] - t2[i][j]);
				totalNum += t1[i][j];
			}
		}
		
		return dist/totalNum;
	}
	
	public static double getL1Dist(int[] t1, double[] t2) {
		double dist = 0;
		int totalNum = 0;
		for (int i=0; i<t1.length; i++) {
			dist += Math.abs(t1[i] - t2[i]);
			totalNum += t1[i];
		}	
		return dist/totalNum;
	}
	
	public static double getL1Dist(int[] t1, int[] t2) {
		double dist = 0;
		int totalNum = 0;
		for (int i=0; i<t1.length; i++) {
			dist += Math.abs(t1[i] - t2[i]);
			totalNum += t1[i];
		}	
		return dist/totalNum;
	}
}
