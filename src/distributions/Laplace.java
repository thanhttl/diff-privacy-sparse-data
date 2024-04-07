package distributions;

import java.util.Random;

// assume mean=0
public class Laplace {
	public static double getPDF(double x, double scale) {
		return Math.exp(-Math.abs(x)/scale)/(2*scale);
	}
	
	public static double getSample(double scale, Random rand) {
		double uniform = rand.nextDouble() - 0.5;
		if (uniform < 0)
			return scale*Math.log(1-2*Math.abs(uniform));
		else
			return -scale*Math.log(1-2*Math.abs(uniform));
	}
	
	public static void main(String[] args) {
		double stdDev = 3;
		Random rand = new Random(123456);
		for (int i=0; i<20; i++) {
			System.out.println(getSample(stdDev, rand));
		}
	}
}
