package anonymization;

import java.util.Random;

import distributions.Geometric;
import distributions.Laplace;

public class Mechanisms {
	
	public static double getLaplace(double input, double scale, Random rand) {
		double noise = Laplace.getSample(scale, rand);
		return (input + noise);
	}
	
	public static int getGeometric(int input, double alpha, Random rand) {
//		int noise = Geometric.getSample(alpha, rand);
//		return (input + noise);
		return 0;
	}
}
