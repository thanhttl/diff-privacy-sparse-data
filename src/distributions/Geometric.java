package distributions;

import java.util.Random;

// two-sided geometric distribution
// Pr(Z=z) = (1-alpha)/(1+alpha) * alpha^|z|
// alpha in [0,1]
public class Geometric {
	public static double alpha;
	public static double ProbAt0;
	public static double ProbGreater0;
	public static double ProbGeq0;
	public static double logAlpha;

	
	public static double getPDF(int z, double alpha) {
		return (1-alpha)/(1+alpha) * Math.pow(alpha, Math.abs(z));
	}
			
	public static void setAlpha(double a) {
		alpha = a;
		ProbAt0 = (1-alpha)/(1+alpha);
		ProbGreater0 = alpha/(1+alpha);
		ProbGeq0 = ProbAt0 + ProbGreater0;
		logAlpha = Math.log(alpha);
	}
	
	public static int getSample(Random rand) {
		double randDouble =  rand.nextDouble();
		
		if (randDouble < ProbAt0)
			return 0;
		int sign;
		if (randDouble < ProbGeq0) {
			sign = 1;
			randDouble -= ProbAt0;
		}
		else { 
			sign = -1;
			randDouble -= (ProbGeq0);
		}
		//System.out.println("alpha/(1+alpha)=" + ProbGreater0 + ", rand=" + randDouble);
		int n = (int) Math.floor(Math.log(alpha - (1+alpha)*randDouble)/logAlpha);
		return sign*n;
	}
	
	public static void main(String[] args) {
	}
	
}
