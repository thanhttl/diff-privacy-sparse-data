package datastructures;

// Psuedo-random number generator
public class PRNG {
	static int MOD=2147483647;
	static int HL=31;
	
	public static long hash31(long a, long b, long x) {
		// return a hash of x using a and b mod (2^31 - 1)
		// may need to do another mod afterwards, or drop high bits
		// depending on d, number of bad guys
		// 2^31 - 1 = 2147483647
		//  result = ((long long) a)*((long long) x)+((long long) b);
		
		long result = a*x + b;
		result = ((result >> HL) + result) & MOD;
		return result;
	}
	
	public static long fourwise(long a, long b, long c, long d, long x) {
		long result = hash31(hash31(hash31(a,b,x),c,x),d,x);
		return result;
	}

	
}
