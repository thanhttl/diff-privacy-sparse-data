package wavelets;

import java.util.HashMap;
import java.util.Random;

import distributions.Laplace;

public class sparse_1d_wavelet {
	int logm;
	HashMap<Integer, Double> coeffMap; // use hashmap instead of array for sparse representation
	double[] coeffs;
	double[] normalize;
	Pair[] transformCoeff; // <<xlevel, offset>, val>
	boolean runInverse;
	
	public sparse_1d_wavelet(int logm, boolean runInverse) {
		this.logm = logm;
		this.runInverse = runInverse;
		
		normalize = new double[logm+1];
		for (int i=0; i<logm; i++) {
			normalize[i] = Math.pow(2, -(i+1));
		}
		normalize[logm] = Math.pow(2, -logm);
		if (!runInverse)
			coeffMap = new HashMap<Integer, Double>(logm*1000000/2);
		else
			coeffs = new double[1 << logm];

		transformCoeff = new Pair[logm+1];
		for (int i=0; i<=logm; i++) {
			transformCoeff[i] = new Pair(new Coeff1(0,0), 0);
		}
	}
	
	public sparse_1d_wavelet(int logm, HashMap<Integer,Double> map, boolean runInverse) {
		this.logm = logm;
		this.runInverse = runInverse;
		
		normalize = new double[logm+1];
		for (int i=0; i<logm; i++) {
			normalize[i] = Math.pow(2, -(i+1));
		}
		normalize[logm] = Math.pow(2, -logm);
		coeffMap = map;

		transformCoeff = new Pair[logm+1];
		for (int i=0; i<=logm; i++) {
			transformCoeff[i] = new Pair(new Coeff1(0,0), 0);
		}
	}
	
	public sparse_1d_wavelet(int logm, double[] coeffArr) {
		this.logm = logm;
		runInverse = true;
		normalize = new double[logm+1];
		for (int i=0; i<logm; i++) {
			normalize[i] = Math.pow(2, -(i+1));
		}
		normalize[logm] = Math.pow(2, -logm);
		coeffs = coeffArr;
		transformCoeff = new Pair[logm+1];
		for (int i=0; i<=logm; i++) {
			transformCoeff[i] = new Pair(new Coeff1(0,0), 0);
		}
	}
	
	// pos is the position in the input array
	public Pair[] haar(int pos) {
		Pair[] result = transformCoeff;
		for (int i=0; i<logm; i++) {
			Coeff1 coeff = (Coeff1) result[i].coeff;
			coeff.level = i;
			coeff.offset = (pos >> (i+1));
			result[i].setValue(1.0- (((pos>>i)&1)<<1));
		}
		Coeff1 coeff = (Coeff1) result[logm].coeff; 
		coeff.level = logm;
		coeff.offset = 0;
		result[logm].setValue(1.0);
		return result;
	}
	

	// this method just generates the Laplace noise. for measuring time to noise adding step
	public void addLaplacian_Time(double scale, Random rand) {
		// weighted as in differential privacy - wavelet paper
		long m = (long) 1 << logm;
		System.out.println(m);
		long currWeight =  m;
		// scale is \lambda in the paper
		Laplace.getSample(scale/currWeight, rand);
		Laplace.getSample(scale/currWeight, rand);
		int currIndex = 2;
		
		for (int i=1; i<=logm; i++) {
			long numEntries = (long) 1 << i;
			currWeight /= 2;
			double thisScale = scale/currWeight;
			for (long j=0; j<numEntries; j++) {
				Laplace.getSample(thisScale, rand);
				currIndex++;
			}	
		}
	}
	
	// this method stores the whole noisy intermediate wavelet. work when it can fit in memory.
	// scale = (1+logn) / epsilon
	 public void addLaplacian(double scale, Random rand) {
		// weighted as in differential privacy - wavelet paper
		int m = (1 << logm);
		int currWeight =  m;
		// scale is \lambda in the paper
		if (!runInverse) {
			coeffMap.put(0, coeffMap.get(0)+Laplace.getSample(scale/currWeight, rand));
			coeffMap.put(1, coeffMap.get(1)+Laplace.getSample(scale/currWeight, rand));
		}
		else {
			coeffs[0] += Laplace.getSample(scale/currWeight, rand);
			coeffs[1] += Laplace.getSample(scale/currWeight, rand);
		}
		int currIndex = 2;
		
		for (int i=1; i<logm; i++) {
			int numEntries = 1 << i;
			currWeight /= 2;
			double thisScale = scale/currWeight;
			for (int j=0; j<numEntries; j++) {
				if (!runInverse)
					coeffMap.put(currIndex, coeffMap.get(currIndex) + Laplace.getSample(thisScale, rand));
				else
					coeffs[currIndex] += Laplace.getSample(thisScale, rand);
				currIndex++;
			}
		}	
	}
	 
	 public static int getWeight(int level, int logm) {
		 return (1<<logm) >> level;
	 }
	 
	// inverse of the transform function
	public void inverse(double[] data) {
		int m = (1 << logm);
//		if (data.length != m)
//			data = new double[m];
		for (int i=0; i<m; i++) {
			Pair[] hx = haar(i);
			double sum = 0;
			for (int j=0; j<hx.length; j++) {
				int offset = ((Coeff1) hx[j].coeff).offset;
				int level = ((Coeff1) hx[j].coeff).level;
				int coeff;
				if (level==logm) coeff = 0;
				else coeff = (1 << (logm-1-level)) + offset;
				sum += (hx[j].value*coeffs[coeff]);
			}
			if (data != null)
				data[i] = sum;
		}
	}
	  
	public void inverse_Time() {
		long m = (long) 1 << logm;
		for (long i=0; i<m; i++) {
			Pair[] hx = haar((int) Math.min(i, Integer.MAX_VALUE));
			double sum = 0;
			for (int j=0; j<hx.length; j++) {
				int offset = ((Coeff1) hx[j].coeff).offset;
				int level = ((Coeff1) hx[j].coeff).level;
				int coeff;
				if (level==logm) coeff = 0;
				else coeff = (1 << (logm-1-level)) + offset;
				sum += (hx[j].value*1.0); // coeff is set to 1.0, just for time measurement
			}
		}
	}
	
	 public void update(int pos, double val) {
		Pair[] tran = transform(pos, val);
		for (Pair p: tran) {
			Coeff1 coeff = (Coeff1) p.coeff;
			double value = p.value;
			int coeffInt;
			if (coeff.level == logm)
				coeffInt = 0;
			else
				coeffInt = (1<<(logm-1-coeff.level)) + coeff.offset;
			if (runInverse)
				coeffs[coeffInt] += value;
			else {
				if (coeffMap.containsKey(coeffInt))
					coeffMap.put(coeffInt, coeffMap.get(coeffInt)+value);
				else
					coeffMap.put(coeffInt, value);
			}
		}
	 }
	
	private Pair[] transform(int pos, double val) {
		Pair[] hx = haar(pos);
		Pair[] result = transformCoeff;
		for (int i=0; i<hx.length; i++) {
	//		int xl = Math.min(i, logn-1);
			int xl = i;
			Pair p = hx[i];
			Pair newP = result[i];
			//Coeff1 c = (Coeff1) newP.coeff; c.setLevel(xl); c.setOffset(p.pos);
			newP.setValue(p.value*val*normalize[xl]); //new Pair(new Coeff1(xl, p.pos), p.value*val*normalize[xl]);
		}
		return result;
	}
	
//	public void print() {
//		Iterator<Integer> iter = coeffs.keySet().iterator();
//		while (iter.hasNext()) {
//			Integer coeff = iter.next();
//			double val = coeffs.get(coeff);
//			System.out.println("Coeff=" + coeff + ", Val=" + val);
//		}
//		
//	}
	
	public static void main(String[] args) {
		sparse_1d_wavelet wavelet = new sparse_1d_wavelet(3, true);
		wavelet.update(0, 9); wavelet.update(1, 3);
		wavelet.update(2, 6); wavelet.update(3, 2);
		wavelet.update(4, 8); wavelet.update(5, 4);
		wavelet.update(6, 5); wavelet.update(7, 7);
		double scale = (1+3)/0.1;
		wavelet.addLaplacian(scale, new Random(24345));
		double[] data = new double[8]; 
		wavelet.inverse(data);
		for (int i=0; i<data.length; i++) {
			System.out.println(data[i]);
		}
		
	}
	
	public static class Pair {
		Coeff coeff;
		double value;
		
		public Pair(Coeff c, double v) {
			coeff = c; value = v;
		}		
		public void setValue(double v) {
			value = v;
		}
		public void setCoeff(Coeff c) {
			coeff = c;
		}
		public String toString() {
			return "(" + coeff + ", " + value + ")";
		}
	}
	
	public static class Coeff {
	}
	
	public static class Coeff1 extends Coeff {
		int level;
		int offset;
		public Coeff1(int xlevel, int offset) {
			this.level = xlevel; this.offset = offset;
		}
		
		public String toString() {
			return "level=" + level + ", offset=" + offset;
		}
		public void setLevel(int l) {
			level = l;
		}
		public void setOffset(int o) {
			offset = o;
		}
	}
	
	public static class Coeff2 extends Coeff {
		int xlevel, ylevel;
		int xoffset, yoffset;
		
		public Coeff2(int xl, int yl, int xoffset, int yoffset) {
			xlevel = xl; ylevel = yl; this.xoffset = xoffset; this.yoffset = yoffset;
		}
		public String toString() {
			return "(xlevel=" + xlevel + ", xoff=" + xoffset + "ylevel=" + ylevel + ", yoff=" + yoffset + ")";
		}
		
		public void reset(int xl, int yl, int xoffset, int yoffset) {
			xlevel = xl; ylevel = yl; this.xoffset = xoffset; this.yoffset = yoffset;
		}
	}
}
