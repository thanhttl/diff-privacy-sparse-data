	// TODO: Add Laplacian to the wavelet coefficients
	// Run accuracy experiments

package wavelets;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

import main.Run;

import datastructures.UtilityObject.NonZeroEntry;
import distributions.Laplace;

import wavelets.sparse_1d_wavelet.Coeff1;
import wavelets.sparse_1d_wavelet.Coeff2;
import wavelets.sparse_1d_wavelet.Pair;

public class sparse_2d_wavelet_array {
	int logx;
	int logy;
	double[] coeffs;
	double[] normalize;
//	double[] power;
	sparse_1d_wavelet x_wavelet;
	sparse_1d_wavelet y_wavelet;
	Pair[] tempPairs;
	
	public sparse_2d_wavelet_array(int logx, int logy) {
		this.logx = logx;
		this.logy = logy;
		int numItems = logx + logy + 1;
		normalize = new double[numItems];
		
		for (int i=0; i<numItems; i++) {
			// this normalization is to get unit vectors of coefficients
			// normalize[i] = Math.pow(2, -0.5*(i+2));
			// according to the wavelet paper (ICDE10)
			normalize[i] = Math.pow(2, -(i+2));
//			power[i] = Math.pow(2, i);
		}
		coeffs = new double[(1<<logx)*(1<<logy)];
		
		x_wavelet = new sparse_1d_wavelet(logx, true);
		y_wavelet = new sparse_1d_wavelet(logy, true);
		tempPairs = new Pair[(logx+1)*(logy+1)];
		for (int i=0; i<tempPairs.length; i++) {
			tempPairs[i] = new Pair(new Coeff2(0,0,0,0), 0);
		}
	}
	
	public void update(int x, int y, double z) {
		Pair[] tran = transform(x, y, z);
		for (Pair p: tran) {
			Coeff2 coeff = (Coeff2) p.coeff;
			
			int coeffX;
			if (coeff.xlevel == logx)
				coeffX = 0;
			else
				coeffX = (1<<(logx-1-coeff.xlevel)) + coeff.xoffset;
			int coeffY;
			if (coeff.ylevel == logy)
				coeffY = 0;
			else 
				coeffY = (1<<(logy-1-coeff.ylevel)) + coeff.yoffset;
			
			int coeffXY = coeffX*(1<<logy)+coeffY;
			
			double val = p.value;
			coeffs[coeffXY] += val;
							
		}
	}
	
	private Pair[] transform(int x, int y, double val) {
		Pair[] hx = x_wavelet.haar(x);
		Pair[] hy = y_wavelet.haar(y);
		Pair[] result = tempPairs;
		
		int index = 0;
		for (int i=0; i<hx.length; i++) {
			for (int j=0; j<hy.length; j++) {
				Pair xpair = hx[i];
				Pair ypair = hy[j];
				Coeff1 xcoeff = (Coeff1) xpair.coeff;
				Coeff1 ycoeff = (Coeff1) ypair.coeff;
	//			Coeff2  coeff = new Coeff2(xcoeff.level, ycoeff.level, xcoeff.offset, ycoeff.offset); //xpair.coeff, ypair.coeff);
				int xl = Math.min(xcoeff.level, logx-1);
				int yl = Math.min(ycoeff.level, logy-1);
				double newVal = xpair.value*ypair.value*val*normalize[xl+yl];
				((Coeff2) result[index].coeff).reset(xcoeff.level, ycoeff.level, xcoeff.offset, ycoeff.offset);
				result[index].value = newVal;
				index++;
				//result[index++] = new Pair(coeff, newVal);
			}
		}
		return result;
	}
	
	public double inverse(int x, int y) {
		Pair[] coeffTran = transform(x, y, 1);
		double result = 0;
		for (Pair p: coeffTran) {
			Coeff2 coeff = (Coeff2) p.coeff;
			int coeffX;
			if (coeff.xlevel == logx)
				coeffX = 0;
			else
				coeffX = (1<<(logx-1-coeff.xlevel)) + coeff.xoffset;
			int coeffY;
			if (coeff.ylevel == logy)
				coeffY = 0;
			else 
				coeffY = (1<<(logy-1-coeff.ylevel)) + coeff.yoffset;
			
			int coeffXY = coeffX*(1<<logy)+coeffY;
			int xl = Math.min(coeff.xlevel, logx-1);
			int yl = Math.min(coeff.ylevel, logy-1);
			
			result += p.value*coeffs[coeffXY]/normalize[xl+yl];
		}
		return result;
	}
	
	public void addLaplacian(double scale, Random rand) {
		for (int i=0; i<=logx; i++) {
			for (int j=0; j<=logx; j++) {
				int xlevel = Math.max(0, i-1);
				int ylevel = Math.max(0,j-1);
				int weight = sparse_1d_wavelet.getWeight(xlevel, logx)*sparse_1d_wavelet.getWeight(ylevel, logy);
				//int weight = (int) (1.0/normalize[xlevel+ylevel]);
				int numXEntries = 1 << xlevel;
				int numYEntries = 1 << ylevel;
				int coeffX, coeffY;
				for (int u=0; u<numXEntries; u++) {
					for (int v=0; v<numYEntries; v++) {
						if (i == 0)
							coeffX = 0;
						else
							coeffX = (1<<xlevel) + u;
						
						if (j == 0)
							coeffY = 0;
						else 
							coeffY = (1<<ylevel) + v;
						
						int coeffXY = coeffX*(1<<logy)+coeffY;
						//System.out.println("i=" + i + ", j=" + j + ", Weight = " + weight + ", Coeff=" + coeffXY);
						
						coeffs[coeffXY] += Laplace.getSample(scale/weight, rand);
						
					}
				}
				
			}
		}
		
	}
	
//	public int getSize() {
//		return coeffs.size();
//	}
	
	public static void main(String[] args) {
		test();
//		int logx = 2; int logy = 2;
//		sparse_2d_wavelet wavelet = new sparse_2d_wavelet(logx,logy);
//		wavelet.update(0,0,2); wavelet.update(0,1,4); wavelet.update(0,2,1); wavelet.update(0,3,0);
//		wavelet.update(1,0,6); wavelet.update(1,1,8); wavelet.update(1,2,2); wavelet.update(1,3,4);
//		wavelet.update(2,0,1); wavelet.update(2,1,3); wavelet.update(2,2,0); wavelet.update(2,3,6);
//		wavelet.update(3,0,2); wavelet.update(3,1,2); wavelet.update(3,2,4); wavelet.update(3,3,4);
//		Iterator<Integer> iter = wavelet.coeffs.keySet().iterator();
//		while (iter.hasNext()) {
//			int coeff = iter.next();
//			System.out.println(coeff + ", " + wavelet.coeffs.get(coeff));
//		}
//		System.out.println(wavelet.inverse(0, 0));
//		System.out.println(wavelet.inverse(0, 1));
//		System.out.println(wavelet.inverse(1, 0));
//		System.out.println(wavelet.inverse(1, 1));
//		
//		wavelet.addLaplacian(0.8, new Random(23355));
	}
	
	private static void test() {
		int logx = 10; int logy = 10;
		sparse_2d_wavelet_array wavelet = new sparse_2d_wavelet_array(logx, logy);
		double epsilon = 0.1;
		
		int xSize = 1 << logx; int ySize = 1<<logy;
		int inputSize = xSize*ySize;
		
		ArrayList<NonZeroEntry> input = new ArrayList<NonZeroEntry>();
		Random rand = new Random(12345);
		for (int i=0; i<xSize; i++) {
			for (int j=0; j<ySize; j++) {
				double val = rand.nextInt(100);
				input.add(new NonZeroEntry(i*ySize+j, val));
				wavelet.update(i, j, val);
			}
		}
		
		double scale = (1+logx+3)*(1+logy+3)/epsilon;
		wavelet.addLaplacian(scale, new Random(23355));
		
		double[] output = new double[inputSize];
		for (int x=0; x<xSize; x++) {
			for (int y=0; y<ySize; y++) {
				output[x*ySize+y] = wavelet.inverse(x, y);
			}
		}
		
		double[] rangePercents = {1.0/inputSize, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
		for (int i=0; i<rangePercents.length; i++) {
			int rangeLength = (int) Math.max(1, rangePercents[i]*inputSize);
			Run.ErrorProfile ep = new Run.ErrorProfile(10000, true);
			Run.computeRangeQueryError(input, output, rangeLength, ep);
			System.out.println("Range=" + rangeLength + ", Error=" + ep.absError/ep.numQueries);
		}
	}
	
}
