package datastructures;

import java.util.Random;

import wavelets.simple_haar;

public class HaarWavelet {
	
	simple_haar wavelet;
	Random rand;
	
	public HaarWavelet(double[] inputArr, int seed) {
		wavelet = new simple_haar();
		wavelet.wavelet_calc(inputArr);
		rand = new Random(seed);
	}
	
	public void addLaplacian(double scale, boolean isWeighted) {
		wavelet.addLaplacian(scale, isWeighted, rand);
		wavelet.inverse();
	}
	
	// call this method after adding laplacian noise and calling inverse function.
	public double rangeQuery(int startIndex, int endIndex) {
		double[] data = wavelet.getData();
		double sum = 0;
		for (int i=startIndex; i<=endIndex; i++)
			sum += data[i];
		return sum;
	}
	
	public double[] getData() {
		return wavelet.getData();
	}
}
