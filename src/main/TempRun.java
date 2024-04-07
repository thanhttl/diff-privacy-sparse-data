package main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Random;

import main.Run.ErrorProfile;

import datastructures.ObjectPool;
import datastructures.UtilityObject.NonZeroEntry;

public class TempRun {
	
	public static void main(String[] args) throws Exception {
		//readTrace();
		computeWaveletError_2d();
	}
	
	public static void computeWaveletError_2d() throws Exception {
		String inputFile = "/state/partition1/ttran/PrivSparseData/data/internet_traffic/combinedTrace.txt";
		String outputFile = "/state/partition1/ttran/PrivSparseData/data/internet_traffic/out/out_wavelet.txt";
		
		//int xSize = 8377; int ySize = 52292;
		int xSize = 5000; int ySize = 50000;
		
		double[][] input = new double[xSize][ySize];
		double[][] output = new double[xSize][ySize];
		
		read_2d(inputFile, input, xSize, ySize);
		read_2d(outputFile, output, xSize, ySize);
		
		double[] rangePercents = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000};
		
		int maxNumQueries = 10000;
		for (int i=0; i<rangePercents.length; i++) {
			//int rangeLength = (int) Math.max(1, rangePercents[i]*xSize);
			int rangeLength = (int) rangePercents[i];
			Run.ErrorProfile ep = new Run.ErrorProfile(maxNumQueries, true);
			computeRangeQueryError_Random(input, output, rangeLength, ep);
			System.out.println("Range=" + rangeLength + ", Error=" + ep.absError/ep.numQueries);
		}
		
		
	}
	
public static void computeRangeQueryError_Random(double[][] input, double[][] output, int rangeLength, ErrorProfile ep) {
		
		
		
		Random rand = new Random(12345);
		
		for (int k=0 ; k<ep.numQueries; k++) {
			double sampleSubsetSum = 0;
			double trueSubsetSum = 0;
			
			int startX = rand.nextInt(input.length - rangeLength);
			int startY = rand.nextInt(input[0].length - rangeLength);
			
			for (int i=startX; i< startX+rangeLength; i++) {
				for (int j=startY; j<startY+rangeLength; j++) {
					sampleSubsetSum += output[i][j];
					trueSubsetSum += input[i][j];
				}
			}
					
			double absError = Math.abs(sampleSubsetSum - trueSubsetSum);
			
			ep.addError(absError, trueSubsetSum);
		}
	}
	
	
public static void computeRangeQueryError(double[][] input, double[][] output, int rangeLength, ErrorProfile ep) {
		
		double sampleSubsetSum = 0;
		double trueSubsetSum = 0;
		
		for (int i=0; i<rangeLength; i++) {
			for (int j=0; j<rangeLength; j++) {
				sampleSubsetSum += output[i][j];
				trueSubsetSum += input[i][j];
			}
		}
				
		double absError = Math.abs(sampleSubsetSum - trueSubsetSum);
		
		ep.addError(absError, trueSubsetSum);		
//		if (pw != null && writeToFile)
//			pw.println(trueSubsetSum + ", " + estSubsetSum);
		
		int numQueries = ep.numQueries; //Math.min(input.length-rangeLength, MAX_NUM_QUERIES);
		for (int i=0; i<numQueries-1; i++) {
			int x;
			if (i+rangeLength < input.length)
				x = i+rangeLength;
			else
				x = input.length-1;
			sampleSubsetSum += (output[x][i+rangeLength] - output[x-rangeLength][i]);
			trueSubsetSum += (input[x][i+rangeLength] - input[x-rangeLength][i]);
			
			for (int j=x-rangeLength+1; j<x; j++) {
				sampleSubsetSum += (output[j][i+rangeLength] - output[j][i]);
				trueSubsetSum += (input[j][i+rangeLength] - input[j][i]);
			}
			for (int j=i+1; j<i+rangeLength; j++) {
				sampleSubsetSum += (output[x][j] - output[x-rangeLength][j]);
				trueSubsetSum += (input[x][j] - input[x-rangeLength][j]);
			}
			
			absError = Math.abs(sampleSubsetSum - trueSubsetSum);
			ep.addError(absError, trueSubsetSum);
//			if (pw != null && writeToFile)
//				pw.println(trueSubsetSum + ", " + estSubsetSum);	
		}	
	}
	
	private static void read_2d(String fileName, double[][] data, int maxX, int maxY) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		//numEntries = Integer.valueOf(reader.readLine());
	
		while (true) {
			String line = reader.readLine();
			if (line == null)
				break;
			else {
				String[] tokens = line.split(",");
				 if (tokens.length == 3) {
					int x = Integer.valueOf(tokens[0]);
					int y = Integer.valueOf(tokens[1]);
					if (x >= maxX || y >= maxY) 
						continue;
					
					double count = Double.valueOf(tokens[2]);
					data[x][y] = count;
				}
				
			}
		}
	}
	
	public static void computeWaveletError_1d() throws Exception {
		String inputFile = "/state/partition1/ttran/PrivSparseData/data/internet_traffic/combinedTrace.txt";
		String outputFile = "/state/partition1/ttran/PrivSparseData/data/internet_traffic/out/out_wavelet.txt";
		
		//int xSize = 8377; int ySize = 52292;
		int xSize = 4000; int ySize = 52292;
		int inputSize = xSize*ySize;
		
		double[] input = new double[inputSize];
		double[] output = new double[inputSize];
		
		read(inputFile, input, xSize, ySize);
		read(outputFile, output, xSize, ySize);
		
		//double[] rangePercents = {1.0/inputSize, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
		double[] rangePercents = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 3000, 4000};
		
		for (int i=0; i<rangePercents.length; i++) {
			int rangeLength = (int) Math.max(1, rangePercents[i]*rangePercents[i]);
			Run.ErrorProfile ep = new Run.ErrorProfile(10000, true);
			Run.computeRangeQueryError(input, output, rangeLength, ep);
			System.out.println("Range=" + rangeLength + ", Error=" + ep.absError/ep.numQueries);
		}
		
		
	}
	
	private static void read(String fileName, double[] data, int maxX, int multiple) throws Exception {
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		//numEntries = Integer.valueOf(reader.readLine());
	
		while (true) {
			String line = reader.readLine();
			if (line == null)
				break;
			else {
				int position=0; double count=0;
				String[] tokens = line.split(",");
				if (tokens.length == 2) {
					position = Integer.valueOf(tokens[0]);
					count = Integer.valueOf(tokens[1]);
				}
				else if (tokens.length == 3) {
					int x = Integer.valueOf(tokens[0]);
					int y = Integer.valueOf(tokens[1]);
					if (x >= maxX || y >= multiple) 
						continue;
					position = x*multiple+ y;
					count = Double.valueOf(tokens[2]);
				}
				data[position] = count;
			}
		}
	}
	
	public static void readTrace() {
		String fileName =  "data/internet_traffic/out/combinedTrace.txt";
		int numEntries = 0;
		try {
			BufferedReader reader=null;
			reader = new BufferedReader(new FileReader(fileName));
			numEntries = Integer.valueOf(reader.readLine());
			
			double sum = 0;
			int numNonZeros = 0;
			while (true) {
				String line = reader.readLine();
				if (line == null)
					break;
				else {
					int position=0, count=0;
					String[] tokens = line.split(",");
					if (tokens.length == 2) {
						position = Integer.valueOf(tokens[0]);
						count = Integer.valueOf(tokens[1]);
					}
					else if (tokens.length == 3) {
						//position = Integer.valueOf(tokens[0])*MULTIPLE + Integer.valueOf(tokens[1]);
						count = Integer.valueOf(tokens[2]);
						sum += count;
						numNonZeros ++;
					}
				}
			}	                     
			System.out.println("Total sum: " + sum);
			System.out.println("Non-zero avg: " + (sum/numNonZeros));
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
}
