package main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;

import wavelets.sparse_1d_wavelet;
import wavelets.sparse_2d_wavelet;
import wavelets.sparse_2d_wavelet_array;

import datahandler.census.CensusTracts;
import datastructures.DyadicRange_Full;
import datastructures.DyadicRange_Sparse;
import datastructures.DyadicRange_Sparse_1d;
import datastructures.DyadicRange_Sparse_2d;
import datastructures.HaarWavelet;
import datastructures.ObjectPool;
import datastructures.PrioritySampling;
import datastructures.Sampling;
import datastructures.Filtering_Absolute;
import datastructures.FilterPrioritySampling;
import datastructures.Filtering_Signed;
import datastructures.ThresholdSampling;
import datastructures.UtilityObject.*;
import distributions.Geometric;
import distributions.Laplace;

public class Run {	
	public enum Mechanism {Filter1, Filter2, Threshold, Priority, Geometric, 
		Uniform, Filter_Priority, Test_TopDown, Wavelet1, Wavelet2};

	static boolean scaledRanges = false; // scaled by sum(M')/sum(M'')
	static int sensitivity = 1;
	static double epsilon = 0.1;
	static Run.Mechanism mechanism = Run.Mechanism.Filter_Priority; //Filter2; //Wavelet2;  //Filter_Priority;
	static boolean useDyadic = false;
	
	static long inputSize = (int) 1E6; //(long) Math.pow(2,30);//1E7; //m
	static int numNonZeros = (int) 1E5; //n
	static boolean fixedInputSize = false; //true;
	static double zeroPercent = (inputSize-numNonZeros)*1.0/inputSize; //0.9;
	static double nonZeroMean = 100;
	static double nonZeroStdDev = 20;
	static boolean measuringTime = true;
	static boolean computeError = true;
	
	static int totalNonZerosSampled = 0;
	static double sampleConst = 1;
	static int numRepeats = 20;
	static boolean skewedData = false;
	static int theta = 40; //500; // 70; //40;
	static double rangePercent = 0;
	static boolean useGivenThreshold = true; // for filtering, use the user-specified theta directly. 
	static DyadicRange_Sparse dyadicRange = null;
	static boolean useMedianRelativeError = true;
//	static boolean writeToFile = false;
	static boolean testingMode = false; // use materialization method
	static boolean useZeroIndices = false; // zero indices are stored as a list {[start, end]}
	
	static boolean useRealData = false;
	static boolean webTrace = false;
	static boolean oneDimensional = true; //false;
	static String fileName = "data/income_census/out/grid1.txt";
	//static String fileName = "data/internet_traffic/out/combinedTrace.txt";
	// for web trace
	//static int xSize = 8377; static int ySize = 52292;
	static int xSize = 8000; static int ySize = 50000;
	static int MULTIPLE = ySize; // use to converse 2d coeffs to 1d, this is the cardinality of the 2nd attribute
	static int MAX_NUM_QUERIES = 100000; // maximum number of queries evaluated for a range percentage
	
	public static void main(String[] args) {			
		int genDataSeed = 12345; int noiseSeed = 35791; int sampleSeed=73937; int binomialSeed=92381;
		Random nextSeedRand = new Random(266771);	
		
		parseCommandLineOptions(args, 0);
					
		long totalAnonTime = 0;
		
		double[] data = null; 
		double[] output = null; 
		HashMap<Integer, Double> coeffMap = null;
		
		if (fixedInputSize)
			numNonZeros = (int) Math.round(inputSize*(1.0-zeroPercent));
		else
			inputSize = (long) Math.round(numNonZeros/(1.0-zeroPercent));
		
		ArrayList<NonZeroEntry> inputNonZeros = new ArrayList<NonZeroEntry>();
		
		// output is described by the non-zeros and zero indices
		ArrayList<ZeroIndexInterval> inputZeroIndices = null;
		if (useZeroIndices)
			inputZeroIndices = new ArrayList<ZeroIndexInterval>();
			
		ArrayList<NonZeroEntry> outputNonZeros = new ArrayList<NonZeroEntry>();
		
		int totalSampleSize = 0;
		
		if (useRealData) {
			inputSize = readFreqMatrixData(fileName, inputNonZeros);
			zeroPercent = (inputSize-inputNonZeros.size())*1.0/inputSize;
			System.out.println("input size=" + inputSize);
			// FOR WEB TRACE
			if (webTrace) {
				inputSize = xSize*ySize;
				System.out.println("n=" + inputNonZeros.size() + ", m=" + inputSize);
				zeroPercent = 1.0*(inputSize-inputNonZeros.size())/inputSize;
			}
		}
		
		//double[] rangePercents = {1.0/inputSize, 1E-5, 5E-5, 1E-4, 5E-5, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
		double[] rangePercents = {1.0/inputSize, 0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
		
		if (useRealData) {
			int numDigits = (int) Math.log10(inputSize);
			int numAdded = (int) (inputSize/(Math.pow(10, numDigits))) -1;
			rangePercents = new double[2*numDigits+1+numAdded];
			rangePercents[0] = 1.0/inputSize; int index =1;
			for (int i=0; i<numDigits; i++) {
				rangePercents[index] = 5*rangePercents[index-1]; index++;
				rangePercents[index] = 2*rangePercents[index-1]; index++;
			}
			double base = rangePercents[index-2];
			for (int i=0; i<numAdded; i++) {
				rangePercents[index] = rangePercents[index-1]+base; index++;
			//	System.out.println(inputSize*rangePercents[index-1]);
			}
		}
		if (rangePercent > 0) { // use a specified range percentage
			rangePercents = new double[1];
			rangePercents[0] = rangePercent;
		}
		
		int targetSampleSize = (int) Math.round(sampleConst*inputSize*(1-zeroPercent));
		
		if (!useDyadic) {
			if  (mechanism == Mechanism.Filter1 || mechanism == Mechanism.Filter2) {
				if (!useGivenThreshold) {
					double alpha = Math.exp(-epsilon/sensitivity);							
					theta = computeFilteringTheta(targetSampleSize, mechanism, alpha, (int)(inputSize*zeroPercent), (int)(inputSize*(1-zeroPercent)), nonZeroMean);
					System.out.println("theta=" + theta);
				}
			}		
		}
//System.out.println(theta);

		ErrorProfile[] rangeXZProfiles = new ErrorProfile[rangePercents.length];
		
		if (computeError) {
			for (int i=0; i<rangePercents.length; i++) {
				int rangeLength = (int) Math.round(rangePercents[i]*inputSize);	
				int numQueries;
				if (inputSize - rangeLength + 1 < MAX_NUM_QUERIES)
					numQueries = (int) inputSize - rangeLength +1;
				else
					numQueries = MAX_NUM_QUERIES;
				
				rangeXZProfiles[i] = new ErrorProfile(numQueries, useMedianRelativeError);		
			}
		}
		
		totalNonZerosSampled = 0;
		double[] absoluteErrors = new double[rangePercents.length];
		double[] relativeErrors = new double[rangePercents.length]; 
		
		int intInputSize;
		if (inputSize < Integer.MAX_VALUE)
			intInputSize = (int) inputSize;
		else
			intInputSize = Integer.MAX_VALUE;
		
		int logm = (int) (Math.log(inputSize)/Math.log(2));
		if (mechanism == Mechanism.Wavelet1 || mechanism == Mechanism.Wavelet2) {
			if (inputSize > (long) 1 << logm) {
				logm += 1; // need to add zero entries to the end of the array
				inputSize = (long) 1 << logm;
			}
		}
		
		double alpha = Math.exp(-epsilon/sensitivity);
		Geometric.setAlpha(alpha);
		
		if (!useRealData && (mechanism == Mechanism.Geometric && computeError)) {
			data = new double[(int) inputSize];
			output = new double[(int) inputSize];
		}
		
		for (int i=0; i<numRepeats; i++) {
			System.out.println(i);
			if (!useRealData) {
				if ((mechanism == Mechanism.Geometric && computeError) || mechanism == Mechanism.Wavelet1) {
					createSynData(data, intInputSize, zeroPercent, nonZeroMean, nonZeroStdDev, true, skewedData, genDataSeed + nextSeedRand.nextInt());
				}
				else {
					createSynData_SparseForm(inputNonZeros, inputZeroIndices, numNonZeros, inputSize, nonZeroMean, nonZeroStdDev, genDataSeed + nextSeedRand.nextInt());
				}
			}		

			if ( !useDyadic || mechanism == Mechanism.Wavelet1 || mechanism == Mechanism.Wavelet2) {
				long start = System.nanoTime();
				int sampleSize=0;
				
				
				if (mechanism == Mechanism.Geometric) {
					if (computeError)
						addGeometricNoise(data, output, alpha, noiseSeed + nextSeedRand.nextInt());
					else
						generateGeometricNoise(inputSize, alpha, noiseSeed + nextSeedRand.nextInt());
				}
				else if (mechanism == Mechanism.Wavelet1) {
					HaarWavelet wavelet = new HaarWavelet(data, noiseSeed + nextSeedRand.nextInt());
					double scale = sensitivity/epsilon; // sensitivity = 1
					wavelet.addLaplacian(scale, true);			
				}
				else if (mechanism == Mechanism.Wavelet2) {
					if (oneDimensional) {
						sparse_1d_wavelet wavelet;
						if (logm <= 29) {
							if (data == null)
								data = new double[(int) inputSize];						
							wavelet = new sparse_1d_wavelet(logm, data);
						}
						else {
							if (coeffMap == null) 
								coeffMap = new HashMap<Integer, Double>(100000000);
							else 
								coeffMap.clear();
							
							wavelet = new sparse_1d_wavelet(logm, coeffMap,  false);
						}
						for (NonZeroEntry entry: inputNonZeros) {
							wavelet.update(entry.getPosition(), entry.getValue());
						}
						double scale = (1+logm)*sensitivity/epsilon;
						Random rand = new Random(noiseSeed + nextSeedRand.nextInt());
	
						if (logm <= 29) {
							wavelet.addLaplacian(scale, rand);
							if (computeError && output == null) 
								output = new double[1 << logm];
							wavelet.inverse(output);
						}
						else {
							wavelet.addLaplacian_Time(scale, rand);
							wavelet.inverse_Time();
						}
					}
					else { // 2d wavelets
						sparse_2d_wavelet_array wavelet;
						int logx = (int) (Math.log(xSize)/Math.log(2)) + 1;
						int logy = (int) (Math.log(ySize)/Math.log(2)) + 1;
						System.out.println("logx=" + logx + ", logy=" + logy);
						wavelet = new sparse_2d_wavelet_array(logx, logy);
						for (NonZeroEntry entry : inputNonZeros) {
							int pos = entry.getPosition();
							int x = pos / MULTIPLE;
							int y = pos % MULTIPLE;
							wavelet.update(x, y, entry.getValue());
						}
						
						System.out.println("Done updating");
						double scale = (1+logx)*(1+logy)*sensitivity/epsilon;
						Random rand = new Random(noiseSeed + nextSeedRand.nextInt());
						wavelet.addLaplacian(scale, rand);
						System.out.println("Done adding noise");
						if (output == null) {
							System.out.println("input size=" + (int) inputSize);
							//output = new double[(int) inputSize];
							output = new double[(int) xSize*ySize];
						}
						
						try {
							//PrintWriter outPW = new PrintWriter(new FileWriter("data/internet_traffic/out/out_wavelet.txt", false));
							for (int x=0; x<xSize; x++) {
								for (int y=0; y<ySize; y++) {
									int pos = x*MULTIPLE + y;
									double val = wavelet.inverse(x, y);
									output[pos] = val;
							//		outPW.println(x + "," + y + "," + val);
								}
							}
							//outPW.close();
						}
						catch (Exception e) {
							e.printStackTrace();
						}
						System.out.println("Done inversing");
					}
				}
				else {							
					if (testingMode) {
						data = new double[(int)inputSize];
						for (NonZeroEntry entry: inputNonZeros)
							data[entry.getPosition()] = entry.getValue();
						output = new double[(int)inputSize];
						//sampleSize = anonymize_basic(data, output, mechanism, theta, targetSampleSize, epsilon, noiseSeed + 51*i, sampleSeed + 101*i);
						sampleSize = anonymize_basic(data, output, mechanism, theta, targetSampleSize, epsilon, noiseSeed + nextSeedRand.nextInt(), sampleSeed + nextSeedRand.nextInt());
						for (int j=0; j<output.length; j++) {
							if (output[j] != 0 && data[j] !=0)
								totalNonZerosSampled ++;
						}
					}
					else {
						sampleSize = anonymize(inputNonZeros, inputZeroIndices, (int) inputSize, outputNonZeros, mechanism, theta, targetSampleSize, epsilon, noiseSeed + nextSeedRand.nextInt(), sampleSeed + nextSeedRand.nextInt(), binomialSeed + nextSeedRand.nextInt());	
					}					
				}
				
				long end = System.nanoTime();
				totalAnonTime += (end-start);
											
				totalSampleSize += sampleSize;	 
				
				if (computeError) {
					if (mechanism != Mechanism.Geometric && ! testingMode) {
						PositionComparator comparator = new PositionComparator();
						Collections.sort(inputNonZeros, comparator);	
						Collections.sort(outputNonZeros, comparator);
					}
//					if (mechanism == Mechanism.Wavelet2) {
//						data = new double[output.length];
//						for (NonZeroEntry entry : inputNonZeros)
//							data[entry.getPosition()] = entry.getValue();
//					}
					
					for (int k=0; k<rangePercents.length; k++) {
						//computeRangeQueryError(data, output, rangePercents[k], 1, rangeXZProfiles[k], errorWriter);
						int rangeSize = (int) Math.max(Math.round(rangePercents[k]*inputSize), 1);
						
//						if (mechanism == Mechanism.Geometric) {
//							for (int l=0; l<rangeXZProfiles[k].numQueries; l++)
//								rangeXZProfiles[k].absError += getGeometricError(rangeSize, new Random(noiseSeed + nextSeedRand.nextInt()));
//						}
					// *******
					//	if (mechanism == Mechanism.Geometric || mechanism == Mechanism.Wavelet2 || testingMode) {
						
						rangeXZProfiles[k].rangeSize = rangeSize;
						
						if (mechanism == Mechanism.Geometric || testingMode) {
							computeRangeQueryError(data, output, rangeSize, rangeXZProfiles[k]);
						}
						else if (mechanism == Mechanism.Wavelet2) {
							computeRangeQueryError(inputNonZeros, output, rangeSize, rangeXZProfiles[k]);
						}
						else {
							computeRangeQueryError(inputNonZeros, outputNonZeros, rangeSize, rangeXZProfiles[k]);
						}
						if (useMedianRelativeError) {
							double relError = rangeXZProfiles[k].getMedianRelError();
//							System.out.println("rel error=" + relError);
							relativeErrors[k] += relError;
						}
						else {
							relativeErrors[k] += (rangeXZProfiles[k].absError/rangeXZProfiles[k].trueSum);
						}
						
						absoluteErrors[k] += rangeXZProfiles[k].absError;
						rangeXZProfiles[k].resetError();
//						System.out.println(k + ", " + absoluteErrors[k]);
					}
					// reset output for wavelets and geometric
					if (mechanism == Mechanism.Wavelet2 || mechanism == Mechanism.Geometric) {
						if (data != null) {
							for (int k=0; k<data.length; k++) {
								data[k] = 0;
							}
						}
						if (output != null) {
							for (int k=0; k<output.length; k++) {	
								output[k] = 0;
							}
						}
					}
				}
				if (!useRealData) {
					ObjectPool.returnNonZeroEntries(inputNonZeros);
					inputNonZeros.clear();
				}
				
				ObjectPool.returnNonZeroEntries(outputNonZeros);
				outputNonZeros.clear();
				if (useZeroIndices) {
					ObjectPool.returnZeroIndexIntevals(inputZeroIndices);
					inputZeroIndices.clear();
				}
			}
			else { // dyadicRange
				if (dyadicRange != null)
					dyadicRange.reset();
				else {
					if (oneDimensional)
						dyadicRange = new DyadicRange_Sparse_1d((int) inputSize);
					else 
						dyadicRange = new DyadicRange_Sparse_2d(xSize, ySize);	
				}
				
				long start=0;
				if (measuringTime)
					start = System.nanoTime();
				System.out.println("start updating");
				if (oneDimensional)
					((DyadicRange_Sparse_1d) dyadicRange).update(inputNonZeros);
				else
					((DyadicRange_Sparse_2d) dyadicRange).update(inputNonZeros);
				System.out.println("done updating");
				int sampleSize = dyadicRange.anonymize(mechanism, theta, targetSampleSize, epsilon, noiseSeed + nextSeedRand.nextInt(), sampleSeed + nextSeedRand.nextInt(), binomialSeed + nextSeedRand.nextInt());
				System.out.println("done anonymized");
				if (measuringTime) {
					long end = System.nanoTime();
					totalAnonTime += (end-start);
					System.out.println(i + "," + (end-start));
				}
				
				totalSampleSize += sampleSize;
				
				if (computeError) {
					int oneDimSize;
					if (oneDimensional) oneDimSize = (int) inputSize;
					else oneDimSize = xSize;
						
					for (int k=0; k<rangePercents.length; k++) {
						int rangeSize = (int) Math.max(Math.round(rangePercents[k]*oneDimSize),1);
						rangeXZProfiles[k].rangeSize = rangeSize;
						computeRangeQueryError_Dyadic(dyadicRange, rangeSize, rangeXZProfiles[k]);
						
						if (useMedianRelativeError) {
							relativeErrors[k] += rangeXZProfiles[k].getMedianRelError();
						}
						else {
							relativeErrors[k] += (rangeXZProfiles[k].absError/rangeXZProfiles[k].trueSum);
						}
						absoluteErrors[k] += rangeXZProfiles[k].absError;
						rangeXZProfiles[k].resetError();
						
					}
				}
			}
		}
	
		int sampleSize = totalSampleSize/numRepeats;
		if (measuringTime) {
			String prefix = "";
			if (useDyadic) prefix = "D-";
			if (mechanism != Mechanism.Filter_Priority)
				System.out.println(prefix + mechanism + ", m=" + inputSize + ", sampleSize=" + sampleSize + ", ZeroPercent=" + zeroPercent + ", Time(ms)= " + totalAnonTime/1E6/(numRepeats)); // + ", NonZerosRetained=" + totalNonZerosSampled/numRepeats);
			else
				System.out.println(prefix + "FilterPriority_" + theta + ", m=" + inputSize +  ", sampleSize=" + sampleSize + ", ZeroPercent=" + zeroPercent + ", Time(ms)=" + totalAnonTime/1E6/(numRepeats)); // + ", NonZerosRetained=" + totalNonZerosSampled/numRepeats);
		}
		
		if (computeError) {
			for (int k=0; k<rangePercents.length; k++) {
				double avgAbsError = absoluteErrors[k]/numRepeats/rangeXZProfiles[k].numQueries;
				double avgRelError = relativeErrors[k]/numRepeats;
				System.out.println("Range=" + rangeXZProfiles[k].rangeSize + ", SampleSize=" + totalSampleSize/numRepeats + ", AbsError=" 
						+ avgAbsError + ", RelError=" + avgRelError);
				//		+ ", " + totalNonZerosSampled/numRepeats + ", DyadicTopDown: " + totalNodesTouched/numRepeats);
			}
		}
	}
	
	// the data is shuffled, non-zeros are placed uniformly (i.e., shuffled) or placed near one end of the array 
	public static void createSynData(double[] counts, int inputSize, double zeroPercent, double mean, double stdDev, boolean isShuffled, boolean isSkewed, int seed) {
		int numNonZeros = (int) Math.round(inputSize*(1.0-zeroPercent));
		Random rand = new Random(seed);
		Random rand2 = new Random(8879 + 3*seed);
		
		for (int i=0; i<numNonZeros; i++) {
			int sample = (int) Math.round(rand.nextGaussian()*stdDev + mean);
			if (sample < 1)
				sample = 1;
			counts[i] = sample;
			
			if (isSkewed) {
				int position = i + (int) Math.abs((rand2.nextGaussian()*numNonZeros/5));
				while (true) {
					if (position > inputSize) position = inputSize - i;	
					if (counts[position] ==0) {
						counts[position] = sample;
						break;
					}
					else {
						position = i + (int) Math.abs((rand2.nextGaussian()*numNonZeros/5));
					}
				}
			}
		}
		
		if (!isSkewed && isShuffled) {
			for (int i=0; i<numNonZeros-1; i++) {
				int nextPos = i + 1 + rand2.nextInt(inputSize-i-1);
				double temp = counts[nextPos];
				counts[nextPos] = counts[i];
				counts[i] = temp;
			}
		}
	}
	
	
	// the data is shuffled, non-zeros are distributed uniformly
	// output ArrayList<NonZeroEntry> nonZeros and int[] array of zero indices
	public static void createSynData_SparseForm(ArrayList<NonZeroEntry> nonZeros, ArrayList<ZeroIndexInterval> zeroIndices, int numNonZeros,  long inputSize, double mean, double stdDev, int seed) {
		Random rand = new Random(seed);
		Random rand2 = new Random(128879 + seed*3);
		
		TreeSet<Integer> nonZeroIndices = new TreeSet<Integer>();
		
		int upperInt = (int) inputSize;
		if (inputSize > Math.pow(2, 31))
			upperInt = Integer.MAX_VALUE;
		
		for (int i=0; i<numNonZeros; i++) {
			int sample = (int) Math.round(rand.nextGaussian()*stdDev + mean);
			if (sample < mean - 3*stdDev) sample = (int) (mean - 3*stdDev);
			else if (sample > mean + 3*stdDev) sample = (int) (mean + 3*stdDev);
			//	sample = 1;
			
			int position = rand2.nextInt(upperInt);
			while (nonZeroIndices.contains(position))
				position = rand2.nextInt(upperInt);
			nonZeroIndices.add(position);
			
			NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
			entry.setPosition(position);
			entry.setValue(sample);
			nonZeros.add(entry);	
		}
		
		if (useZeroIndices) {
			int start = 0; 
			Iterator<Integer> nonZeroIter = nonZeroIndices.iterator();
			while (nonZeroIter.hasNext()) {
				int end = nonZeroIter.next()-1;
				if (end < start) {
					start = end + 2;
					continue;
				}
				
				ZeroIndexInterval zeroIndex = ObjectPool.allocateZeroIndexInterval();
				zeroIndex.setStartIndex(start);
				zeroIndex.setEndIndex(end); // end index is inclusive
				zeroIndices.add(zeroIndex);
				start = end + 2;
			}
			
			if (start <= inputSize - 1) {
				ZeroIndexInterval zeroIndex = ObjectPool.allocateZeroIndexInterval();
				zeroIndex.setStartIndex(start);
				zeroIndex.setEndIndex((int) inputSize-1); 
				zeroIndices.add(zeroIndex);
			}
			
			int numZeros = (int) inputSize - numNonZeros;
			int totalZerosSofar = 0;
			for (ZeroIndexInterval interval :  zeroIndices) {
				totalZerosSofar += (interval.getEndIndex() - interval.getStartIndex() + 1);
				double prob = 1.0*totalZerosSofar/numZeros;
				interval.setCumProb(prob);
			}
			
		}
	}

	

	public static void addLaplaceNoise(double[] input, double[] output, double epsilon, int seed) {
		double scale = sensitivity/epsilon;
		Random rand = new Random(seed);
		
		for (int i=0; i<input.length; i++) {
			double noise = Laplace.getSample(scale, rand);
			output[i] = input[i] + noise;		
		}
	}
	
	public static void addGeometricNoise(double[] input, double[] output, double alpha, int seed) {
		Random rand = new Random(seed);
		for (int i=0; i<input.length; i++) {
			// noise can be positive or negative
			int noise = Geometric.getSample(rand); 
			output[i] = input[i] + noise;
			//System.out.println("In=" + input[i] + ", Output=" + output[i]);
		}
	}

	public static void generateGeometricNoise(long num, double alpha, int seed) {
		Random rand = new Random(seed);
		for (long i=0; i<num; i++)
			Geometric.getSample(rand); 
	}
	
	public static double getGeometricError(int rangeSize, Random rand) {
		double sum = 0;
		for (int i=0; i<rangeSize; i++) {
			sum += Geometric.getSample(rand);
		}
		return Math.abs(sum);
	}

	public static int anonymize(ArrayList<NonZeroEntry> nonZeros, ArrayList<ZeroIndexInterval> zeroIndices, int inputSize, ArrayList<NonZeroEntry> output, Mechanism mechanism, int theta, int K, double epsilon, int noiseSeed, int sampleSeed, int zeroPosSeed)  {	
		double alpha = Math.exp(-epsilon/sensitivity);
		
		int sampleSize;
		Sampling sampling=null;
		if (mechanism == Mechanism.Filter1) {			
			sampling = new Filtering_Signed(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);		
		}
		else if (mechanism == Mechanism.Filter2) {	
			sampling = new Filtering_Absolute(alpha, theta, noiseSeed, sampleSeed, zeroPosSeed);
		}	
		else if (mechanism == Mechanism.Threshold) {
			double absValueSum = 0;
			for (int i=0; i<nonZeros.size(); i++)
				absValueSum += Math.abs(nonZeros.get(i).getValue());
			
			absValueSum += 2*alpha/(1-alpha*alpha)*(inputSize - nonZeros.size());
			int D = (int) (absValueSum*(1.0)/K);
			System.out.println("Threshold=" + D);
			
			sampling = new ThresholdSampling(alpha, D, noiseSeed, sampleSeed, zeroPosSeed);		
		}				
		else if (mechanism == Mechanism.Priority) {
			sampling = new PrioritySampling(alpha, K, noiseSeed, sampleSeed, zeroPosSeed);			
		}
		else if (mechanism == Mechanism.Filter_Priority) {
			sampling = new FilterPrioritySampling(alpha, theta, K, noiseSeed, sampleSeed, zeroPosSeed);		
		}
		
		sampleSize = sampling.anonymize(nonZeros, zeroIndices, inputSize, output);
		return sampleSize;
	}

	
	public static int anonymize_basic(double[] input, double[] output, Mechanism mechanism, int theta, int K, double epsilon, int noiseSeed, int sampleSeed)  {	
		if (mechanism == Mechanism.Uniform) {
			int inputSize = input.length;
			int sampleSize = (int) K;
			int[] indices = new int[inputSize];
			for (int i=0; i< inputSize; i++)
				indices[i] = i;
			// shuffling
			Random rand = new Random(sampleSeed);
			for (int i=0; i< inputSize; i++) {
			    int randomPosition = rand.nextInt(inputSize);
			    double temp = input[i];
			    input[i] = input[randomPosition];
			    input[randomPosition] = temp;
			}
			Random rand2 = new Random(noiseSeed);
			double scale = sensitivity/epsilon;
			for (int i=0; i<sampleSize; i++) {
				output[indices[i]] = (int) Math.round(input[indices[i]] + Laplace.getSample(scale, rand2));  
			}
			return sampleSize;	
		}
		else { 
			double alpha = Math.exp(-epsilon/sensitivity);
			int sampleSize = 0;
			if (mechanism == Mechanism.Geometric) {
				addGeometricNoise(input, output, alpha, noiseSeed);
				return input.length;
			}
			else if (mechanism == Mechanism.Filter1) {			
				//System.out.println("param=" + param);
				//theta = solveForTheta(alpha, stats, param);
			
				Filtering_Signed ts = new Filtering_Signed(alpha, theta, noiseSeed, sampleSeed, 0);
				sampleSize = ts.anonymize_basic(input, output);
			}
			else if (mechanism == Mechanism.Filter2) {	
				Filtering_Absolute ts = new Filtering_Absolute(alpha, theta, noiseSeed, sampleSeed, 0);
				sampleSize = ts.anonymize_basic(input, output);
			}
			
			else if (mechanism == Mechanism.Threshold) {
				double absValueSum = 0;
				for (int i=0; i<input.length; i++)
					absValueSum += Math.abs(input[i]);
				
				absValueSum += 2*alpha/(1-alpha*alpha)*input.length*zeroPercent;
				int D = (int) (absValueSum*(1.0)/K);
				// System.out.println("D=" + D);
				Sampling wrs = new ThresholdSampling(alpha, D, noiseSeed, sampleSeed, 0);
				sampleSize = wrs.anonymize_basic(input, output);
			}				
			else if (mechanism == Mechanism.Priority) {
				// set K to be the same as the number of non-zero entries
				// int K = (int) Math.round(counts.length/(stats.zeroNonZeroRatio + 1));	
				PrioritySampling ps = new PrioritySampling(alpha, K, noiseSeed, sampleSeed, 0);				
				sampleSize = ps.anonymize_basic(input, output);	
			}
			else if (mechanism == Mechanism.Filter_Priority) {
				// set K to be the same as the number of non-zero entries
				// int K = (int) Math.round(counts.length/(stats.zeroNonZeroRatio + 1));			
				FilterPrioritySampling ps = new FilterPrioritySampling(alpha, theta, K, noiseSeed, sampleSeed, 0);	
				sampleSize = ps.anonymize_basic(input, output);
			}
			return sampleSize;
		}
	}
	
//	public static DyadicRange_Sparse_1d anonymizeDyadicRanges(ArrayList<NonZeroEntry> nonZeros, int inputSize, Mechanism mechanism, int theta, int sampleSize, double epsilon, int noiseSeed, int sampleSeed, int zeroPosSeed)  {
//		DyadicRange_Sparse_1d dr = new DyadicRange_Sparse_1d(nonZeros, inputSize);
//		
//		dr.anonymize(mechanism, theta, sampleSize, epsilon, noiseSeed, sampleSeed, zeroPosSeed);
//		
//		return dr;
//	}
	
	
	public static void computeRangeQueryError(double[] input, double[] output, int rangeLength, ErrorProfile ep) {
		
		double sampleSubsetSum = 0;
		double trueSubsetSum = 0;
		for (int i=0; i<rangeLength; i++) {
			sampleSubsetSum += output[i];
			trueSubsetSum += input[i];	
		}
		
		double estSubsetSum;
		estSubsetSum = sampleSubsetSum;
			
		double absError = Math.abs(estSubsetSum - trueSubsetSum);
		
		ep.addError(absError, trueSubsetSum);		
//		if (pw != null && writeToFile)
//			pw.println(trueSubsetSum + ", " + estSubsetSum);
		
		int numQueries = ep.numQueries; //Math.min(input.length-rangeLength, MAX_NUM_QUERIES);
		for (int i=0; i<numQueries-1; i++) {
			sampleSubsetSum += (output[i+rangeLength] - output[i]);
			trueSubsetSum += (input[i+rangeLength] - input[i]);
			
			estSubsetSum = sampleSubsetSum;
			
			absError = Math.abs(estSubsetSum - trueSubsetSum);
			ep.addError(absError, trueSubsetSum);
//			if (pw != null && writeToFile)
//				pw.println(trueSubsetSum + ", " + estSubsetSum);	
		
		}	
	}
	
	public static void computeRangeQueryError(ArrayList<NonZeroEntry> input, ArrayList<NonZeroEntry> output, int rangeSize, ErrorProfile ep) {
		int startInputIndex = 0;
		int endInputIndex = 0;
		for (; endInputIndex < input.size(); endInputIndex++) {
			if (input.get(endInputIndex).getPosition() >= rangeSize) {
				break;
			}
		}
		endInputIndex--;
		
		int startOutputIndex = 0;
		int endOutputIndex = 0;
		for (; endOutputIndex < output.size(); endOutputIndex++) {
			if (output.get(endOutputIndex).getPosition() >= rangeSize) break;
		}
		endOutputIndex--;
		
		double trueSubsetSum = 0;
		for (int i=startInputIndex; i<=endInputIndex; i++)
			trueSubsetSum += input.get(i).getValue();
		
		double estSubsetSum = 0;
		for (int i=startOutputIndex; i<=endOutputIndex; i++)
			estSubsetSum += output.get(i).getValue();
			
		double absError = Math.abs(estSubsetSum - trueSubsetSum);
		ep.addError(absError, trueSubsetSum);
		
		int numQueries = ep.numQueries;
		for (int i=1; i<numQueries; i++) {
			if (startInputIndex<=input.size()-1 && input.get(startInputIndex).getPosition() < i) {
				trueSubsetSum -= input.get(startInputIndex).getValue();
				startInputIndex++;
			}
			if (endInputIndex<input.size()-1 && input.get(endInputIndex+1).getPosition() < i+rangeSize) {
				endInputIndex++;
				trueSubsetSum += input.get(endInputIndex).getValue();
			}
			if (startOutputIndex<=output.size()-1 && output.get(startOutputIndex).getPosition() < i) {
				estSubsetSum -= output.get(startOutputIndex).getValue();
				startOutputIndex++;
			}
			if (endOutputIndex<output.size()-1 && output.get(endOutputIndex+1).getPosition() < i+rangeSize) {
				endOutputIndex++;
				estSubsetSum += output.get(endOutputIndex).getValue();
			}
			
			absError = Math.abs(estSubsetSum - trueSubsetSum);
			ep.addError(absError, trueSubsetSum);
		}
	}
	
	public static void computeRangeQueryError(ArrayList<NonZeroEntry> input, double[] output, int rangeSize, ErrorProfile ep) {
		int startInputIndex = 0;
		int endInputIndex = 0;
		for (; endInputIndex < input.size(); endInputIndex++) {
			if (input.get(endInputIndex).getPosition() >= rangeSize) {
				break;
			}
		}
		endInputIndex--;
		
		
		
		double trueSubsetSum = 0;
		for (int i=startInputIndex; i<=endInputIndex; i++)
			trueSubsetSum += input.get(i).getValue();
		
		double estSubsetSum = 0;
		for (int i=0; i<rangeSize; i++)
			estSubsetSum += output[i];
			
		double absError = Math.abs(estSubsetSum - trueSubsetSum);
		ep.addError(absError, trueSubsetSum);
		
		int numQueries = ep.numQueries;
		for (int i=1; i<numQueries; i++) {
			if (startInputIndex<=input.size()-1 && input.get(startInputIndex).getPosition() < i) {
				trueSubsetSum -= input.get(startInputIndex).getValue();
				startInputIndex++;
			}
			if (endInputIndex<input.size()-1 && input.get(endInputIndex+1).getPosition() < i+rangeSize) {
				endInputIndex++;
				trueSubsetSum += input.get(endInputIndex).getValue();
			}
			
			estSubsetSum += (output[i+rangeSize-1] - output[i-1]);
			
			absError = Math.abs(estSubsetSum - trueSubsetSum);
			ep.addError(absError, trueSubsetSum);
		}
	}
	
	public static void computeRangeQueryError_Dyadic(DyadicRange_Sparse dyadic, int rangeLength, ErrorProfile ep) {
		//int orgInputSize = dyadic.getOrgInputSize();
		int numQueries = ep.numQueries;
		if (dyadic instanceof DyadicRange_Sparse_1d) {
			DyadicRange_Sparse_1d dyadic1d = (DyadicRange_Sparse_1d) dyadic;
			for (int i=0; i<numQueries; i++) {
				double estSubsetSum = dyadic1d.rangeQuery(i, i+rangeLength-1, true);
				double trueSubsetSum = dyadic1d.rangeQuery(i, i+rangeLength-1, false);
	
				double absError = Math.abs(estSubsetSum - trueSubsetSum);
				ep.addError(absError, trueSubsetSum);
	//			if (pw != null && writeToFile)
	//				pw.println(trueSubsetSum + ", " + sampleSubsetSum);
				
			}
		}
		else {
			DyadicRange_Sparse_2d dyadic2d = (DyadicRange_Sparse_2d) dyadic;
			for (int i=0; i<numQueries; i++) {
				int endX;
				if (i+rangeLength > xSize)
					endX = xSize-1;
				else
					endX = i+rangeLength-1;
				double estSubsetSum = dyadic2d.rangeQuery(endX+1-rangeLength, endX, i, i+rangeLength-1, true);
				double trueSubsetSum = dyadic2d.rangeQuery(endX+1-rangeLength, endX, i, i+rangeLength-1, false);
	
				double absError = Math.abs(estSubsetSum - trueSubsetSum);
				ep.addError(absError, trueSubsetSum);
	//			if (pw != null && writeToFile)
	//				pw.println(trueSubsetSum + ", " + sampleSubsetSum);
				
			}
		}
	}
	
	
	public static int readFreqMatrixData(String fileName, ArrayList<NonZeroEntry> data) {
		int numEntries = 0;
		try {
			BufferedReader reader=null;
			reader = new BufferedReader(new FileReader(fileName));
			numEntries = Integer.valueOf(reader.readLine());
			
			while (true) {
				String line = reader.readLine();
				if (line == null)
					break;
				else {
//					int commaPos = line.indexOf(",");
//					int position = Integer.valueOf(line.substring(0, commaPos));
//					if (position > numEntries) {
//						System.err.println("The position of a non-zero has to be less than the data size, m");
//						System.exit(1);
//					}
//					int count = Integer.valueOf(line.substring(commaPos+1));
//					if (count <= 0) {
//						System.err.println("The value of a non-zero has to be greater than 0.");
//						System.exit(1);
//					}
					int position=0, count=0;
					String[] tokens = line.split(",");
					if (tokens.length == 2) {
						position = Integer.valueOf(tokens[0]);
						count = Integer.valueOf(tokens[1]);
					}
					else if (tokens.length == 3) {
						int x = Integer.valueOf(tokens[0]);
						int y = Integer.valueOf(tokens[1]);
						if (x >= xSize || y >= ySize)
							continue;
						position = x*MULTIPLE + y;
						count = Integer.valueOf(tokens[2]);
					}
					NonZeroEntry entry = ObjectPool.allocateNonZeroEntry();
					entry.setPosition(position);
					entry.setValue(count);
					data.add(entry);
				}
			}	                           
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		return numEntries;
		
	}
	
	// e.g., stateCode: 34, stateName: nj
	public static double[] readOnTheMapData1(String stateCode, String stateName, int numTractsDigits, int maxNumLines) {
		String[] fileNames = {"data/onthemap/al_mt_tracts.txt", "data/onthemap/ne_wy_tracts.txt"};
		HashMap<String, ArrayList<String>> tracts = CensusTracts.readLocationList(fileNames, numTractsDigits);
		System.out.println("Number of counties: " + tracts.size());
		
		HashMap<Long, Integer> freqMatrix = CensusTracts.readFreqMatrix_old(stateName, stateCode, tracts, maxNumLines);
		
		double[] counts = new double[freqMatrix.size()];
		Iterator<Integer> countIter = freqMatrix.values().iterator();
		int i=0;
		int totalNum = 0;
		int numZeroEntries = 0;
		while (countIter.hasNext()) {
			counts[i] = countIter.next();
			totalNum += counts[i];
			if (counts[i] == 0)
				numZeroEntries ++;
			i++;
		}
		
		System.out.println("Frequency matrix size: " + freqMatrix.size());
		System.out.println("Total num of people: " + totalNum);
		System.out.println("Num zero entries: " + numZeroEntries);
		System.out.println("Percentage of sparsity: " + ((double) numZeroEntries)/freqMatrix.size());
		System.out.println("Avg non-zero entry: " + ((double) totalNum)/(freqMatrix.size()-numZeroEntries));
		return counts;
	}
	
	static public int computeFilteringTheta(int sampleSize, Mechanism mechanism, double alpha, int numZeros, int numNonZeros, double nonZeroMean) {
		double l=  (1.0)*sampleSize/numNonZeros;;
		double ratio = (1.0)*numZeros/numNonZeros;
		int theta;
		
		if (mechanism == Mechanism.Filter1) {
			if ((ratio*Math.pow(alpha, nonZeroMean) + 1)/(1+alpha) >= l )
				theta = (int) (Math.log(l*(1+alpha)/(ratio + Math.pow(alpha, -nonZeroMean)))/Math.log(alpha));
			else {
				double delta = Math.pow((l-1)*(1+alpha), 2) + 4*ratio*Math.pow(alpha, nonZeroMean+1);
				theta = (int) (Math.log(((l-1)*(1+alpha) + Math.sqrt(delta))/(2*ratio))/Math.log(alpha));
			}
			return theta;
			
		}	
		else if (mechanism == Mechanism.Filter2) {	
			if (2*ratio*Math.pow(alpha,nonZeroMean) + Math.pow(alpha, 2*nonZeroMean) + 1 >= l*(1+alpha) ) {
				theta = (int) Math.round(Math.log(l*(1+alpha)/(2*ratio + Math.pow(alpha, -nonZeroMean+1) + 1))/Math.log(alpha));
			}
			else {
				double delta = Math.pow((l-1)*(1+alpha), 2) + 4*(2*ratio+Math.pow(alpha, nonZeroMean))*Math.pow(alpha, nonZeroMean+1);
				theta = (int) Math.round(Math.log(((l-1)*(1+alpha) + Math.sqrt(delta))/(2*(2*ratio+Math.pow(alpha, nonZeroMean))))/Math.log(alpha));
			}
			return theta;
		}
		return 0;
	}
	
	static public void parseCommandLineOptions(String[] args, int start) {
		for (int i=start; i<args.length; i++) {
			if (args[i].startsWith("--epsilon=")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				epsilon = Double.valueOf(arg);
//				System.out.println("Epsilon=" + epsilon);
			}
			else if (args[i].startsWith("--sensitivity=")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				sensitivity = Integer.valueOf(arg);
//				System.out.println("Sensitivity=" + sensitivity);
			} 
			else if (args[i].startsWith("--mechanism=")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				if (arg.equals("filter1"))
					mechanism = Run.Mechanism.Filter1;
				else if (arg.endsWith("filter2"))
					mechanism = Run.Mechanism.Filter2;
				else if (arg.equals("threshold"))
					mechanism = Run.Mechanism.Threshold;
				else if (arg.equals("priority"))
					mechanism = Run.Mechanism.Priority;
				else if (arg.equals("filter_priority"))
					mechanism = Run.Mechanism.Filter_Priority;
				else if (arg.equals("geometric") || arg.equals("laplace"))
					mechanism = Run.Mechanism.Geometric;
				else if (arg.equals("test_topdown"))
					mechanism = Run.Mechanism.Test_TopDown;
				else if (arg.equals("wavelet1")) 
                    mechanism = Run.Mechanism.Wavelet1;
				else if (arg.equals("wavelet2")) 
                    mechanism = Run.Mechanism.Wavelet2;

			}
			else if (args[i].startsWith("--materialize")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				testingMode = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--useDyadic")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				useDyadic = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--useRealData")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				useRealData = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--fileName")) {
				fileName = args[i].substring(args[i].indexOf("=")+1);
			}
			else if (args[i].startsWith("--m=")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				inputSize = Integer.valueOf(arg);
				fixedInputSize = true;
			}
			else if (args[i].startsWith("--n=")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				numNonZeros = Integer.valueOf(arg);
				fixedInputSize = false;
			}
			else if (args[i].startsWith("--mean")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				nonZeroMean = Double.valueOf(arg);
			}
			else if (args[i].startsWith("--stddev")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				nonZeroStdDev = Double.valueOf(arg);
			}
			else if (args[i].startsWith("--zeroPercent")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				zeroPercent = Double.valueOf(arg);
			}
			else if (args[i].startsWith("--measureTime")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				measuringTime = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--computeError")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				computeError = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--sampleConst")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				sampleConst = Double.valueOf(arg);
			}
			else if (args[i].startsWith("--numRepeats")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				numRepeats= Integer.valueOf(arg);
			}
			else if (args[i].startsWith("--skewedData")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				skewedData = Boolean.valueOf(arg);
			}
			else if (args[i].startsWith("--theta")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				theta = Integer.valueOf(arg);
			}
			else if (args[i].startsWith("--rangePercent")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				rangePercent= Double.valueOf(arg);
			}
			else if (args[i].startsWith("--useGivenThreshold")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				useGivenThreshold = Boolean.valueOf(arg);
			}
			// for consistency check
			else if (args[i].startsWith("--consistencyCheck")) {
				String arg = args[i].substring(args[i].indexOf("=")+1);
				DyadicRange_Full.consistencyCheck = Boolean.valueOf(arg);
			}
//			else if (args[i].startsWith("--errorFileName")) {
//				String arg = args[i].substring(args[i].indexOf("=")+1);
//				errorFileName = arg;
//			}
		
		}
	}

	public static class ErrorProfile {
		public double absError; // absolute error
		double trueSum; // to compute relative error
		//int numRelUpdates;
		//int numAbsUpdates;
		// for using median of relative errors
		double[] relErrors;
		int queryCount;
		public int numQueries;
//		int numUpdates = 0;
		boolean useMedianRelativeErrors;
		int rangeSize;
		
		public ErrorProfile() {
			this(0, false);// this is for not using median of relative errors
		}
		
		public ErrorProfile(int numQueries, boolean useMedian) {
			absError = 0;
			trueSum = 0;
			relErrors = new double[numQueries];
			queryCount = 0;
			this.numQueries = numQueries;
			useMedianRelativeErrors = true;
		}
	
		
		public void addError(double error, double trueVal) {
			absError += error;
			trueSum += trueVal;
			
			double relError;
			if (trueVal == 0 && error == 0) {
				relError = 1;
			}
			else if (trueVal == 0) {
				relError = Double.MAX_VALUE;
			}
			else {
				relError = Math.abs(error/trueVal);
			}
			if (useMedianRelativeErrors)
				relErrors[queryCount++] = relError;			
		}
		
		public double getMedianRelError() {
			Arrays.sort(relErrors);
			if (relErrors.length % 2 != 0)
				return relErrors[relErrors.length/2];
			else
				return (relErrors[relErrors.length/2-1] + relErrors[relErrors.length/2])/2;
		}
		
		public void resetError() {
			queryCount = 0;
			absError = 0;
			trueSum = 0;
			
		}
	}
	
	
}


