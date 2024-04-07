package main;

import java.util.Random;

import main.Run.ErrorProfile;
import main.Run.Mechanism;

import distributions.Geometric;

public class RunGeometric {
	
	static int inputSize = 400000000; //8377*52292;
	static double epsilon = 0.1;
	static double sensitivity = 1;
	static double nonZeroMean = 139043;
	
	public static void main(String[] args) {
		
		int numRepeats = 5;
		int noiseSeed = 35791;
		Random nextSeedRand = new Random(266771);	
		
		double[] data = new double[inputSize];
		
		double alpha = Math.exp(-epsilon/sensitivity);
		Geometric.setAlpha(alpha);
		
		//double[] rangePercents = {1.0/inputSize, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
		int numDigits = (int) Math.log10(inputSize);
		int numAdded = (int) (inputSize/(Math.pow(10, numDigits)));
		double[] rangePercents = new double[2*numDigits+1+numAdded];
		rangePercents[0] = 1.0/inputSize; int index =1;
		for (int i=0; i<numDigits; i++) {
			rangePercents[index] = 5*rangePercents[index-1]; index++;
			rangePercents[index] = 2*rangePercents[index-1]; index++;
		}
		double base = rangePercents[index-2];
		for (int i=0; i<numAdded; i++) {
			rangePercents[index] = rangePercents[index-1]+base; index++;
		}
		
		
		ErrorProfile[] rangeXZProfiles = new ErrorProfile[rangePercents.length];
		int maxNumQueries = 1000000;
		for (int i=0; i<rangePercents.length; i++) {
			int rangeSize = (int) Math.max(Math.round(rangePercents[i]*inputSize), 1);
			int numQueries = Math.min(maxNumQueries, inputSize-rangeSize+1); 
			rangeXZProfiles[i] = new ErrorProfile(numQueries, true);
			rangeXZProfiles[i].rangeSize = rangeSize;
		}
		
		double[] absErrors = new double[rangePercents.length];
		double[] relErrors = new double[rangePercents.length];
		
		for (int k=0; k<numRepeats; k++) {
			System.out.println(k);
			Random rand = new Random(noiseSeed + nextSeedRand.nextInt());
			for (int i=0; i<inputSize; i++) {
				data[i] = Geometric.getSample(rand);
			}
					
		
			for (int i=0; i<rangePercents.length; i++) {
				//computeRangeQueryError(data, output, rangePercents[k], 1, rangeXZProfiles[k], errorWriter);
				
				computeRangeQueryError(data, rangeXZProfiles[i].rangeSize, rangeXZProfiles[i]);
				absErrors[i] += (rangeXZProfiles[i].absError/rangeXZProfiles[i].numQueries);
				relErrors[i] += (rangeXZProfiles[i].getMedianRelError());
				rangeXZProfiles[i].resetError();
			}
		}
		
		for (int i=0; i<rangePercents.length; i++) {
			int rangeSize = (int) Math.max(Math.round(rangePercents[i]*inputSize), 1);
			System.out.println(rangeSize + ", " + absErrors[i]/numRepeats + ", " + relErrors[i]/numRepeats);
		}
			
	}
	
	public static void computeRangeQueryError(double[] input,  int rangeLength, ErrorProfile ep) {
		
		double subsetSum = 0;
		for (int i=0; i<rangeLength; i++) {
			subsetSum += input[i];	
		}
		
		double absError = Math.abs(subsetSum);
		
		ep.addError(absError, subsetSum + nonZeroMean*rangeLength);		
//		if (pw != null && writeToFile)
//			pw.println(trueSubsetSum + ", " + estSubsetSum);
		
		int numQueries = ep.numQueries; //Math.min(input.length-rangeLength, MAX_NUM_QUERIES);
		for (int i=0; i<numQueries-1; i++) {
			subsetSum += (input[i+rangeLength] - input[i]);
			
			absError = Math.abs(subsetSum);
			ep.addError(absError, subsetSum + nonZeroMean*rangeLength);
//			if (pw != null && writeToFile)
//				pw.println(trueSubsetSum + ", " + estSubsetSum);	
		
		}	
	}
}
