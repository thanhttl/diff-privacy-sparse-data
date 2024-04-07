package datahandler;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import datastructures.UtilityObject.NonZeroEntry;

public class CensusIncomeReader {
	public static String fileName = "data/census_income/usa.dat"; // num people = 11343120
	public static String outputFileName = "data/census_income/out/grid3.txt";
	
	public static int ageMultiple = 5; //5;
	public static int incomeMultiple = 5000; //5000;
	public static int birthplaceNumDigits = 5; //2; // maximum 5 digits
	public static int multiple = 20; // from 5% sample
	static long ageSetSize = Math.round(91/ageMultiple);
	static long incomeSetSize = 1010/(incomeMultiple/1000);
	static long birthplaceSetSize = 344;
	static long occupationSetSize = 506;
	
	public static void readCensusIncome() throws Exception {
		HashMap<Long,NonZeroEntry> dataMap = new HashMap<Long,NonZeroEntry>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		
		String line;

		HashSet<Integer> ageSet = new HashSet<Integer>();
		HashSet<Integer> birthplaceSet = new HashSet<Integer>();
		HashSet<Integer> occSet = new HashSet<Integer>();
		HashSet<Integer> incomeSet = new HashSet<Integer>();
		//int numPeople = 0;
		
		while ((line = reader.readLine()) != null) {
			int age = Integer.valueOf(line.substring(2,5));
			age = age/ageMultiple;
			int birthplace = Integer.valueOf(line.substring(5,5+birthplaceNumDigits)); // zip code?
			int occupation = Integer.valueOf(line.substring(10,13));
			int income = Integer.valueOf(line.substring(13, 19));	
			if (income >= 0)
				income = income/incomeMultiple;
			else
				income = income/incomeMultiple - 1;
//			System.out.println(income);
			// lowest-income = -10000, to make income>=0
			income += (10000/incomeMultiple);
	//		income += 2;
	//		if (income < 0) {
	//			System.out.println(line + ", income=" + income);
	//			System.exit(1);
	//		}
			// System.out.println("Age=" + age + ", Birthplace=" + birthplace + ", Occp=" + occupation + ", Income=" + income);
			ageSet.add(age);
			birthplaceSet.add(birthplace);
			occSet.add(occupation);
			incomeSet.add(income);
			
			long position = (income+ occupation*incomeSetSize + birthplace*incomeSetSize*occupationSetSize 
					+ age*incomeSetSize*occupationSetSize*birthplaceSetSize);
			if (!dataMap.containsKey(position)) {
				NonZeroEntry entry = new NonZeroEntry(0, 1);
				dataMap.put(position, entry);
			}
			else {
				NonZeroEntry entry = dataMap.get(position);
				entry.setValue(entry.getValue()+1);
			}
//			numPeople ++;
//			if (numPeople % 10000 == 0) {
//				System.out.println(numPeople);
//			}
		}
//		System.out.println("age: " + ageSet.size() + ", birthplace: " + birthplaceSet.size() + ", occupation: " + occSet.size() + ",income: " + incomeSet.size());
		double n = dataMap.size();
		//double m = ageSet.size()*birthplaceSet.size()*occSet.size()*incomeSet.size();
		double m = ageSetSize*incomeSetSize*birthplaceSetSize*occupationSetSize;
		System.out.println("Size, n: " + n);
		System.out.println("Domain size, m: " + m);
		System.out.println("Density: " + n/m);
		
		int sum = 0;
		for (NonZeroEntry entry : dataMap.values()) {
			entry.setValue(entry.getValue()*multiple);
			sum += entry.getValue();
		}
		System.out.println("Avg: " + sum*1.0/dataMap.size());
		
		ArrayList<Integer> sortedAgeSet = new ArrayList<Integer>(ageSet);
		Collections.sort(sortedAgeSet);
		ArrayList<Integer> sortedBirthplaceSet = new ArrayList<Integer>(birthplaceSet);
		Collections.sort(sortedBirthplaceSet);
		ArrayList<Integer> sortedOccSet = new ArrayList<Integer>(occSet);
		Collections.sort(sortedOccSet);
		ArrayList<Integer> sortedIncomeSet = new ArrayList<Integer>(incomeSet);
		Collections.sort(sortedIncomeSet);
		
		System.out.println("Size of age set: " + sortedAgeSet.size() + ", Smallest: " + sortedAgeSet.get(0) + ", Largest: " + sortedAgeSet.get(ageSet.size()-1));
		System.out.println("Size of birthplace set: " + sortedBirthplaceSet.size() + ", Smallest: " + sortedBirthplaceSet.get(0) + ", Largest: " + sortedBirthplaceSet.get(birthplaceSet.size()-1));
		System.out.println("Size of occupation set: " + sortedOccSet.size() + ", Smallest: " + sortedOccSet.get(0) + ", Largest: " + sortedOccSet.get(occSet.size()-1));
		System.out.println("Size of income set: " + sortedIncomeSet.size() + ", Smallest: " + sortedIncomeSet.get(0) + ", Largest: " + sortedIncomeSet.get(incomeSet.size()-1));
//		System.exit(0);
		
		PrintWriter outWriter= new PrintWriter(new FileWriter(outputFileName, false));
		outWriter.println((long) m);
		Iterator<Long> keyIter = dataMap.keySet().iterator();
		while (keyIter.hasNext()) {
			Long pos = keyIter.next();
			int val = (int) (dataMap.get(pos).getValue());
			outWriter.println(pos + "," + val);
		}
		
//		for (int i=0; i<ageSet.size(); i++) {
//			for (int j=0; j<birthplaceSet.size(); j++) {
//				for (int k=0; k<occSet.size(); k++) {
//					for (int l=0; l<incomeSet.size(); l++) {
//						int oldPos = (int) (sortedIncomeSet.get(l)+sortedOccSet.get(k)*1E2
//								+sortedBirthplaceSet.get(j)*1E5 + sortedAgeSet.get(i)*1E7);
//						if (dataMap.containsKey(oldPos)) {
//							int val = (int) (dataMap.get(oldPos).getValue());
//							int newPos = l + k*incomeSet.size() + j*incomeSet.size()*occSet.size()
//								+ i*incomeSet.size()*occSet.size()*birthplaceSet.size();
//							outWriter.println(newPos + "," + val);
//						}
//					}
//				}
//			}
//		}
		outWriter.close();
		return;
	}
	
	public static void  main(String[] args) throws Exception{
		readCensusIncome();
	}
	
}
