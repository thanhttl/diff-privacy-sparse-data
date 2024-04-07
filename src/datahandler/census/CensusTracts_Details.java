package datahandler.census;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.StringTokenizer;



public class CensusTracts_Details {
	static int numTractDigits=2;
	
	static HashSet<String> excludedStateCodes;
	
	public static void main(String[] args) {
		HashMap<String,Integer> states = readStateCodes("data/map_info/state_codes.txt");
		System.out.println(states.size());
		ArrayList<String> excludedStates = new ArrayList<String>();
		excludedStates.add("ct"); excludedStates.add("dc"); excludedStates.add("ma"); excludedStates.add("nh");
		excludedStates.add("pr"); excludedStates.add("vi");
		
		excludedStateCodes = new HashSet<String>(10);
		excludedStateCodes.add("09"); excludedStateCodes.add("11"); excludedStateCodes.add("25");
		excludedStateCodes.add("33"); excludedStateCodes.add("72");
		
		String[] fileNames = {"data/map_info/al_mt_tracts.txt", "data/map_info/ne_wy_tracts.txt"};
		HashMap<String, ArrayList<String>> tractMap = CensusTracts.readLocationList(fileNames, numTractDigits);
		
		ArrayList<String> allTracts = new ArrayList<String>();
 		Iterator<String> countyIter = tractMap.keySet().iterator();
 		HashSet<String> statesIncluded = new HashSet<String>();
 		
		while (countyIter.hasNext()) {
			String county = countyIter.next();
			String state = county.substring(0,2).toLowerCase();
			if (states.containsKey(state) && !excludedStates.contains(state)) {
				statesIncluded.add(state);
				ArrayList<String> tracts = tractMap.get(county);
				for (String tractID : tracts) {
						allTracts.add(county.substring(2) + tractID);
				}
			}
		}
		// write the list of tracts to file
		String tractFile = "data/map_freq_matrix/validTracts_" + numTractDigits + "_details.txt";
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(tractFile));
			for (String tract: allTracts)
				pw.println(tract);
			pw.close();
		}
		catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		System.out.println("Num states included: " + statesIncluded.size());
		System.out.println("Num of coarsened tracts: " + allTracts.size());
		
		HashMap<Long, Integer> freqMatrix = new HashMap<Long, Integer>();
		Iterator<String> stateIter = states.keySet().iterator();
		HashSet<String> validTractSet = new HashSet<String>(allTracts);
		int index = 1;
		
		while (stateIter.hasNext()) {
			String state = stateIter.next();
			System.out.println(index + ". " + state);
			//if (state.equals("mi"))
				readFreqMatrix(state, states.get(state).toString(), freqMatrix, validTractSet, numTractDigits, -1);
			System.out.println("-- Num non zeros: " + freqMatrix.size());
			index++;
		}
		
		
		ArrayList<Long> countArray = new ArrayList<Long>(freqMatrix.keySet());
		Collections.sort(countArray);
		String outFile = "data/map_freq_matrix/wholeUS_" + numTractDigits + ".txt";
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
			int i=0;
			int totalNum = 0;
			int numZeroEntries = 0;
			for (Long loc: countArray) {
				int count = freqMatrix.get(loc);
				pw.println(loc + ", " + count);
				totalNum += count;
				if (count == 0)
					numZeroEntries ++;
				i++;
			}
			pw.close();
			
			System.out.println("Total num of people: " + totalNum);
			System.out.println("Num non-zeros: " + freqMatrix.size());
			System.out.println("Percentage of sparsity: " + (1- ((double) freqMatrix.size())/Math.pow(allTracts.size(), 2)));
			System.out.println("Avg non-zero entry: " + ((double) totalNum)/(freqMatrix.size()));
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
//		String outFile = "data/map_freq_matrix/wholeUS_" + numTractDigits + ".txt";
//		try {
//			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
//			Iterator<String> keyIter = freqMatrix.keySet().iterator();
//			int i=0;
//			int totalNum = 0;
//			int numZeroEntries = 0;
//			while (keyIter.hasNext()) {
//				String key = keyIter.next();
//				int count = freqMatrix.get(key);
//				pw.println(key + ", " + count);
//				totalNum += count;
//				if (count == 0)
//					numZeroEntries ++;
//				i++;
//			}
//			pw.close();
//			
//			System.out.println("Total num of people: " + totalNum);
//			System.out.println("Num non-zeros: " + freqMatrix.size());
//			System.out.println("Percentage of sparsity: " + (1- ((double) freqMatrix.size())/Math.pow(allTracts.size(), 2)));
//			System.out.println("Avg non-zero entry: " + ((double) totalNum)/(freqMatrix.size()));
//		}
//		catch(Exception e) {
//			e.printStackTrace();
//			System.exit(1);
//		}
//		
		
	}
	
	public static void readState(String stateName, String stateCode) {
		
		//String stateCode = "34"; String stateName = "nj"; // Texas
		HashMap<String, ArrayList<String>> tractMap = null; // TODO
		HashMap<Long, Integer> freqMatrix = readFreqMatrix_old(stateName, stateCode, tractMap, -1);
		
		long start = System.currentTimeMillis();
		String outFile = "data/map_freq_matrix/" + stateName + "_" + numTractDigits + ".txt";
		try {
			PrintWriter pw = new PrintWriter(new FileWriter(outFile));
			Iterator<Long> keyIter = freqMatrix.keySet().iterator();
			int i=0;
			int totalNum = 0;
			int numZeroEntries = 0;
			while (keyIter.hasNext()) {
				long key = keyIter.next();
				int count = freqMatrix.get(key);
				pw.println(key + ", " + count);
				totalNum += count;
				if (count == 0)
					numZeroEntries ++;
				i++;
			}
			pw.close();
			
			System.out.println("Frequency matrix size: " + freqMatrix.size());
			System.out.println("Total num of people: " + totalNum);
			System.out.println("Num zero entries: " + numZeroEntries);
			System.out.println("Percentage of sparsity: " + ((double) numZeroEntries)/freqMatrix.size());
			System.out.println("Avg non-zero entry: " + ((double) totalNum)/(freqMatrix.size()-numZeroEntries));
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		long end = System.currentTimeMillis();
		System.out.println("Read duration:" + (end-start)/1000/60 + " mins");
	}
	
	// read the census tract file, generate the list of available census tracts for each
	// (state, county).
	// create the 2-d matrix where source and destination correspond to a possible pair
	// of census tracts.
	// key: state-county string (5-digit), 
	// value: array list of tract numbers (6-digit)
	// numDigits < 6 indicates coarsened tracts -- just use up to numDigits
	public static HashMap<String, ArrayList<String>> readLocationList(String[] fileNames, int numDigits) {
		HashMap<String, ArrayList<String>> censusTracts;
		censusTracts = new HashMap<String, ArrayList<String>>();
		// AB12 StateName
		// AB12XYZCountyName (the 2nd token is "County")
		// AB12XYZ TractNumber
		for (String fileName: fileNames) {
			// read the file
			BufferedReader reader=null;
			try {
				reader = new BufferedReader(new FileReader(fileName));
				reader.readLine(); // first line contains attribute names
			}
			catch(Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			
			String countyID = null;
			ArrayList<String> tractList = null;
			while (true) {
				try {
					String line = reader.readLine();
					if (line == null)
						break;
					StringTokenizer st = new StringTokenizer(line, " ");
					// each line should have at least two tokens
					String firstToken = st.nextToken();
					if (firstToken.length() == 4) {
						// skip
					}
					else if (firstToken.length() > 7) {
						if (tractList == null) {
							tractList = new ArrayList<String>();
							countyID = firstToken.substring(0, 7);
						}
						else {
							censusTracts.put(countyID, tractList);
							tractList = new ArrayList<String>();
							countyID = firstToken.substring(0, 7);
						}
					}
					else { // census tract info
						String tractID = st.nextToken();
						int index = tractID.indexOf(".");
						if (index == -1) {
							// add leading zeros to get 6 digits
							int numZeros = 4 - tractID.length();
							for (int i=0; i<numZeros; i++) {
								tractID = "0" + tractID;
							}
							tractID = tractID + "00";
						}
						else {
							int numZeros = 7 - tractID.length();
							for (int i=0; i<numZeros; i++) {
								tractID = "0" + tractID;
							}
							tractID = tractID.substring(0,4) + tractID.substring(5,7);	
						}
						
						if (numDigits < 6)
							tractID = tractID.substring(0, numDigits);
						
						//System.out.println(tractID);
						String candidate = tractID;
						if (!tractList.contains(candidate))
							tractList.add(candidate);
					}
				}
				catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
		}
		return censusTracts;
	}
	
	public static HashMap<String,Integer> readStateCodes(String fileName) {
		HashMap<String,Integer> stateMap = new HashMap<String,Integer>();
		BufferedReader reader=null;
		try {
			reader = new BufferedReader(new FileReader(fileName));
			while (true) {
				String line = reader.readLine();
				if (line == null)
					break;
				StringTokenizer st = new StringTokenizer(line, ", ");
				String stateAbbr = st.nextToken();
				Integer code = Integer.valueOf(st.nextToken());
				stateMap.put(stateAbbr, code);
			}	
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		return stateMap;
	}
	
	
	// the key is a String of source_census_tract + dest_census_tract  
	public static void readFreqMatrix(String stateName, String stateCode, HashMap<Long, Integer> fm, HashSet<String> tracts, int numDigits, int maxNumLines) {
		String mainFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_main_ja_2008_1.csv";
		String auxFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_aux_ja_2008_1.csv";
		if (stateName.equalsIgnoreCase("nc")) {
			mainFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_main_ja_2007_1.csv";
			auxFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_aux_ja_2007_1.csv";
		}
		String[] files = {mainFile, auxFile};
		for (String fileName : files) {
			BufferedReader reader=null;
			try {
				reader = new BufferedReader(new FileReader(fileName));
				reader.readLine(); // first line contains attribute names
			}
			catch(Exception e) {
				e.printStackTrace();
				continue;
				//System.exit(1);
			}
			System.out.println("Reading file: " + fileName);
			int numLines = 0;
			//int tractIDLength = allTracts.get(0).length();
			//System.out.println("Tract id length=" + tractIDLength + "; " + allTracts.get(0));
			while (true) {
				try {
					String line = reader.readLine();
					if (line == null)
						break;
					StringTokenizer st = new StringTokenizer(line, ", ");
					
					String work = st.nextToken().substring(0, 5+numDigits); //source
					String home = st.nextToken().substring(0, 5+numDigits); //home
					if (excludedStateCodes.contains(work.substring(0,2)) || excludedStateCodes.contains(work.substring(0,2)) || !tracts.contains(work) || !tracts.contains(home))
						continue;
					
					int count = Integer.valueOf(st.nextToken()); // total
					int age1 = Integer.valueOf(st.nextToken());
					int age2 = Integer.valueOf(st.nextToken());
					int age3 = Integer.valueOf(st.nextToken());
					int earn1 = Integer.valueOf(st.nextToken());
					int earn2 = Integer.valueOf(st.nextToken());
					int earn3 = Integer.valueOf(st.nextToken());
					int ind1 = Integer.valueOf(st.nextToken());
					int ind2 = Integer.valueOf(st.nextToken());
					int ind3 = Integer.valueOf(st.nextToken());
					//int age = 1000 + 100*age1 + 10*age2 + age3;
					//int earn = 1000 + 100*earn1 + 10*earn2 + earn3;
					//int ind = 1000 + 100*ind1 + 10*ind2 + ind3;
					int age=0;
					if (age1 == 1) age=2;
					else if (age2 == 1) age=1;
					else if (age3 == 1) age=0;
					int earn=0;
					if (earn1 == 1) earn=2;
					else if (earn2 == 1) earn=1;
					else if (earn3 == 1) earn=0;
					int ind=0;
					if (ind1 == 1) ind=2;
					else if (ind2 == 1) ind=1;
					else if (ind3 == 1) ind=0;
					
					int code = age*9 + earn*3 + ind;
					//String sourceDestStr = work + home + Integer.toString(age) + Integer.toString(earn) + Integer.toString(ind);
					String sourceDestStr = "1" + work + home + Integer.toString(code);
					Long sourceDest = Long.valueOf(sourceDestStr);
					if (fm.containsKey(sourceDest)) {
						fm.put(sourceDest, fm.get(sourceDest) + count);
						numLines ++;
					}
					else {
						fm.put(sourceDest, 1);
						
						numLines ++;
	//					System.out.println("Line: " + line);
	//					System.out.println("---" + sourceDest);
	//					System.out.println("Invalid source-dest pair");
						//System.exit(1);
					}
					
					if (maxNumLines > 0 && numLines == maxNumLines)
						break;
				}
				catch (Exception e) {
					e.printStackTrace();
//					System.exit(1);
				}
			}
			System.out.println("Read " + numLines + " lines. Size of freq matrix= " + fm.size());
		}
	}
	
	
	// the key is a String of source_census_tract + dest_census_tract  
	public static HashMap<Long, Integer> readFreqMatrix_old(String stateName, String stateCode, HashMap<String, ArrayList<String>> tractMap, int maxNumLines) {
		String mainFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_main_ja_2008_1.csv";
		String auxFile = "data/onthemap/" + stateName + "/od/" + stateName + "_od_aux_ja_2008_1.csv";
		HashMap<Long, Integer> fm = new HashMap<Long, Integer>();
		
		ArrayList<String> allTracts = new ArrayList<String>();
 		Iterator<String> countyIter = tractMap.keySet().iterator();
		while (countyIter.hasNext()) {
			String county = countyIter.next();
			ArrayList<String> tracts = tractMap.get(county);
			for (String tractID : tracts) 
				allTracts.add(county + tractID);
		}
		
		int numEntries = 0;
		System.out.println("Number of census tracts: " + allTracts.size());
		for (String sourceTract: allTracts) {
			if (sourceTract.startsWith(stateCode)) {
				for (String destTract: allTracts) {
					if (destTract.startsWith(stateCode)) {	
						// put a dummy "1" at front to avoid the truncation of leading zeros
						fm.put(Long.valueOf("1" + sourceTract.substring(2) + destTract.substring(2)), 0);
						numEntries ++;
						// if (numEntries % 1000 == 0)
						//	System.out.println(numEntries);
					}
				}
			}
		}
		
		String[] files = {mainFile, auxFile};
		for (String fileName : files) {
			BufferedReader reader=null;
			try {
				reader = new BufferedReader(new FileReader(fileName));
				reader.readLine(); // first line contains attribute names
			}
			catch(Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
			System.out.println("Reading file: " + fileName);
			int numLines = 0;
			int tractIDLength = allTracts.get(0).length();
			//System.out.println("Tract id length=" + tractIDLength + "; " + allTracts.get(0));
			while (true) {
				try {
					String line = reader.readLine();
					if (line == null)
						break;
					StringTokenizer st = new StringTokenizer(line, ", ");
					
					String work = st.nextToken(); //source
					String home = st.nextToken(); //home
					int count = Integer.valueOf(st.nextToken());
					String sourceDestStr = "1" + work.substring(2, tractIDLength) + home.substring(2, tractIDLength);
					Long sourceDest = Long.valueOf(sourceDestStr);
					if (fm.containsKey(sourceDest)) {
						fm.put(sourceDest, fm.get(sourceDest) + count);
						numLines ++;
					}
					else {
						fm.put(sourceDest, 1);
						numLines ++;
	//					System.out.println("Line: " + line);
	//					System.out.println("---" + sourceDest);
	//					System.out.println("Invalid source-dest pair");
						//System.exit(1);
					}
					
					if (maxNumLines > 0 && numLines == maxNumLines)
						break;
				}
				catch (Exception e) {
					e.printStackTrace();
					//System.exit(1);
				}
			}
			System.out.println("Read " + numLines + " lines. Size of freq matrix= " + fm.size());
		}
		return fm;
	}
	
	public static class LocCount {
		public String id;
		public int count;
		public LocCount(String s, int c) {
			id = s; count = c;
		}
	}
}
