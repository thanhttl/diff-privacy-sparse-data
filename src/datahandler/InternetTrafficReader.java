package datahandler;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;


public class InternetTrafficReader {
	public static String fileName = "data/internet_traffic/"; //trace.txt";
	public static String outputFileName = "data/internet_traffic/out/"; //trace1.txt";
	
	public static void readTrafficData() throws Exception {
		String[] fileIDs = {"trace1.txt", "trace2.txt", "trace3.txt", "trace4.txt"};
		for (String fileID : fileIDs) {
			System.out.println(fileID);
			String fullFileName = fileName + fileID;
			
			ArrayList<Request> requestList = new ArrayList<Request>();
			TreeSet<Long> clientIPs = new TreeSet<Long>();
			TreeSet<Long> serverIPs = new TreeSet<Long>();
			
			long power3 = 256*256*256; long power2 = 256*256;
			String line;
			int numBadRequests = 0;
			
			
				
			BufferedReader reader = new BufferedReader(new FileReader(fullFileName));
			while ((line = reader.readLine()) != null) {
				String[] tokens = line.split(" ");
				String client = tokens[3].substring(0, tokens[3].indexOf(":"));
				String[] cs = client.split("[.]");
				long clientInt = Integer.valueOf(cs[0])*power3 + Integer.valueOf(cs[1])*power2 + Integer.valueOf(cs[2])*256 + Integer.valueOf(cs[3]);
				String server = tokens[4].substring(0, tokens[4].indexOf(":"));
				String[] ss = server.split("[.]");
				long serverInt = Integer.valueOf(ss[0])*power3 + Integer.valueOf(ss[1])*power2 + Integer.valueOf(ss[2])*256 + Integer.valueOf(ss[3]);
				String responseLength = tokens[11];
				if (responseLength.length() > 9)
					numBadRequests ++;
				else {
					Request request = new Request(clientInt, serverInt, Integer.valueOf(responseLength));
					requestList.add(request);			
					clientIPs.add(clientInt);
					serverIPs.add(serverInt);
				}
			}	
			
			System.out.println("Num of clients: " + clientIPs.size());
			System.out.println("Num of servers: " + serverIPs.size());
			System.out.println("Num bad requests: " + numBadRequests);
			
			Collections.sort(requestList);
			
			HashMap<Long, Integer> clientMap = new HashMap<Long, Integer>();
			HashMap<Long, Integer> serverMap = new HashMap<Long, Integer>();
			Iterator<Long> clientIter = clientIPs.iterator();
			int index = 0;
			while (clientIter.hasNext()) {
				long ip = clientIter.next();
				clientMap.put(ip, index++);
			}
			Iterator<Long> serverIter = serverIPs.iterator();
			index = 0;
			while (serverIter.hasNext()) {
				long ip = serverIter.next();
				serverMap.put(ip, index++);
			}
			
			String fullOutputFileName = outputFileName + fileID;
			PrintWriter outWriter= new PrintWriter(new FileWriter(fullOutputFileName, false));
			outWriter.println(clientIPs.size()*serverIPs.size());
			
			Request oldRequest = requestList.get(0);
			int currLength = oldRequest.numPackets;
			long totalLength = 0;
			for (int i=1; i< requestList.size(); i++) {
				Request currRequest = requestList.get(i);
				if (currRequest.clientIP==oldRequest.clientIP && currRequest.serverIP==oldRequest.serverIP)
					currLength += currRequest.numPackets;
				else {
					//outWriter.println(clientMap.get(oldRequest.clientIP) + "," + serverMap.get(oldRequest.serverIP) + "," + currLength);
					outWriter.println(oldRequest.clientIP + "," + oldRequest.serverIP + "," + currLength);
					totalLength += currLength;
					oldRequest = currRequest;
					currLength = currRequest.numPackets;
				}
			}
			//outWriter.println(clientMap.get(oldRequest.clientIP) + "," + serverMap.get(oldRequest.serverIP) + "," + currLength);
			outWriter.println(oldRequest.clientIP + "," + oldRequest.serverIP + "," + currLength);
			totalLength += currLength;
			outWriter.close();
			System.out.println("Total length=" + totalLength);
		}
	}
	
	public static void mergeFiles() throws Exception {
		ArrayList<Request> requestList = new ArrayList<Request>();
		TreeSet<Long> clientIPs = new TreeSet<Long>();
		TreeSet<Long> serverIPs = new TreeSet<Long>();
		
		String[] fileIDs = {"trace1.txt", "trace2.txt", "trace3.txt", "trace4.txt"};
		
		for (String fileID : fileIDs) {
			String fileName = outputFileName + fileID;
			BufferedReader reader = new BufferedReader(new FileReader(fileName));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] tokens = line.split(",");
				if (tokens.length != 3) continue;
				
				long clientInt = Long.valueOf(tokens[0]);
				long serverInt = Long.valueOf(tokens[1]);
				int length = Integer.valueOf(tokens[2]);
				
				Request request = new Request(clientInt, serverInt, length);
				requestList.add(request);			
				clientIPs.add(clientInt);
				serverIPs.add(serverInt);
			}	
		}	
		System.out.println("Num of clients: " + clientIPs.size());
		System.out.println("Num of servers: " + serverIPs.size());
		
		Collections.sort(requestList);
		
		HashMap<Long, Integer> clientMap = new HashMap<Long, Integer>();
		HashMap<Long, Integer> serverMap = new HashMap<Long, Integer>();
		Iterator<Long> clientIter = clientIPs.iterator();
		int index = 0;
		while (clientIter.hasNext()) {
			long ip = clientIter.next();
			clientMap.put(ip, index++);
		}
		Iterator<Long> serverIter = serverIPs.iterator();
		index = 0;
		while (serverIter.hasNext()) {
			long ip = serverIter.next();
			serverMap.put(ip, index++);
		}
		
		String fullOutputFileName = outputFileName + "combinedTrace.txt";
		PrintWriter outWriter= new PrintWriter(new FileWriter(fullOutputFileName, false));
		outWriter.println(clientIPs.size()*serverIPs.size());
		
		Request oldRequest = requestList.get(0);
		int currLength = oldRequest.numPackets;
		long totalLength = 0;
		for (int i=1; i< requestList.size(); i++) {
			Request currRequest = requestList.get(i);
			if (currRequest.clientIP==oldRequest.clientIP && currRequest.serverIP==oldRequest.serverIP)
				currLength += currRequest.numPackets;
			else {
				outWriter.println(clientMap.get(oldRequest.clientIP) + "," + serverMap.get(oldRequest.serverIP) + "," + currLength);
				//outWriter.println(oldRequest.clientIP + "," + oldRequest.serverIP + "," + currLength);
				totalLength += currLength;
				oldRequest = currRequest;
				currLength = currRequest.numPackets;
			}
		}
		outWriter.println(clientMap.get(oldRequest.clientIP) + "," + serverMap.get(oldRequest.serverIP) + "," + currLength);
		//outWriter.println(oldRequest.clientIP + "," + oldRequest.serverIP + "," + currLength);
		totalLength += currLength;
		outWriter.close();
		System.out.println("Total length=" + totalLength);
	}
	
	public static void main(String[] args) throws Exception {
		//readTrafficData();
		mergeFiles();
	}
	
	public static class Request implements Comparable {
		long clientIP;
		long serverIP;
		int numPackets;
		
		public Request(long c, long s, int n) {
			clientIP = c; serverIP = s; numPackets = n;
		}
		public int compareTo(Object o) {
			Request other = (Request) o;
			if (clientIP > other.clientIP) {
				return 1;
			}
			else if (clientIP < other.clientIP) {
				return -1;
			}
			else {
				if (serverIP > other.serverIP)
					return 1;
				else if (serverIP < other.serverIP)
					return -1;
				else 
					return 0;
			}
		}
 	}
}
