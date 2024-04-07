package datastructures;

import java.util.Comparator;


public class UtilityObject {
	public static class NonZeroEntry implements Comparable {
		int position;
		double value;
		Double rand = null;
//		double priority = -1;
		
		
		public NonZeroEntry(int p, double v) {
			position = p; value = v; 
			//priority = Math.abs(value)/rand;
		}

		
		public int compareTo(Object other) {
			if (other instanceof NonZeroEntry) {
				NonZeroEntry otherEntry = (NonZeroEntry) other;
				if (this.getPriority() > otherEntry.getPriority())
					return 1;
				else if (this.getPriority() == otherEntry.getPriority())
					return 0;
				else 
					return -1;
			}
			else {
				System.err.println("Comparing objects of different types");
				System.exit(1);
				return 0;
			}
		}
	
		public void setPosition(int pos) {
			position = pos;
		}
		public void setValue(double val) {
			value = val;
		}
		public void setRand(double r) {
			rand = r;
		}
		public double getPriority() {
			return Math.abs(value)/rand;
		}
		
		public int getPosition() {
			return (int)position;
		}
		public double getValue() {
			return value;
		}
		public double getRand() {
			return rand;
		}
		public String toString() {
			return "{Val=" + value + ", Pos=" + position + ", Pri=" + Math.abs(value)/rand +"}";
		}
	}
	
	public static class NonZeroEntry_Long implements Comparable {
		long position;
		double value;
		double rand;
		
		
		public NonZeroEntry_Long(long p, double v) {
			position = p; value = v; 
			//priority = Math.abs(value)/rand;
		}

		
		public int compareTo(Object other) {
			if (other instanceof NonZeroEntry_Long) {
				NonZeroEntry_Long otherEntry = (NonZeroEntry_Long) other;
				if (this.getPriority() > otherEntry.getPriority())
					return 1;
				else if (this.getPriority() == otherEntry.getPriority())
					return 0;
				else 
					return -1;
			}
			else {
				System.err.println("Comparing objects of different types");
				System.exit(1);
				return 0;
			}
		}
	
		public void setPosition(long pos) {
			position = pos;
		}
		public void setValue(double val) {
			value = val;
		}
		public void setRand(double r) {
			rand = r;
		}
		public double getPriority() {
			return Math.abs(value)/rand;
		}
		
		public long getPosition() {
			return position;
		}
		public double getValue() {
			return value;
		}
		public double getRand() {
			return rand;
		}
		public String toString() {
			return "{Val=" + value + ", Pos=" + position + ", Pri=" + Math.abs(value)/rand +"}";
		}
	}
	
	public static class NonZeroEntry_2D implements Comparable {
		int x;
		int y;
		double value;
		double rand;
//		double priority = -1;
		//private double priority;
		
		public NonZeroEntry_2D(int x, int y, double v,  double r) {
			this.x = x; this.y = y; value = v; rand = r;
		}
		
//		public NonZeroEntry(int p, int p2, double v,  double r) {
//			position = p; position2 = p2; value = v; rand = r;
//		}
		
		public int compareTo(Object other) {
			if (other instanceof NonZeroEntry) {
				NonZeroEntry otherEntry = (NonZeroEntry) other;
				if (this.getPriority() > otherEntry.getPriority())
					return 1;
				else if (this.getPriority() == otherEntry.getPriority())
					return 0;
				else 
					return -1;
			}
			else {
				System.err.println("Comparing objects of different types");
				System.exit(1);
				return 0;
			}
		}
	
		public void setXY(int x, int y) {
			this.x = x; this.y = y;
		}
		public void setValue(double val) {
			value = val;
		}
		public void setRand(double r) {
			rand = r;
		}
		public double getPriority() {
			return Math.abs(value)/rand;
		}
		
		public int getX() {
			return x;
		}
		public int getY() {
			return y;
		}
		public double getValue() {
			return value;
		}
		public double getRand() {
			return rand;
		}
		public String toString() {
			return "{Val=" + value + ", X=" + x + ", Y=" + y + ", Pri=" + Math.abs(value)/rand +"}";
		}
	}
	
	public static class ZeroIndexInterval {
		int startIndex;
		int endIndex;
		double cumProb; // prob of zero entries less than this endIndex
		
		public ZeroIndexInterval(int s, int e) {
			startIndex = s;
			endIndex = e;
		}
		
		public void setStartIndex(int s) {
			startIndex = s;
		}
		
		public void setEndIndex(int e) {
			endIndex = e;
		}
		
		public void setCumProb(double prob) {
			cumProb = prob;
		}
		
		public int getStartIndex() {
			return startIndex;
		}
		public int getEndIndex() {
			return endIndex;
		}
	}
	
	public static class PositionComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			if (o1 instanceof NonZeroEntry && o2 instanceof NonZeroEntry) {
				NonZeroEntry e1 = (NonZeroEntry) o1;
				NonZeroEntry e2 = (NonZeroEntry) o2;
				if (e1.getPosition() >= e2.getPosition())
					return 1;
				else 
					return -1;
			}
			return 0;
		}
	}
	
	public static class PositionComparator_Long implements Comparator {
		public int compare(Object o1, Object o2) {
			if (o1 instanceof NonZeroEntry_Long && o2 instanceof NonZeroEntry_Long) {
				NonZeroEntry_Long e1 = (NonZeroEntry_Long) o1;
				NonZeroEntry_Long e2 = (NonZeroEntry_Long) o2;
				if (e1.getPosition() >= e2.getPosition())
					return 1;
				else 
					return -1;
			}
			return 0;
		}
	}
	
	
}
