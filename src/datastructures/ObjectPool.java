package datastructures;

import java.util.ArrayList;
import java.util.LinkedList;

import datastructures.DyadicRange_Full.BinaryTree;
import datastructures.UtilityObject.*;

public class ObjectPool {
	static LinkedList<BinaryTree> treeList = new LinkedList<BinaryTree>();
	static LinkedList<NonZeroEntry> nonZeroEntryList = new LinkedList<NonZeroEntry>();
	static LinkedList<ZeroIndexInterval> zeroIndexIntervalList = new LinkedList<ZeroIndexInterval>();
	
	static LinkedList<NonZeroEntry_2D> nonZeroEntry_2D_List = new LinkedList<NonZeroEntry_2D>();
	static LinkedList<NonZeroEntry_Long> longNonZeroEntryList = new LinkedList<NonZeroEntry_Long>();
	
	public static BinaryTree allocalteTree() {
		if (!treeList.isEmpty())
			return treeList.remove();
		else
			return new BinaryTree(0, null, null);
	}
	
	public static void returnTree(BinaryTree tree) {
		if (tree != null) {
			BinaryTree left = tree.getLeft();
			BinaryTree right = tree.getRight();
			tree.setLeft(null);
			tree.setRight(null);
			treeList.add(tree);
			returnTree(left);
			returnTree(right);
		}
	}
	
	public static NonZeroEntry allocateNonZeroEntry() {
		if (!nonZeroEntryList.isEmpty())
			return nonZeroEntryList.remove();
		else
			return new NonZeroEntry(0,0);
	}
	
	public static void returnNonZeroEntry(NonZeroEntry entry) {
		nonZeroEntryList.add(entry);
	}
	
	public static void returnNonZeroEntries(ArrayList<NonZeroEntry> entries) {
		for (NonZeroEntry entry: entries)
			nonZeroEntryList.add(entry);
	}
	
	public static NonZeroEntry_Long allocateNonZeroEntry_Long() {
		if (!longNonZeroEntryList.isEmpty())
			return longNonZeroEntryList.remove();
		else
			return new NonZeroEntry_Long(0,0);
	}
	
	public static void returnNonZeroEntry_Long(NonZeroEntry_Long entry) {
		longNonZeroEntryList.add(entry);
	}
	
	public static void returnNonZeroEntries_Long(ArrayList<NonZeroEntry_Long> entries) {
		for (NonZeroEntry_Long entry: entries)
			longNonZeroEntryList.add(entry);
	}
	
	public static NonZeroEntry_2D allocateNonZeroEntry_2D() {
		if (!nonZeroEntry_2D_List.isEmpty())
			return nonZeroEntry_2D_List.remove();
		else
			return new NonZeroEntry_2D(0,0,0,0);
	}
	
	public static void returnNonZeroEntry_2D(NonZeroEntry_2D entry) {
		nonZeroEntry_2D_List.add(entry);
	}
	
	public static void returnNonZeroEntries_2D(ArrayList<NonZeroEntry_2D> entries) {
		for (NonZeroEntry_2D entry: entries)
			nonZeroEntry_2D_List.add(entry);
	}
	
	public static ZeroIndexInterval allocateZeroIndexInterval() {
		if (!zeroIndexIntervalList.isEmpty())
			return zeroIndexIntervalList.remove();
		else
			return new ZeroIndexInterval(0,0);
	}
	
	public static void returnZeroIndexInterval(ZeroIndexInterval interval) {
		zeroIndexIntervalList.add(interval);
	}
	
	public static void returnZeroIndexIntevals(ArrayList<ZeroIndexInterval> intervals) {
		for (ZeroIndexInterval interval : intervals)
			zeroIndexIntervalList.add(interval);
	}
	
}
