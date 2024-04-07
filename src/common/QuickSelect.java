package common;

import java.util.*;

/* Quick selection algorithm.
Ê Ê Ê* Places the kth smallest item in a[k-1].
Ê Ê Ê* @param <C>
Ê Ê Ê* @param a an array of Comparable items.
Ê Ê Ê* @param k the desired rank (1 is minimum) in the entire array.
Ê Ê Ê*/
public class QuickSelect {

	public static void quickselect( Comparable [ ] a, int k ) {
        quickselect( a, 0, a.length - 1, k );
    }

	private static final int CUTOFF = 10;

    private static void quickselect( Comparable [ ] a, int low, int high, 
									 int k ) {
        if( low + CUTOFF > high )
            insertionSort( a, low, high );
        else {
            // Sort low, middle, high
            int middle = ( low + high ) / 2;
            if( a[ middle ].compareTo( a[ low ] ) < 0 )
                swapReferences( a, low, middle );
            if( a[ high ].compareTo( a[ low ] ) < 0 )
                swapReferences( a, low, high );
            if( a[ high ].compareTo( a[ middle ] ) < 0 )
                swapReferences( a, middle, high );

               // Place pivot at position high - 1
            swapReferences( a, middle, high - 1 );
            Comparable pivot = a[ high - 1 ];

                // Begin partitioning
            int i, j;
            for( i = low, j = high - 1; ; )
            {
                while( a[ ++i ].compareTo( pivot ) < 0 )
                    ;
                while( pivot.compareTo( a[ --j ] ) < 0 )
                    ;
                if( i >= j )
                    break;
                swapReferences( a, i, j );
            }
 
                // Restore pivot
            swapReferences( a, i, high - 1 );
 
				// Recurse on the relevant sub-array
            if( k <= i )
                quickselect( a, low, i - 1, k );
            else if( k > i + 1 )
                quickselect( a, i + 1, high, k );
        }
    }

    private static void swapReferences( Object [ ] a, int index1, int index2 ) {
        Object tmp = a[ index1 ];
        a[ index1 ] = a[ index2 ];
        a[ index2 ] = tmp;
    }
 
	private static void insertionSort( Comparable [ ] a, int low, int high ) {
        for( int p = low + 1; p <= high; p++ ) {
            Comparable tmp = a[ p ];
            int j;
 
            for( j = p; j > low && tmp.compareTo( a[ j - 1 ] ) < 0; j-- )
                a[ j ] = a[ j - 1 ];
            a[ j ] = tmp;
        }
    }
 
	public static void quickselect( ArrayList<Comparable>  a, int k ) {
        quickselect( a, 0, a.size() - 1, k );
    }


    private static void quickselect( ArrayList<Comparable>  a, int low, int high, int k ) {
        if( low + CUTOFF > high )
            insertionSort( a, low, high );
        else {
            // Sort low, middle, high
            int middle = ( low + high ) / 2;
            if( a.get( middle ).compareTo( a.get( low ) ) < 0 )
                swapReferences( a, low, middle );
            if( a.get( high ).compareTo( a.get( low ) ) < 0 )
                swapReferences( a, low, high );
            if( a.get( high ).compareTo( a.get( middle ) ) < 0 )
                swapReferences( a, middle, high );

               // Place pivot at position high - 1
            swapReferences( a, middle, high - 1 );
            Comparable pivot = a.get( high - 1 );

                // Begin partitioning
            int i, j;
            for( i = low, j = high - 1; ; )
            {
                while( a.get( ++i ).compareTo( pivot ) < 0 )
                    ;
                while( pivot.compareTo( a.get( --j ) ) < 0 )
                    ;
                if( i >= j )
                    break;
                swapReferences( a, i, j );
            }
 
                // Restore pivot
            swapReferences( a, i, high - 1 );
 
				// Recurse on the relevant sub-array
            if( k <= i )
                quickselect( a, low, i - 1, k );
            else if( k > i + 1 )
                quickselect( a, i + 1, high, k );
        }
    }

    private static void swapReferences( ArrayList<Comparable>  a, int index1, int index2 ) {
        Comparable tmp = a.get( index1 );
        a.set(index1, a.get( index2 ));
        a.set(index2, tmp);
    }
 
	private static void insertionSort( ArrayList<Comparable> a, int low, int high ) {
        for( int p = low + 1; p <= high; p++ ) {
            Comparable tmp = a.get( p );
            int j;
 
            for( j = p; j > low && tmp.compareTo( a.get( j - 1 ) ) < 0; j-- )
                a.set( j , a.get( j - 1 ));
            a.set( j , tmp);
        }
    }
 
	
	public static void main(String[] args) {
		for (int SIZE = 100; SIZE < 1000000; SIZE *= 2) {
	    	long start, end;
	    	long elapsed1 = 0, elapsed2 = 0, elapsed3 = 0;
	    	Comparable [ ] a = new Comparable [ SIZE ];
 
			// sorted input
	    	for( int i = 0; i < SIZE; i++ ) {
				a[ i ] = new Integer( i );
			}
	    	start = System.currentTimeMillis();
	    	quickselect( a, SIZE / 2 );
	    	end = System.currentTimeMillis();
	    	elapsed1 = end - start;
 
			// reverse sorted input
		    for( int i = 0; i < SIZE; i++ ) {
				a[ i ] = new Integer( SIZE - i );
			}
	    	start = System.currentTimeMillis();
	    	quickselect( a, SIZE / 2);
	    	end = System.currentTimeMillis();
	    	elapsed2 = end - start;
 
			// random input
		    Random r = new Random();
		    for( int i = 0; i < SIZE; i++ ) {
				a[ i ] = new Integer( r.nextInt(SIZE) );
			}
		    start = System.currentTimeMillis();
	    	quickselect( a, SIZE / 2);
	    	end = System.currentTimeMillis();
		    elapsed3 = end - start;
 
		    System.out.println("size: " + SIZE + 
							   "\tsorted: " + elapsed1 + 
						       "\treverse: " + elapsed2 + 
							   "\trandom: " + elapsed3);
		}
	}
 
}
 