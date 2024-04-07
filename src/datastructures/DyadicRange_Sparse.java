package datastructures;

import main.Run.Mechanism;

public interface DyadicRange_Sparse {
	public void reset();
	public int anonymize(Mechanism mechanism, int theta, int K, double epsilon, int noiseSeed, int sampleSeed, int zeroPosSeed);
}
