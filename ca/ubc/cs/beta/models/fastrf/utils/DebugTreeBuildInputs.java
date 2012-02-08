package ca.ubc.cs.beta.models.fastrf.utils;

import java.io.Serializable;
import ca.ubc.cs.beta.models.fastrf.*;

public final class DebugTreeBuildInputs implements Serializable {
	private static final long serialVersionUID = -3226427544555346376L;
	public double[][] allTheta; 
	public double[][] allX;
	public int[][] dataIdxs;
	public double[] y;
	public boolean[] cens;
	public RegtreeBuildParams params;
}
