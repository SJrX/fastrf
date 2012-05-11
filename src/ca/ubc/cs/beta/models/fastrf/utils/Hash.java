package ca.ubc.cs.beta.models.fastrf.utils;

import java.util.Arrays;
/**
 * Utility class to handle generating lossy float based hashes on objects
 * (If not in strictfp mode, the last few digits on 
 * @author seramage
 *
 */
public strictfp class Hash {
	
	public static int hash(int[] obj)
	{
		return Arrays.hashCode(obj);
	}

	public static int hash(int[][] obj)
	{
		return Arrays.deepHashCode(obj);
	}
	
	protected static float[] cFloat(double[] value)
	{
		if(value == null) return null;
		float[] arr = new float[value.length];
		for(int i=0; i < value.length; i++)
		{
			arr[i] = (float) value[i];
		}
		return arr;
	}
	
	protected static float[][] cFloat(double[][] value)
	{
		if(value == null) return null;
		float[][] arr = new float[value.length][];
		for(int i=0; i < value.length; i++)
		{
			arr[i] = cFloat(value[i]);
		}
		return arr;
	}

	public static int hash(double[] obj)
	{
		
		return Arrays.hashCode(cFloat(obj));
	}
	public static int hash(double[][] obj)
	{
		
		return Arrays.deepHashCode(cFloat(obj));
	}
}
