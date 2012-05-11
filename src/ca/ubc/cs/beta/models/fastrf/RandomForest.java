package ca.ubc.cs.beta.models.fastrf;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.*;

public strictfp class RandomForest implements java.io.Serializable {
    private static final long serialVersionUID = 5204746081208095705L;
    
    public int numTrees;
    public Regtree[] Trees;
    public int logModel;
    public static final double MIN_VARIANCE_RESULT = -1 * Math.pow(10,-6);
    public double minVariance;

	private RegtreeBuildParams buildParams;
    
    public RandomForest(int numtrees, RegtreeBuildParams buildParams) {
        if (numtrees <= 0) {
            throw new RuntimeException("Invalid number of regression trees in forest: " + numtrees);
        }
        this.logModel = buildParams.logModel;
        numTrees = numtrees;
        Trees = new Regtree[numtrees];
        this.minVariance = buildParams.minVariance;
        this.buildParams = buildParams;
    }
    
    
    
    public boolean equals(Object o)
    {
    	
    	if(o instanceof RandomForest)
    	{
    		RandomForest rf = (RandomForest) o;
    		
    		if(numTrees != rf.numTrees) return false;
    		if(logModel != rf.logModel) return false;
    		return Arrays.equals(Trees, rf.Trees);    		
    		
    	} 
    	return false;
    }
    
    
    public int hashCode()
    {
    	return logModel ^ 2*numTrees ^ Arrays.deepHashCode(Trees);
    }
    
    public int matlabHashCode()
    {
    	return Math.abs(hashCode()) % 32452867;
    }
    /**
     *
     *
     * N - Number of configurations 
     * K - Number of parameters for a configuration 
     * M - Number of instances
     * L - Number of features for an instance
     * P - Number of Runs preformed
     *
     * @param numTrees - Number of trees in the Random Forest
     * @param allTheta - N x K matrix of parameter values [ Each row is a configuration, each entry in a row represents the value for that parameter (in that configuration)]
     * @param allX - M x L matrix of instance features [ Each row is all the features for a single instance, each entry in a row represents the value of that feature].
     * @param theta_inst_idxs - P x 2 - Pairs of indexes ( A,B) where A points to a row in allTheta, and B points to a row in allX [ So an configuration, instance pair]. 
     * @param y - P x 1 - Response values for the corresponding pairs in theta_inst_idxes.
     * @param cens - P x 1 - (If it's false the run was not capped, if true, then y is a lower bound of the run time).
     * @param params - 
     *
     */
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, RegtreeBuildParams params) {
        Random r = params.random;
        if (r == null) {
            r = new Random();
            if (params.seed != -1) {
                r.setSeed(params.seed);
            }
        }        
        
        
        
        int N = y.length;
        
        // Do bootstrap sampling for data for each tree.
        int[][] dataIdxs = new int[numTrees][N];
        for (int i = 0; i < numTrees; i++) {
            for (int j = 0; j < N; j++) {
                dataIdxs[i][j] = r.nextInt(N);
            }
        }
        return learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, dataIdxs, params);
    }
    
    public static void fixInputs(double[] input)
    {
    	for(int j=0; j < input.length; j++)
		{
			
			double val = input[j];
			long raw = Double.doubleToLongBits(val);
			val = Double.longBitsToDouble(raw);
			input[j] = (float) val;
		}
    }
    
    public static void fixInputs(double[][] input)
    {
    	for(int i=0; i < input.length; i++)
    	{
    		fixInputs(input[i]);
    	}
    }
    
    
    public static String print(Object o)
    {
    	
		
    	if(o instanceof Integer) 
    	{
    		return o.toString();
    	}
    	
    	if(o instanceof double[])
    	{
    		return Arrays.toString((double[])o);
    	}
    	
    	if(o instanceof int[][])
    	{
    		return Arrays.deepToString((int[][])o);
    	}
    	
    	if(o instanceof double[][])
    	{
    		return Arrays.deepToString((double[][]) o);
    	}
    	
    	if(o instanceof Object)
    	{
    		return o.toString(); 
    	}
    	
    	System.out.println(o);
    	throw new IllegalStateException();
    }
    
    public static boolean equalTest(Object o, Object o2)
    {
    	
		
    	if(o instanceof Integer) 
    	{
    		return o.equals(o2);
    	}
    	
    	if(o instanceof double[])
    	{
    		return Arrays.equals((double[])o,(double[])o2);
    	}
    	
    	if(o instanceof int[][])
    	{
    		return Arrays.deepEquals((int[][])o,(int[][])o2);
    	}
    	
    	if(o instanceof double[][])
    	{
    		return Arrays.deepEquals((double[][]) o,(double[][]) o2);
    	}
    	
    	if(o instanceof Object)
    	{
    		return o.equals(o2); 
    	}
    	
    	System.out.println(o);
    	throw new IllegalStateException();
    }
    
    public static void main3(String[] args)
    {
    	RandomForest matlabForest = fromForestFile("/tmp/RandomForest4433899701602217560Build");
    	RandomForest javaForest = fromForestFile("/tmp/RandomForest8044841660237959103Build");
    	
		
    	
    	
	double[][] configs = {{2.0, 3.0, 2.0, 2.0, 4.0, 12.0, 5.0, 3.0, 5.0, 3.0, 19.0, 6.0, 4.0, 5.0, 2.0, 1.0, 1.0, 3.0, 16.0, 3.0, 5.0, 2.0, 1.0, 2.0, 7.0, 2.0}};
	int[] treesUsed = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	
	System.out.println(matlabForest.equals(javaForest));
	System.out.println(Arrays.deepToString(RandomForest.applyRoundedMarginal(matlabForest, treesUsed, configs)));
	System.out.println(Arrays.deepToString(RandomForest.applyRoundedMarginal(javaForest, treesUsed, configs)));
	
	
		
	}
    	
    
    public static void main2(String[] args)
	{
    	RandomForest f1 = fromFile("/tmp/RandomForest7400514883271338896Build");
    	RandomForest f2 = fromFile("/tmp/RandomForest7510701929783488816Build");
    	
    	System.out.println(f1.matlabHashCode());
    	System.out.println("---------");
    	System.out.println(f2.matlabHashCode());
    	
    	Object[][] arr2 = arr;
    	System.out.println("f1var:"+Arrays.toString(f1.Trees[2].var));
    	System.out.println("f2var:"+Arrays.toString(f2.Trees[2].var));
    	System.out.println(Arrays.equals(f1.Trees[2].var, f2.Trees[2].var));
    	System.out.println("f1catsplit:"+Arrays.deepToString(f1.Trees[2].catsplit));
    	System.out.println("f2catsplit:"+Arrays.deepToString(f2.Trees[2].catsplit));
    	System.out.println(Arrays.deepEquals(f1.Trees[2].catsplit, f2.Trees[2].catsplit));
    	
    	
    	System.out.println("f1cut:"+Arrays.toString(f1.Trees[2].cut));
    	System.out.println("f2cut:"+Arrays.toString(f2.Trees[2].cut));
    	System.out.println(Arrays.equals(f1.Trees[2].cut, f2.Trees[2].cut));
    	
    	
    	for(int i=0; i < arr[1].length; i++)
    	{
    		System.out.println("a1:(" + i+")" + print(arr[0][i]));
    		System.out.println("a2:(" + i+")" + print(arr[1][i]));
    		System.out.println("=:(" + i + ")" + equalTest(arr[0][i],arr[1][i]));
    		System.out.println("");
    		
    	}
    	System.out.println(f1.equals(f2));
    	
		
    	
    	
	double[][] configs = {{2.0, 3.0, 2.0, 2.0, 4.0, 12.0, 5.0, 3.0, 5.0, 3.0, 19.0, 6.0, 4.0, 5.0, 2.0, 1.0, 1.0, 3.0, 16.0, 3.0, 5.0, 2.0, 1.0, 2.0, 7.0, 2.0}};
	int[] treesUsed = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	
	System.out.println(Arrays.deepToString(RandomForest.applyMarginal(f1, treesUsed, configs)));
	System.out.println(Arrays.deepToString(RandomForest.applyMarginal(f2, treesUsed, configs)));
	
	
		
	}
    private static Object[][] arr = new Object[2][];
    
    
    public static RandomForest fromFile(String s)
    {
    	File f = new File(s);
		ObjectInputStream in;
		try {
		in = new ObjectInputStream(new FileInputStream(f));
		int numTrees = (int) in.readInt();
		double[][] allTheta  = (double[][]) in.readObject();
		double[][] allX = (double[][]) in.readObject();
		int[][] theta_inst_idxs = (int [][]) in.readObject();
		double[] y = (double[]) in.readObject();
		int[][] dataIdxs = (int [][]) in.readObject();
		RegtreeBuildParams params = (RegtreeBuildParams) in.readObject();
		
		in.close();
	
		fixInputs(allTheta);
		fixInputs(y);
		fixInputs(allX);
		Object[] obj = {numTrees, allTheta,allX, theta_inst_idxs, y, dataIdxs, params};
		
		
		if(arr[0] == null)
		{
			arr[0] = obj;
		} else
		{
			arr[1] = obj;
		}
		
		
			
		return learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, dataIdxs, params);
		}catch(Exception e)
		{
			throw new RuntimeException(e);
		}
    }
    
    public static RandomForest fromForestFile(String s)
    {
    	File f = new File(s);
		ObjectInputStream in;
		try {
		in = new ObjectInputStream(new FileInputStream(f));
		RandomForest forest= (RandomForest) in.readObject();
				
		in.close();
		
		return forest;

		}catch(Exception e)
		{
			throw new RuntimeException(e);
		}
    }
    
    
    public static void main(String[] args)
    {
    	String badParams = "/tmp/RandomForestParams3027553244461234091Build";
    	String goodParams = "/tmp/RandomForestParams4017574993756479929Build";
    	
    	File f1 = new File(goodParams);
    	File f2 = new File(badParams);
		ObjectInputStream in,in2;
		try {
		in = new ObjectInputStream(new FileInputStream(f1));
		in2 = new ObjectInputStream(new FileInputStream(f2));
		
		int numTrees = (int) in.readInt();
		int numTrees2 = (int) in2.readInt();
		System.out.println(numTrees - numTrees2);
		
		double[][] allTheta  = (double[][]) in.readObject();
		
		double[][] allTheta2  = (double[][]) in2.readObject();
		
		for(int i=0; i < allTheta.length; i++)
		{
			for(int j=0; j < allTheta[i].length; j++)
			{
				if(allTheta[i][j] - allTheta2[i][j] != 0)
				{
				System.out.println(i + "," + j + ":" + allTheta[i][j] + " " + allTheta2[i][j]);
				
				
				}
			}
		}
		
		System.out.println(Arrays.deepEquals(allTheta, allTheta2));
		System.out.println("a");
		
		double[][] allX = (double[][]) in.readObject();
		double[][] allX2 = (double[][]) in2.readObject();
		System.out.println(Arrays.deepEquals(allX, allX2));
		
		System.out.println("b");
		
		int[][] theta_inst_idxs = (int [][]) in.readObject();
		int[][] theta_inst_idxs2 = (int [][]) in2.readObject();
		System.out.println(Arrays.deepEquals(theta_inst_idxs, theta_inst_idxs2));
		System.out.println("c");
		
		
		double[] y = (double[]) in.readObject();
		double[] y2 = (double[]) in2.readObject();
		System.out.println(Arrays.equals(y, y2));
		
	
		System.out.println("d");
		
		int[][] dataIdxs = (int [][]) in.readObject();
		int[][] dataIdxs2 = (int [][]) in2.readObject();
		System.out.println(Arrays.deepEquals(dataIdxs, dataIdxs2));
		System.out.println("e");
		
		RegtreeBuildParams params = (RegtreeBuildParams) in.readObject();
		
		in.close();
		in2.close();
		} catch(Exception e)
		{
			
		}
		
		
    }
    public static void save(RandomForest forest)
    {
    	
    		
    		
    		try {
    		File f = File.createTempFile("RandomForest", "Build");
    		ObjectOutputStream o = new ObjectOutputStream(new FileOutputStream(f));
    		
    		
    		o.writeObject(forest);

    		
    		
    		System.out.println("Forest Saved To:" + f.getAbsolutePath() );
    		o.close();
    		} catch(IOException e)
    		{
    			System.err.println(e);
    		}
    }
    
    /**
     * Learns a random forest putting the specified data points into each tree.
     * @see RegtreeFit.fit
     */
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, int[][] dataIdxs, RegtreeBuildParams params) {
    	/*
    	fixInputs(allTheta);
		fixInputs(y);
		fixInputs(allX);
        System.out.println("New Hash Codes");
        System.out.println(RandStateHash.randHash(params.random));
    	System.out.println(numTrees);
        System.out.println(hash(allTheta));
        System.out.println(hash(allX));
        System.out.println(hash(theta_inst_idxs));
        System.out.println(hash(y));
        System.out.println(hash(dataIdxs));
        System.out.println(hash(allTheta));
        */
//        System.out.println(Arrays.deepToString(allX));
       //System.out.println(Arrays.deepToString(theta_inst_idxs));
//        System.out.println(Arrays.toString(y));
//        System.out.println(Arrays.deepToString(dataIdxs));
       /**
        * TODO Add validaiton for index errors
        */
    	
   	boolean writeOutput = false;
	if(writeOutput)
	{
		
		
		try {
		File f = File.createTempFile("RandomForestParams", "Build");
		ObjectOutputStream o = new ObjectOutputStream(new FileOutputStream(f));
		
		o.writeInt(numTrees);
		o.writeObject(allTheta);
		o.writeObject(allX);
		o.writeObject(theta_inst_idxs);
		o.writeObject(y);
		o.writeObject(dataIdxs);
		o.writeObject(params);
		
		
		System.out.println("Calls written & deleted to:" + f.getAbsolutePath() );
		o.close();
		} catch(IOException e)
		{
			System.err.println(e);
		}
		
	}
	
        if (dataIdxs.length != numTrees) {
            throw new RuntimeException("length(dataIdxs) must be equal to numtrees.");
        }
    
        
        
        RandomForest rf = new RandomForest(numTrees, params);
        for (int i = 0; i < numTrees; i++) {
            int N = dataIdxs[i].length;
            int[][] this_theta_inst_idxs = new int[N][];
            double[] thisy = new double[N];
            for (int j=0; j<N; j++) {
                int idx = dataIdxs[i][j];
                this_theta_inst_idxs[j] = theta_inst_idxs[idx];
                thisy[j] = y[idx];
            }
            rf.Trees[i] = RegtreeFit.fit(allTheta, allX, this_theta_inst_idxs, thisy, params);
        }
        return rf;
    }
    

    
    /**
     * Propogates data points down the regtree and returns a numtrees*X.length vector of node #s
     * specifying which node each data point falls into.
     * @params X a numdatapoints*numvars matrix
     */
    public static int[][] fwd(RandomForest forest, double[][] X) {
        int[][] retn = new int[forest.numTrees][X.length];
        for (int i=0; i < forest.numTrees; i++) {
            int[] result = RegtreeFwd.fwd(forest.Trees[i], X);
            System.arraycopy(result, 0, retn[i], 0, result.length);
        }
        return retn;
    }
    
    /**
     * Gets a prediction for the given instantiations of features
     * @returns a matrix of size X.length*2 where index (i,0) is the prediction for X[i] 
     * and (i,1) is the variance of that prediction. See Matlab code for how var is calculated.
     */
    public static double[][] apply(RandomForest forest, double[][] X) {
		double[][] retn = new double[X.length][2]; // mean, var
        for (int i=0; i < forest.numTrees; i++) {
            int[] result = RegtreeFwd.fwd(forest.Trees[i], X);
            for (int j=0; j < X.length; j++) {
				double pred = forest.Trees[i].nodepred[result[j]];
				double var = forest.Trees[i].nodevar[result[j]];
				/*
                if (forest.logModel>0) {
                    pred = Math.log10(pred);
                    System.out.println("Variance for Tree[  " + i + "] is " + var);
                }*/
                retn[j][0] += pred;
                retn[j][1] += var+pred*pred;
            }
        }
        for (int i=0; i < X.length; i++) {
            retn[i][0] /= forest.numTrees;
            retn[i][1] /= forest.numTrees;
            retn[i][1] -= retn[i][0]*retn[i][0];
            retn[i][1] = retn[i][1] * ((forest.numTrees+0.0)/Math.max(1, forest.numTrees-1));
            
            if(forest.logModel > 0)
            {
	            /**
	             * Take mean and variance in non-log space, transform them into ln space, then Linearly Transform them to log-10 space. 
	             * Get the parameters of the log normal distribution, with that mean and variance (in non-log space).
	             */
	            double test_mu_n = retn[i][0];
	            double test_var_n = retn[i][1];
	            
	            double var_ln = Math.log(test_var_n/(test_mu_n*test_mu_n) + 1);
	            double	mu_ln = Math.log(test_mu_n) - var_ln/2;
	            
	            double var_l10 = var_ln / Math.log(10) / Math.log(10);
	            double mu_l10 = mu_ln / Math.log(10); 
	            
	            retn[i][0] = mu_l10;
	            retn[i][1] = var_l10;
            }
            
            if(retn[i][1] < MIN_VARIANCE_RESULT)
            {
            	System.err.println("[WARN]: Variance is less than " + MIN_VARIANCE_RESULT + " > " + retn[i][1]);
            	assert(retn[i][1] > MIN_VARIANCE_RESULT); //Assert negative variance only comes from numerical issues (and they shouldn't make it too small)
            }
            retn[i][1] = Math.max(forest.minVariance, retn[i][1]);
        
            
            
        }
        
        
        return retn;
    }
    
    public static double round(double val)
    {
    	float fval = (float) val;
    	
    	int bits = Float.floatToRawIntBits(fval);
    	bits = bits & 0xFFFFF800;
    	float nval = Float.intBitsToFloat(bits);
    	
    	//System.out.println(fval + "=>" + nval);
    	return nval;
    }
    

    public static double[][] applyRoundedMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta) {
    	double[][] results = applyMarginal(forest, tree_idxs_used, Theta);
    	for(int i=0; i < results.length; i++)
    		for(int j=0; j < results[i].length; j++)
    		{
    			double result = results[i][j];
    			
    			results[i][j] = round(result);
    		}
        return results;
    }
    
    public static double[][] applyRoundedMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta, double[][] X) {
    	double[][] results = applyMarginal(forest, tree_idxs_used, Theta, X);
    	for(int i=0; i < results.length; i++)
    		for(int j=0; j < results[i].length; j++)
    		{
    			double result = results[i][j];
    			
    			results[i][j] = round(result);
    		}
        return results;
    }
    
    
    
    public static double[][] applyMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta) {
    	return applyMarginal(forest, tree_idxs_used, Theta, null);
    }
    
    /**
     * Gets a prediction for the given configurations and instances
     * @returns a matrix of size Theta.length*2 where index (i,0) is the prediction for Theta[i] 
     * and (i,1) is the variance of that prediction. See Matlab code for how var is calculated.
     * @see RegtreeFwd.marginalFwd
     */
    public static double[][] applyMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta, double[][] X) {
        int nTheta = Theta.length, nTrees = tree_idxs_used.length;
		double[][] retn = new double[nTheta][2]; // mean, var
        
        for (int i=0; i < nTrees; i++) {
            Object[] result = RegtreeFwd.marginalFwd(forest.Trees[tree_idxs_used[i]], Theta, X);
            double[] preds = (double[])result[0];
            double[] vars = (double[])result[1];

            for (int j=0; j < nTheta; j++) {
                double pred = preds[j];
                if (forest.logModel>0) {
                    pred = Math.log10(pred);
                }
                retn[j][0] += pred;
                retn[j][1] += vars[j]+pred*pred;
            }
        }
        
        for (int i=0; i < nTheta; i++) {
            retn[i][0] /= nTrees;
            retn[i][1] /= nTrees;
			retn[i][1] -= retn[i][0]*retn[i][0];
            retn[i][1] = retn[i][1] * ((forest.numTrees+0.0)/Math.max(1, forest.numTrees-1));
            
            if(retn[i][1] < MIN_VARIANCE_RESULT) 
            {
            	System.err.println("[WARN]: Variance is less than " + MIN_VARIANCE_RESULT + " > " + retn[i][1]);
            	assert(retn[i][1] > MIN_VARIANCE_RESULT); //Assert negative variance only comes from numerical issues (and they shouldn't make it too small)
            }
            retn[i][1] = Math.max(forest.minVariance, retn[i][1]);
        }
        return retn;
    }
    
    /** 
     * Prepares the random forest for marginal predictions.
     * @see RegtreeFwd.preprocess_inst_splits
     */
    public static RandomForest preprocessForest(RandomForest forest, double[][] X) {
        RandomForest prepared = new RandomForest(forest.numTrees,forest.buildParams);
        for (int i=0; i < forest.numTrees; i++) {
            prepared.Trees[i] = RegtreeFwd.preprocess_inst_splits(forest.Trees[i], X);
        }
        return prepared;
    }
}
