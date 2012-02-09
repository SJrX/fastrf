package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import ca.ubc.cs.beta.models.fastrf.utils.Utils;

public class RandomForest implements java.io.Serializable {
    private static final long serialVersionUID = 5204746081208095705L;
    
    public int numTrees;
    public Regtree[] Trees;
    public int logModel;
    
    public RandomForest(int numtrees, int logModel) {
        if (numtrees <= 0) {
            throw new RuntimeException("Invalid number of regression trees in forest: " + numtrees);
        }
        this.logModel = logModel;
        numTrees = numtrees;
        Trees = new Regtree[numtrees];
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
    /**
     * Learns a random forest putting the specified data points into each tree.
     * @see RegtreeFit.fit
     */
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, int[][] dataIdxs, RegtreeBuildParams params) {       
//        System.out.println(numTrees);
//        System.out.println(Arrays.deepToString(allTheta));
//        System.out.println(Arrays.deepToString(allX));
//        System.out.println(Arrays.deepToString(theta_inst_idxs));
//        System.out.println(Arrays.toString(y));
//        System.out.println(Arrays.deepToString(dataIdxs));
        if (dataIdxs.length != numTrees) {
            throw new RuntimeException("length(dataIdxs) must be equal to numtrees.");
        }
        
        RandomForest rf = new RandomForest(numTrees, params.logModel);
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
                if (forest.logModel>0) {
                    pred = Math.log10(pred);
                }
                retn[j][0] += pred;
                retn[j][1] += forest.Trees[i].nodevar[result[j]]+pred*pred;
            }
        }
        for (int i=0; i < X.length; i++) {
            retn[i][0] /= forest.numTrees;
            retn[i][1] /= forest.numTrees;
            retn[i][1] -= retn[i][0]*retn[i][0];
            retn[i][1] = retn[i][1] * ((forest.numTrees+0.0)/Math.max(1, forest.numTrees-1));
        }
        return retn;
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
        }
        return retn;
    }
    
    /** 
     * Prepares the random forest for marginal predictions.
     * @see RegtreeFwd.preprocess_inst_splits
     */
    public static RandomForest preprocessForest(RandomForest forest, double[][] X) {
        RandomForest prepared = new RandomForest(forest.numTrees, forest.logModel);
        for (int i=0; i < forest.numTrees; i++) {
            prepared.Trees[i] = RegtreeFwd.preprocess_inst_splits(forest.Trees[i], X);
        }
        return prepared;
    }
}
