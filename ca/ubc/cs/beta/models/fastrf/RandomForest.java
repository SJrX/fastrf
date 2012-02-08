package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import ca.ubc.cs.beta.models.fastrf.utils.Utils;

public class RandomForest implements java.io.Serializable {
    private static final long serialVersionUID = 5204746081208095703L;
    
    public int numTrees;
    public Regtree[] Trees;
    
    public RandomForest(int numtrees) {
        if (numtrees <= 0) {
            throw new RuntimeException("Invalid number of regression trees in forest: " + numtrees);
        }
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
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, boolean[] cens, RegtreeBuildParams params) {
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
        return learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, cens, dataIdxs, params);
    }
    /**
     * Learns a random forest putting the specified data points into each tree.
     * @see RegtreeFit.fit
     */
    public static RandomForest learnModel(int numTrees, double[][] allTheta, double[][] allX, int[][] theta_inst_idxs, double[] y, boolean[] cens, int[][] dataIdxs, RegtreeBuildParams params) {       

        System.out.println(numTrees);
        System.out.println(Arrays.deepToString(allTheta));
        System.out.println(Arrays.deepToString(allX));
        System.out.println(Arrays.deepToString(theta_inst_idxs));
        System.out.println(Arrays.toString(y));
        System.out.println(Arrays.deepToString(dataIdxs));
        if (true) throw new RuntimeException("YAY!");
        
        
        if (dataIdxs.length != numTrees) {
            throw new RuntimeException("length(dataIdxs) must be equal to numtrees.");
        }
        
		for (int i = 0; i < cens.length; i++) {
			if (cens[i]) throw new RuntimeException("Censored data is not supported!");
		}
          
        RandomForest rf = new RandomForest(numTrees);
        for (int i = 0; i < numTrees; i++) {
            int N = dataIdxs[i].length;
            int[][] this_theta_inst_idxs = new int[N][];
            double[] thisy = new double[N];
            boolean[] thiscens = new boolean[N];
            for (int j=0; j<N; j++) {
                int idx = dataIdxs[i][j];
                this_theta_inst_idxs[j] = theta_inst_idxs[idx];
                thisy[j] = y[idx];
                thiscens[j] = cens[idx];
            }
            rf.Trees[i] = RegtreeFit.fit(allTheta, allX, this_theta_inst_idxs, thisy, thiscens, params);
        }
		
        return rf;
    }
    
    /**
     * Propogates data points down the regtree and returns a numtrees*X.length vector of node #s
     * specifying which node each data point falls into.
     * @params X a numdatapoints*numvars [numvars is concatenation of K parameter values and L feature values] matrix
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
     * Propogates configurations down the regtree and returns a numtrees*2 matrix containing 
     * leafIdxs and thetaIdxs.
     * @see RegtreeFwd.fwdThetas
     */
    public static Object[][] fwdThetas(RandomForest forest, double[][] Theta) {
        Object[][] retn = new Object[forest.numTrees][2];
        for (int i=0; i < forest.numTrees; i++) {
            Object[] result = RegtreeFwd.fwdThetas(forest.Trees[i], Theta);
            int[] leafIdxs = (int[]) result[0];
            int[][] ThetaIdxs = (int[][]) result[1];
            
            retn[i][0] = leafIdxs;
            retn[i][1] = ThetaIdxs;
        }
        return retn;
    }
    
    /**
     * Gets a prediction for the given instantiations of features
     * @returns a matrix of size X.length*2 where index (i,0) is the prediction for X[i] 
     * and (i,1) is the variance of that prediction.
     */
    public static double[][] apply(RandomForest forest, double[][] X) {
		double[][] retn = new double[X.length][2]; // mean, var
        for (int i=0; i < forest.numTrees; i++) {
            int[] result = RegtreeFwd.fwd(forest.Trees[i], X);
            for (int j=0; j < X.length; j++) {
				double pred = forest.Trees[i].nodepred[result[j]];
                retn[j][0] += pred;
                retn[j][1] += forest.Trees[i].nodevar[result[j]]+pred*pred;
            }
        }
        for (int i=0; i < X.length; i++) {
            retn[i][0] /= forest.numTrees;
            retn[i][1] /= forest.numTrees;
            retn[i][1] -= retn[i][0]*retn[i][0];
        }
        return retn;
    }
    
    public static double[][] applyMarginal(RandomForest forest, int[] tree_idxs_used, double[][] Theta) {
        return applyMarginal(forest, tree_idxs_used, Theta, null);
    }
    
    /**
     * Gets a prediction for the given configurations and instances
     * @returns a matrix of size Theta.length*2 where index (i,0) is the prediction for Theta[i] 
     * and (i,1) is the variance of that prediction.
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
                retn[j][0] += preds[j];
                retn[j][1] += vars[j]+preds[j]*preds[j];
            }
        }
        
        for (int i=0; i < nTheta; i++) {
            retn[i][0] /= nTrees;
            retn[i][1] /= nTrees;
			retn[i][1] -= retn[i][0]*retn[i][0];
        }
        return retn;
    }
    
    /**
     * Marginalizes over toBeMarginalized and returns a marginal prediction for each row of X.
     * Different from applyMarginalModel in that toBeMarginalized is not a list of data points but we use the cartesian product of the rows of toBeMarginalized.
     * @params X a N*d matrix (N d-variable sets) to get the marginal for
     * @params toBeMarginalized a cell array with p rows where each row specifies the values for a variable to marginalize over.
     * @params variableIndicesForColumnsOfX a d*1 vector of the variable numbers each column of X represents. 1-indexed.
       d+p = numvars. The columns of toBeMarginalized represent the variables {0:numvars-1}-variableIndicesForColumnsOfX, in sorted increasing order.
     * @returns a N*1 vector of marginalized predictions.
     */
    public static double[] getMarginal(RandomForest forest, double[][] X, double[][] toBeMarginalized, int[] variableIndicesForColumnsOfX) {  
		// TODO: lognormal        
		if (X == null) {
            X = new double[1][0];
        }
        if (toBeMarginalized == null) {
            toBeMarginalized = new double[0][];
        }
        if (variableIndicesForColumnsOfX == null) {
            variableIndicesForColumnsOfX = new int[0];
        }
        if (X[0].length + toBeMarginalized.length != forest.Trees[0].npred) {
            throw new RuntimeException("d+p must be equal to numvars.");
        }
        if (X[0].length != variableIndicesForColumnsOfX.length) {
            throw new RuntimeException("The number of columns of X and the length of variableIndicesForColumnsOfX must match up.");
        }
        
        // Fill in auxillary arrays
        int[] toBeMarginalizedVariables = new int[forest.Trees[0].npred-variableIndicesForColumnsOfX.length];
        int[] isXvar = new int[forest.Trees[0].npred+1];
        int counter = 0, next = 1;
        for (int i=0; i < variableIndicesForColumnsOfX.length; i++) {
            while (variableIndicesForColumnsOfX[i] != next) {
                toBeMarginalizedVariables[counter++] = next++;
            }
            isXvar[variableIndicesForColumnsOfX[i]] = i+1;
            next++;
        }
        while(next <= forest.Trees[0].npred) {
            toBeMarginalizedVariables[counter++] = next++;
        }
        
        // marginalized means for each X in each tree.
        double[][] treeMeans = new double[forest.numTrees][X == null ? 1 : X.length];
        
        LinkedList<Integer> queue = new LinkedList<Integer>();
        
        // Loop through each tree, for each tree calculate treeMeans
        for (int m=0; m < forest.numTrees; m++) {
            Regtree tree = forest.Trees[m];
        
            // element i of result is # of entries of cartesian product of toBeMarginalized that fall into node i.
            double[] results = new double[tree.numNodes];
            Arrays.fill(results, 1);
            
            // For each row of toBeMarginalized, forward it and for each leaf store # of elements that fall into it.
            for (int i=0; i < toBeMarginalized.length; i++) {
                int nextvar = toBeMarginalizedVariables[i];

                // counts[j] is # of times some element from row i falls into node j.
                int[] counts = new int[tree.numNodes];
                int numValues = toBeMarginalized[i].length;
                
                for (int j=0; j < numValues; j++) {
                    queue.add(0);
                    while(!queue.isEmpty()) {
                        int thisnode = queue.poll();
                        while(true) {
                            int splitvar = tree.var[thisnode];
                            double cutoff = tree.cut[thisnode];
                            int left_kid = tree.children[thisnode][0];
                            int right_kid = tree.children[thisnode][1];

                            if (splitvar == 0) {
                                // We are in leaf node. store results.
                                counts[thisnode]++;
                                break;
                            } else if (Math.abs(splitvar) != nextvar) {
                                // Splitting on different var - pass this down both children
                                queue.add(right_kid);
                                thisnode = left_kid;
                            } else {
                                if (splitvar > 0) { // continuous
                                    thisnode = (toBeMarginalized[i][j] <= cutoff ? left_kid : right_kid);
                                } else { // categorical
                                    int x = (int)toBeMarginalized[i][j];
                                    int split = tree.catsplit[(int)cutoff][x-1];
                                    if (split == 0) thisnode = left_kid;
                                    else if (split == 1) thisnode = right_kid;
                                    else throw new RuntimeException("Missing value -- not allowed in this implementation.");
                                }
                            }
                        }
                    }
                }
                
                // Multiply with existing results (the cardinality of the cartesian product of the actual values so far)
                for (int j=0; j < tree.numNodes; j++) {
                    results[j] *= counts[j]*1.0/numValues;
                }
            }
            
            //=== Calculate marginal means for X.
            if (X[0].length==0) {
                // Special case if X is empty.
                double sum = 0;
                for (int i=0; i < results.length; i++) sum += results[i];
                if (sum - 1 > 1e-6) {
                    throw new RuntimeException("Something is wrong. Sum: " + sum + "(" + Arrays.toString(results) + ")");
                }
                for (int i=0; i < results.length; i++) {
                    treeMeans[m][0] += results[i]*tree.nodepred[i];
                }
            } else {
                // For each row of X, get the marginal mean from this tree.
                for (int i=0; i < X.length; i++) {
                    queue.add(0);
                    while(!queue.isEmpty()) {
                        int thisnode = queue.poll();
                        while(true) {
                            int splitvar = tree.var[thisnode];
                            double cutoff = tree.cut[thisnode];
                            int left_kid = tree.children[thisnode][0];
                            int right_kid = tree.children[thisnode][1];
                            int varidx = isXvar[Math.abs(splitvar)];

                            if (splitvar == 0) {
                                // We are in leaf node. store results.
                                treeMeans[m][i] += results[thisnode]*tree.nodepred[thisnode];
                                break;
                            } else if (varidx == 0) {
                                // Splitting on different var - pass this down both children
                                queue.add(right_kid);
                                thisnode = left_kid;
                            } else {
                                if (splitvar > 0) { // continuous
                                    thisnode = (X[i][varidx-1] <= cutoff ? left_kid : right_kid);
                                } else { // categorical
                                    int x = (int)X[i][varidx-1];
                                    int split = tree.catsplit[(int)cutoff][x-1];
                                    if (split == 0) thisnode = left_kid;
                                    else if (split == 1) thisnode = right_kid;
                                    else throw new RuntimeException("Missing value -- not allowed in this implementation.");
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // calculate the mean of treeMeans
        double[] retn = new double[X.length];
        for (int i=0; i < retn.length; i++) {
            double sum = 0;
            for (int m=0; m < forest.numTrees; m++) sum += treeMeans[m][i];
            retn[i] = sum * 1.0 / forest.numTrees;
        }
        return retn;
    }
    
    /** 
     * Prepares the random forest for marginal predictions.
     * @see RegtreeFwd.preprocess_inst_splits
     */
    public static RandomForest preprocessForest(RandomForest forest, double[][] X) {
        RandomForest prepared = new RandomForest(forest.numTrees);
        for (int i=0; i < forest.numTrees; i++) {
            prepared.Trees[i] = RegtreeFwd.preprocess_inst_splits(forest.Trees[i], X);
        }
        return prepared;
    }
    
    /**
     * Collect response values from leaves for the given data points.
     * @returns a 3*X.length*totalNumberOfDataPointsInLeaves matrix of y_pred, cens_pred, and weights.
     */
    public static Object[] collectData(RandomForest forest, double[][] X) {
        int[][] leafNodes = fwd(forest, X);
        return collectData(forest, leafNodes);
    }
    
    /**
     * Collect response values from the given leaves.
     * @paramss leafNodes a numTrees*X.length matrix of the nodes of the data points. (eg. from fwd)
     * @returns a 3*X.length*totalNumberOfDataPointsInLeaves matrix of y_pred, cens_pred, and weights.
     */
    public static Object[] collectData(RandomForest forest, int[][] leafNodes) {
        if (!forest.Trees[0].resultsStoredInLeaves) {
            throw new RuntimeException("Cannot collect data if they were not stored.");
        }
        Object[] retn = new Object[3];

        int numdata = leafNodes[0].length;
        
        double[][] y_pred_all = new double[numdata][];
        boolean[][] cens_pred_all = new boolean[numdata][];
        double[][] weights_all = new double[numdata][];
        
        for (int i=0; i < numdata; i++) {
            int size = 0;
            for (int j=0; j < forest.numTrees; j++) {
                size += forest.Trees[j].nodesize[leafNodes[j][i]];
            }
            
            double[] y_pred = new double[size];
            boolean[] cens_pred = new boolean[size];
            double[] weights = new double[size];
            
            int counter = 0;
            for (int j=0; j < forest.numTrees; j++) {
                Regtree tree = forest.Trees[j];
                int node = leafNodes[j][i];
                int nodesize = tree.nodesize[node];
                for (int k=0; k < nodesize; k++) {
                    y_pred[counter+k] = tree.ysub[node][k];
                    cens_pred[counter+k] = tree.is_censored[node][k];
                    if (!cens_pred[counter+k]) y_pred[counter+k] += 1e-10;
                }
                Arrays.fill(weights, counter, counter+nodesize, 1.0/nodesize);
                counter += nodesize;
            }
            
            y_pred_all[i] = y_pred;
            cens_pred_all[i] = cens_pred;
            weights_all[i] = weights;
        }
        retn[0] = y_pred_all;
        retn[1] = cens_pred_all;
        retn[2] = weights_all;
        return retn;
    }
}
