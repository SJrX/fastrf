package ca.ubc.cs.beta.models.fastrf;

import java.util.*;
import ca.ubc.cs.beta.models.fastrf.utils.*;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;

public class RegtreeFit {
    
    private static Random r;
    private static long seed;
    //*
    private static final int RAND_MAX = Integer.MAX_VALUE - 1;
    private static int rand() {
        int retn = r.nextInt(Integer.MAX_VALUE);
        return retn;
    }
    /*/
    private static final int RAND_MAX = 2147483646;
    private static int rand() {
        return (int)(seed = (seed*22695477+1)%(RAND_MAX+1));
    }
    //*/
    
    private static final double INVALID_CRITVAL = -1e13;
    
    public static Regtree fit(double[][] allTheta, double[][] allX, double[] y, boolean[] cens, RegtreeBuildParams params) {    
        Random r = params.random;
        if (r == null) {
            r = new Random();
            if (params.seed != -1) {
                r.setSeed(params.seed);
            }
        }        
        
        int N = y.length;
        int numTheta = allTheta.length;
        int numX = (allX == null ? 0 : allX.length);
        
        // Do bootstrap sampling for data for each tree.
        int[][] dataIdxs = new int[N][2];
        for (int i = 0; i < N; i++) {
            dataIdxs[i][0] = (numTheta == 0 ? 0 : r.nextInt(numTheta));
            dataIdxs[i][1] = (numX == 0 ? 0 : r.nextInt(numX));
        }
        return fit(allTheta, allX, dataIdxs, y, cens, params);
    }
    
    
    private static int[][] dataIdxs;
    private static double[] y;
    private static boolean[] cens;
    private static double kappa;
    

	private static double ybar;
	private static double[] catmeans;
	private static int[] catcounts;
	private static double[] ycum;
	private static int[] ycountcum;
	private static int[] uniqueIdxs;
	private static int[] dataRowsHere;
    
    private static double[] km_mean;
	private static int[] numUptoCategory;
	private static int[] sortedByCategory;
	private static double[] y_in;
	private static double[] c_in;
	private static int[] sorder;
	private static int[] maxlocs;
	
	private static int numleft;
	private static int numright;
	private static int[] leftside;
	private static int[] rightside;

	private static int[] period;
	private static double[] period_time;
	private static int[] N_all;
	private static int[] O_all;
	private static int[] C_all;
	private static int[] N_1;
	private static int[] O_1;
	private static int[] N_add;
	private static int num_periods;
    /**
     * Fits a regression tree.
     * @params allTheta, allX: matrices of all of the configurations/instances
     * @params dataIdxs row i is data point i with X=[allTheta[dataIdxs[i][1]], allX[dataIdxs[i][2]]] and y=y[i].
     * @params y a vector of size X.length of response values
     * @params cens a vector of size X.length indicating whether the corresponding response is right-censored.
     * @params params see RegtreeBuildParams
     */
    public static Regtree fit(double[][] allTheta, double[][] allX, int[][] dataIdxs, double[] y, boolean[] cens, RegtreeBuildParams params) {
    	boolean printDebug = false;
    	long startTime = new Date().getTime();
    	long currentTime = startTime;
    	
        
        if (dataIdxs == null || dataIdxs.length == 0) throw new RuntimeException("Cannot build a tree with no data.");
        int N = dataIdxs.length;
        if (y.length != N) throw new RuntimeException("The number of data points and the number of responses must be the same.");
        
        r = params.random;
        if (r == null) {
            r = new Random();
            if (params.seed != -1) {
                r.setSeed(params.seed);
            }
        }
        seed = params.seed;
        
        // Calculate input data dimensions
        int numTheta = (allTheta == null ? 0 : allTheta.length);
        int numX = (allX == null ? 0 : allX.length);
        int numThetavars = (allTheta == null ? 0 : allTheta[0].length);
        int numXvars = (allX == null ? 0 : allX[0].length);
        int nvars = numThetavars + numXvars;
        
        //=== Start: drop rows of allTheta and allX that we don't have data for.
        boolean[] hasThetaIdx = new boolean[numTheta+1];
        boolean[] hasXIdx = new boolean[numX+1];
        
        if (numTheta != 0) {
            for (int i=0; i < N; i++) {
                hasThetaIdx[dataIdxs[i][0]] = true;
            }
        }
        if (numX != 0) {
            for (int i=0; i < N; i++) {
                hasXIdx[dataIdxs[i][1]] = true;
            }
        }
        
        int[] numMissingThetaIdxsBeforeThis = new int[numTheta+1];
        int[] numMissingXIdxsBeforeThis = new int[numX+1];
        int numMissing = 0;
        for (int i=0; i <= numTheta; i++) {
            numMissingThetaIdxsBeforeThis[i] = numMissing;
            if (!hasThetaIdx[i]) numMissing++;            
        }
        numMissing = 0;
        for (int i=0; i <= numX; i++) {
            numMissingXIdxsBeforeThis[i] = numMissing;
            if (!hasXIdx[i]) numMissing++;
        }
        
        double[][] actualAllTheta = new double[numTheta - numMissingThetaIdxsBeforeThis[numTheta]][];
        for (int i=0, counter=0; i < numTheta; i++) {
            if (hasThetaIdx[i]) actualAllTheta[counter++] = allTheta[i];
        }
        allTheta = actualAllTheta;
        numTheta = allTheta.length;
        
        double[][] actualAllX = new double[numX - numMissingXIdxsBeforeThis[numX]][];
        for (int i=0, counter=0; i < numX; i++) {
            if (hasXIdx[i]) actualAllX[counter++] = allX[i];
        }
        allX = actualAllX;
        numX = allX.length;
        
        int[][] newDataIdxs = new int[N][2];
        for (int i=0; i < N; i++) {
            int thetaIdx = dataIdxs[i][0], xIdx = dataIdxs[i][1];
            if (numTheta != 0) {
                newDataIdxs[i][0] = thetaIdx - numMissingThetaIdxsBeforeThis[thetaIdx];
            } else {
                newDataIdxs[i][0] = 0;
            }
            if (numX != 0) {
                newDataIdxs[i][1] = xIdx - numMissingXIdxsBeforeThis[xIdx];            
            } else {
                newDataIdxs[i][1] = 0;
            }
        }
        dataIdxs = newDataIdxs;
        //=== End: drop rows of allTheta and allX that we don't have data for.
        
        RegtreeFit.dataIdxs = dataIdxs;
    	RegtreeFit.y = y;
    	RegtreeFit.cens = cens;
        
        //=== Extract data from the input params.
        int[] catDomainSizes = params.catDomainSizes;
        System.out.println("Params:" + params);
        System.out.println("Nvars" + nvars);
        System.out.println("AllX" + Arrays.deepToString(allX));
        System.out.println("AllTheta" + Arrays.deepToString(allTheta));
        System.out.println("dataIdxs" + Arrays.deepToString(dataIdxs));
        System.out.println("y" + Arrays.toString(y));
        System.out.println("cens:" + Arrays.toString(cens));
        
        if(true) throw new RuntimeException("Error");
        
        if (catDomainSizes != null && catDomainSizes.length != nvars) {
            throw new RuntimeException("catDomainSizes must be of the same length as size(X, 2)");
        }
        int maxDomSize = 0;
        for (int i=0; i < nvars; i++) {
            if (catDomainSizes[i] > maxDomSize) maxDomSize = catDomainSizes[i];
        }
        
        int[][] condParents = params.condParents;
        int[][][] condParentVals = params.condParentVals; 
        
        kappa = params.kappa;

        //=== Extract tuning parameters.
        double ratioFeatures = params.ratioFeatures;
        int splitMin = params.splitMin;

        boolean hascens = false;
        for (int i=0; i < N; i++) {
            if (cens[i]) {
                hascens = true;
                break;
			}
		}
        if (hascens) throw new RuntimeException("Censoring Unimplemented!");
    
        //========== Initialize stuffs. Special codes for root node. ====
        int[] nodenumber = new int[2*N];
        int[] nodesize = new int[2*N];
        
        int[] cutvar = new int[2*N];
        double[] cutpoint = new double[2*N];
        int[] leftchildren = new int[2*N];
        int[] rightchildren = new int[2*N];
        int[] parent = new int[2*N];
        double[][] ysub = new double[2*N][];
        boolean[][] censsub = new boolean[2*N][];
        
        int ncatsplit = 0;
        int[][] catsplit = new int[2*N][];
        
        int[] randomPermutation = new int[nvars];
        double[] variableValuesHere = new double[N];
        dataRowsHere = new int[Math.max(N, maxDomSize)];
        sortedByCategory = dataRowsHere;
        uniqueIdxs = new int[N];
        rows_km = uniqueIdxs; // share this memory with 2 variables
        catmeans = new double[maxDomSize];
        km_mean = catmeans; // share this memory with 2 variables
        catcounts = new int[maxDomSize];
        
        // For categorical splits
        int numBestLeft = 0;
        int numBestRight = 0;
        int[] bestLeft = new int[maxDomSize];
        int[] bestRight = new int[maxDomSize];
        leftside = new int[maxDomSize];
        rightside = new int[maxDomSize];     
        
        maxlocs = new int[Math.max(N, maxDomSize)-1];
        
        period = null;
        period_time = null;
        N_all = null;
		O_all = null;
		C_all = null;
		N_1 = null;
		O_1 = null;
		N_add = null;
        if (hascens) {
            period = new int[N];
            period_time = new double[N+1];
            N_all = new int[N];
            O_all = new int[N];
            C_all = new int[N];
            N_1 = new int[N];
            O_1 = new int[N];
            N_add = new int[N+1];
        } 
        
        ycum = new double[Math.max(N, maxDomSize)+1];
        y_in = ycum;
        c_in = variableValuesHere;
        
        ycountcum = new int[Math.max(N, maxDomSize)+1];
        numUptoCategory = ycountcum;
        
        // For passing data to children
        boolean[] yGoesLeft = new boolean[N];
        boolean[] primaryGoesLeft = new boolean[Math.max(numTheta, numX)];
        
        double ystd = Utils.var(y);        
        
        //=== Start: pre-sort each variable
        // The entries of sortedTheta and sortedX are indices into allTheta/allX.
        int[][] sortedTheta = new int[numThetavars][];
        int[][] sortedX = new int[numXvars][];
        int[] index_into_dataIdxs_here = new int[N];
        
        double[] temp = new double[Math.max(numTheta, numX)];
        for (int i=0; i < numThetavars; i++) {
            if (catDomainSizes[i] != 0) continue;
            for (int j=0; j < numTheta; j++) {
                temp[j] = allTheta[j][i];
            }
            sortedTheta[i] = new int[numTheta];
            rankSort(temp, numTheta, sortedTheta[i]);
        }
        for (int i=0; i < numXvars; i++) {
            if (catDomainSizes[i+numThetavars] != 0) continue;
            for (int j=0; j < numX; j++) {
                temp[j] = allX[j][i];
            }
            sortedX[i] = new int[numX];
            rankSort(temp, numX, sortedX[i]);
        }
        rankSort(y, N, index_into_dataIdxs_here); 
        temp = null;
        //=== End: pre-sort each variable

        //=== Start: initialize ynodeTheta and ynodeX for the root node.
        //=== For each node, y_node holds the node's indices into sortedY,  
        //=== ynodeTheta holds one array per row of allTheta, each of them holding indices into y_node (most of them empty)
        //=== e.g. index_into_dataIdxs_here has values for configurations [3,5,3,7]; then ynodeTheta[3-1] is [0,2], ynodeTheta[5-1] is [1], etc
        //=== Analogously for ynodeX and allX
        //===
        //=== The temporary Thetacount[i] holds how many times configuration i appears. E.g., in the example above, Thetacount[3-1]=2. Similarly for Xcount.  
        int[][] y_node = new int[2*N][];
        int[][][] y_Theta = new int[2*N][][]; 
        int[][][] y_X = new int[2*N][][];
        
        int[][] ynodeTheta;
        int[][] ynodeX;
        
        y_node[0] = index_into_dataIdxs_here;
        
        if (N * Math.log10(N) < numTheta || numTheta == 0) {
        	y_Theta[0] = null;
        } else {
        	ynodeTheta = new int[numTheta][];
        	int[] Thetacount = new int[numTheta];
        	for (int i=0; i < N; i++) {
                Thetacount[dataIdxs[i][0]]++;
        	}
        	for (int i=0; i < numTheta; i++) {
                ynodeTheta[i] = new int[Thetacount[i]];
            }
        	for (int i=0; i < N; i++) {
                int dataIdx = index_into_dataIdxs_here[i];
                int idx = dataIdxs[dataIdx][0];
                ynodeTheta[idx][--Thetacount[idx]] = i;
        	}
        	y_Theta[0] = ynodeTheta;
        	Thetacount = null;
        }
        
        if (N * Math.log10(N) < numX || numX == 0) {
        	y_X[0] = null;
        } else {
        	ynodeX = new int[numX][];
        
	        int[] Xcount = new int[numX];
	        for (int i=0; i < N; i++) {
	            Xcount[dataIdxs[i][1]]++;
	        }
	        for (int i=0; i < numX; i++) {
	            ynodeX[i] = new int[Xcount[i]];
	        }
        
	        for (int i=0; i < N; i++) {
	            int idx = dataIdxs[index_into_dataIdxs_here[i]][1];
	            ynodeX[idx][--Xcount[idx]] = i;
	        }
	        y_X[0] = ynodeX;
	        Xcount = null;
        }
        //=== End: initialize ynodeTheta and ynodeX for the root node.

        
        sorder = new int[Math.max(N, maxDomSize)];
        nodesize[0] = N;
        
        //========== Gather data for building the tree (we only build the actual tree afterwards using that data). ==========
        int[] stack = new int[N]; // Stack for DFS
        stack[0] = 0;
        int stacktop = 0; // Top of the stack
        int numNodes = 1; // Number of nodes in the tree so far
        
        
        if (printDebug)
        	System.out.println("Setup took " + (-currentTime + (currentTime = new Date().getTime())) + " milliseconds.");
        
        while (stacktop >= 0) {
        	//== Get the data for this node.
            int tnode = stack[stacktop--];
            
            index_into_dataIdxs_here = y_node[tnode];
            ynodeTheta = y_Theta[tnode];
            ynodeX = y_X[tnode];
            y_node[tnode] = null; // only needed up to here, then freed for garbage collection
            y_Theta[tnode] = null; // ditto
            y_X[tnode] = null; // ditto
            
        	//== Compute some basic stats for this node.
            int Nnode = nodesize[tnode];
            if (Nnode == 0) throw new RuntimeException("ERROR! Nnode is 0 (split gave zero data points to this node!?)");
            int NnodeUncens = 0;
            double ysum = 0, ysumOfSq = 0;
            double ymax = -1e13, ymin = 1e13;
            for (int i=0; i < Nnode; i++) {
                int idx = index_into_dataIdxs_here[i];
                ysum += y[idx];
                ysumOfSq += y[idx]*y[idx];
                
                if (y[idx] > ymax) ymax = y[idx];
                if (y[idx] < ymin) ymin = y[idx];
                
                if (!cens[idx]) NnodeUncens++;
            }
            ybar = ysum / Nnode;
            double mincost = (Nnode == 1 ? 0 : (ysumOfSq - ysum*ysum/Nnode) / (Nnode-1));
            boolean impure = (mincost > 1e-20 * ystd);
            impure = (ymax - ymin > 1e-10);

            cutvar[tnode] = 0; // this marks the current node as a leaf for now until we decide to split it
            
            //System.out.println("in node " + tnode + " theta sort: " + (ynodeTheta == null ? "naive" : "pre") + ", x sort: " + (ynodeX == null ? "naive" : "pre"));
            /*
            System.out.print("In node " + tnode + " parent " + parent[tnode] + " size " + nodesize[tnode] + " y [");
            for (int i=0; i < Nnode; i++) System.out.print(y[index_into_dataIdxs_here[i]] + " ");
            System.out.println("]");//*/

            if (impure && NnodeUncens >= splitMin) { // split only impure nodes with more than a threshold of uncensored values
                //=== Start: handle conditional parameters. 
                int nvarsenabled = 0; // #variables that are active for sure given the variable instantiations up to the root 
                if (condParents == null) {
                    nvarsenabled = nvars;
                    for (int i=0; i < nvars; i++) {
                        randomPermutation[i] = i;
                    }
                } else {
                	//=== Determine which variables are active in this node. This is currently O(#vars + #cond. vars * average #parents of cond. vars * domain size^2).
                	//=== For large domain sizes, this could be slow; domain size^2 could be easily sped up to domain size * log(domain size) by using sets
                    for (int i=0; i < nvars; i++) {
                        boolean isenabled = true;
                        if (condParents[i] != null) {
                            for (int j=0; j < condParents[i].length; j++) {
                                int[] compatibleValues = getCompatibleValues(tnode, condParents[i][j], N, parent, cutvar, cutpoint, leftchildren, rightchildren, catsplit, catDomainSizes);

                                for (int k=0; k < compatibleValues.length; k++) {
                                    boolean isokvalue = false;
                                    for (int l=0; l < condParentVals[i][j].length; l++) {
                                        if (compatibleValues[k] == condParentVals[i][j][l]) {
                                            isokvalue = true;
                                            break;
                                        }
                                    }
                                    if (!isokvalue) {
                                        isenabled = false;
                                        break;
                                    }
                                }
                                if (!isenabled) break;
                            }
                        }
                        if (isenabled) {
                            randomPermutation[nvarsenabled++] = i;
                        }
                    }
                }
                //=== End: handle conditional parameters
                shuffle(randomPermutation, nvarsenabled);
                
                num_periods = 0;
                if (NnodeUncens != Nnode) {
                    // First, compute period for uncensored data. (Censored data with same y value could be before or after the first uncensored one in the ordering, thus we need 2 passes.)                                     
                    Arrays.fill(period, 0, Nnode, 0);
                    
                    double last_time = -1e10;
                    for (int i=0; i < Nnode; i++) {
                        int idx = index_into_dataIdxs_here[i];
                        if (!cens[idx]) {
                            if (y[idx] > last_time + 1e-6) {
                                num_periods++;
                                period_time[num_periods] = y[idx];
                                last_time = y[idx];
                            }
                            period[i] = num_periods;
                        }
                    }
                    
                    Arrays.fill(O_all, 0, num_periods, 0);
                    Arrays.fill(C_all, 0, num_periods, 0);
                    // Then, fill in period for censored data
                    int curr_period = 0;
                    double this_time = 0;
                    for (int i=0; i < Nnode; i++) {
                        int idx = index_into_dataIdxs_here[i];
                        if (cens[idx]) {
                            if (curr_period == num_periods) {
                                this_time = 1e9;
                            } else {
                                this_time = period_time[curr_period+1];
                            }
                            while ( y[idx] > this_time - 1e-6 ) {
                                curr_period++;
                                if (curr_period == num_periods) {
                                    this_time = 1e9;
                                } else {
                                    this_time = period_time[curr_period+1];
                                }
                            }
                            period[i] = curr_period;
                            if (curr_period > 0) {
                                C_all[curr_period-1]++;
                            }
                        } else {
                            O_all[period[i]-1]++;
                        }
                    }
                    
                    // Go through periods, and collect N.
                    N_all[0] = Nnode;
                    for (int p = 1; p < num_periods; p++) {
                        N_all[p] = N_all[p-1] - O_all[p-1] - C_all[p-1];
                    }
                }
                //=== END censoring stuff
                
                //=== Keep track of the best cut/var found.
                int bestvar = -1;
                double bestcut = 0;
                double bestcrit = -1e12;
                
                //=== Try splitting each variable at every split point and pick best split ===
                for (int i=0; i < nvarsenabled; i++) {
                    int nextvar = randomPermutation[i];
                    boolean is_nextvar_cat = (catDomainSizes[nextvar] != 0);
                    
                    int varIdx, numData;
                    int[][] sortedData, ynodeData;
                    double[][] allData;
                    int is_X;
                    if (nextvar < numThetavars) {
                        varIdx = nextvar;
                        numData = numTheta;
                        sortedData = sortedTheta;
                        ynodeData = ynodeTheta;
                        allData = allTheta;
                        is_X = 0;
                    } else {
                        varIdx = nextvar - numThetavars;
                        numData = numX;
                        sortedData = sortedX;
                        ynodeData = ynodeX;
                        allData = allX;
                        is_X = 1;
                    }
                    
                    double critval=INVALID_CRITVAL, cutval=INVALID_CRITVAL;
                    if (is_nextvar_cat) { // Categorical variable
                        int domSize = catDomainSizes[nextvar];

                        // compute critval
                        if (hascens) { // Censored Categorical
                            critval = critval_cat_logrank(varIdx, is_X, allData, Nnode, index_into_dataIdxs_here, domSize);
                        } else { // Uncensored Categorical
                        	critval = critval_cat(varIdx, is_X, allData, Nnode, index_into_dataIdxs_here, domSize);
                        }
                        //=== Change best split if this one is best so far.
                        if (critval > bestcrit + 1e-10) {
                            bestcrit = critval;
                            bestvar = nextvar;
                            numBestLeft = numleft;
                            numBestRight = numright;
                            for (int j=0; j < numBestLeft; j++) {
                                bestLeft[j] = leftside[j];
                            }
                            for (int j=0; j < numBestRight; j++) {
                                bestRight[j] = rightside[j];
                            }
                        }
                    } else { // Continuous variable  
                    	// Get the values of y that we have in this node, in order corresponding to sorted variable values
                        int[] results = prepare_for_cont_critval(varIdx, is_X, numData, sortedData, allData, Nnode, index_into_dataIdxs_here, ynodeData, variableValuesHere);
                        if (results == null) continue;
                        int numUniqData = results[0];
                        int numUniqValues = results[1];
                        
                        // Compute critval
                        if (hascens) { // Censored Continuous
                            double[] result = critval_cont_logrank(numUniqValues, uniqueIdxs, index_into_dataIdxs_here, ynodeData, dataRowsHere, variableValuesHere);
                            critval = result[0];
                            cutval = result[1];
                        } else { // Noncensored Continuous
                        	double[] result = critval_cont(numUniqData, numUniqValues, uniqueIdxs, index_into_dataIdxs_here, ynodeData, dataRowsHere, variableValuesHere);
                            critval = result[0];
                            cutval = result[1];
                        }
                        //=== Change best split if this one is best so far.
                        if (critval > bestcrit + 1e-10) {
                            bestcrit = critval;
                            bestvar = nextvar;
                            bestcut = cutval;
                        }
                    }
                    if (critval == INVALID_CRITVAL) {
                        continue;
                    }
                    /*
                    System.out.print("var " + nextvar + " critval " + critval + (!is_nextvar_cat ? " cutval " + cutval : "") + " x [");
                    for (int j=0; j < Nnode; j++) {
                        int idx = index_into_dataIdxs_here[j];
                        System.out.print(allData[dataIdxs[idx][is_X]][varIdx] + " ");
                    }
                    System.out.println(" ]");//*/
                    
                    if (i >= Math.max(1, (int)(ratioFeatures*nvarsenabled)) - 1 && bestcrit > -1e11) {
                        //System.out.println("Breaking: i " + i + " num " + (Math.max(1, (int)((ratioFeatures*nvarsenabled)))-1) + " bestcrit " + bestcrit);
                        break;
                    }
                }
                //System.out.println("Bestvar " + bestvar + " critval " + bestcrit);
                //=== Best split point has been found. Split this node using the best rule found.
                if (bestvar != -1) {
                    int numPrimary;
                    int[][] ynodePrimary;
                    double[][] allPrimary;
                    int is_X;
                    int varIdx;
                    if (bestvar < numThetavars) {
                        varIdx = bestvar;
                        numPrimary = numTheta;         ynodePrimary = ynodeTheta;
                        allPrimary = allTheta;
                        is_X = 0;
                    } else {
                        varIdx = bestvar - numThetavars;
                        numPrimary = numX;             ynodePrimary = ynodeX;
                        allPrimary = allX;
                        is_X = 1;
                    }
                    
                    int nleft = 0, nright = 0;
                    if (catDomainSizes[bestvar]!=0) { 
                        cutvar[tnode] = -(bestvar+1); // negative indicates cat. var. split
                        cutpoint[tnode] = ncatsplit; // index into catsplit cell array. 0-indexed!!!
                        
                        /* 1: To get all compatible values, walk up the tree, looking
                         * for a split on the same parameter. If none is found
                         * take the initial domain of that parameter.
                         */
                        int[] compatibleValues = getCompatibleValues(tnode, bestvar, N, parent, cutvar, cutpoint, leftchildren, rightchildren, catsplit, catDomainSizes);
                        
                        int[] missing_values_for_left = new int[compatibleValues.length];
                        int[] missing_values_for_right = new int[compatibleValues.length];
                        
                        // 2: For each compatible but missing value choose a side u.a.r.
                        // TODO: If compatibleValues and bestLeft/bestRight were all sorted, this can be faster by a lot
						int num_missing_to_left = 0, num_missing_to_right = 0;
                        for (int i=0; i < compatibleValues.length; i++) {
                            int nextValue = compatibleValues[i];
                            int j;
                            for (j=0; j < numBestLeft; j++) {
                                if (nextValue == bestLeft[j]) break;
                            }
                            if (j == numBestLeft) {
                                for (j=0; j < numBestRight; j++) {
                                    if (nextValue == bestRight[j]) break;
                                }
                                if (j == numBestRight) {
                                    // Missing but compatible value: choose side u.a.r.
                                    if (rand() % 2 == 0) {
                                        missing_values_for_left[num_missing_to_left++] = nextValue;
                                    } else {
                                        missing_values_for_right[num_missing_to_right++] = nextValue;
                                    }
                                }
                            }
                        }
                        
                        // 3: Store the information                        
						catsplit[ncatsplit] = new int[num_missing_to_left + numBestLeft];
                        for (int i=0; i<num_missing_to_left; i++) {
                            int nextval = missing_values_for_left[i];
                            leftside[nextval-1] = 1; // because catsplit is 1-indexed
                            catsplit[ncatsplit][i] = nextval;
                        }
						for (int i=0; i<numBestLeft; i++) {
                            int nextval = bestLeft[i];
                            leftside[nextval-1] = 1; // because catsplit is 1-indexed
                            catsplit[ncatsplit][num_missing_to_left+i] = nextval;
                        }
						Arrays.sort(catsplit[ncatsplit]);
						
                        catsplit[ncatsplit+N] = new int[num_missing_to_right + numBestRight];
                        for (int i=0; i<num_missing_to_right; i++) {
                            int nextval = missing_values_for_right[i];
                            leftside[nextval-1] = 0; // because catsplit is 1-indexed
                            catsplit[ncatsplit+N][i] = nextval;
                        }
						for (int i=0; i<numBestRight; i++) {
                            int nextval = bestRight[i];
                            leftside[nextval-1] = 0; // because catsplit is 1-indexed
                            catsplit[ncatsplit+N][num_missing_to_right+i] = nextval;
                        }
						Arrays.sort(catsplit[ncatsplit+N]);
						ncatsplit++;
                        
                        if (ynodePrimary == null) {
                            for (int i=0; i < Nnode; i++) {
                                int idx = index_into_dataIdxs_here[i];
                                double xVal = allPrimary[dataIdxs[idx][is_X]][varIdx];
                                boolean onleft = false;
                                for (int j=0; j < numBestLeft; j++) {
                                    if ((int)(xVal+0.5) == bestLeft[j]) {
                                        onleft = true;
                                        break;
                                    }
                                }
                                if (onleft) {
                                    nleft++;
                                } else {
                                    nright++;
                                }
                                yGoesLeft[i] = onleft;
                            }
                        } else {
                            for (int i=0; i < numPrimary; i++) {
                                if (leftside[(int)(allPrimary[i][varIdx]-0.5)] == 1) {
                                    primaryGoesLeft[i] = true;
                                    if (ynodePrimary[i] != null) {
                                        for (int j=0; j < ynodePrimary[i].length; j++) {
                                            nleft++;
                                            yGoesLeft[ynodePrimary[i][j]] = true;
                                        }
                                    }
                                } else {
                                    primaryGoesLeft[i] = false;
                                    if (ynodePrimary[i] != null) {
                                        for (int j=0; j < ynodePrimary[i].length; j++) {
                                            nright++;
                                            yGoesLeft[ynodePrimary[i][j]] = false;
                                        }
                                    }
                                }
                            }
                        }
                    } else { 
                        cutvar[tnode] = bestvar + 1; // splitting on cont. var
                        cutpoint[tnode] = bestcut;
                        
                        if (ynodePrimary == null) {
                            for (int i=0; i < Nnode; i++) {
                                int idx = index_into_dataIdxs_here[i];
                                double xVal = allPrimary[dataIdxs[idx][is_X]][varIdx];
                                if (xVal <= bestcut) {
                                    nleft++;
                                    yGoesLeft[i] = true;
                                } else {
                                    nright++;
                                    yGoesLeft[i] = false;
                                }
                            }
                        } else {
                            for (int i=0; i < numPrimary; i++) {
                                if (allPrimary[i][varIdx] <= bestcut) {
                                    primaryGoesLeft[i] = true;
                                    if (ynodePrimary[i] != null) {
                                        for (int j=0; j < ynodePrimary[i].length; j++) {
                                            nleft++;
                                            yGoesLeft[ynodePrimary[i][j]] = true;
                                        }
                                    }
                                } else {
                                    primaryGoesLeft[i] = false;
                                    if (ynodePrimary[i] != null) {
                                        for (int j=0; j < ynodePrimary[i].length; j++) {
                                            nright++;
                                            yGoesLeft[ynodePrimary[i][j]] = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (nleft == 0 || nright == 0) {
                        throw new RuntimeException("Empty side after splitting!");
                    }
                    
                    // Create y_node, y_Theta and y_X for children.
                    int numSecondary;
                    int[][] ynodeSecondary;
                    int[][][] y_Primary, y_Secondary;
                    if (bestvar < numThetavars) {
                        numSecondary = numX;
                        ynodeSecondary = ynodeX;
                        y_Primary = y_Theta;
                        y_Secondary = y_X;
                    } else {
                        numSecondary = numTheta;
                        ynodeSecondary = ynodeTheta;
                        y_Primary = y_X;
                        y_Secondary = y_Theta;
                    }
                    
                    int[] ynodeLeft = new int[nleft];
                    int[] ynodeRight = new int[nright];
                    for (int i=0, leftCounter=0, rightCounter=0; i < Nnode; i++) {
                        if (yGoesLeft[i]) {
                            ynodeLeft[leftCounter] = index_into_dataIdxs_here[i];
                            index_into_dataIdxs_here[i] = leftCounter++; // Keep track of what the new index is so we can update yPrimaryLeft/Right, which indexes into ynodeLeft/Right
                        } else {
                            ynodeRight[rightCounter] = index_into_dataIdxs_here[i];
                            index_into_dataIdxs_here[i] = rightCounter++;
                        }
                    }
                    y_node[numNodes] = ynodeLeft;
                    y_node[numNodes+1] = ynodeRight;
                    
                    boolean naiveSortPrimaryLeft = (ynodePrimary == null || (nleft * Math.log10(nleft) < numPrimary));
                    boolean naiveSortPrimaryRight = (ynodePrimary == null || (nright * Math.log10(nright) < numPrimary));
                    
                    int[][] yPrimaryLeft = naiveSortPrimaryLeft ? null : new int[numPrimary][];
                    int[][] yPrimaryRight = naiveSortPrimaryRight ? null : new int[numPrimary][];
                    if (!naiveSortPrimaryLeft || !naiveSortPrimaryRight) {
	                    for (int i=0; i < numPrimary; i++) {
	                        if (primaryGoesLeft[i]) {
	                        	if (!naiveSortPrimaryLeft) {
		                            yPrimaryLeft[i] = ynodePrimary[i];
		                            if (yPrimaryLeft[i] != null) {
		                                for (int j=0; j < yPrimaryLeft[i].length; j++) {
		                                    yPrimaryLeft[i][j] = index_into_dataIdxs_here[yPrimaryLeft[i][j]]; // Update to the new idxs.
		                                }
		                            }
	                        	}
	                        }
	                        else if (!naiveSortPrimaryRight) {
	                            yPrimaryRight[i] = ynodePrimary[i];
	                            if (yPrimaryRight[i] != null) {
	                                for (int j=0; j < yPrimaryRight[i].length; j++) {
	                                    yPrimaryRight[i][j] = index_into_dataIdxs_here[yPrimaryRight[i][j]];
	                                }
	                            }
	                        }
	                    }
                    }
                    y_Primary[numNodes] = yPrimaryLeft;
                    y_Primary[numNodes+1] = yPrimaryRight;
                    
                    boolean naiveSortSecondaryLeft = (ynodeSecondary == null || (nleft * Math.log10(nleft) < numSecondary));
                    boolean naiveSortSecondaryRight = (ynodeSecondary == null || (nright * Math.log10(nright) < numSecondary));

                    int[][] ySecondaryLeft = naiveSortSecondaryLeft ? null : new int[numSecondary][];
                    int[][] ySecondaryRight = naiveSortSecondaryRight ? null : new int[numSecondary][];
                    if (!naiveSortSecondaryLeft || !naiveSortSecondaryRight) {
	                    for (int i=0; i < numSecondary; i++) {
	                        if (ynodeSecondary[i] != null) {
	                            int[] thisynodeSecondary = ynodeSecondary[i];
	                            int numySecondaryLeft = 0, numySecondaryRight = 0;
	                            for (int j=0; j < thisynodeSecondary.length; j++) {
	                                if (yGoesLeft[thisynodeSecondary[j]]) numySecondaryLeft++;
	                                else numySecondaryRight++;
	                            }
	                            
	                            int[] ySecondaryLeft_i = naiveSortSecondaryLeft ? null : new int[numySecondaryLeft];
	                            int[] ySecondaryRight_i = naiveSortSecondaryRight ? null : new int[numySecondaryRight];
	                            numySecondaryLeft = 0;
	                            numySecondaryRight = 0;
	                            for (int j=0; j < thisynodeSecondary.length; j++) {
	                                if (yGoesLeft[thisynodeSecondary[j]]) {
	                                	if (!naiveSortSecondaryLeft) {
	                                		ySecondaryLeft_i[numySecondaryLeft++] = index_into_dataIdxs_here[thisynodeSecondary[j]];
	                                	}
	                                } else if (!naiveSortSecondaryRight) {
	                                	ySecondaryRight_i[numySecondaryRight++] = index_into_dataIdxs_here[thisynodeSecondary[j]];
	                                }
	                            }
	
	                            if (numySecondaryLeft != 0) {
	                                ySecondaryLeft[i] = ySecondaryLeft_i;
	                            }
	                            if (numySecondaryRight != 0) {
	                                ySecondaryRight[i] = ySecondaryRight_i;
	                            }
	                        }
	                    }
                    }
                    y_Secondary[numNodes] = ySecondaryLeft;
                    y_Secondary[numNodes+1] = ySecondaryRight;
                    
                    leftchildren[tnode] = numNodes;
                    rightchildren[tnode] = numNodes+1;
                    nodenumber[numNodes] = numNodes;
                    nodenumber[numNodes+1] = numNodes+1;
                    parent[numNodes] = tnode;
                    parent[numNodes+1] = tnode;
                    
                    nodesize[numNodes] = nleft;
                    nodesize[numNodes+1] = nright;
                    
                    stack[++stacktop] = numNodes;
                    stack[++stacktop] = numNodes+1;
                    numNodes += 2; 
                }
            }
            // Leaf => store results falling here (don't store them everywhere to avoid O(N^2) storage)
            // Save *runtimes*, not losses. 
            if (cutvar[tnode] == 0) {
                ysub[tnode] = new double[Nnode];
                censsub[tnode] = new boolean[Nnode];
                for (int i=0; i < Nnode; i++) {
                    int idx = index_into_dataIdxs_here[i];
                    ysub[tnode][i] = y[idx];
                    censsub[tnode][i] = cens[idx];
                }
            }
            if (printDebug) {
            	long diff = (-currentTime + (currentTime = new Date().getTime()));
            	if (diff > 1000)
            		System.out.println("Node " + tnode + " had " + Nnode + " data points and took " + diff + " milliseconds.");
            }
        }
        
        //==================== Build the actual tree ===========================
        Regtree tree = new Regtree(numNodes, ncatsplit, params.storeResponses, params.logModel);
        
        System.arraycopy(nodenumber, 0, tree.node, 0, numNodes);
        System.arraycopy(parent, 0, tree.parent, 0, numNodes);
        System.arraycopy(cutvar, 0, tree.var, 0, numNodes);
        System.arraycopy(cutpoint, 0, tree.cut, 0, numNodes);
        System.arraycopy(nodesize, 0, tree.nodesize, 0, numNodes);
        tree.npred = nvars;
        
        int nextnode=-1;
        for (int i=0; i < ncatsplit; i++) {
            while(cutvar[++nextnode] >= 0);
            int[] tmp = new int[catDomainSizes[-cutvar[nextnode]-1]];
            Arrays.fill(tmp, -1);
			int cs_idx = (int)cutpoint[nextnode];
            int[] cs = catsplit[cs_idx];
            for (int j=0; j < cs.length; j++) {
                tmp[cs[j]-1] = 0;
            }
            cs = catsplit[cs_idx + N];
            for (int j=0; j < cs.length; j++) {
                tmp[cs[j]-1] = 1;
            }
            tree.catsplit[cs_idx] = tmp;
	   }
        
        for (int i=0; i < numNodes; i++) {
            tree.children[i][0] = leftchildren[i];
            tree.children[i][1] = rightchildren[i];
            
            int Nnode = cutvar[i] == 0 ? nodesize[i] : 0;
            if (Nnode != 0) {
                if (params.storeResponses) {
                    tree.ysub[i] = new double[Nnode];
                    tree.is_censored[i] = new boolean[Nnode]; 
                    System.arraycopy(ysub[i], 0, tree.ysub[i], 0, Nnode);
                    System.arraycopy(censsub[i], 0, tree.is_censored[i], 0, Nnode);
                } else {
                    double sum = 0, sumOfSq = 0;
                    for (int j=0; j < Nnode; j++) {
                        double next = ysub[i][j];
                        sum += next;
                        sumOfSq += next * next;
                    }
                    tree.ysub[i][0] = sum;
                    tree.ysub[i][1] = sumOfSq;
                }
            }
        }
        tree.recalculateStats();
        
        
        // Free up static fields for GC
        rows_km = null;
        km_mean = null;
        numUptoCategory = null;
        sortedByCategory = null;
        sorder = null;
        maxlocs = null;
        RegtreeFit.dataIdxs = null;
        RegtreeFit.y = null;
        RegtreeFit.cens = null;
        
        catmeans = null;
        catcounts = null;
        ycum = null;
        ycountcum = null;
        uniqueIdxs = null;
        dataRowsHere = null;
        
        leftside = null;
        rightside = null;
        
        period = null;
    	period_time = null;
    	N_all = null;
    	O_all = null;
    	C_all = null;
    	N_1 = null;
    	O_1 = null;
    	N_add = null;
    	
    	
    	if (printDebug)
    		System.out.println("Building the tree took a total of " + (new Date().getTime() - startTime) + " milliseconds.");   	
        return tree;
    }
    
    private static int[] prepare_for_cont_critval(int varIdx, int is_X, int numData, int[][] sortedData, double[][] allData, int Nnode, int[] index_into_dataIdxs_here, int[][] ynodeData, double[] variableValuesHere) {
    	int numUniqData = 0;
        int numUniqValues = 0;
    	if (ynodeData == null) { // do Nnode log Nnode sorting
            for (int j=0; j < Nnode; j++) {
                int idx = index_into_dataIdxs_here[j];
                variableValuesHere[j] = allData[dataIdxs[idx][is_X]][varIdx];
            }
            rankSort(variableValuesHere, Nnode, sorder);
            
            double prevValue = variableValuesHere[sorder[0]];
            uniqueIdxs[numUniqValues++] = 0;
            for (int j=1; j < Nnode; j++) {
                double nextValue = variableValuesHere[sorder[j]];
                if (prevValue + 1e-10 < nextValue) {
                    uniqueIdxs[numUniqValues++] = j;
                }
                prevValue = nextValue;
            }
            numUniqData = Nnode;
        } else { 
            double prevValue = 0;
            for (int j=0; j < numData; j++) {
                int nextIdx = sortedData[varIdx][j];
                int[] yhere = ynodeData[nextIdx];
                if (yhere != null) {
                    double nextValue = allData[nextIdx][varIdx];
                    if (numUniqValues == 0 || prevValue + 1e-10 < nextValue) {
                        uniqueIdxs[numUniqValues] = numUniqData; // the start of a new value.
                        variableValuesHere[numUniqValues++] = nextValue;
                    }
                    prevValue = nextValue;
                    dataRowsHere[numUniqData++] = nextIdx;
                }
            }
        }
        if (numUniqValues <= 1) return null;
        return new int[]{numUniqData, numUniqValues};
	}

	private static double[] critval_cont(int numUniqData, int numUniqValues, int[] uniqueIdxs, int[] index_into_dataIdxs_here, int[][] ynodeData, int[] dataRowsHere, double[] variableValuesHere) {
        double critval = INVALID_CRITVAL;

        ycum[0] = 0;
        ycountcum[0] = 0;
        if (ynodeData == null) { // did Nnode log Nnode sorting
            for (int j=1; j <= numUniqData; j++) {
                ycum[j] = ycum[j-1] + y[index_into_dataIdxs_here[sorder[j-1]]] - ybar;
                ycountcum[j] = ycountcum[j-1] + 1;
            }
        } else {
            for (int j=1; j <= numUniqData; j++) {
                int[] ynodeIdxsHere = ynodeData[dataRowsHere[j-1]];
                int numYValuesHere = ynodeIdxsHere.length;
                double sumYValuesHere = 0;
                for (int k : ynodeIdxsHere) {
                    sumYValuesHere += y[index_into_dataIdxs_here[k]];
                }
                ycum[j] = ycum[j-1] + sumYValuesHere - numYValuesHere * ybar; // centered cumulative sum
                ycountcum[j] = ycountcum[j-1] + numYValuesHere;
            }
        }
        double ytotal = ycum[numUniqData];
        int numytotal = ycountcum[numUniqData];

        int numlocs_with_max_crit = 0;
        for (int j=1; j < numUniqValues; j++) {
            int idx = uniqueIdxs[j];
            double yc = ycum[idx];
            double ssx = yc*yc/ycountcum[idx] + (ytotal-yc)*(ytotal-yc)/(numytotal-ycountcum[idx]);
            if (ssx > critval - 1e-10) {
                if (ssx > critval + 1e-10) {
                    critval = ssx;
                    numlocs_with_max_crit = 0;
                }
                maxlocs[numlocs_with_max_crit++] = j-1;
            }
        }
        int maxloc = maxlocs[rand() % numlocs_with_max_crit];

        //=== Get cutval.
        double u = rand() * 1.0 / RAND_MAX;
        // if points are close the just take average. If points are farther sample randomly from lerp
        double prev, next;
        if (ynodeData == null) {
            prev = variableValuesHere[sorder[uniqueIdxs[maxloc]]];
            next = variableValuesHere[sorder[uniqueIdxs[maxloc+1]]];
        } else {
            prev = variableValuesHere[maxloc];
            next = variableValuesHere[maxloc+1];
        }
        
        //System.out.println("*** u: " + u + " prev: " + prev + " next: " + next + " sort: " + (ynodeData == null ? "naive" : "pre"));
        
        double cutval = 0;
        if (next - prev < 1.9*1e-6) {
            cutval = (next + prev) / 2;
        } else {
            cutval = (1-u)*(prev + 1e-6) + u*(next-1e-6);
            if (cutval < prev + 1e-8 || cutval > next - 1e-8) {
                throw new RuntimeException("random splitpoint has to lie in between the upper and lower limit");
            }
        }
        return new double[]{critval, cutval};
	}
	private static double[] critval_cont_logrank(int numUniqValues, int[] uniqueIdxs, int[] index_into_dataIdxs_here, int[][] ynodeData, int[] dataRowsHere, double[] variableValuesHere) {
    	Arrays.fill(N_1, 0, num_periods, 0);
        Arrays.fill(O_1, 0, num_periods, 0);
        
        //=== For adding one group at a time, compute resulting logrank statistic.
        double critval = INVALID_CRITVAL;
        int numlocs_with_max_crit = 0;
        for (int j=1; j < numUniqValues; j++) {
            int low = uniqueIdxs[j-1];
            int high = uniqueIdxs[j];
            
            int maxperiod = -1;
            if (ynodeData == null) { // did Nnode log Nnode sorting
                for (int k=low; k < high; k++) {
                    maxperiod = Math.max(maxperiod, period[sorder[k]]);
                }
            } else {
                for (int k=low; k < high; k++) {
                    int[] ynodeIdxsHere = ynodeData[dataRowsHere[k]];
                    for (int l : ynodeIdxsHere) {
                        maxperiod = Math.max(maxperiod, period[l]);
                    }
                }
            }
            Arrays.fill(N_add, 0, maxperiod, 0);
            
            if (ynodeData == null) {
                for (int k = low; k < high; k++) {
                    int idx = sorder[k];
                    int p = period[idx]-1;
                    if (p >= 0) {
                        if (!cens[index_into_dataIdxs_here[idx]]) {
                            O_1[p]++;
                        }
                        N_add[p]++;
                    }
                }
            } else {
                for (int k = low; k < high; k++) {
                    int[] ynodeIdxsHere = ynodeData[dataRowsHere[k]];
                    for (int l : ynodeIdxsHere) {
                        int p = period[l]-1;
                        if (p >= 0) {
                            if (!cens[index_into_dataIdxs_here[l]]) {
                                O_1[p]++;
                            }
                            N_add[p]++;
                        }
                    }
                }
            }
            
            int N_cumul = 0;
            for (int p = maxperiod-1; p >= 0; p--) {
                N_cumul += N_add[p];
                N_1[p] += N_cumul;
            }
            
            double logrank = Math.abs(logrank_statistic(num_periods, N_all, O_all, N_1, O_1));
            if (logrank > critval - 1e-10) {
                if (logrank > critval + 1e-10) {
                    critval = logrank;
                    numlocs_with_max_crit = 0;
                }
                maxlocs[numlocs_with_max_crit++] = j-1;
            }
        }
        int maxloc = maxlocs[rand() % numlocs_with_max_crit];

        //=== Get cutval.
        double u = rand() * 1.0 / RAND_MAX;
        // if points are close the just take average. If points are farther sample randomly from lerp
        double prev, next;
        if (ynodeData == null) {
            prev = variableValuesHere[sorder[maxloc]];
            next = variableValuesHere[sorder[maxloc+1]];
        } else {
            prev = variableValuesHere[maxloc];
            next = variableValuesHere[maxloc+1];
        }

        double cutval = 0;
        if (next - prev < 1.9*1e-6) {
            cutval = (next + prev) / 2;
        } else {
            cutval = (1-u)*(prev + 1e-6) + u*(next-1e-6);
            if (cutval < prev + 1e-8 || cutval > next - 1e-8) {
                throw new RuntimeException("random splitpoint has to lie in between the upper and lower limit");
            }
        }
        return new double[]{critval, cutval};
	}
	
	private static double critval_cat(int varIdx, int var_is_X, double[][] allData, int Nnode, int[] index_into_dataIdxs_here, int domSize) {
		double critval = INVALID_CRITVAL;
		
		// Sort by category means
        Arrays.fill(catmeans, 0, domSize, 0);
        Arrays.fill(catcounts, 0, domSize, 0);
        
        for (int j=0; j < Nnode; j++) {
        	// Calculate categorical sums and # of data points in each category
            int idx = index_into_dataIdxs_here[j];            
            int category = (int)(allData[dataIdxs[idx][var_is_X]][varIdx] - 0.5);
            catmeans[category] += y[idx];
            catcounts[category]++;
        }
        
        int numtotal = 0;
        // calculate categorical means
        for (int j=0; j < domSize; j++) {
            if (catcounts[j] != 0) {
                numtotal++;
                catmeans[j] /= catcounts[j];
            }
            else catmeans[j] = Double.POSITIVE_INFINITY; // We don't have a data point in this category
        }
        if (numtotal <= 1) return critval;

        rankSort(catmeans, domSize, sorder);
        
        // Calculuate cumulative sums and counts.
        ycum[0] = 0;
        ycountcum[0] = 0;
        for (int j=1; j <= numtotal; j++) {
            int idx = sorder[j-1];
            ycum[j] = ycum[j-1] + catcounts[idx] * (catmeans[idx] - ybar);
            ycountcum[j] = ycountcum[j-1] + catcounts[idx];
        }
        double ytotal = ycum[numtotal];
        int numytotal = ycountcum[numtotal];

    	int numlocs_with_max_crit = 0;
        for (int j=1; j < numtotal; j++) {
            double ssx = ycum[j]*ycum[j]/ycountcum[j] + (ytotal-ycum[j])*(ytotal-ycum[j])/(numytotal-ycountcum[j]);
            if (ssx > critval - 1e-10) {
                if (ssx > critval + 1e-10) {
                    critval = ssx;
                    numlocs_with_max_crit = 0;
                }
                maxlocs[numlocs_with_max_crit++] = j;
            }
        }
        
        int maxloc = maxlocs[rand() % numlocs_with_max_crit];
        
        numleft = maxloc;
        numright = numtotal - numleft;
        for (int j=0; j < numleft; j++) {
            leftside[j] = sorder[j] + 1;
        }
        for (int j=0; j < numright; j++) {
            rightside[j] = sorder[j+numleft] + 1;
        }
		return critval;
	}
	private static double critval_cat_logrank(int varIdx, int var_is_X, double[][] allData, int Nnode, int[] index_into_dataIdxs_here, int domSize) {
		double critval = INVALID_CRITVAL;
		
		// Get Kaplan-Meier estimates of the mean for each category       
        Arrays.fill(numUptoCategory, 0, domSize+1, 0); 
        for (int j=0; j < Nnode; j++) {
            int idx = index_into_dataIdxs_here[j];
            int category = (int)(allData[dataIdxs[idx][var_is_X]][varIdx] - 0.5);
            numUptoCategory[category]++;
        }
        
        for (int j=1; j <= domSize; j++) {
            numUptoCategory[j] += numUptoCategory[j-1];
        }
        // index i of numUptoCategory contains the number of data points with category numbers <= i
        
        for (int j=0; j < Nnode; j++) {
        	//=== Sort all data points in this node by category and place result in dataRowsHere (counting sort).
            int idx = index_into_dataIdxs_here[j];
            int category = (int)(allData[dataIdxs[idx][var_is_X]][varIdx] - 0.5);
            sortedByCategory[--numUptoCategory[category]] = j;
        }
        int[] categoryStartIdx = numUptoCategory; // the array now contains the start idx of each category.
        int numtotal = 0;
        for (int j=0; j < domSize; j++) {
            int low = categoryStartIdx[j];
            int high = categoryStartIdx[j+1];
            
            int numUncens = 0, numCens = 0;
            for (int k=high-1; k >= low; k--) { // Iterate in reverse order so c_in and y_in are sorted
                int idx = index_into_dataIdxs_here[sortedByCategory[k]];
                if (cens[idx]) {
                    c_in[numCens++] = y[idx];
                } else {
                    y_in[numUncens++] = y[idx];
                }
            }
            if (high-low == 0) {
                km_mean[j] = Double.POSITIVE_INFINITY;
            } else {
                km_mean[j] = kaplan_meier_mean(numUncens, numCens, y_in, c_in, kappa);
                numtotal++;
            }
        }
        if (numtotal <= 1) return critval;
        
        rankSort(km_mean, domSize, sorder);
        
        Arrays.fill(N_1, 0, num_periods, 0);
        Arrays.fill(O_1, 0, num_periods, 0);
        //=== For adding one group at a time, compute resulting logrank statistic.
        int numlocs_with_max_crit = 0;
        for (int j=1; j < numtotal; j++) {
            int category = sorder[j-1];
            int low = categoryStartIdx[category];
            int high = categoryStartIdx[category+1];
            
            int maxperiod = -1;
            for (int k=low; k < high; k++) {
                int idx = sortedByCategory[k];
                maxperiod = Math.max(maxperiod, period[idx]);
            }
            Arrays.fill(N_add, 0, maxperiod, 0);
            
            for (int k = low; k < high; k++) {
                int idx = sortedByCategory[k];
                int p = period[idx]-1;
                if (p >= 0) {
                    if (!cens[index_into_dataIdxs_here[idx]]) {
                        O_1[p]++;
                    }
                    N_add[p]++;
                }
            }
            
            int N_cumul = 0;
            for (int p = maxperiod-1; p >= 0; p--) {
                N_cumul += N_add[p];
                N_1[p] += N_cumul;
            }
            
            double logrank = Math.abs(logrank_statistic(num_periods, N_all, O_all, N_1, O_1));
            
            if (logrank > critval - 1e-10) {
                if (logrank > critval + 1e-10) {
                    critval = logrank;
                    numlocs_with_max_crit = 0;
                }
                maxlocs[numlocs_with_max_crit++] = j;
            }
        }
        
        int maxloc = maxlocs[rand() % numlocs_with_max_crit];
        
        numleft = maxloc;
        numright = numtotal - numleft;
        for (int j=0; j < numleft; j++) {
            leftside[j] = sorder[j] + 1;
        }
        for (int j=0; j < numright; j++) {
            rightside[j] = sorder[j+numleft] + 1;
        }
        return critval;
	}

	// NOTE: the y_km and c_km array you pass in WILL be modified!
    // Assumes y_km and c_km are sorted already.
    private static int[] rows_km;
    private static double kaplan_meier_mean(int numUncens, int numCens, double[] y_km, double[] c_km, double upper_bound) {
        int N_i = numUncens + numCens;
        if (N_i == 0) {
            throw new RuntimeException("Kaplan-Meier estimator undefined for empty population.");
        }
        double km_mean = 0;
        
        if (numUncens == 0) {
            km_mean = upper_bound;
        } else {
            int num_rows = 0;
            for (int i=0; i < numUncens-1; i++) {
                if (y_km[i] + 1e-10 < y_km[i+1]) {
                    rows_km[num_rows++] = i;
                }
            }
            rows_km[num_rows++] = numUncens-1;
            
            double prod = 1;
            int c_idx = 0;
            for (int i=0; i < num_rows; i++) {
                double t_i = y_km[rows_km[i]];
                while(c_idx < numCens && c_km[c_idx] < t_i) {
                    N_i--;
                    c_idx++;
                }
                int d_i;
                if (i==0) {
                    d_i = rows_km[i]+1;
                    km_mean += prod * y_km[rows_km[i]];
                } else {
                    d_i = rows_km[i] - rows_km[i-1];
                    km_mean += prod * (y_km[rows_km[i]] - y_km[rows_km[i-1]]);
                }
                prod *= (N_i-d_i)*1.0/N_i;
                N_i -= d_i;
            }
            /* Deal with remaining censored values (if the highest value is
             * uncensored, then prod=0, so nothing happens then)
             */
            km_mean += prod * (upper_bound-y_km[numUncens-1])/2;
        }
        return km_mean;
    }
    
    private static double logrank_statistic(int num_periods, int[] N, int[] O, int[] N_1, int[] O_1) {
        /* Compute logrank statistic for the specified "numbers at risk" (N),
         * observed events (O), "numbers at risk" in group 1 (N_1), and observed
         * events in group 1 (O_1)
         *
         * If N(i) <= 1, we don't include it in the sum (undefined variance, 0/0)
         */

        double sum_V=0, numerator=0;
        for (int p=0; p<num_periods; p++) {
            int Np = N[p];
            if (Np > 1) {
                int N1p = N_1[p];
                double E = O[p] * ((N1p+0.0)/Np);
                double V = E * (1 - ((N1p+0.0)/Np)) * (Np-O[p]) / (Np-1.0);
                sum_V += V;
                numerator += (O_1[p]-E);
            }
        }

        double denominator = Math.sqrt(sum_V);
        
        if (Math.abs(denominator) < 1e-6) {
            if (Math.abs(numerator) < 1e-6) {
                return 0;
            } else {
                throw new RuntimeException("Division by zero in function logrank_statistic.");
            }
        }
        
        return numerator/denominator;
    }
    
    
    
    
    
    //=== Get the values of variable var's domain that it could potentially take at this node (the values that are *compatible* with the splits going to this node)
    private static int[] getCompatibleValues(int currnode, int var, int N, int[] parent, int[] cutvar, double[] cutpoint, int[] leftchildren, int[] rightchildren, int[][] catsplit, int[] catDomainSizes) {
        int[] compatibleValues = null;
        while (currnode > 0) {
            int parent_node = parent[currnode];
            if (-cutvar[parent_node]-1 == var) { // cutvar is negative for categorical splits
                int catsplit_index = (int)cutpoint[parent_node];

                if (leftchildren[parent_node] == currnode) {
                    compatibleValues = catsplit[catsplit_index];
                } else if (rightchildren[parent_node] == currnode) {
                    compatibleValues = catsplit[catsplit_index+N];
                } else {
                    throw new RuntimeException("currnode must be either left or right child of its parent.");
                }
                break;
            }
            currnode = parent_node;
        }
        if (currnode == 0) {
            compatibleValues = new int[catDomainSizes[var]];
            for (int i=0; i < compatibleValues.length; i++) compatibleValues[i] = i+1;
        }
        return compatibleValues;
    }
    
    
    //======================================================================\\
    //                        BEGIN HELPER FUNCTIONS                        \\
    //======================================================================\\
    private static void rankSort(double[] arr, int len, int[] sorder) {
        for (int i=0; i<len; i++) {
            sorder[i] = i;
        }
        dp_quick(arr, sorder, 0, len-1);
    }
    
    private static void shuffle(int[] arr, int n) {
        for (int i=0; i < n-1; i++) {
            int j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = arr[j];
            arr[j] = arr[i];
            arr[i] = t;
        }
    }
    
    private static void dp_quick(double[] input, int[] sorder, int min, int max) {
        while (max - min > 0) {
            int i = min, j = max;
            double pivot = input[sorder[(i+j) >> 1]];
            do {
                while(input[sorder[i]] < pivot) i++;
                while(input[sorder[j]] > pivot) j--;
                if (i >= j) break;
                int t = sorder[i];
                sorder[i] = sorder[j];
                sorder[j] = t;
            } while(++i < --j);

            while (min < j && input[sorder[j]] == pivot) j--;
            if (min < j) dp_quick(input, sorder, min, j);

            while (i < max && input[sorder[i]] == pivot) i++;
            min = i; // dp_quick(input, sorder, i, max);
        }
    }
    
    public static void serializeDataForDebug(double[][] allTheta, double[][] allX, int[][] dataIdxs, double[] y, boolean[] cens, RegtreeBuildParams params, String filename) {   
    	DebugTreeBuildInputs t = new DebugTreeBuildInputs();
    	t.allTheta = allTheta;
    	t.allX = allX;
    	t.dataIdxs = dataIdxs;
    	t.y = y;
    	t.cens = cens;
    	t.params = params;

	    try {
			FileOutputStream fileOut = new FileOutputStream(filename + ".ser");
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(t);
			out.close();
			fileOut.close();
		} catch (IOException i) {
			i.printStackTrace();
		}
    }
}