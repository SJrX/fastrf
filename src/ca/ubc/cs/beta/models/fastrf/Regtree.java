package ca.ubc.cs.beta.models.fastrf;

import java.util.Arrays;

import ca.ubc.cs.beta.models.fastrf.utils.Utils;
import static ca.ubc.cs.beta.models.fastrf.utils.Hash.*;

public strictfp class Regtree implements java.io.Serializable {    
    private static final long serialVersionUID = -7861532246973394125L;
    
    public int numNodes;
    public int[] node;
    public int[] parent;
    public double[][] ysub;
    public int[] var;
    public double[] cut;
    public int[][] children;
    public int[] nodesize;
    public int npred;
    public int[][] catsplit;
    
    public double[] nodepred;
    public double[] nodevar;
    
    public boolean resultsStoredInLeaves;
    
    public boolean preprocessed;
    public double[] weightedpred;
    public double[] weightedvar;
    public double[] weights;
    
    public int logModel;

    public Regtree(int numNodes, int logModel) {
        this.numNodes = numNodes;
        this.logModel = logModel;
        preprocessed = false;
    }
    
    
    public boolean equals(Object o)
    {
    	if(o instanceof Regtree)
    	{
    		Regtree rt = (Regtree) o;
    		if (npred != rt.npred) return false;
    		if (resultsStoredInLeaves != rt.resultsStoredInLeaves) return false;
    		if (preprocessed != rt.preprocessed) return false;
    		if (logModel != rt.logModel) return false;
    		
    		if (!Arrays.equals(var,rt.var)) return false;
    		if (!Arrays.equals(cut,rt.cut)) return false;
    		
    		if (!Arrays.deepEquals(catsplit,rt.catsplit)) return false;
    		if (!Arrays.deepEquals(ysub,rt.ysub)) return false;
    		if (!Arrays.deepEquals(children,rt.children)) return false;
    		
    		if (!Arrays.equals(parent,rt.parent)) return false;
    		
    		
    		
    		if (!Arrays.equals(nodesize,rt.nodesize)) return false;
    		
    		

    		if (!Arrays.equals(nodepred,rt.nodepred)) return false;
    		if (!Arrays.equals(nodevar,rt.nodevar)) return false;

    		

    		
    		if (!Arrays.equals(weightedpred,rt.weightedpred)) return false;
    		if (!Arrays.equals(weightedvar,rt.weightedvar)) return false;
    		if (!Arrays.equals(weights,rt.weights)) return false;

    		return true;

    		
    	}
    	return false;
    	
    }
    
    
    
    public int hashCode()
    {
    	System.out.println("=======");
    	
    	int hash = 0;
    	
    	hash = hash(node) + 31*hash;
    	System.out.println("HASH1:"+hash(node));
    	Arrays.toString(node);
    	
    	hash = hash(parent) + 31*hash;
    	System.out.println("HASH2:"+hash(parent));
    	Arrays.toString(parent);
    	
    	hash = hash(var) + 31*hash;
    	System.out.println("HASH3:"+hash(var));
    	
    	
    	
    	hash = hash(cut) + 31*hash;
    	System.out.println("HASH4:"+hash(cut));
    	
    	hash = hash(nodesize) + 31*hash;
    	System.out.println("HASH5:"+hash(nodesize));
    	
    	hash = npred + 31*hash;
    	System.out.println("HASH6:"+npred);
    	
    	hash = hash(children) + 31*hash;
    	System.out.println("HASH7:"+hash(children));
    	
    	hash = hash(catsplit) + 31*hash;
    	System.out.println("HASH8:"+hash(catsplit));
    	
    	hash = hash(ysub) + 31*hash;
    	System.out.println("HASH9:"+hash(ysub));
    	
    	
    	
    	hash = hash(nodepred) + 31*hash;
    	System.out.println("HASH0:"+hash(nodepred));
    	
    	hash = hash(nodevar) + 31*hash;
    	System.out.println("HASHA:"+hash(nodevar));
    	
    	hash = (resultsStoredInLeaves?1:0) + 31*hash;
    	System.out.println("HASHB:"+resultsStoredInLeaves);
    	
    	hash = preprocessed?1:0 + 31*hash;
    	System.out.println("HASHC:"+preprocessed);
    	
    	hash = hash(weightedpred) + 31*hash;
    	System.out.println("HASHD:"+hash(weightedpred));
    	
    	hash = hash(weightedvar) + 31*hash;
    	System.out.println("HASHE:"+hash(weightedvar));
    	
    	hash = hash(weights) + 31*hash;
    	System.out.println("HASHF:"+hash(weights));
    	
    	hash = logModel + 31*hash;
    	System.out.println("HASHG:"+logModel);
    	
    	return hash;
    	
    	
    }
    
    
    public Regtree(int numNodes, int ncatsplit, boolean storeResultsInLeaves, int logModel) {
        this(numNodes, logModel);
        
        resultsStoredInLeaves = storeResultsInLeaves;
        
        node = new int[numNodes];
        parent = new int[numNodes];
        var = new int[numNodes];
        cut = new double[numNodes];
        children = new int[numNodes][2];
        nodesize = new int[numNodes];
        catsplit = new int[ncatsplit][];
        
        if (resultsStoredInLeaves) {
            ysub = new double[numNodes][];
        } else {
            ysub = new double[numNodes][2]; // sum, sumofsq
        }
    }
    
    public Regtree(Regtree t) {
        this(t.numNodes, t.catsplit.length, t.resultsStoredInLeaves, t.logModel);
        
        npred = t.npred;     

        System.arraycopy(t.node, 0, node, 0, numNodes);
        System.arraycopy(t.parent, 0, parent, 0, numNodes);
        System.arraycopy(t.var, 0, var, 0, numNodes);
        System.arraycopy(t.cut, 0, cut, 0, numNodes);
        System.arraycopy(t.nodesize, 0, nodesize, 0, numNodes);
        
        for (int i=0; i < t.catsplit.length; i++) {
            catsplit[i] = new int[t.catsplit[i].length];
            catsplit[i] = new int[t.catsplit[i].length];
            System.arraycopy(t.catsplit[i], 0, catsplit[i], 0, t.catsplit[i].length);
        }
        
        for (int i=0; i < numNodes; i++) {
            children[i][0] = t.children[i][0];
            children[i][1] = t.children[i][1];
            
            if (resultsStoredInLeaves) {
                int Nnode = (var[i] == 0 ? nodesize[i] : 0);
                if (Nnode != 0) {
                    ysub[i] = new double[Nnode];
                    System.arraycopy(t.ysub[i], 0, ysub[i], 0, Nnode);
                }
            } else {
                System.arraycopy(t.ysub[i], 0, ysub[i], 0, t.ysub[i].length);
            }
        }
        
        preprocessed = t.preprocessed;
        if (preprocessed) {
            weights = new double[numNodes];
            weightedpred = new double[numNodes];
			weightedvar = new double[numNodes];
            System.arraycopy(t.weights, 0, weights, 0, numNodes);
            System.arraycopy(t.weightedpred, 0, weightedpred, 0, numNodes);
			System.arraycopy(t.weightedvar, 0, weightedvar, 0, numNodes);
        }
        
        recalculateStats();
    }
    
    /**
     * Gets a prediction for the given instantiations of features
     * @returns a matrix of size X.length*2 where index (i,0) is the prediction for X[i] 
     * and (i,1) is the variance of that prediction.
     */
    public static double[][] apply(Regtree tree, double[][] X) {
        int[] nodes = RegtreeFwd.fwd(tree, X);
        
        double[][] retn = new double[X.length][2]; // mean, var
        for (int i=0; i < X.length; i++) {
            retn[i][0] = tree.nodepred[nodes[i]];
            retn[i][1] = tree.nodevar[nodes[i]];
        }
        return retn;
    }
    
    /**
     * Gets a prediction for the given configurations and instances
     * @see RegtreeFwd.marginalFwd
     */
    public static Object[] applyMarginal(Regtree tree, double[][] Theta, double[][] X) {
        return RegtreeFwd.marginalFwd(tree, Theta, X);
    }
    
    public static void update(Regtree tree, double[][] newx, double[] newy) {
		if (tree.preprocessed) {
            throw new RuntimeException("Cannot update preprocessed forests.");
        }
        if (null == newx || null == newy) {
            throw new RuntimeException("Input newx or newy to update is null.");
        }
        if (newx.length != newy.length) {
            throw new RuntimeException("Argument sizes mismatch.");
        }
        if (tree.logModel > 0){
            for (int i = 0; i < newy.length; i++){
                newy[i] = Math.pow(10,newy[i]);
            }
        }
        
        int[] nodes = RegtreeFwd.fwd(tree, newx);
        boolean[] nodeChanged = new boolean[tree.node.length];
        
        for (int i = 0; i < newx.length; i++) {
            int node = nodes[i];
            nodeChanged[node] = true;

            int Nnode = tree.nodesize[node];

            if (tree.resultsStoredInLeaves) {
                double[] newysub = new double[Nnode+1];
                if (Nnode != 0) {
                    System.arraycopy(tree.ysub[node], 0, newysub, 0, Nnode);
                }
                newysub[Nnode] = newy[i];
                tree.ysub[node] = newysub;
            } else {
                tree.ysub[node][0] += newy[i]; // sum
                tree.ysub[node][1] += newy[i]*newy[i]; // sum of squares
            }

            tree.nodesize[node]++;
            while(node != 0) {
                node = tree.parent[node];
                tree.nodesize[node]++;
            }
        }
        
        for (int i=0; i < nodeChanged.length; i++) {
            if (nodeChanged[i]) {
                tree.recalculateStats(i);
            }
        }
    }

    /**
     * Recalculate statistic (mean, var) of the entire tree
     */
    public void recalculateStats() {
        nodepred = new double[numNodes];
        nodevar = new double[numNodes];
        
        for (int i=0; i < numNodes; i++) {
            recalculateStats(i);
        }
    }
    
    /**
     * Recalculate statistic (mean, var) of the specified node
     */
    public void recalculateStats(int node) {
        if (var[node] != 0) return;
        
        if (resultsStoredInLeaves) {
            nodepred[node] = Utils.mean(ysub[node]);
            nodevar[node] = Utils.var(ysub[node]);
        } else {
            double sum = ysub[node][0], sumOfSq = ysub[node][1];
            int N = nodesize[node];
            
            nodepred[node] = sum/N;
            nodevar[node] = (sumOfSq - sum*sum/N) / Math.max(N-1, 1);
        }
    }
}
