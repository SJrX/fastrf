package ca.ubc.cs.beta.models.fastrf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Arrays;

/**
 * @param catDomainSizes a vector of size X[0].length indicating the size of the domain for the corresponding categorical feature, or 0 if the feature is continuous.
 *                       catDomainSizes should be the length of the number of Parameters + Number of Features.
 * @param condParents the conditional parents of the given variable (column index into X). Empty if no parents. NOT SORTED
 * @param condParentVals the OK values of all the conditional parents of the given variable (each element is a matrix). Indexed wrt condParents.
 * @param splitMin the minimum number of data points in each node.
 * @param ratioFeatures the percentage (0, 1] of the features to use in deciding the split at each node.
 * @param logModel whether to build a model in log space
 * @param storeResponses whether to store the responses themselves in the leaves or just to store some statistic (sum, sumofsquares, and leaf size)
 * @param seed -1 means don't use a seed (i.e. create a new Random but don't call setSeed). 
 * @param random will be used instead of seed if it's not null.
 * @param minVariance - Minimum Variance value that will ever be returned on apply call
 */
public strictfp class RegtreeBuildParams implements java.io.Serializable {    
	public RegtreeBuildParams(int numVars, boolean doBootstrapping) {
		this(numVars, doBootstrapping, 10);
	}
	
	public RegtreeBuildParams(int numVars, boolean doBootstrapping, int splitMin) {
		this(numVars, doBootstrapping, splitMin, 2.0/3);
	}
	
	public RegtreeBuildParams(boolean doBootstrapping, int splitMin, int[] catDomainSizes) {
		this( doBootstrapping, splitMin, 2.0/3, catDomainSizes);
	}
	
	public RegtreeBuildParams(int numVars, boolean doBootstrapping, int splitMin, double ratioFeatures) {
		this( doBootstrapping, splitMin, ratioFeatures, new int[numVars]); //initializes catDomainSizes to all zero
	}
	
	public RegtreeBuildParams(boolean doBootstrapping, int splitMin, double ratioFeatures, int[] catDomainSizes) {
		super();
		this.catDomainSizes = catDomainSizes; 
		this.doBootstrapping = doBootstrapping;
		this.splitMin = splitMin;
		this.ratioFeatures = ratioFeatures;
		this.random = new Random(3);
	}
	

	private static final long serialVersionUID = -6803645785543626390L;
	public int[] catDomainSizes;
    
	//DEPRECATED
	//public int[][] condParents = null;
    //public int[][][] condParentVals = null;
    
    /*
     * mapping from variable index to disjunctions of condition clauses (i.e., conjunctions) with parent variable index
     */
    public Map<Integer, int[][]> nameConditionsMapParentsArray;
    /*
     * mapping from variable index to disjunctions of condition clauses (i.e., conjunctions) with parent variable value*s* (based on conditional operator, only one value (>,<,!=,==) or several values ("in") 
     */
    public Map<Integer, double[][][]> nameConditionsMapParentsValues;
    /*
     * mapping from variable index to dijunctions of condition clauses (i.e., conjunctions) with conditional operators (<,>,!=,==,"in") 
     */
    public Map<Integer, int[][]> nameConditionsMapOp;
    
    public double minVariance;
    public String toString()
    {
        try{
        StringBuffer sb = new StringBuffer();
        sb.append("catDomain: (" + ((catDomainSizes != null) ? catDomainSizes.length : "null") +")");
        sb.append(Arrays.toString(catDomainSizes));
        sb.append("\n");
        sb.append("\n");
        //sb.append("condParent: (" + ((condParents != null) ? condParents.length : "null") +")");
        //sb.append(Arrays.deepToString(condParents));
        sb.append("\n");
        //sb.append("condParentVals: (" + ((condParentVals != null) ? condParentVals.length : "null") + ")");
        //sb.append(Arrays.deepToString(condParentVals));
        sb.append("\nSplitMin:" + splitMin);
        sb.append("\nRatioFeatures:" + ratioFeatures);
        //sb.append("\nCutOffPenaltyFactor:" + cutoffPenaltyFactor);
        sb.append("\nLogModel:" + logModel);
        sb.append("\nStoreResponses:" + storeResponses);
       
        return sb.toString();
        } catch(RuntimeException e)
        {
            e.printStackTrace();
            return "RuntimeException occured building RegtreeBuildPalams";
        }
        
        
    }
    
    public int splitMin;
    public void setSplitMin(int splitMin) {
		this.splitMin = splitMin;
	}

    public boolean doBootstrapping = true;

	public double ratioFeatures;
    public int logModel;
    public boolean storeResponses;
    
    public long seed = -1;
    public Random random = null;
	public boolean brokenVarianceCalculation = true;

	public void setLogModel(int logModel) {
		this.logModel = logModel;
	}

	public static RegtreeBuildParams copy(RegtreeBuildParams bp, int size)
	{
			RegtreeBuildParams bpNew = new RegtreeBuildParams(size, bp.doBootstrapping);
			
			bpNew.logModel = bp.logModel;
			bpNew.brokenVarianceCalculation = bp.brokenVarianceCalculation;
			
			bpNew.catDomainSizes = new int[size];
			//bpNew.condParents = new int[size][];
			//bpNew.condParentVals = new int[size][][];
			System.arraycopy(bp.catDomainSizes, 0,bpNew.catDomainSizes, 0, size);
			//System.arraycopy(bp.condParents, 0,bpNew.condParents, 0, size);
			//System.arraycopy(bp.condParentVals, 0,bpNew.condParentVals, 0, size);
			bpNew.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>(bp.nameConditionsMapParentsArray);
			bpNew.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>(bp.nameConditionsMapParentsValues);
			bpNew.nameConditionsMapOp = new HashMap<Integer, int[][]>(bp.nameConditionsMapOp);
			System.arraycopy(bp.topolical_order_vars, 0, bpNew.topolical_order_vars, 0, size);
			
			bpNew.random = bp.random;
			bpNew.ratioFeatures = bp.ratioFeatures;
			bpNew.minVariance = bp.minVariance;
			bpNew.splitMin = bp.splitMin;
			bpNew.storeResponses = bp.storeResponses;
			
			bpNew.seed = bp.seed;
			return bpNew;
	}
	
	public int[] topolical_order_vars;
	
	/*
	 * generates once a topological order of the variables, such that
	 * X is listed before Y, if there is a conditional of the form "Y | X"
	 */
	public void generate_topological_order(){
		if (this.topolical_order_vars != null) { // do this only once
			return;
		}
		
        if (this.nameConditionsMapParentsArray == null || this.nameConditionsMapParentsArray.isEmpty()) {
        	int nvars = this.catDomainSizes.length;
        	this.topolical_order_vars = new int[nvars];
            for (int i=0; i <nvars ; i++) {
            	topolical_order_vars[i] = i;
            }
        }
        else {
        	int nvars = this.catDomainSizes.length;
        	int running_index = 0;
        	int iterations = 0;
        	this.topolical_order_vars = new int[nvars];
        	List<Integer> remaining = new ArrayList<Integer>();
            for (int i=0; i <nvars ; i++) {
            	remaining.add(i);
            }
        	List<Integer> assigned = new ArrayList<Integer>();
        	while (running_index != nvars && iterations < 1000) {
        		iterations++;
        		List<Integer> to_be_removed = new ArrayList<Integer>();
        		for (int var : remaining){
        			if (this.nameConditionsMapParentsArray.get(var) == null) { // no parents
						this.topolical_order_vars[running_index++] = var;
        				assigned.add(var);
        				to_be_removed.add(var);
            		}
        			else {
        				// maybe more restrictive than necessary
        				// potentially, only the parents of one clause have to active in the end
        				// but this would change the topo order based on the assigned parents
        				// and we would need to recompute this order for each split which is ineffective
        				boolean all_parents = true;
                		for(int i=0; i < this.nameConditionsMapParentsArray.get(var).length; i++) { //iterate over disjunctions
                			for(int j=0; j < this.nameConditionsMapParentsArray.get(var)[i].length; j++) { //iterate over conjunctions
                				Integer parent_idx = this.nameConditionsMapParentsArray.get(var)[i][j];
                				if (! assigned.contains(parent_idx)) {
                					all_parents = false;
                					break;
                				}
                			}
                			if (! all_parents){
                				break;
                			}
                		}
                		if (all_parents){
                			this.topolical_order_vars[running_index++] = var;
            				assigned.add(var);
            				to_be_removed.add(var);
                		}
        			}
        		}
        		for (Integer var_rm: to_be_removed){
        			remaining.remove(var_rm);
        		}
        	}
        	if (iterations == 1000){
        		throw new IllegalArgumentException("Could not find topological order of variables -- probably there is a cylce in the parameters' conditionals.");
        	}
        	
        }
        
     }
	
}
