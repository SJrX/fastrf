package fastrf;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import ca.ubc.cs.beta.models.fastrf.RandomForest;
import ca.ubc.cs.beta.models.fastrf.RegtreeBuildParams;

public class fastrf {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	public enum A {
		a, b
	};
	
	public enum B {
		a, b
	};

	@Test
	public void test() {
		
		
		//RegtreeBuildParams(boolean doBootstrapping, int splitMin, double ratioFeatures, int[] catDomainSizes)
		

		//buildparams.splitMin = 1;
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.random = new Random(3);
		buildparams.catDomainSizes = catDomainSizes;
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "1" | "0" == 1.0
		/*
		int[][] op = {{0}}; //EQ
		buildparams.nameConditionsMapOp.put(1, op);
		int[][] par = {{0}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(1, par);
		double[][][] val = {{{1}}}; //value 1
		buildparams.nameConditionsMapParentsValues.put(1, val);
		 */
		
		// "0" | "1" < 2 // should be always active - no influence on tree
		/*
		int[][] op2 = {{2}}; //LE
		buildparams.nameConditionsMapOp.put(0, op2);
		int[][] par2 = {{1}}; //Parent second theta
		buildparams.nameConditionsMapParentsArray.put(0, par2);
		double[][][] val2 = {{{2.}}}; //value 2
		buildparams.nameConditionsMapParentsValues.put(0, val2);
		*/
		
		// "0" | "1" < -1  -> never split on "0"
		/*int[][] op2 = {{2}}; //LE
		buildparams.nameConditionsMapOp.put(0, op2);
		int[][] par2 = {{1}}; //Parent second theta
		buildparams.nameConditionsMapParentsArray.put(0, par2);
		double[][][] val2 = {{{-1}}}; //value 2
		buildparams.nameConditionsMapParentsValues.put(0, val2);
		*/
		
		// "1" | "1" < -1  -> never split on "1"
		/*int[][] op3 = {{2}}; //LE
		buildparams.nameConditionsMapOp.put(1, op3);
		int[][] par3 = {{1}}; //Parent second theta
		buildparams.nameConditionsMapParentsArray.put(1, par3);
		double[][][] val3 = {{{-1}}}; //value 2
		buildparams.nameConditionsMapParentsValues.put(1, val3);
		*/
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {1}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
	}
	

	@Test
	/*
	 * test EQ by adding a "==" condition that disables parameter "0" always;
	 * since "0" is the only informative variable, always the same value should be predicted 
	 */
	public void testEQContNonActive() {

		
		
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.random = new Random(3);
		//buildparams.catDomainSizes = catDomainSizes;
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" == 1.0
		int[][] op = {{0}}; //EQ
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{1.}}}; //value 1
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("EQ NonActive Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] == pred[1][0]);
		assertTrue(pred[2][0] == pred[1][0]);
	}
	
	@Test
	/*
	 * test NEQ by adding a "!=" condition that does nothing;
	 */
	public void testNEQCont() {
	
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.catDomainSizes = catDomainSizes;
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" != 2.0
		int[][] op = {{1}}; //NEQ
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{2.}}}; //value 0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("NEQ Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] ==  pred[2][0]);
		assertTrue(pred[1][0] != pred[2][0]);
	}
	
	@Test
	/*
	 * test NEQ by adding a "!=" condition that disables parameter "0" always;
	 * since "0" is the only informative variable, always the same value should be predicted 
	 */
	public void testNEQContNonActive() {
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" != 0.0
		int[][] op = {{1}}; //NEQ
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{0.}}}; //value 0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("NEQ NonActive Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] == pred[1][0]);
		assertTrue(pred[2][0] == pred[1][0]);
	}
	
	@Test
	/*
	 * test LE by adding a "<" condition that does nothing;
	 */
	public void testLECont() {
		
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.catDomainSizes = catDomainSizes;
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" < 1.0
		int[][] op = {{2}}; //LE
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{1.}}}; //value 0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("LE Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] ==  pred[2][0]);
		assertTrue(pred[1][0] != pred[2][0]);
	}
	
	@Test
	/*
	 * test LE by adding a "<" condition that disables parameter "0" always;
	 * since "0" is the only informative variable, always the same value should be predicted 
	 */
	public void testLEContNonActive() {
		
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" < 0.0
		int[][] op = {{2}}; //LE
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{0.}}}; //value 0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("LE NonActive Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] == pred[1][0]);
		assertTrue(pred[2][0] == pred[1][0]);
	}
	
	@Test
	/*
	 * test GR by adding a ">" condition that does nothing
	 */
	public void testGRCont() {
		
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" > -0.5
		int[][] op = {{3}}; //GR
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{-0.5}}}; //value 1
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("GR Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] ==  pred[2][0]);
		assertTrue(pred[1][0] != pred[2][0]);
	}
	
	@Test
	/*
	 * test GR by adding a ">" condition that disables parameter "0" always;
	 * since "0" is the only informative variable, always the same value should be predicted 
	 */
	public void testGRContNonActive() {
		
		int[] catDomainSizes = {2,0,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" > 1.0
		int[][] op = {{3}}; //GR
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{1.}}}; //value 1
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,0.}, {2.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("GR NonActive Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] == pred[1][0]);
		assertTrue(pred[2][0] == pred[1][0]);
	}
	
	@Test
	/*
	 * test IN by adding a "in" condition that disables parameter "0" always;
	 * since "0" is the only informative variable, always the same value should be predicted 
	 */
	public void testINContNonActive() {
		
		
		int[] catDomainSizes = {2,1,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" in {2.0,3.0}
		int[][] op = {{4}}; //IN
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{2.0,3.0}}}; //value 2.0,3.0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,1.}, {2.,1.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("IN NonActive Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] == pred[1][0]);
		assertTrue(pred[2][0] == pred[1][0]);
	}
	
	@Test
	/*
	 * test IN by adding a "in" condition that does nothing;
	 */
	public void testINCont() {
		
		//theta 0 has two categorical values: 1,2; all others are continuous
		int[] catDomainSizes = {2,1,0};
		RegtreeBuildParams buildparams = new RegtreeBuildParams(false, 1, 1, catDomainSizes);
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		
		buildparams.nameConditionsMapOp = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsArray = new HashMap<Integer, int[][]>();
		buildparams.nameConditionsMapParentsValues = new HashMap<Integer, double[][][]>();
		
		// "0" | "1" in {1.0}
		int[][] op = {{4}}; //IN
		buildparams.nameConditionsMapOp.put(0, op);
		int[][] par = {{1}}; //Parent first theta
		buildparams.nameConditionsMapParentsArray.put(0, par);
		double[][][] val = {{{1.}}}; //value 0.0
		buildparams.nameConditionsMapParentsValues.put(0, val);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,1.}, {2.,1.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1},{0,1},{1,0}};
		double[] y = {4., 1., 5., 2.};
				
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {2.,0.,0.}, {1.,1.,1.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println("IN Test:");
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
		assertTrue(pred[0][0] ==  pred[2][0]);
		assertTrue(pred[1][0] != pred[2][0]);
	}

}
