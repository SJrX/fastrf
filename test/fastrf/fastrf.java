package fastrf;

import static org.junit.Assert.*;

import java.util.Arrays;
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
		RegtreeBuildParams buildparams = new RegtreeBuildParams();
		buildparams.random = new Random(3);
		buildparams.splitMin = 1;
		int[] catDomainSizes = {0,0,0};
		buildparams.catDomainSizes = catDomainSizes;
		buildparams.minVariance = Math.pow(10,-14);
		buildparams.ratioFeatures = 1;
		//RandomForest forest = new RandomForest(1, buildparams);
		
		int numTrees = 10;
		double[][] allTheta = {{1.,1.}, {0.,0.} };
		
		//no features
		double[][] allX = {{0}, {0}};
		int[][] theta_inst_idxs = {{0,0},{1,1}};
		double[] y = {4., 2.};
				
	
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildparams);
	
		double[][] x = {{1.,1.,0.,}, {0.,0.,0.}};
		double[][] pred = RandomForest.apply(rf, x);
		System.out.println(Arrays.deepToString(pred));
		//fail("Not yet implemented");
		
	}

}
