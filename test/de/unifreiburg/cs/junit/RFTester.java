package de.unifreiburg.cs.junit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ca.ubc.cs.beta.models.fastrf.RandomForest;
import ca.ubc.cs.beta.models.fastrf.RegtreeBuildParams;

import com.opencsv.CSVReader;

import org.junit.*;

import com.sun.xml.internal.ws.client.ClientSchemaValidationTube;

import static org.junit.Assert.*;

public class RFTester {
	public static String rfDeployment = null;
	public static Long seedOffset = Math.abs((new Random()).nextLong());

	//Boiler plate copied from SMAC, but does not work here...
	/*
  	@BeforeClass
	public static void initDeploymentToTest()
	{

		String s = ClassLoader.getSystemClassLoader().getResource("lastbuild-deploy.txt").getFile();
		//System.out.println("CL Name:"+ ClassLoader.getSystemClassLoader().getClass().getName());
		if(s == null || s.trim().length() == 0)
		{
			throw new AssertionError("Could not find deployment file lastbuild-deploy.txt on classpath");
		}
		File f = new File(s);
		try {
			BufferedReader r = new BufferedReader(new FileReader(f));
			rfDeployment = r.readLine();
			System.out.println("RF DEPLOYMENT: " + rfDeployment);
			r.close();
		} catch (FileNotFoundException e) {
			throw new AssertionError("Could open the deployment file lastbuild-deploy.txt");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new AssertionError(e);
		}

		//System.out.println(s);

		//System.out.println(f.getAbsolutePath())
	}
*/

	
	@Test
	public void testRfFromCSV() throws IOException{
		CSVReader reader = new CSVReader(new FileReader("test_files/mini_csv.csv"));
		List<String[]> lines = reader.readAll();
		reader.close();

		String[] header = lines.get(0);
		
		
		double y[] = new double[lines.size()];
		double Theta[][] = new double[lines.size()][];
		for (String[] line : lines) {
			for (int i = 0; i < line.length; i++) {
				System.out.print(line[i] + " ");
			}
			System.out.println();			
		}
	}
	
	@Test
	public void testRfPredictionForAllThetaNull(){
		double[][] allX = new double[10][1];
		for (int i = 0; i < allX.length; i++) {
			allX[i][0] = i;
		}
		
		double[] y = new double[allX.length];
		for (int i = 0; i < allX.length; i++) {
			y[i] = 10*i;
		}
		
		int[][] theta_inst_idxs = new int[allX.length][2];
		for (int i = 0; i < allX.length; i++) {
			theta_inst_idxs[i][0] = 0;
			theta_inst_idxs[i][1] = i;
		}
		
		RegtreeBuildParams buildParams = new RegtreeBuildParams(1, false, 1);
				
		int numTrees = 10;
		RandomForest rf = RandomForest.learnModel(numTrees, null, allX, theta_inst_idxs, y, buildParams);
				
		double[][] meanvar = RandomForest.apply(rf, allX);
		for (int i = 0; i < meanvar.length; i++) {
			assertEquals(meanvar[i][0], y[i], 1e-10);
		}
	}
	
	/* In the following test, there are two occurrences for each element of X,  
	 * which the RF cannot put into separate leaves,  
	 * so it has to predict the mean across the two corresponding y values.
	 */
	@Test
	public void testRfPredictionNonUniqueX(){
		double[][] allTheta = new double[1][1];
		allTheta[0][0] = 1;
		int modFactor = 500;
	
		double[][] allX = new double[1000][1];
		for (int i = 0; i < allX.length; i++) {
			allX[i][0] = i%modFactor;
		}
		
		double[] y = new double[allX.length];
		for (int i = 0; i < allX.length; i++) {
			y[i] = 10*i;
		}
		
		int[][] theta_inst_idxs = new int[allX.length][2];
		for (int i = 0; i < allX.length; i++) {
			theta_inst_idxs[i][0] = 0;
			theta_inst_idxs[i][1] = i;
		}
		
		RegtreeBuildParams buildParams = new RegtreeBuildParams(2, false, 1);
				
		int numTrees = 1;
		RandomForest rf = RandomForest.learnModel(numTrees, allTheta, allX, theta_inst_idxs, y, buildParams);
		
		
		double[][] combinedMatrixForTest = new double[theta_inst_idxs.length][2];
		for (int i = 0; i < combinedMatrixForTest.length; i++) {
			combinedMatrixForTest[i][0] = allTheta[theta_inst_idxs[i][0]][0];
			combinedMatrixForTest[i][1] = allX[theta_inst_idxs[i][1]][0];
		}
		
		double[][] meanvar = RandomForest.apply(rf, combinedMatrixForTest);
		for (int i = 0; i < meanvar.length; i++) {
//			System.out.println( meanvar[i][0] + " " + (y[i%modFactor] + y[modFactor + i%modFactor])/2);
			assertEquals(meanvar[i][0], (y[i%modFactor] + y[modFactor + i%modFactor])/2, 1e-10);
		}
	}
	
	
}

