package de.unifreiburg.cs.junit;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Random;

import org.junit.Test;
import static org.junit.Assert.assertTrue;

import ca.ubc.cs.beta.models.fastrf.RandomForest;
import ca.ubc.cs.beta.models.fastrf.RegtreeBuildParams;
import ca.ubc.cs.beta.models.fastrf.utils.CsvToDataConverter;
import ca.ubc.cs.beta.models.fastrf.utils.RfData;

public class DataTester {
	public static String rfDeployment = null;
	public static Long seedOffset = Math.abs((new Random()).nextLong());

	@Test
	public void testMakeUnique() throws IOException{
		int[] thetaColIdxs = {0};
		int[] xColIdxs = {1};
		int yColIdx= 2;
		int[] catColIdxs = {};
		String filename = "test_files/mini_with_repetitions.csv";
		
		CsvToDataConverter converter = new CsvToDataConverter(filename, thetaColIdxs, xColIdxs, yColIdx, catColIdxs);
		RfData data = converter.readDataFromCsvFile(filename);
		assertEquals(data.getTheta().length, 5);
		assertEquals(data.getX().length, 5);
		
		data.makeUnique(true);
		assertEquals(data.getTheta().length, 2);
		assertEquals(data.getX().length, 5);
		
		data.makeUnique(false);
		assertEquals(data.getTheta().length, 2);
		assertEquals(data.getX().length, 3);

		data.makeUnique(true);
		assertEquals(data.getTheta().length, 2);
		assertEquals(data.getX().length, 3);

		data.makeUnique(false);
		assertEquals(data.getTheta().length, 2);
		assertEquals(data.getX().length, 3);		
	}
	
	@Test
	public void testRfFromCSV() throws IOException{
		int[] thetaColIdxs = {0};
		int[] xColIdxs = {1};
		int yColIdx= 2;
		int[] catColIdxs = {};
		String filename = "test_files/mini.csv";
		
		CsvToDataConverter converter = new CsvToDataConverter(filename, thetaColIdxs, xColIdxs, yColIdx, catColIdxs);
		RfData data = converter.readDataFromCsvFile(filename);
		
		RandomForest rf = RandomForest.learnModel(1, data.getTheta(), data.getX(), data.getTheta_inst_idxs(), data.getY(), new RegtreeBuildParams(2, false, 1));

		//=== Predictions on the training data must be exact.
		double[][] newX = data.buildMatrixForApply();
		double[][] meanvar = RandomForest.apply(rf, newX);
		for (int i = 0; i < data.getY().length; i++) {
			assertEquals(meanvar[i][0], data.getY()[i], 1e-10);
		}
	}

	
	/* 
	 * Assert that predictions of a deterministic tree with splitMin 1 are exact on a small set of unique entries. 
	 */
	@Test
	public void testRfFromMiniLDOF() throws IOException{
		int[] thetaColIdxs = {2,3,4,5};
		int[] xColIdxs = new int[94];
		for (int i = 0; i < xColIdxs.length; i++) {
			xColIdxs[i] = i+6;
		}
		int yColIdx= 101;
		int[] catColIdxs = {};
		//frame,smac configuration id,parameter1,parameter2,parameter3,parameter4,feature1,feature2,feature3,feature4,feature5,feature6,feature7,feature8,feature9,feature10,feature11,feature12,feature13,feature14,feature15,feature16,feature17,feature18,feature19,feature20,feature21,feature22,feature23,feature24,feature25,feature26,feature27,feature28,feature29,feature30,feature31,feature32,feature33,feature34,feature35,feature36,feature37,feature38,feature39,feature40,feature41,feature42,feature43,feature44,feature45,feature46,feature47,feature48,feature49,feature50,feature51,feature52,feature53,feature54,feature55,feature56,feature57,feature58,feature59,feature60,feature61,feature62,feature63,feature64,feature65,feature66,feature67,feature68,feature69,feature70,feature71,feature72,feature73,feature74,feature75,feature76,feature77,feature78,feature79,feature80,feature81,feature82,feature83,feature84,feature85,feature86,feature87,feature88,feature89,feature90,feature91,feature92,feature93,feature94,feature95,performance

		String trainFilename = "test_files/mini_ldof.csv";
		CsvToDataConverter converter = new CsvToDataConverter(trainFilename, thetaColIdxs, xColIdxs, yColIdx, catColIdxs);
		RfData trainData = converter.readDataFromCsvFile(trainFilename);
		//		trainData.makeUnique(true);
		//		trainData.makeUnique(false);

		RegtreeBuildParams regTreeBuildParams = new RegtreeBuildParams(false, 1, converter.getCatDomainSizes());  // , true, 10);
		RandomForest rf = RandomForest.learnModel(1, trainData.getTheta(), trainData.getX(), trainData.getTheta_inst_idxs(), trainData.getY(), regTreeBuildParams);

		//=== Predictions on the training data must be exact.
		double[][] meanvar = RandomForest.apply(rf, trainData.buildMatrixForApply());
		for (int i = 0; i < trainData.getY().length; i++) {
			assertEquals(meanvar[i][0], trainData.getY()[i], 1e-10);
		}
	}
	
	/* 
	 * Assert that predictions of an RF are better than the mean predictor. 
	 */
	@Test
	public void testBetterThanMeanPredictions() throws IOException{
		int[] thetaColIdxs = {2,3,4,5};
		int[] xColIdxs = new int[94];
		for (int i = 0; i < xColIdxs.length; i++) {
			xColIdxs[i] = i+6;
		}
		int yColIdx= 101;
		int[] catColIdxs = {};
		//frame,smac configuration id,parameter1,...,parameter4,feature1,...,feature95,performance
		//String trainFilename = "test_files/LDOF_traindata_subset.csv";
		String trainFilename = "test_files/train_ldof_data_shuffled_first1000.csv";
		CsvToDataConverter converter = new CsvToDataConverter(trainFilename, thetaColIdxs, xColIdxs, yColIdx, catColIdxs);
		RfData trainData = converter.readDataFromCsvFile(trainFilename);
	
		//=== Learn RF.
		RegtreeBuildParams regTreeBuildParams = new RegtreeBuildParams(true, 10, converter.getCatDomainSizes());
		RandomForest rf = RandomForest.learnModel(10, trainData.getTheta(), trainData.getX(), trainData.getTheta_inst_idxs(), trainData.getY(), regTreeBuildParams);
	
		//=== Determine mean predictor.
		double rmse = 0;
		double rmseOfMeanPred = 0;
		double meanPred = 0; 
		for (int i = 0; i < trainData.getY().length; i++) {
			meanPred += trainData.getY()[i]; 
		}
		meanPred /= trainData.getY().length;
				
		//=== Load test set.
		String testFilename = "test_files/test_ldof_data_shuffled_first10000.csv";
		RfData testData = converter.readDataFromCsvFile(testFilename);
		testData.makeUnique(true);
		testData.makeUnique(false);
	
		//=== Assert that our model is better than the mean predictor.
		double[][] meanvar = RandomForest.apply(rf, testData.buildMatrixForApply());
		for (int i = 0; i < testData.getY().length; i++) {
			rmse += Math.pow( (testData.getY()[i]-meanvar[i][0]), 2); 
			rmseOfMeanPred += Math.pow( (testData.getY()[i] - meanPred), 2);
		}
		rmse = Math.sqrt(rmse/testData.getY().length);
		rmseOfMeanPred = Math.sqrt(rmseOfMeanPred/testData.getY().length);
		System.out.println("RMSE = " + rmse + "; rmse of mean pred: " + rmseOfMeanPred);
		assertTrue(rmse < rmseOfMeanPred);
	}
	
}