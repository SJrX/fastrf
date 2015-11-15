package ca.ubc.cs.beta.models.fastrf;

import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import static org.junit.Assert.*;


/**
 * Created by Steve Ramage <seramage@cs.ubc.ca> on 11/14/15
 */
public class RegtreeFitTestCategoricalDomain {


    public static RegtreeBuildParams getRegtreeBuildParamsSplitMinThreeSingleThetaSingleFeature()
    {


        return getRegtreeBuildParamsSingleThetaSingleFeature(3, 0.83333333);

    }


    public static RegtreeBuildParams getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature()
    {


        return getRegtreeBuildParamsSingleThetaSingleFeature(1, 0.83333333);

    }

    public static RegtreeBuildParams getRegtreeBuildParamsSingleThetaSingleFeature(int splitMin, double ratioFeatures)
    {


        RegtreeBuildParams rbp = new RegtreeBuildParams(2,false);

        rbp.splitMin = splitMin;
        rbp.seed = new Random().nextInt();

        System.err.println("Seed for run was: " + rbp.seed);
        rbp.ratioFeatures = ratioFeatures;

        rbp.catDomainSizes = new int[]{10, 0};

        rbp.nameConditionsMapOp = new HashMap<>();
        rbp.nameConditionsMapParentsArray = new HashMap<>();
        rbp.nameConditionsMapParentsValues = new HashMap<>();


        rbp.logModel = 0;
        rbp.storeResponses = false;

        return rbp;

    }

    public static RegtreeBuildParams getRegtreeBuildParamsDoubleThetaFirstParamBinaryConditionalSingleFeature()
    {


        RegtreeBuildParams rbp = new RegtreeBuildParams(3,false);

        rbp.splitMin = 1;
        rbp.seed = new Random().nextInt();

        System.err.println("Seed for run was: " + rbp.seed);
        rbp.ratioFeatures = 0.83333333333;

        rbp.catDomainSizes = new int[]{2,10, 0};

        //rbp.condParents = new int[][]{{},{0},{}};
        //rbp.condParentVals = new int[][][]{{},{{2}},{}};


        rbp.nameConditionsMapOp = new HashMap<>();
        rbp.nameConditionsMapOp.put(1,new int[][]{{4}});

        rbp.nameConditionsMapParentsArray = new HashMap<>();
        rbp.nameConditionsMapParentsArray.put(1,new int[][]{{0}});

        rbp.nameConditionsMapParentsValues = new HashMap<>();
        rbp.nameConditionsMapParentsValues.put(1,new double[][][]{{{1.0}}});


        /*
        		System.out.println(Arrays.deepToString(pcs.getNameConditionsMapOp().get(1)));
		System.out.println(Arrays.deepToString(pcs.getNameConditionsMapParentsArray().get(1)));
		System.out.println(Arrays.deepToString(pcs.getNameConditionsMapParentsValues().get(1)));


         */
        rbp.logModel = 0;
        rbp.storeResponses = false;

        return rbp;

    }



    @Test
    /**
     * The expected tree here should be a single constant point.
     */
    public void testRegTreeFitSingleCategoricalPoint()
    {
        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 2 } };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}};
        double[] y = { 0.5};


        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 1 node", 1,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node prediction should be ~0.5", 0.5, tree.nodepred[0], 0.000001);
        assertEquals("Tree root node variance should be ~0", 0, tree.nodevar[0], 0.000001 );
        assertEquals("Node size of root should be 1", 1, tree.nodesize[0]);
    }

    @Test
    /**
     * The expected tree here should be a single point with some variance.
     */
    public void testRegTreeFitTwoObservationsAtSamePoint()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 2} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {0,0}};
        double[] y = { 0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 1 node", 1,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node prediction should be ~0.5", 0.5, tree.nodepred[0], 0.000001);

        /**
         * Sample Variance is ((2/3-1/2)^2+(1/3-1/2)^2)*1/(2-1)
         */
        assertEquals("Tree root node variance should be ~0", 0.0555555555, tree.nodevar[0], 0.000001 );
    }

    @Test
    /**
     * The expected tree here should have depth two with one node in each leaf.
     */
    public void testRegTreeFitTwoDistinctSamplesSameFeature()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 2}, {9} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0}};
        double[] y = { 0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 3 nodes", 3,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node should have split on first variable",-1, tree.var[0]);

        assertEquals("Second value should go to the first child", tree.catsplit[0][1] , 0);

        assertEquals("Ninth value should go to the second child", tree.catsplit[0][8] , 1);


        //Not sure what cut has here
        // assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.20) );
        assertEquals("Node size of root should be 2", 2, tree.nodesize[0]);


        assertEquals("Tree root node should have first child as node 1", tree.children[0][0], 1);
        assertEquals("Node size of first child should be 1", 1, tree.nodesize[1]);
        assertEquals("Tree first child node prediction should be ~0.3333333333", 0.3333333333, tree.nodepred[1], 0.000001);

        assertEquals("Tree root node should have second child as node 2", tree.children[0][1], 2);
        assertEquals("Node size of second child should be 1", 1, tree.nodesize[2]);
        assertEquals("Tree second child node prediction should be ~0.666666666", 0.6666666666, tree.nodepred[2], 0.000001);

    }


    @Test
    /**
     * The expected tree here should have depth two with one node in each leaf.
     *
     * Even values should go one way and odd to another.
     */
    public void testRegTreeFitTenDistinctSamplesSameFeature()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 1}, {2}, {3}, {4} ,{5},{6},{7},{8},{9},{10} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0} , {2,0},{3,0},{4,0},{5,0},{6,0},{7,0},{8,0},{9,0}};
        double[] y = { 0.3333333333, 0.6666666666, 0.3333333333, 0.6666666666,0.3333333333, 0.6666666666,0.3333333333, 0.6666666666,0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 3 nodes", 3,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node should have split on first variable",-1, tree.var[0]);


        // I only test that one valid case of the invariant, I am not sure why in the code this seems to work.
        // (First/Odd categorical value, goes to first child and Second/Even goes to second), the other would be valid
        // but it never came up. If this failed, well then also check for {1,0...}.
        assertTrue("Even nodes should go one way and odd nodes the other way",Arrays.equals(tree.catsplit[0] , new int[]{0,1,0,1,0,1,0,1,0,1}));


        //Not sure what cut has here
        // assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.20) );
        assertEquals("Node size of root should be 2", 10, tree.nodesize[0]);


        assertEquals("Tree root node should have first child as node 1", tree.children[0][0], 1);
        assertEquals("Node size of first child should be 1", 5, tree.nodesize[1]);
        assertEquals("Tree first child node prediction should be ~0.3333333333", 0.3333333333, tree.nodepred[1], 0.000001);

        assertEquals("Tree root node should have second child as node 2", tree.children[0][1], 2);
        assertEquals("Node size of second child should be 1", 5, tree.nodesize[2]);
        assertEquals("Tree second child node prediction should be ~0.666666666", 0.6666666666, tree.nodepred[2], 0.000001);

    }





    @Test
    /**
     * The expected tree here should have depth three with a structure:
     *        0
     *       / \
     *      1  2
     *        / \
     *       3  4
     *
     *
     *
     */
    public void testRegTreeFitThreeSampleOneCategoricalSameFeature()
    {
        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 1}, {2}, {3}};

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0}, {2,0}};

        /**
         * We are essentially modelling y = 10/3*x
         */
        double[] y = { 0.3333333333, 0.633333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 5 nodes", 5,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);

        assertTrue("Tree root node should have split on first variable in node 0, and 2 and not split anywhere else", Arrays.equals(new int[]{-1, 0, -1, 0, 0},tree.var));

        assertTrue("Node sizes of nodes should be 3 for root, 2 for node #2 and 1 for everything else", Arrays.equals(new int[] { 3,1,2,1,1}, tree.nodesize));

        assertEquals("Tree root should be a categorical split and have a cut of 0", 0, tree.cut[0], 0.000001);

        assertEquals("Node #2 should be a categorical split and have index one", 1, tree.cut[2], 0.000001);

        assertTrue("All nodevar should be 0", Arrays.equals(new double[] { 0,0,0,0,0}, tree.nodevar));

        assertEquals("Tree root node should have first child as node 1", tree.children[0][0], 1);
        assertEquals("Tree root node should have second child as node 2", tree.children[0][1], 2);

        assertEquals("Tree first child node prediction should be ~0.3333333333", 0.3333333333, tree.nodepred[1], 0.000001);


        assertEquals("Node size of second child should be 2", 2, tree.nodesize[2]);

        assertEquals("Tree node #2 should have first child as node 3", tree.children[2][0], 3);
        assertEquals("Tree node #2 should have second child as node 4", tree.children[2][1], 4);

        assertEquals("Tree node #3 prediction should be ~0.666666666", 0.6333333333, tree.nodepred[3], 0.000001);
        assertEquals("Tree node #4 prediction should be ~0.666666666", 0.6666666666, tree.nodepred[4], 0.000001);

    }


    @Test

    /**
     * The expected tree here is a single node since split min is too high
     */
    public void testRegTreeFitTwoDistinctSamplesSameFeatureSplitMinThree()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinThreeSingleThetaSingleFeature();

        double[][] allTheta = { { 0.10}, {0.20} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0}};
        double[] y = { 0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 1 nodes", 1,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);


        assertEquals("Node size of root should be 2", 2, tree.nodesize[0]);


        assertEquals("Tree root node prediction should be ~0.5", 0.5, tree.nodepred[0], 0.000001);

        /**
         * Sample Variance is ((2/3-1/2)^2+(1/3-1/2)^2)*1/(2-1)
         */
        assertEquals("Tree root node variance should be ~0", 0.0555555555, tree.nodevar[0], 0.000001 );

    }

    @Test
    /**
     * The expected tree here should have depth three as follows:
     *
     *                       0
     *                      / \
     *                     1   2
     *                        / \
     *                       3  4
     *
     * First split should be the always active variable in the first domain.
     * Next split should break into even and odds.
     *
     *
     *
     * Even values should go one way and odd to another.
     */
    public void testRegTreeFitTenDistinctSamplesSameParameterActive()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsDoubleThetaFirstParamBinaryConditionalSingleFeature();

        double[][] allTheta = { {2, 1}, {2,2}, {2,3}, {2,4} ,{2,5},{2,6},{2,7},{2,8},{2,9},{2,10} , {1,1} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0} , {2,0},{3,0},{4,0},{5,0},{6,0},{7,0},{8,0},{9,0}, {10,0}};
        double[] y = { 0.5, 0.75,0.5, 0.75,0.5, 0.75,0.5, 0.75,0.5, 0.75, 0.1};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 5 nodes", 5,tree.numNodes);
        assertEquals("Tree should have 3 predictors", 3, tree.npred);


        assertEquals("Tree root node should have split on second variable",-1, tree.var[0]);
        //I only tested one of the valid splits not sure why it always works but {1,0} would also be valid here.
        assertTrue("Even nodes should go one way and odd nodes the other way",Arrays.equals(tree.catsplit[0] , new int[]{0,1}));


        //Not sure what cut has here
        // assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.20) );
        assertTrue("Tree should match the structure specified in the file", Arrays.equals(tree.nodesize,new int[]{11,1,10,5,5}));

        assertEquals("Tree node #1 node prediction should be ~0.1", 0.1, tree.nodepred[1], 0.000001);

        assertEquals("Tree node #3 node prediction should be ~0.5", 0.5, tree.nodepred[3], 0.000001);

        assertEquals("Tree node #4 prediction should be ~0.75", 0.75, tree.nodepred[4], 0.000001);

    }


    @Test
    /**
     * The expected tree here should have depth three as follows:
     *
     *                       0
     *                      / \
     *                     1   2
     *
     *
     * First split should be the always active variable in the first domain.
     * Next split should break into even and odds.
     *
     *
     *
     * Even values should go one way and odd to another.
     */
    public void testRegTreeFitTenDistinctSamplesSameParameterInActive()
    {

        RegtreeBuildParams rbp = getRegtreeBuildParamsDoubleThetaFirstParamBinaryConditionalSingleFeature();

        double[][] allTheta = { {1, 1}, {1,2}, {1,3}, {1,4} ,{1,5},{1,6},{1,7},{1,8},{1,9},{1,10} , {2,1} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0} , {2,0},{3,0},{4,0},{5,0},{6,0},{7,0},{8,0},{9,0}, {10,0}};
        double[] y = {0.5, 0.75, 0.5, 0.75, 0.5, 0.75, 0.5, 0.75, 0.5, 0.75, 0.1};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 3 nodes", 3,tree.numNodes);
        assertEquals("Tree should have 3 predictors", 3, tree.npred);


        assertEquals("Tree root node should have split on second variable",-1, tree.var[0]);
        //I only tested one of the valid splits not sure why it always works but {0,1} would also be valid here.
        assertTrue("Even nodes should go one way and odd nodes the other way",Arrays.equals(tree.catsplit[0] , new int[]{1,0}));


        //Not sure what cut has here
        // assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.20) );
        assertTrue("Tree should match the structure specified in the file", Arrays.equals(tree.nodesize,new int[]{11,1,10}));

        assertEquals("Tree node #1 node prediction should be ~0.1", 0.1, tree.nodepred[1], 0.000001);

        assertEquals("Tree node #3 node prediction should be ~0.625", 0.625, tree.nodepred[2], 0.000001);
        assertEquals("Tree node #3 node predicted variance should 1/(10-9)*[5*(0.5-0.625)^2+5*(0.75-0.625)^2] should be ~0.0173611111", 0.0173611111, tree.nodevar[2], 0.000001);



    }
}