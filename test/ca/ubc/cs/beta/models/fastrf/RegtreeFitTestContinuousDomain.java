package ca.ubc.cs.beta.models.fastrf;

import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;


/**
 * Created by Steve Ramage <seramage@cs.ubc.ca> on 11/14/15
 */
public class RegtreeFitTestContinuousDomain {

    public void testEmptyDataIndexes()
    {
        fail();
    }


    public void testResponseDataIndexesMismatchThrowsException()
    {
        fail();
    }



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


        RegtreeBuildParams rbp = new RegtreeBuildParams();

        rbp.splitMin = splitMin;
        rbp.seed = new Random().nextInt();

        System.err.println("Seed for run was: " + rbp.seed);
        rbp.ratioFeatures = ratioFeatures;

        rbp.catDomainSizes = new int[]{0, 0};
        rbp.condParents = null;
        rbp.logModel = 0;
        rbp.storeResponses = false;

        return rbp;

    }


    @Test
    /**
     * The expected tree here should be a single constant point.
     */
    public void testRegTreeFitSingleContinuousPoint()
    {
        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 0.25} };

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

        double[][] allTheta = { { 0.25} };

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

        double[][] allTheta = { { 0.10}, {0.20} };

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0}};
        double[] y = { 0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 3 nodes", 3,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node should have split on first variable",1, tree.var[0]);

        assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.20) );
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
     */
    public void testRegTreeFitTwoDistinctSamplesSameTheta()
    {
        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 0.25}};

        double[][] allX = { { 0.8} , {0.9} };
        int[][] dataIndxs = { { 0, 0}, {0,1}};
        double[] y = { 0.3333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 3 nodes", 3,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertEquals("Tree root node should have split on second variable (feature)",2, tree.var[0]);

        assertTrue("Tree root node should have a cut between 0.1 and 0.2", (tree.cut[0] > 0.80) && (tree.cut[0] < 0.90) );
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
    public void testRegTreeFitThreeSampleOneContinuousSameFeature()
    {
        RegtreeBuildParams rbp = getRegtreeBuildParamsSplitMinOneSingleThetaSingleFeature();

        double[][] allTheta = { { 0.10}, {0.19}, {0.20}};

        double[][] allX = { { 0.75 } };
        int[][] dataIndxs = { { 0, 0}, {1,0}, {2,0}};

        /**
         * We are essentially modelling y = 10/3*x
         */
        double[] y = { 0.3333333333, 0.633333333333, 0.6666666666};

        Regtree tree = RegtreeFit.fit(allTheta,allX, dataIndxs,y,rbp);

        assertEquals("Tree should have 5 nodes", 5,tree.numNodes);
        assertEquals("Tree should have 2 predictors", 2, tree.npred);
        assertTrue("Tree root node should have split on first variable in node 0, and 2 and not split anywhere else", Arrays.equals(new int[]{1, 0, 1, 0, 0},tree.var));

        assertTrue("Tree root node should have a cut between 0.1 and 0.19", (tree.cut[0] > 0.10) && (tree.cut[0] < 0.19) );
        assertTrue("Tree node 2 should have a cut between 0.19 and 0.20", (tree.cut[2] > 0.19) && (tree.cut[2] < 0.20) );

        assertTrue("Node sizes of nodes should be 3 for root, 2 for node #2 and 1 for everything else", Arrays.equals(new int[] { 3,1,2,1,1}, tree.nodesize));

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

}
