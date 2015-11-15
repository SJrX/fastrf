package ca.ubc.cs.beta.models.fastrf;

import de.unifreiburg.cs.junit.DataTester;
import de.unifreiburg.cs.junit.RFTester;
import fastrf.fastrf;
import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 * Created by Steve Ramage <seramage@cs.ubc.ca> on 11/14/15
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({
        RegtreeFitTestContinuousDomain.class,
        RegtreeFitTestCategoricalDomain.class,
        DataTester.class,
        RFTester.class,
        fastrf.class
})
public class RandomForestTestSuite {
}