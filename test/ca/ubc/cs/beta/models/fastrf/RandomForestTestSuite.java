package ca.ubc.cs.beta.models.fastrf;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 * Created by Steve Ramage <seramage@cs.ubc.ca> on 11/14/15
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({
        RegtreeFitTestContinuousDomain.class,
        RegtreeFitTestCategoricalDomain.class
})
public class RandomForestTestSuite {
}
