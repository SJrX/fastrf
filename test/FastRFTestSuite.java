import de.unifreiburg.cs.junit.DataTester;
import de.unifreiburg.cs.junit.RFTester;
import fastrf.fastrf;
import junit.framework.TestSuite;
import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 * Test Suite for Fast RF Code
 *
 * @author Steve Ramage <seramage@cs.ubc.ca>
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({
        DataTester.class,
        RFTester.class,
        fastrf.class
})

public class FastRFTestSuite {

}
