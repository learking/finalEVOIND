package test.beast.evolution.substmodel;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.TwoStructs;
import junit.framework.TestCase;

public class TwoStructsRateMatrixTest extends TestCase {

	@Test
	public void testSetupRateMatrix() throws Exception {
		RealParameter f = new RealParameter(new Double[]{0.2, 0.1, 0.3, 0.4});
		Frequencies nucleoFrequencies = new Frequencies();
		nucleoFrequencies.initByName("frequencies", f, "estimate", false);
		
		RealParameter codonProb_struct1 = new RealParameter(new Double[]{
				0.031350262,0.016678685,0.024086503,0.013927975,0.025667266,0.019513173,0.001036515,
				0.019661785,0.004234123,0.021709737,0.023116887,0.008893605,0.023613895,0.016321916,
				0.021248437,0.034586834,0.013469791,0.013974397,0.019040362,0.023127154,0.012276834,
				0.016818336,0.013005605,0.016724293,0.031102112,0.025201710,0.002783288,0.028634640,
				0.011845951,0.035306332,0.007706424,0.016484928,0.001266775,0.015199627,0.024995155,
				0.027195620,0.022472747,0.008187028,0.029989578,0.003142386,0.027191536,0.010992749,
				0.010295569,0.005584570,0.029299030,0.006562422,0.004537246,0.009528992,0.002336504,
				0.002458360,0.004808568,0.014593801,0.006875164,0.009996637,0.019877895,0.012218437,
				0.026721705,0.015533992,0.020767555,0.021984981,0.012235612});
		
		RealParameter codonProb_struct2 = new RealParameter(new Double[]{
				0.027881210,0.022293149,0.026883246,0.023120719,0.016919530,0.019623787,
				0.006413637,0.014555090,0.011296171,0.018025663,0.010311862,0.012645179,
				0.007781835,0.018769459,0.018094488,0.017920131,0.023053429,0.015897817,
				0.022228268,0.016487978,0.015231865,0.015117623,0.003877001,0.015223928,
				0.007608210,0.009564091,0.009518513,0.004670664,0.006888680,0.017072271,
				0.032558763,0.013889674,0.032962994,0.028136097,0.031783136,0.029180570,
				0.016575473,0.026630714,0.005404942,0.021145367,0.021327834,0.022109745,
				0.015010434,0.012191939,0.008375515,0.013699401,0.027602561,0.012560591,
				0.018228453,0.018905132,0.011479767,0.015714268,0.003529814,0.014752038,
				0.009355456,0.013942310,0.009702750,0.008341407,0.018867836,0.011491272,
				0.019568250});
		
		TwoStructs twoStructs = new TwoStructs();       
		twoStructs.initByName("percent", "0.3", "kappa", "2", "omega", "1.3", "nucleoFrequencies", nucleoFrequencies, "codonProbStruct1", codonProb_struct1, "codonProbStruct2", codonProb_struct2);  
		
		twoStructs.prepareMatricesForTest();
		
		//write out codonProbResults
		System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/TwoStructs_rootDist.txt")));  
		System.out.println(Arrays.toString(twoStructs.getDiagMatrix()));
		
		//and rateM
		System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/TwoStructs_rateM.txt")));  
		System.out.println(Arrays.deepToString(twoStructs.getRateMatrix()));
		
		/*Python verify
		 * 
		 * import numpy as np
from numpy.testing import (run_module_suite, TestCase,
        decorators, assert_, assert_allclose)
import re

def read2dArray(inputFile):
    f = open(inputFile, "r")
    inputData = f.readlines()
    f.close()
    resultArray = []
    if len(inputData) != 1:
        print("there should be only 1 line for 2dArray")
    else:
        #get 2d arrays using regular expression                                                        
        allRows = re.split('\],\s\[', inputData[0])
        for r in allRows:
            #extract numbers in each row                                                               
            tmpNumbers = [float(x) for x in r.replace(',','').replace('[','').replace(']','').split()]
            resultArray.append(tmpNumbers)
    return resultArray

rootDist = [0.0289219256, 0.0206088098, 0.0260442231, 0.0203628958, 0.0195438508, 0.0195906028, 0.004800500399999999, 0.0160870985, 0.009177556600000001, 0.0191308852, 0.014153369499999999, 0.0115197068, 0.012531453, 0.018035196099999998, 0.019040672699999997, 0.022920141899999996, 0.0201783376, 0.015320791, 0.021271896199999998, 0.0184797308, 0.0143453557, 0.015627836899999997, 0.006615582199999999, 0.0156740375, 0.014656380600000001, 0.014255376699999999, 0.0074979454999999995, 0.011859856799999999, 0.008375861299999999, 0.022542489300000003, 0.025103061299999996, 0.014668250199999998, 0.023454128299999998, 0.024255155999999996, 0.029746741699999996, 0.028585084999999996, 0.018344655199999997, 0.0210976082, 0.012780332799999999, 0.0157444727, 0.0230869446, 0.018774646199999998, 0.0135959745, 0.0102097283, 0.0146525695, 0.011558307299999999, 0.0206829665, 0.0116511113, 0.013460868299999998, 0.0139711004, 0.009478407299999998, 0.015378127899999998, 0.004533419, 0.013325417699999998, 0.0125121877, 0.0134251481, 0.014808436499999997, 0.010499182499999999, 0.0194377517, 0.0146393847, 0.017368458599999997]

root_prior_distn1d = np.array(rootDist, float)

# Define the rate matrix
rateM = np.array(read2dArray("/home/kuangyu/Desktop/TwoStructs_rateM.txt"),float)

# Check that the stationary distribution is correct.
equilibrium_rates = np.dot(root_prior_distn1d, rateM)

# check that the expected rate is 1 (or -1 if we use diagonal)
sum(np.diag(rateM) * root_prior_distn1d)


		 * 
		 */
	}
	
}
