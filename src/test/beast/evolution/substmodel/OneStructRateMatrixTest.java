package test.beast.evolution.substmodel;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.OneStruct;
import junit.framework.TestCase;

public class OneStructRateMatrixTest extends TestCase {

	public void testSetupRateMatrix() throws Exception {
		
		RealParameter f = new RealParameter(new Double[]{0.2, 0.1, 0.3, 0.4});
		Frequencies nucleoFrequencies = new Frequencies();
		nucleoFrequencies.initByName("frequencies", f, "estimate", false);

		RealParameter codonProb = new RealParameter(new Double[]{
				0.031350262,0.016678685,0.024086503,0.013927975,0.025667266,0.019513173,0.001036515,
				0.019661785,0.004234123,0.021709737,0.023116887,0.008893605,0.023613895,0.016321916,
				0.021248437,0.034586834,0.013469791,0.013974397,0.019040362,0.023127154,0.012276834,
				0.016818336,0.013005605,0.016724293,0.031102112,0.025201710,0.002783288,0.028634640,
				0.011845951,0.035306332,0.007706424,0.016484928,0.001266775,0.015199627,0.024995155,
				0.027195620,0.022472747,0.008187028,0.029989578,0.003142386,0.027191536,0.010992749,
				0.010295569,0.005584570,0.029299030,0.006562422,0.004537246,0.009528992,0.002336504,
				0.002458360,0.004808568,0.014593801,0.006875164,0.009996637,0.019877895,0.012218437,
				0.026721705,0.015533992,0.020767555,0.021984981,0.012235612});

		OneStruct oneStruct = new OneStruct();       
		oneStruct.initByName("kappa", "2", "omega", "1.3", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb);       

		oneStruct.prepareMatricesForTest();
		double[][] currRateM = oneStruct.getRateMatrix();
		
		double sumValue = 0;
		for (int i=0; i < currRateM.length; i++){
			for (int j=0; j < currRateM.length; j++){
				sumValue += currRateM[i][j];
			}
		}
		System.out.println("sumValue:" + sumValue);
		
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/curr_OneStruct_rateM.txt")));        
		//System.out.println(Arrays.deepToString(currRateM));
	}
}