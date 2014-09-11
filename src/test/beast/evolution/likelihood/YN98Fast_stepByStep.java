package test.beast.evolution.likelihood;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98Fast;
import junit.framework.TestCase;

public class YN98Fast_stepByStep extends TestCase {

		@Test
		public void testRateMatrix_allEqual() throws Exception {			
	        RealParameter f = new RealParameter(new Double[]{0.25,0.25,0.25,0.25});
	        Frequencies nucleoFrequencies = new Frequencies();
	        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
	        YN98Fast yn98 = new YN98Fast();
	        yn98.initByName("kappa", "1", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);
	        
	        yn98.prepareMatricesForTest();

	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_symmM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_symmM.txt")));
	        System.out.println(Arrays.deepToString(yn98.getSymmMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_diagM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_diagM.txt")));
	        System.out.println(Arrays.toString(yn98.getDiagMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_rateM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_rateM.txt")));
	        System.out.println(Arrays.deepToString(yn98.getRateMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_freq_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_freq.txt")));
	        System.out.println(Arrays.toString(yn98.getFrequencies()));
		}
		
		@Test
		public void testRateMatrix_allDiff() throws Exception {			
	        RealParameter f = new RealParameter(new Double[]{0.308,0.228,0.23,0.234});
	        Frequencies nucleoFrequencies = new Frequencies();
	        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
	        YN98Fast yn98 = new YN98Fast();
	        yn98.initByName("kappa", "1", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);
	        
	        yn98.prepareMatricesForTest();

	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_symmM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_symmM.txt")));
	        System.out.println(Arrays.deepToString(yn98.getSymmMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_diagM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_diagM.txt")));
	        System.out.println(Arrays.toString(yn98.getDiagMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_rateM_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_rateM.txt")));
	        System.out.println(Arrays.deepToString(yn98.getRateMatrix()));
	        
	        System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/YN98_freq_updated.txt")));
	        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98_freq.txt")));
	        System.out.println(Arrays.toString(yn98.getFrequencies()));
		}
		
}
