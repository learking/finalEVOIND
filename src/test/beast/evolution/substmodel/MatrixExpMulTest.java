package test.beast.evolution.substmodel;

import java.util.Arrays;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import junit.framework.TestCase;

public class MatrixExpMulTest extends TestCase {
	
	public void testJBLAS3() throws Exception{
		double distance = 0.3;
		DoubleMatrix doubleRateMatrix;
		double[][] rateMatrix = new double[][]{{1,2,3,4},{4,3,2,1},{1,1,2,2},{2,2,1,1}}; 
		doubleRateMatrix = new DoubleMatrix(rateMatrix);
		
		double[][] tmpM1 = MatrixFunctions.expm(doubleRateMatrix.mul(distance)).toArray2();
		double[][] tmpM2 = MatrixFunctions.expm(doubleRateMatrix.mul(distance)).toArray2();

		System.out.println(Arrays.deepToString(tmpM1));
		System.out.println(Arrays.deepToString(tmpM2));		
		
	}
}
