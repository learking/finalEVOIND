package beast.evolution.substitutionmodel;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Node;

public class OneStruct extends GeneralSubstitutionModel {
    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "kappa parameter in YN98 model", Validate.REQUIRED);
    //frequencies for pi_A, pi_C, pi_G and pi_T
    public Input<Frequencies> nucleoFreqInput =
            new Input<Frequencies>("nucleoFrequencies", "substitution model equilibrium state frequencies", Validate.REQUIRED);
    public Input<RealParameter> codonProbInput = new Input<RealParameter>("codonProb", "probabilities for each codon in OneStruct model", Validate.REQUIRED);    
    double[] codonProb;
    
    Frequencies nucleoFrequencies;
    double[][] symmMatrix;
    double[] diagMatrix;
    
    @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
    	double distance = (fStartTime - fEndTime) * fRate;

        int i, j, k;
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
            	System.out.println("rateMatrix was updated");
                setupRelativeRates();
                setupRateMatrix();
                
/*        		double sumValue = 0;
        		for (int rowNr=0; rowNr < rateMatrix.length; rowNr++){
        			for (int colNr=0; colNr < rateMatrix.length; colNr++){
        				sumValue += rateMatrix[rowNr][colNr];
        			}
        		}
        		System.out.println("sumValue:" + sumValue);*/
        		
                double[][] copyRateMatrix = new double[nrOfStates][nrOfStates];
        		for (int rowNr=0; rowNr < rateMatrix.length; rowNr++){
        			for (int colNr=0; colNr < rateMatrix.length; colNr++){
        				copyRateMatrix[rowNr][colNr] = rateMatrix[rowNr][colNr];
        			}
        		}
                eigenDecomposition = eigenSystem.decomposeMatrix(copyRateMatrix);
                updateMatrix = false;
            }
        }
        //System.out.println(node.getNr() + " " + Arrays.deepToString(rateMatrix));
/*		double sumValue = 0;
		for (int rowNr=0; rowNr < rateMatrix.length; rowNr++){
			for (int colNr=0; colNr < rateMatrix.length; colNr++){
				sumValue += rateMatrix[rowNr][colNr];
			}
		}
		System.out.println(node.getNr() + " " + "sumValue:" + sumValue);*/
        
        // is the following really necessary?
        // implemented a pool of iexp matrices to support multiple threads
        // without creating a new matrix each call. - AJD
        // a quick timing experiment shows no difference - RRB
        double[] iexp = new double[nrOfStates * nrOfStates];
        // Eigen vectors
        double[] Evec = eigenDecomposition.getEigenVectors();
        // inverse Eigen vectors
        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
        // Eigen values
        double[] Eval = eigenDecomposition.getEigenValues();
        for (i = 0; i < nrOfStates; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < nrOfStates; j++) {
                iexp[i * nrOfStates + j] = Ievc[i * nrOfStates + j] * temp;
            }
        }

        int u = 0;
        for (i = 0; i < nrOfStates; i++) {
            for (j = 0; j < nrOfStates; j++) {
                temp = 0.0;
                for (k = 0; k < nrOfStates; k++) {
                    temp += Evec[i * nrOfStates + k] * iexp[k * nrOfStates + j];
                }

                matrix[u] = Math.abs(temp);
                u++;
            }
        }
        
        //System.out.println(node.getNr() + " probabilities:" + Arrays.toString(matrix));
    } // getTransitionProbabilities
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public OneStruct() {
    	frequenciesInput.setRule(Validate.OPTIONAL);
        ratesInput.setRule(Validate.OPTIONAL);
        try {
        	frequenciesInput.setValue(null, this);
        	ratesInput.setValue(null, this);
        } catch (Exception e) {
        	e.printStackTrace();
			// TODO: handle exception
		}
    }
    
    @Override
    public void initAndValidate() throws Exception {
        if (frequenciesInput.get() != null) {
            throw new Exception("the frequencies attribute should not be used. Use the nucleoFrequencies instead");
        }
        if (ratesInput.get() != null) {
            throw new Exception("the rates attribute should not be used. Use the individual rates rateAC, rateCG, etc, instead.");
        }

        // set codonProb here (codonProb will not be changed after initialization)
        codonProb = new double[codonProbInput.get().getDimension()];
        for (int i = 0; i < codonProb.length; i++) {
        	codonProb[i] = codonProbInput.get().getValue(i);
        }
        // sanity check
        double totalCodonP = 0.0;
        for(int i=0;i<codonProb.length;i++)
        {
        	totalCodonP = totalCodonP + codonProb[i];
        }
        if (Math.abs(totalCodonP - 1.0) > 1e-4) {
        	System.out.format("totalCodonP value is:" + "%f%n", totalCodonP); 
            throw new Exception("Codon probabilities should sum up to 1.");
        }
        
        nucleoFrequencies = nucleoFreqInput.get();
        updateMatrix = true;
        nrOfStates = 61;

        eigenSystem = createEigenSystem();
        
        symmMatrix = new double[nrOfStates][nrOfStates];
        diagMatrix = new double[nrOfStates];
        rateMatrix = new double[nrOfStates][nrOfStates];

    }
    
    void setupDiagMatrix(){
    	diagMatrix = getCodonProb();
    }
    
    double[] getCodonProb(){
    	Double[] codonProb = codonProbInput.get().getValues();
    	double[] codonProbResult = new double[codonProb.length];
    	for (int i = 0; i < codonProb.length; i++) {
    		codonProbResult[i] = codonProb[i];
    	}
    	return codonProbResult;
    }
    /**
     * sets up rate matrix *
     */
    void setupSymmMatrix(){
    	
    	final double k = kappaInput.get().getValue();
    	final double omega = omegaInput.get().getValue();
    	
    	//0:A 1:C 2:G 3:T
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	//P(C|GE)P(AA|GE, RSA)
    	double[] codonP = getCodonProb();
    	//compute P(c|pi): 61 values
    	double[] p0 = new double[nrOfStates];
    	//ratio between codonP and p0
    	double[] pRatio = new double[nrOfStates];
    	
    	p0[0] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[0]; // AAA
    	p0[1] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[1]; // AAC
    	p0[2] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[2]; // AAG
    	p0[3] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[3]; // AAT
    	p0[4] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[0]; // ACA
    	p0[5] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[1]; // ACC
    	p0[6] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[2]; // ACG
    	p0[7] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[3]; // ACT
    	p0[8] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[0]; // AGA
    	p0[9] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[1]; // AGC
    	p0[10] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[2]; // AGG
    	p0[11] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[3]; // AGT
    	p0[12] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[0]; // ATA
    	p0[13] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[1]; // ATC
    	p0[14] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[2]; // ATG
    	p0[15] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[3]; // ATT
    	p0[16] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[0]; // CAA
    	p0[17] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[1]; // CAC
    	p0[18] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[2]; // CAG
    	p0[19] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[3]; // CAT
    	p0[20] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[0]; // CCA
    	p0[21] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[1]; // CCC
    	p0[22] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[2]; // CCG
    	p0[23] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[3]; // CCT
    	p0[24] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[0]; // CGA
    	p0[25] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[1]; // CGC
    	p0[26] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[2]; // CGG
    	p0[27] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[3]; // CGT
    	p0[28] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[0]; // CTA
    	p0[29] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[1]; // CTC
    	p0[30] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[2]; // CTG
    	p0[31] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[3]; // CTT
    	p0[32] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[0]; // GAA
    	p0[33] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[1]; // GAC
    	p0[34] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[2]; // GAG
    	p0[35] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[3]; // GAT
    	p0[36] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[0]; // GCA
    	p0[37] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[1]; // GCC
    	p0[38] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[2]; // GCG
    	p0[39] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[3]; // GCT
    	p0[40] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[0]; // GGA
    	p0[41] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[1]; // GGC
    	p0[42] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[2]; // GGG
    	p0[43] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[3]; // GGT
    	p0[44] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[0]; // GTA
    	p0[45] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[1]; // GTC
    	p0[46] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[2]; // GTG
    	p0[47] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[3]; // GTT
    	p0[48] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[1]; // TAC
    	p0[49] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[3]; // TAT
    	p0[50] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[0]; // TCA
    	p0[51] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[1]; // TCC
    	p0[52] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[2]; // TCG
    	p0[53] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[3]; // TCT
    	p0[54] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[1]; // TGC
    	p0[55] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[2]; // TGG
    	p0[56] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[3]; // TGT
    	p0[57] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[0]; // TTA
    	p0[58] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[1]; // TTC
    	p0[59] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[2]; // TTG
    	p0[60] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[3]; // TTT
    	
    	//get array that contains codonP and p0 ratio before proceeding to construct symm matrix (save computation)
    	for (int i = 0; i < pRatio.length; i++) {
    		pRatio[i] = codonP[i] / p0[i];
    	}
    	
    	//System.out.print("codonP:");
    	//System.out.println(Arrays.toString(codonP));
    	//System.out.print("p0:");
    	//System.out.println(Arrays.toString(p0));
    	//System.out.print("pRatio:");
    	//System.out.println(Arrays.toString(pRatio));
    	
    	//construct symmetric codon substitution matrix (61X61):
    	//0: diagonal positions, codons differ at more than one position
    	symmMatrix[0][1] = Math.log(pRatio[1] / pRatio[0]) / (pRatio[1] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->AAC: IS NonSynonymous;
    	symmMatrix[0][2] = Math.log(pRatio[2] / pRatio[0]) / (pRatio[2] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * k ; 	//AAA->AAG: IS Transition;
    	symmMatrix[0][3] = Math.log(pRatio[3] / pRatio[0]) / (pRatio[3] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->AAT: IS NonSynonymous;
    	symmMatrix[0][4] = Math.log(pRatio[4] / pRatio[0]) / (pRatio[4] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->ACA: IS NonSynonymous;
    	symmMatrix[0][8] = Math.log(pRatio[8] / pRatio[0]) / (pRatio[8] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//AAA->AGA: IS Transition;IS NonSynonymous;
    	symmMatrix[0][12] = Math.log(pRatio[12] / pRatio[0]) / (pRatio[12] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->ATA: IS NonSynonymous;
    	symmMatrix[0][16] = Math.log(pRatio[16] / pRatio[0]) / (pRatio[16] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->CAA: IS NonSynonymous;
    	symmMatrix[0][32] = Math.log(pRatio[32] / pRatio[0]) / (pRatio[32] - pRatio[0]) / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//AAA->GAA: IS Transition;IS NonSynonymous;
    	symmMatrix[1][2] = Math.log(pRatio[2] / pRatio[1]) / (pRatio[2] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAC->AAG: IS NonSynonymous;
    	symmMatrix[1][3] = Math.log(pRatio[3] / pRatio[1]) / (pRatio[3] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[0]) * k ; 	//AAC->AAT: IS Transition;
    	symmMatrix[1][5] = Math.log(pRatio[5] / pRatio[1]) / (pRatio[5] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->ACC: IS NonSynonymous;
    	symmMatrix[1][9] = Math.log(pRatio[9] / pRatio[1]) / (pRatio[9] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//AAC->AGC: IS Transition;IS NonSynonymous;
    	symmMatrix[1][13] = Math.log(pRatio[13] / pRatio[1]) / (pRatio[13] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->ATC: IS NonSynonymous;
    	symmMatrix[1][17] = Math.log(pRatio[17] / pRatio[1]) / (pRatio[17] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->CAC: IS NonSynonymous;
    	symmMatrix[1][33] = Math.log(pRatio[33] / pRatio[1]) / (pRatio[33] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//AAC->GAC: IS Transition;IS NonSynonymous;
    	symmMatrix[1][48] = Math.log(pRatio[48] / pRatio[1]) / (pRatio[48] - pRatio[1]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->TAC: IS NonSynonymous;
    	symmMatrix[2][3] = Math.log(pRatio[3] / pRatio[2]) / (pRatio[3] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAG->AAT: IS NonSynonymous;
    	symmMatrix[2][6] = Math.log(pRatio[6] / pRatio[2]) / (pRatio[6] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->ACG: IS NonSynonymous;
    	symmMatrix[2][10] = Math.log(pRatio[10] / pRatio[2]) / (pRatio[10] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//AAG->AGG: IS Transition;IS NonSynonymous;
    	symmMatrix[2][14] = Math.log(pRatio[14] / pRatio[2]) / (pRatio[14] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->ATG: IS NonSynonymous;
    	symmMatrix[2][18] = Math.log(pRatio[18] / pRatio[2]) / (pRatio[18] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->CAG: IS NonSynonymous;
    	symmMatrix[2][34] = Math.log(pRatio[34] / pRatio[2]) / (pRatio[34] - pRatio[2]) / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//AAG->GAG: IS Transition;IS NonSynonymous;
    	symmMatrix[3][7] = Math.log(pRatio[7] / pRatio[3]) / (pRatio[7] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->ACT: IS NonSynonymous;
    	symmMatrix[3][11] = Math.log(pRatio[11] / pRatio[3]) / (pRatio[11] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//AAT->AGT: IS Transition;IS NonSynonymous;
    	symmMatrix[3][15] = Math.log(pRatio[15] / pRatio[3]) / (pRatio[15] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->ATT: IS NonSynonymous;
    	symmMatrix[3][19] = Math.log(pRatio[19] / pRatio[3]) / (pRatio[19] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->CAT: IS NonSynonymous;
    	symmMatrix[3][35] = Math.log(pRatio[35] / pRatio[3]) / (pRatio[35] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//AAT->GAT: IS Transition;IS NonSynonymous;
    	symmMatrix[3][49] = Math.log(pRatio[49] / pRatio[3]) / (pRatio[49] - pRatio[3]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->TAT: IS NonSynonymous;
    	symmMatrix[4][5] = Math.log(pRatio[5] / pRatio[4]) / (pRatio[5] - pRatio[4]) / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACA->ACC: 
    	symmMatrix[4][6] = Math.log(pRatio[6] / pRatio[4]) / (pRatio[6] - pRatio[4]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k ; 	//ACA->ACG: IS Transition;
    	symmMatrix[4][7] = Math.log(pRatio[7] / pRatio[4]) / (pRatio[7] - pRatio[4]) / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACA->ACT: 
    	symmMatrix[4][8] = Math.log(pRatio[8] / pRatio[4]) / (pRatio[8] - pRatio[4]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//ACA->AGA: IS NonSynonymous;
    	symmMatrix[4][12] = Math.log(pRatio[12] / pRatio[4]) / (pRatio[12] - pRatio[4]) / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//ACA->ATA: IS Transition;IS NonSynonymous;
    	symmMatrix[4][20] = Math.log(pRatio[20] / pRatio[4]) / (pRatio[20] - pRatio[4]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//ACA->CCA: IS NonSynonymous;
    	symmMatrix[4][36] = Math.log(pRatio[36] / pRatio[4]) / (pRatio[36] - pRatio[4]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//ACA->GCA: IS Transition;IS NonSynonymous;
    	symmMatrix[4][50] = Math.log(pRatio[50] / pRatio[4]) / (pRatio[50] - pRatio[4]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//ACA->TCA: IS NonSynonymous;
    	symmMatrix[5][6] = Math.log(pRatio[6] / pRatio[5]) / (pRatio[6] - pRatio[5]) / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACC->ACG: 
    	symmMatrix[5][7] = Math.log(pRatio[7] / pRatio[5]) / (pRatio[7] - pRatio[5]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k ; 	//ACC->ACT: IS Transition;
    	symmMatrix[5][9] = Math.log(pRatio[9] / pRatio[5]) / (pRatio[9] - pRatio[5]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//ACC->AGC: IS NonSynonymous;
    	symmMatrix[5][13] = Math.log(pRatio[13] / pRatio[5]) / (pRatio[13] - pRatio[5]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//ACC->ATC: IS Transition;IS NonSynonymous;
    	symmMatrix[5][21] = Math.log(pRatio[21] / pRatio[5]) / (pRatio[21] - pRatio[5]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//ACC->CCC: IS NonSynonymous;
    	symmMatrix[5][37] = Math.log(pRatio[37] / pRatio[5]) / (pRatio[37] - pRatio[5]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//ACC->GCC: IS Transition;IS NonSynonymous;
    	symmMatrix[5][51] = Math.log(pRatio[51] / pRatio[5]) / (pRatio[51] - pRatio[5]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//ACC->TCC: IS NonSynonymous;
    	symmMatrix[6][7] = Math.log(pRatio[7] / pRatio[6]) / (pRatio[7] - pRatio[6]) / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACG->ACT: 
    	symmMatrix[6][10] = Math.log(pRatio[10] / pRatio[6]) / (pRatio[10] - pRatio[6]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//ACG->AGG: IS NonSynonymous;
    	symmMatrix[6][14] = Math.log(pRatio[14] / pRatio[6]) / (pRatio[14] - pRatio[6]) / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//ACG->ATG: IS Transition;IS NonSynonymous;
    	symmMatrix[6][22] = Math.log(pRatio[22] / pRatio[6]) / (pRatio[22] - pRatio[6]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//ACG->CCG: IS NonSynonymous;
    	symmMatrix[6][38] = Math.log(pRatio[38] / pRatio[6]) / (pRatio[38] - pRatio[6]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//ACG->GCG: IS Transition;IS NonSynonymous;
    	symmMatrix[6][52] = Math.log(pRatio[52] / pRatio[6]) / (pRatio[52] - pRatio[6]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//ACG->TCG: IS NonSynonymous;
    	symmMatrix[7][11] = Math.log(pRatio[11] / pRatio[7]) / (pRatio[11] - pRatio[7]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ACT->AGT: IS NonSynonymous;
    	symmMatrix[7][15] = Math.log(pRatio[15] / pRatio[7]) / (pRatio[15] - pRatio[7]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//ACT->ATT: IS Transition;IS NonSynonymous;
    	symmMatrix[7][23] = Math.log(pRatio[23] / pRatio[7]) / (pRatio[23] - pRatio[7]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//ACT->CCT: IS NonSynonymous;
    	symmMatrix[7][39] = Math.log(pRatio[39] / pRatio[7]) / (pRatio[39] - pRatio[7]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//ACT->GCT: IS Transition;IS NonSynonymous;
    	symmMatrix[7][53] = Math.log(pRatio[53] / pRatio[7]) / (pRatio[53] - pRatio[7]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//ACT->TCT: IS NonSynonymous;
    	symmMatrix[8][9] = Math.log(pRatio[9] / pRatio[8]) / (pRatio[9] - pRatio[8]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGA->AGC: IS NonSynonymous;
    	symmMatrix[8][10] = Math.log(pRatio[10] / pRatio[8]) / (pRatio[10] - pRatio[8]) / (nucleoFreqs[0] * nucleoFreqs[2]) * k ; 	//AGA->AGG: IS Transition;
    	symmMatrix[8][11] = Math.log(pRatio[11] / pRatio[8]) / (pRatio[11] - pRatio[8]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGA->AGT: IS NonSynonymous;
    	symmMatrix[8][12] = Math.log(pRatio[12] / pRatio[8]) / (pRatio[12] - pRatio[8]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AGA->ATA: IS NonSynonymous;
    	symmMatrix[8][24] = Math.log(pRatio[24] / pRatio[8]) / (pRatio[24] - pRatio[8]) / (nucleoFreqs[2] * nucleoFreqs[0]) ; 	//AGA->CGA: 
    	symmMatrix[8][40] = Math.log(pRatio[40] / pRatio[8]) / (pRatio[40] - pRatio[8]) / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//AGA->GGA: IS Transition;IS NonSynonymous;
    	symmMatrix[9][10] = Math.log(pRatio[10] / pRatio[9]) / (pRatio[10] - pRatio[9]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGC->AGG: IS NonSynonymous;
    	symmMatrix[9][11] = Math.log(pRatio[11] / pRatio[9]) / (pRatio[11] - pRatio[9]) / (nucleoFreqs[0] * nucleoFreqs[2]) * k ; 	//AGC->AGT: IS Transition;
    	symmMatrix[9][13] = Math.log(pRatio[13] / pRatio[9]) / (pRatio[13] - pRatio[9]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AGC->ATC: IS NonSynonymous;
    	symmMatrix[9][25] = Math.log(pRatio[25] / pRatio[9]) / (pRatio[25] - pRatio[9]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//AGC->CGC: IS NonSynonymous;
    	symmMatrix[9][41] = Math.log(pRatio[41] / pRatio[9]) / (pRatio[41] - pRatio[9]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//AGC->GGC: IS Transition;IS NonSynonymous;
    	symmMatrix[9][54] = Math.log(pRatio[54] / pRatio[9]) / (pRatio[54] - pRatio[9]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//AGC->TGC: IS NonSynonymous;
    	symmMatrix[10][11] = Math.log(pRatio[11] / pRatio[10]) / (pRatio[11] - pRatio[10]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGG->AGT: IS NonSynonymous;
    	symmMatrix[10][14] = Math.log(pRatio[14] / pRatio[10]) / (pRatio[14] - pRatio[10]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGG->ATG: IS NonSynonymous;
    	symmMatrix[10][26] = Math.log(pRatio[26] / pRatio[10]) / (pRatio[26] - pRatio[10]) / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//AGG->CGG: 
    	symmMatrix[10][42] = Math.log(pRatio[42] / pRatio[10]) / (pRatio[42] - pRatio[10]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//AGG->GGG: IS Transition;IS NonSynonymous;
    	symmMatrix[10][55] = Math.log(pRatio[55] / pRatio[10]) / (pRatio[55] - pRatio[10]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//AGG->TGG: IS NonSynonymous;
    	symmMatrix[11][15] = Math.log(pRatio[15] / pRatio[11]) / (pRatio[15] - pRatio[11]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AGT->ATT: IS NonSynonymous;
    	symmMatrix[11][27] = Math.log(pRatio[27] / pRatio[11]) / (pRatio[27] - pRatio[11]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//AGT->CGT: IS NonSynonymous;
    	symmMatrix[11][43] = Math.log(pRatio[43] / pRatio[11]) / (pRatio[43] - pRatio[11]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//AGT->GGT: IS Transition;IS NonSynonymous;
    	symmMatrix[11][56] = Math.log(pRatio[56] / pRatio[11]) / (pRatio[56] - pRatio[11]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//AGT->TGT: IS NonSynonymous;
    	symmMatrix[12][13] = Math.log(pRatio[13] / pRatio[12]) / (pRatio[13] - pRatio[12]) / (nucleoFreqs[0] * nucleoFreqs[3]) ; 	//ATA->ATC: 
    	symmMatrix[12][14] = Math.log(pRatio[14] / pRatio[12]) / (pRatio[14] - pRatio[12]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//ATA->ATG: IS Transition;IS NonSynonymous;
    	symmMatrix[12][15] = Math.log(pRatio[15] / pRatio[12]) / (pRatio[15] - pRatio[12]) / (nucleoFreqs[0] * nucleoFreqs[3]) ; 	//ATA->ATT: 
    	symmMatrix[12][28] = Math.log(pRatio[28] / pRatio[12]) / (pRatio[28] - pRatio[12]) / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//ATA->CTA: IS NonSynonymous;
    	symmMatrix[12][44] = Math.log(pRatio[44] / pRatio[12]) / (pRatio[44] - pRatio[12]) / (nucleoFreqs[3] * nucleoFreqs[0]) * k * omega ; 	//ATA->GTA: IS Transition;IS NonSynonymous;
    	symmMatrix[12][57] = Math.log(pRatio[57] / pRatio[12]) / (pRatio[57] - pRatio[12]) / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//ATA->TTA: IS NonSynonymous;
    	symmMatrix[13][14] = Math.log(pRatio[14] / pRatio[13]) / (pRatio[14] - pRatio[13]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ATC->ATG: IS NonSynonymous;
    	symmMatrix[13][15] = Math.log(pRatio[15] / pRatio[13]) / (pRatio[15] - pRatio[13]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k ; 	//ATC->ATT: IS Transition;
    	symmMatrix[13][29] = Math.log(pRatio[29] / pRatio[13]) / (pRatio[29] - pRatio[13]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//ATC->CTC: IS NonSynonymous;
    	symmMatrix[13][45] = Math.log(pRatio[45] / pRatio[13]) / (pRatio[45] - pRatio[13]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//ATC->GTC: IS Transition;IS NonSynonymous;
    	symmMatrix[13][58] = Math.log(pRatio[58] / pRatio[13]) / (pRatio[58] - pRatio[13]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//ATC->TTC: IS NonSynonymous;
    	symmMatrix[14][15] = Math.log(pRatio[15] / pRatio[14]) / (pRatio[15] - pRatio[14]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ATG->ATT: IS NonSynonymous;
    	symmMatrix[14][30] = Math.log(pRatio[30] / pRatio[14]) / (pRatio[30] - pRatio[14]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//ATG->CTG: IS NonSynonymous;
    	symmMatrix[14][46] = Math.log(pRatio[46] / pRatio[14]) / (pRatio[46] - pRatio[14]) / (nucleoFreqs[3] * nucleoFreqs[2]) * k * omega ; 	//ATG->GTG: IS Transition;IS NonSynonymous;
    	symmMatrix[14][59] = Math.log(pRatio[59] / pRatio[14]) / (pRatio[59] - pRatio[14]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//ATG->TTG: IS NonSynonymous;
    	symmMatrix[15][31] = Math.log(pRatio[31] / pRatio[15]) / (pRatio[31] - pRatio[15]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//ATT->CTT: IS NonSynonymous;
    	symmMatrix[15][47] = Math.log(pRatio[47] / pRatio[15]) / (pRatio[47] - pRatio[15]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//ATT->GTT: IS Transition;IS NonSynonymous;
    	symmMatrix[15][60] = Math.log(pRatio[60] / pRatio[15]) / (pRatio[60] - pRatio[15]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//ATT->TTT: IS NonSynonymous;
    	symmMatrix[16][17] = Math.log(pRatio[17] / pRatio[16]) / (pRatio[17] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CAC: IS NonSynonymous;
    	symmMatrix[16][18] = Math.log(pRatio[18] / pRatio[16]) / (pRatio[18] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k ; 	//CAA->CAG: IS Transition;
    	symmMatrix[16][19] = Math.log(pRatio[19] / pRatio[16]) / (pRatio[19] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CAT: IS NonSynonymous;
    	symmMatrix[16][20] = Math.log(pRatio[20] / pRatio[16]) / (pRatio[20] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CCA: IS NonSynonymous;
    	symmMatrix[16][24] = Math.log(pRatio[24] / pRatio[16]) / (pRatio[24] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CAA->CGA: IS Transition;IS NonSynonymous;
    	symmMatrix[16][28] = Math.log(pRatio[28] / pRatio[16]) / (pRatio[28] - pRatio[16]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CTA: IS NonSynonymous;
    	symmMatrix[16][32] = Math.log(pRatio[32] / pRatio[16]) / (pRatio[32] - pRatio[16]) / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//CAA->GAA: IS NonSynonymous;
    	symmMatrix[17][18] = Math.log(pRatio[18] / pRatio[17]) / (pRatio[18] - pRatio[17]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAC->CAG: IS NonSynonymous;
    	symmMatrix[17][19] = Math.log(pRatio[19] / pRatio[17]) / (pRatio[19] - pRatio[17]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k ; 	//CAC->CAT: IS Transition;
    	symmMatrix[17][21] = Math.log(pRatio[21] / pRatio[17]) / (pRatio[21] - pRatio[17]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CAC->CCC: IS NonSynonymous;
    	symmMatrix[17][25] = Math.log(pRatio[25] / pRatio[17]) / (pRatio[25] - pRatio[17]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CAC->CGC: IS Transition;IS NonSynonymous;
    	symmMatrix[17][29] = Math.log(pRatio[29] / pRatio[17]) / (pRatio[29] - pRatio[17]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CAC->CTC: IS NonSynonymous;
    	symmMatrix[17][33] = Math.log(pRatio[33] / pRatio[17]) / (pRatio[33] - pRatio[17]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//CAC->GAC: IS NonSynonymous;
    	symmMatrix[17][48] = Math.log(pRatio[48] / pRatio[17]) / (pRatio[48] - pRatio[17]) / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//CAC->TAC: IS Transition;IS NonSynonymous;
    	symmMatrix[18][19] = Math.log(pRatio[19] / pRatio[18]) / (pRatio[19] - pRatio[18]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAG->CAT: IS NonSynonymous;
    	symmMatrix[18][22] = Math.log(pRatio[22] / pRatio[18]) / (pRatio[22] - pRatio[18]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CAG->CCG: IS NonSynonymous;
    	symmMatrix[18][26] = Math.log(pRatio[26] / pRatio[18]) / (pRatio[26] - pRatio[18]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CAG->CGG: IS Transition;IS NonSynonymous;
    	symmMatrix[18][30] = Math.log(pRatio[30] / pRatio[18]) / (pRatio[30] - pRatio[18]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CAG->CTG: IS NonSynonymous;
    	symmMatrix[18][34] = Math.log(pRatio[34] / pRatio[18]) / (pRatio[34] - pRatio[18]) / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//CAG->GAG: IS NonSynonymous;
    	symmMatrix[19][23] = Math.log(pRatio[23] / pRatio[19]) / (pRatio[23] - pRatio[19]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CAT->CCT: IS NonSynonymous;
    	symmMatrix[19][27] = Math.log(pRatio[27] / pRatio[19]) / (pRatio[27] - pRatio[19]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CAT->CGT: IS Transition;IS NonSynonymous;
    	symmMatrix[19][31] = Math.log(pRatio[31] / pRatio[19]) / (pRatio[31] - pRatio[19]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CAT->CTT: IS NonSynonymous;
    	symmMatrix[19][35] = Math.log(pRatio[35] / pRatio[19]) / (pRatio[35] - pRatio[19]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//CAT->GAT: IS NonSynonymous;
    	symmMatrix[19][49] = Math.log(pRatio[49] / pRatio[19]) / (pRatio[49] - pRatio[19]) / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//CAT->TAT: IS Transition;IS NonSynonymous;
    	symmMatrix[20][21] = Math.log(pRatio[21] / pRatio[20]) / (pRatio[21] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCA->CCC: 
    	symmMatrix[20][22] = Math.log(pRatio[22] / pRatio[20]) / (pRatio[22] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k ; 	//CCA->CCG: IS Transition;
    	symmMatrix[20][23] = Math.log(pRatio[23] / pRatio[20]) / (pRatio[23] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCA->CCT: 
    	symmMatrix[20][24] = Math.log(pRatio[24] / pRatio[20]) / (pRatio[24] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CCA->CGA: IS NonSynonymous;
    	symmMatrix[20][28] = Math.log(pRatio[28] / pRatio[20]) / (pRatio[28] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CCA->CTA: IS Transition;IS NonSynonymous;
    	symmMatrix[20][36] = Math.log(pRatio[36] / pRatio[20]) / (pRatio[36] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CCA->GCA: IS NonSynonymous;
    	symmMatrix[20][50] = Math.log(pRatio[50] / pRatio[20]) / (pRatio[50] - pRatio[20]) / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CCA->TCA: IS Transition;IS NonSynonymous;
    	symmMatrix[21][22] = Math.log(pRatio[22] / pRatio[21]) / (pRatio[22] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCC->CCG: 
    	symmMatrix[21][23] = Math.log(pRatio[23] / pRatio[21]) / (pRatio[23] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k ; 	//CCC->CCT: IS Transition;
    	symmMatrix[21][25] = Math.log(pRatio[25] / pRatio[21]) / (pRatio[25] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CCC->CGC: IS NonSynonymous;
    	symmMatrix[21][29] = Math.log(pRatio[29] / pRatio[21]) / (pRatio[29] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CCC->CTC: IS Transition;IS NonSynonymous;
    	symmMatrix[21][37] = Math.log(pRatio[37] / pRatio[21]) / (pRatio[37] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CCC->GCC: IS NonSynonymous;
    	symmMatrix[21][51] = Math.log(pRatio[51] / pRatio[21]) / (pRatio[51] - pRatio[21]) / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CCC->TCC: IS Transition;IS NonSynonymous;
    	symmMatrix[22][23] = Math.log(pRatio[23] / pRatio[22]) / (pRatio[23] - pRatio[22]) / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCG->CCT: 
    	symmMatrix[22][26] = Math.log(pRatio[26] / pRatio[22]) / (pRatio[26] - pRatio[22]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CCG->CGG: IS NonSynonymous;
    	symmMatrix[22][30] = Math.log(pRatio[30] / pRatio[22]) / (pRatio[30] - pRatio[22]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CCG->CTG: IS Transition;IS NonSynonymous;
    	symmMatrix[22][38] = Math.log(pRatio[38] / pRatio[22]) / (pRatio[38] - pRatio[22]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CCG->GCG: IS NonSynonymous;
    	symmMatrix[22][52] = Math.log(pRatio[52] / pRatio[22]) / (pRatio[52] - pRatio[22]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CCG->TCG: IS Transition;IS NonSynonymous;
    	symmMatrix[23][27] = Math.log(pRatio[27] / pRatio[23]) / (pRatio[27] - pRatio[23]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CCT->CGT: IS NonSynonymous;
    	symmMatrix[23][31] = Math.log(pRatio[31] / pRatio[23]) / (pRatio[31] - pRatio[23]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CCT->CTT: IS Transition;IS NonSynonymous;
    	symmMatrix[23][39] = Math.log(pRatio[39] / pRatio[23]) / (pRatio[39] - pRatio[23]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CCT->GCT: IS NonSynonymous;
    	symmMatrix[23][53] = Math.log(pRatio[53] / pRatio[23]) / (pRatio[53] - pRatio[23]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CCT->TCT: IS Transition;IS NonSynonymous;
    	symmMatrix[24][25] = Math.log(pRatio[25] / pRatio[24]) / (pRatio[25] - pRatio[24]) / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGA->CGC: 
    	symmMatrix[24][26] = Math.log(pRatio[26] / pRatio[24]) / (pRatio[26] - pRatio[24]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k ; 	//CGA->CGG: IS Transition;
    	symmMatrix[24][27] = Math.log(pRatio[27] / pRatio[24]) / (pRatio[27] - pRatio[24]) / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGA->CGT: 
    	symmMatrix[24][28] = Math.log(pRatio[28] / pRatio[24]) / (pRatio[28] - pRatio[24]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CGA->CTA: IS NonSynonymous;
    	symmMatrix[24][40] = Math.log(pRatio[40] / pRatio[24]) / (pRatio[40] - pRatio[24]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//CGA->GGA: IS NonSynonymous;
    	symmMatrix[25][26] = Math.log(pRatio[26] / pRatio[25]) / (pRatio[26] - pRatio[25]) / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGC->CGG: 
    	symmMatrix[25][27] = Math.log(pRatio[27] / pRatio[25]) / (pRatio[27] - pRatio[25]) / (nucleoFreqs[1] * nucleoFreqs[2]) * k ; 	//CGC->CGT: IS Transition;
    	symmMatrix[25][29] = Math.log(pRatio[29] / pRatio[25]) / (pRatio[29] - pRatio[25]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CGC->CTC: IS NonSynonymous;
    	symmMatrix[25][41] = Math.log(pRatio[41] / pRatio[25]) / (pRatio[41] - pRatio[25]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//CGC->GGC: IS NonSynonymous;
    	symmMatrix[25][54] = Math.log(pRatio[54] / pRatio[25]) / (pRatio[54] - pRatio[25]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//CGC->TGC: IS Transition;IS NonSynonymous;
    	symmMatrix[26][27] = Math.log(pRatio[27] / pRatio[26]) / (pRatio[27] - pRatio[26]) / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGG->CGT: 
    	symmMatrix[26][30] = Math.log(pRatio[30] / pRatio[26]) / (pRatio[30] - pRatio[26]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CGG->CTG: IS NonSynonymous;
    	symmMatrix[26][42] = Math.log(pRatio[42] / pRatio[26]) / (pRatio[42] - pRatio[26]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//CGG->GGG: IS NonSynonymous;
    	symmMatrix[26][55] = Math.log(pRatio[55] / pRatio[26]) / (pRatio[55] - pRatio[26]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//CGG->TGG: IS Transition;IS NonSynonymous;
    	symmMatrix[27][31] = Math.log(pRatio[31] / pRatio[27]) / (pRatio[31] - pRatio[27]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CGT->CTT: IS NonSynonymous;
    	symmMatrix[27][43] = Math.log(pRatio[43] / pRatio[27]) / (pRatio[43] - pRatio[27]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//CGT->GGT: IS NonSynonymous;
    	symmMatrix[27][56] = Math.log(pRatio[56] / pRatio[27]) / (pRatio[56] - pRatio[27]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//CGT->TGT: IS Transition;IS NonSynonymous;
    	symmMatrix[28][29] = Math.log(pRatio[29] / pRatio[28]) / (pRatio[29] - pRatio[28]) / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTA->CTC: 
    	symmMatrix[28][30] = Math.log(pRatio[30] / pRatio[28]) / (pRatio[30] - pRatio[28]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k ; 	//CTA->CTG: IS Transition;
    	symmMatrix[28][31] = Math.log(pRatio[31] / pRatio[28]) / (pRatio[31] - pRatio[28]) / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTA->CTT: 
    	symmMatrix[28][44] = Math.log(pRatio[44] / pRatio[28]) / (pRatio[44] - pRatio[28]) / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//CTA->GTA: IS NonSynonymous;
    	symmMatrix[28][57] = Math.log(pRatio[57] / pRatio[28]) / (pRatio[57] - pRatio[28]) / (nucleoFreqs[3] * nucleoFreqs[0]) * k ; 	//CTA->TTA: IS Transition;
    	symmMatrix[29][30] = Math.log(pRatio[30] / pRatio[29]) / (pRatio[30] - pRatio[29]) / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTC->CTG: 
    	symmMatrix[29][31] = Math.log(pRatio[31] / pRatio[29]) / (pRatio[31] - pRatio[29]) / (nucleoFreqs[1] * nucleoFreqs[3]) * k ; 	//CTC->CTT: IS Transition;
    	symmMatrix[29][45] = Math.log(pRatio[45] / pRatio[29]) / (pRatio[45] - pRatio[29]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//CTC->GTC: IS NonSynonymous;
    	symmMatrix[29][58] = Math.log(pRatio[58] / pRatio[29]) / (pRatio[58] - pRatio[29]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//CTC->TTC: IS Transition;IS NonSynonymous;
    	symmMatrix[30][31] = Math.log(pRatio[31] / pRatio[30]) / (pRatio[31] - pRatio[30]) / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTG->CTT: 
    	symmMatrix[30][46] = Math.log(pRatio[46] / pRatio[30]) / (pRatio[46] - pRatio[30]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//CTG->GTG: IS NonSynonymous;
    	symmMatrix[30][59] = Math.log(pRatio[59] / pRatio[30]) / (pRatio[59] - pRatio[30]) / (nucleoFreqs[3] * nucleoFreqs[2]) * k ; 	//CTG->TTG: IS Transition;
    	symmMatrix[31][47] = Math.log(pRatio[47] / pRatio[31]) / (pRatio[47] - pRatio[31]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//CTT->GTT: IS NonSynonymous;
    	symmMatrix[31][60] = Math.log(pRatio[60] / pRatio[31]) / (pRatio[60] - pRatio[31]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//CTT->TTT: IS Transition;IS NonSynonymous;
    	symmMatrix[32][33] = Math.log(pRatio[33] / pRatio[32]) / (pRatio[33] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GAC: IS NonSynonymous;
    	symmMatrix[32][34] = Math.log(pRatio[34] / pRatio[32]) / (pRatio[34] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * k ; 	//GAA->GAG: IS Transition;
    	symmMatrix[32][35] = Math.log(pRatio[35] / pRatio[32]) / (pRatio[35] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GAT: IS NonSynonymous;
    	symmMatrix[32][36] = Math.log(pRatio[36] / pRatio[32]) / (pRatio[36] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GCA: IS NonSynonymous;
    	symmMatrix[32][40] = Math.log(pRatio[40] / pRatio[32]) / (pRatio[40] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//GAA->GGA: IS Transition;IS NonSynonymous;
    	symmMatrix[32][44] = Math.log(pRatio[44] / pRatio[32]) / (pRatio[44] - pRatio[32]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GTA: IS NonSynonymous;
    	symmMatrix[33][34] = Math.log(pRatio[34] / pRatio[33]) / (pRatio[34] - pRatio[33]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAC->GAG: IS NonSynonymous;
    	symmMatrix[33][35] = Math.log(pRatio[35] / pRatio[33]) / (pRatio[35] - pRatio[33]) / (nucleoFreqs[2] * nucleoFreqs[0]) * k ; 	//GAC->GAT: IS Transition;
    	symmMatrix[33][37] = Math.log(pRatio[37] / pRatio[33]) / (pRatio[37] - pRatio[33]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GAC->GCC: IS NonSynonymous;
    	symmMatrix[33][41] = Math.log(pRatio[41] / pRatio[33]) / (pRatio[41] - pRatio[33]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//GAC->GGC: IS Transition;IS NonSynonymous;
    	symmMatrix[33][45] = Math.log(pRatio[45] / pRatio[33]) / (pRatio[45] - pRatio[33]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GAC->GTC: IS NonSynonymous;
    	symmMatrix[33][48] = Math.log(pRatio[48] / pRatio[33]) / (pRatio[48] - pRatio[33]) / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//GAC->TAC: IS NonSynonymous;
    	symmMatrix[34][35] = Math.log(pRatio[35] / pRatio[34]) / (pRatio[35] - pRatio[34]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAG->GAT: IS NonSynonymous;
    	symmMatrix[34][38] = Math.log(pRatio[38] / pRatio[34]) / (pRatio[38] - pRatio[34]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GAG->GCG: IS NonSynonymous;
    	symmMatrix[34][42] = Math.log(pRatio[42] / pRatio[34]) / (pRatio[42] - pRatio[34]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//GAG->GGG: IS Transition;IS NonSynonymous;
    	symmMatrix[34][46] = Math.log(pRatio[46] / pRatio[34]) / (pRatio[46] - pRatio[34]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GAG->GTG: IS NonSynonymous;
    	symmMatrix[35][39] = Math.log(pRatio[39] / pRatio[35]) / (pRatio[39] - pRatio[35]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GAT->GCT: IS NonSynonymous;
    	symmMatrix[35][43] = Math.log(pRatio[43] / pRatio[35]) / (pRatio[43] - pRatio[35]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//GAT->GGT: IS Transition;IS NonSynonymous;
    	symmMatrix[35][47] = Math.log(pRatio[47] / pRatio[35]) / (pRatio[47] - pRatio[35]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GAT->GTT: IS NonSynonymous;
    	symmMatrix[35][49] = Math.log(pRatio[49] / pRatio[35]) / (pRatio[49] - pRatio[35]) / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//GAT->TAT: IS NonSynonymous;
    	symmMatrix[36][37] = Math.log(pRatio[37] / pRatio[36]) / (pRatio[37] - pRatio[36]) / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCA->GCC: 
    	symmMatrix[36][38] = Math.log(pRatio[38] / pRatio[36]) / (pRatio[38] - pRatio[36]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k ; 	//GCA->GCG: IS Transition;
    	symmMatrix[36][39] = Math.log(pRatio[39] / pRatio[36]) / (pRatio[39] - pRatio[36]) / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCA->GCT: 
    	symmMatrix[36][40] = Math.log(pRatio[40] / pRatio[36]) / (pRatio[40] - pRatio[36]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GCA->GGA: IS NonSynonymous;
    	symmMatrix[36][44] = Math.log(pRatio[44] / pRatio[36]) / (pRatio[44] - pRatio[36]) / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//GCA->GTA: IS Transition;IS NonSynonymous;
    	symmMatrix[36][50] = Math.log(pRatio[50] / pRatio[36]) / (pRatio[50] - pRatio[36]) / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//GCA->TCA: IS NonSynonymous;
    	symmMatrix[37][38] = Math.log(pRatio[38] / pRatio[37]) / (pRatio[38] - pRatio[37]) / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCC->GCG: 
    	symmMatrix[37][39] = Math.log(pRatio[39] / pRatio[37]) / (pRatio[39] - pRatio[37]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k ; 	//GCC->GCT: IS Transition;
    	symmMatrix[37][41] = Math.log(pRatio[41] / pRatio[37]) / (pRatio[41] - pRatio[37]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GCC->GGC: IS NonSynonymous;
    	symmMatrix[37][45] = Math.log(pRatio[45] / pRatio[37]) / (pRatio[45] - pRatio[37]) / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//GCC->GTC: IS Transition;IS NonSynonymous;
    	symmMatrix[37][51] = Math.log(pRatio[51] / pRatio[37]) / (pRatio[51] - pRatio[37]) / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//GCC->TCC: IS NonSynonymous;
    	symmMatrix[38][39] = Math.log(pRatio[39] / pRatio[38]) / (pRatio[39] - pRatio[38]) / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCG->GCT: 
    	symmMatrix[38][42] = Math.log(pRatio[42] / pRatio[38]) / (pRatio[42] - pRatio[38]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GCG->GGG: IS NonSynonymous;
    	symmMatrix[38][46] = Math.log(pRatio[46] / pRatio[38]) / (pRatio[46] - pRatio[38]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//GCG->GTG: IS Transition;IS NonSynonymous;
    	symmMatrix[38][52] = Math.log(pRatio[52] / pRatio[38]) / (pRatio[52] - pRatio[38]) / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//GCG->TCG: IS NonSynonymous;
    	symmMatrix[39][43] = Math.log(pRatio[43] / pRatio[39]) / (pRatio[43] - pRatio[39]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GCT->GGT: IS NonSynonymous;
    	symmMatrix[39][47] = Math.log(pRatio[47] / pRatio[39]) / (pRatio[47] - pRatio[39]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//GCT->GTT: IS Transition;IS NonSynonymous;
    	symmMatrix[39][53] = Math.log(pRatio[53] / pRatio[39]) / (pRatio[53] - pRatio[39]) / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//GCT->TCT: IS NonSynonymous;
    	symmMatrix[40][41] = Math.log(pRatio[41] / pRatio[40]) / (pRatio[41] - pRatio[40]) / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGA->GGC: 
    	symmMatrix[40][42] = Math.log(pRatio[42] / pRatio[40]) / (pRatio[42] - pRatio[40]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k ; 	//GGA->GGG: IS Transition;
    	symmMatrix[40][43] = Math.log(pRatio[43] / pRatio[40]) / (pRatio[43] - pRatio[40]) / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGA->GGT: 
    	symmMatrix[40][44] = Math.log(pRatio[44] / pRatio[40]) / (pRatio[44] - pRatio[40]) / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GGA->GTA: IS NonSynonymous;
    	symmMatrix[41][42] = Math.log(pRatio[42] / pRatio[41]) / (pRatio[42] - pRatio[41]) / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGC->GGG: 
    	symmMatrix[41][43] = Math.log(pRatio[43] / pRatio[41]) / (pRatio[43] - pRatio[41]) / (nucleoFreqs[2] * nucleoFreqs[2]) * k ; 	//GGC->GGT: IS Transition;
    	symmMatrix[41][45] = Math.log(pRatio[45] / pRatio[41]) / (pRatio[45] - pRatio[41]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GGC->GTC: IS NonSynonymous;
    	symmMatrix[41][54] = Math.log(pRatio[54] / pRatio[41]) / (pRatio[54] - pRatio[41]) / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GGC->TGC: IS NonSynonymous;
    	symmMatrix[42][43] = Math.log(pRatio[43] / pRatio[42]) / (pRatio[43] - pRatio[42]) / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGG->GGT: 
    	symmMatrix[42][46] = Math.log(pRatio[46] / pRatio[42]) / (pRatio[46] - pRatio[42]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GGG->GTG: IS NonSynonymous;
    	symmMatrix[42][55] = Math.log(pRatio[55] / pRatio[42]) / (pRatio[55] - pRatio[42]) / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GGG->TGG: IS NonSynonymous;
    	symmMatrix[43][47] = Math.log(pRatio[47] / pRatio[43]) / (pRatio[47] - pRatio[43]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GGT->GTT: IS NonSynonymous;
    	symmMatrix[43][56] = Math.log(pRatio[56] / pRatio[43]) / (pRatio[56] - pRatio[43]) / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GGT->TGT: IS NonSynonymous;
    	symmMatrix[44][45] = Math.log(pRatio[45] / pRatio[44]) / (pRatio[45] - pRatio[44]) / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTA->GTC: 
    	symmMatrix[44][46] = Math.log(pRatio[46] / pRatio[44]) / (pRatio[46] - pRatio[44]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k ; 	//GTA->GTG: IS Transition;
    	symmMatrix[44][47] = Math.log(pRatio[47] / pRatio[44]) / (pRatio[47] - pRatio[44]) / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTA->GTT: 
    	symmMatrix[44][57] = Math.log(pRatio[57] / pRatio[44]) / (pRatio[57] - pRatio[44]) / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//GTA->TTA: IS NonSynonymous;
    	symmMatrix[45][46] = Math.log(pRatio[46] / pRatio[45]) / (pRatio[46] - pRatio[45]) / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTC->GTG: 
    	symmMatrix[45][47] = Math.log(pRatio[47] / pRatio[45]) / (pRatio[47] - pRatio[45]) / (nucleoFreqs[2] * nucleoFreqs[3]) * k ; 	//GTC->GTT: IS Transition;
    	symmMatrix[45][58] = Math.log(pRatio[58] / pRatio[45]) / (pRatio[58] - pRatio[45]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//GTC->TTC: IS NonSynonymous;
    	symmMatrix[46][47] = Math.log(pRatio[47] / pRatio[46]) / (pRatio[47] - pRatio[46]) / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTG->GTT: 
    	symmMatrix[46][59] = Math.log(pRatio[59] / pRatio[46]) / (pRatio[59] - pRatio[46]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//GTG->TTG: IS NonSynonymous;
    	symmMatrix[47][60] = Math.log(pRatio[60] / pRatio[47]) / (pRatio[60] - pRatio[47]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//GTT->TTT: IS NonSynonymous;
    	symmMatrix[48][49] = Math.log(pRatio[49] / pRatio[48]) / (pRatio[49] - pRatio[48]) / (nucleoFreqs[3] * nucleoFreqs[0]) * k ; 	//TAC->TAT: IS Transition;
    	symmMatrix[48][51] = Math.log(pRatio[51] / pRatio[48]) / (pRatio[51] - pRatio[48]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TAC->TCC: IS NonSynonymous;
    	symmMatrix[48][54] = Math.log(pRatio[54] / pRatio[48]) / (pRatio[54] - pRatio[48]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//TAC->TGC: IS Transition;IS NonSynonymous;
    	symmMatrix[48][58] = Math.log(pRatio[58] / pRatio[48]) / (pRatio[58] - pRatio[48]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TAC->TTC: IS NonSynonymous;
    	symmMatrix[49][53] = Math.log(pRatio[53] / pRatio[49]) / (pRatio[53] - pRatio[49]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TAT->TCT: IS NonSynonymous;
    	symmMatrix[49][56] = Math.log(pRatio[56] / pRatio[49]) / (pRatio[56] - pRatio[49]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//TAT->TGT: IS Transition;IS NonSynonymous;
    	symmMatrix[49][60] = Math.log(pRatio[60] / pRatio[49]) / (pRatio[60] - pRatio[49]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TAT->TTT: IS NonSynonymous;
    	symmMatrix[50][51] = Math.log(pRatio[51] / pRatio[50]) / (pRatio[51] - pRatio[50]) / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCA->TCC: 
    	symmMatrix[50][52] = Math.log(pRatio[52] / pRatio[50]) / (pRatio[52] - pRatio[50]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k ; 	//TCA->TCG: IS Transition;
    	symmMatrix[50][53] = Math.log(pRatio[53] / pRatio[50]) / (pRatio[53] - pRatio[50]) / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCA->TCT: 
    	symmMatrix[50][57] = Math.log(pRatio[57] / pRatio[50]) / (pRatio[57] - pRatio[50]) / (nucleoFreqs[3] * nucleoFreqs[0]) * k * omega ; 	//TCA->TTA: IS Transition;IS NonSynonymous;
    	symmMatrix[51][52] = Math.log(pRatio[52] / pRatio[51]) / (pRatio[52] - pRatio[51]) / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCC->TCG: 
    	symmMatrix[51][53] = Math.log(pRatio[53] / pRatio[51]) / (pRatio[53] - pRatio[51]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k ; 	//TCC->TCT: IS Transition;
    	symmMatrix[51][54] = Math.log(pRatio[54] / pRatio[51]) / (pRatio[54] - pRatio[51]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TCC->TGC: IS NonSynonymous;
    	symmMatrix[51][58] = Math.log(pRatio[58] / pRatio[51]) / (pRatio[58] - pRatio[51]) / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//TCC->TTC: IS Transition;IS NonSynonymous;
    	symmMatrix[52][53] = Math.log(pRatio[53] / pRatio[52]) / (pRatio[53] - pRatio[52]) / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCG->TCT: 
    	symmMatrix[52][55] = Math.log(pRatio[55] / pRatio[52]) / (pRatio[55] - pRatio[52]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TCG->TGG: IS NonSynonymous;
    	symmMatrix[52][59] = Math.log(pRatio[59] / pRatio[52]) / (pRatio[59] - pRatio[52]) / (nucleoFreqs[3] * nucleoFreqs[2]) * k * omega ; 	//TCG->TTG: IS Transition;IS NonSynonymous;
    	symmMatrix[53][56] = Math.log(pRatio[56] / pRatio[53]) / (pRatio[56] - pRatio[53]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TCT->TGT: IS NonSynonymous;
    	symmMatrix[53][60] = Math.log(pRatio[60] / pRatio[53]) / (pRatio[60] - pRatio[53]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//TCT->TTT: IS Transition;IS NonSynonymous;
    	symmMatrix[54][55] = Math.log(pRatio[55] / pRatio[54]) / (pRatio[55] - pRatio[54]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGC->TGG: IS NonSynonymous;
    	symmMatrix[54][56] = Math.log(pRatio[56] / pRatio[54]) / (pRatio[56] - pRatio[54]) / (nucleoFreqs[3] * nucleoFreqs[2]) * k ; 	//TGC->TGT: IS Transition;
    	symmMatrix[54][58] = Math.log(pRatio[58] / pRatio[54]) / (pRatio[58] - pRatio[54]) / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TGC->TTC: IS NonSynonymous;
    	symmMatrix[55][56] = Math.log(pRatio[56] / pRatio[55]) / (pRatio[56] - pRatio[55]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGG->TGT: IS NonSynonymous;
    	symmMatrix[55][59] = Math.log(pRatio[59] / pRatio[55]) / (pRatio[59] - pRatio[55]) / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGG->TTG: IS NonSynonymous;
    	symmMatrix[56][60] = Math.log(pRatio[60] / pRatio[56]) / (pRatio[60] - pRatio[56]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TGT->TTT: IS NonSynonymous;
    	symmMatrix[57][58] = Math.log(pRatio[58] / pRatio[57]) / (pRatio[58] - pRatio[57]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTA->TTC: IS NonSynonymous;
    	symmMatrix[57][59] = Math.log(pRatio[59] / pRatio[57]) / (pRatio[59] - pRatio[57]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k ; 	//TTA->TTG: IS Transition;
    	symmMatrix[57][60] = Math.log(pRatio[60] / pRatio[57]) / (pRatio[60] - pRatio[57]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTA->TTT: IS NonSynonymous;
    	symmMatrix[58][59] = Math.log(pRatio[59] / pRatio[58]) / (pRatio[59] - pRatio[58]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTC->TTG: IS NonSynonymous;
    	symmMatrix[58][60] = Math.log(pRatio[60] / pRatio[58]) / (pRatio[60] - pRatio[58]) / (nucleoFreqs[3] * nucleoFreqs[3]) * k ; 	//TTC->TTT: IS Transition;
    	symmMatrix[59][60] = Math.log(pRatio[60] / pRatio[59]) / (pRatio[60] - pRatio[59]) / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTG->TTT: IS NonSynonymous;
    	
    	//let lower triangle equal to upper
        for (int i = 0; i < nrOfStates; i++) {
        	for (int j = 0; j < nrOfStates; j++) {
    			if (i < j){
    				symmMatrix[j][i] = symmMatrix[i][j];
    			}
        	}
        }
        
    	//System.out.print("symmM after setup:");
        //System.out.println(Arrays.deepToString(symmMatrix));
    }
    
    @Override
    protected void setupRelativeRates() {}
    
    @Override
    protected void setupRateMatrix() {
    	setupSymmMatrix();
    	setupDiagMatrix();
        
    	//multiply symmMatrix and diagMatrix together to get rateMatrix
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    rateMatrix[i][j] = symmMatrix[i][j] * diagMatrix[j];
            }
        }
    	
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }

        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
        	//please see Yang and Nielsen (2008) for explanation (symmMatrix times diagonal, while diagonal matrix represents stationary for codons)
            fSubst += -rateMatrix[i][i] * diagMatrix[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
    }
    
    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        synchronized (this) {
            if (updateMatrix) {
                setupRelativeRates();
                setupRateMatrix();
                double[][] copyRateMatrix = new double[nrOfStates][nrOfStates];
        		for (int rowNr=0; rowNr < rateMatrix.length; rowNr++){
        			for (int colNr=0; colNr < rateMatrix.length; colNr++){
        				copyRateMatrix[rowNr][colNr] = rateMatrix[rowNr][colNr];
        			}
        		}
                eigenDecomposition = eigenSystem.decomposeMatrix(copyRateMatrix);
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }
    
    @Override
    public double[] getFrequencies() {
    	double[] codonP = getCodonProb();
    	return codonP;
    }
    
    @Override
    public boolean canHandleDataType(DataType dataType) throws Exception {
        if (dataType instanceof Codon) {
            return true;
        }
        if (dataType instanceof UserDataType) {
        	if (dataType.getStateCount() == 61){
        		return true;
        	}
        }
        throw new Exception("Can only handle codon data");
    }
    
    /***************************************************************************************/
    //for testing purpose, need to be disabled after testing
    /***************************************************************************************/    
    /**
     * access to (copy of) rate matrix *
     */
    @Override
    public double[][] getRateMatrix() {
        return rateMatrix.clone();
    }
    
    public double[][] getSymmMatrix() {
        return symmMatrix.clone();
    }
    
    public void prepareMatricesForTest(){
        setupRelativeRates();
        setupRateMatrix();
    }
    
    public double[] getDiagMatrix() {
        return diagMatrix.clone();
    }
    
}