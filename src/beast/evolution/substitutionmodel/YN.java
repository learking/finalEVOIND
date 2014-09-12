package beast.evolution.substitutionmodel;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;

public class YN extends SubstitutionModel.Base {
	
    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "omega parameter in YN98 model", Validate.REQUIRED);
    public Input<Frequencies> nucleoFreqInput =
            new Input<Frequencies>("nucleoFrequencies", "substitution model equilibrium state frequencies", Validate.REQUIRED);
    
    double[][] rateMatrix;    
    Frequencies nucleoFrequencies;
    
    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;
    
    protected EigenSystem eigenSystem;
    
    protected EigenDecomposition eigenDecomposition;
    private EigenDecomposition storedEigenDecomposition;
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public YN() {
    	frequenciesInput.setRule(Validate.OPTIONAL);
        try {
        	frequenciesInput.setValue(null, this);
        } catch (Exception e) {
        	e.printStackTrace();
			// TODO: handle exception
		}
    }
    
    @Override
    public void initAndValidate() throws Exception {
        nucleoFrequencies = nucleoFreqInput.get();
    	updateMatrix = true;
        nrOfStates = 61;

        eigenSystem = new DefaultEigenSystem(nrOfStates);
        
        rateMatrix = new double[nrOfStates][nrOfStates];        
    } // initAndValidate
    
	@Override
	public void getTransitionProbabilities(Node node, double fStartTime,
			double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;

        int i, j, k;
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupRateMatrix();
                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }

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
	}

    /**
     * sets up rate matrix *
     */
    protected void setupRateMatrix() {
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	final double kappa = kappaInput.get().getValue();
    	final double omega = omegaInput.get().getValue();
    	
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
        }
    	
    	//syn, transversion
        rateMatrix[4][5] = nucleoFreqs[1]; /*ACA->ACC*/ rateMatrix[5][4] = nucleoFreqs[0]; /*ACC->ACA*/
        rateMatrix[4][7] = nucleoFreqs[3]; /*ACA->ACT*/ rateMatrix[7][4] = nucleoFreqs[0]; /*ACT->ACA*/
        rateMatrix[5][6] = nucleoFreqs[2]; /*ACC->ACG*/ rateMatrix[6][5] = nucleoFreqs[1]; /*ACG->ACC*/
        rateMatrix[6][7] = nucleoFreqs[3]; /*ACG->ACT*/ rateMatrix[7][6] = nucleoFreqs[2]; /*ACT->ACG*/
        rateMatrix[8][24] = nucleoFreqs[1]; /*AGA->CGA*/ rateMatrix[24][8] = nucleoFreqs[0]; /*CGA->AGA*/
        rateMatrix[10][26] = nucleoFreqs[1]; /*AGG->CGG*/ rateMatrix[26][10] = nucleoFreqs[0]; /*CGG->AGG*/
        rateMatrix[12][13] = nucleoFreqs[1]; /*ATA->ATC*/ rateMatrix[13][12] = nucleoFreqs[0]; /*ATC->ATA*/
        rateMatrix[12][15] = nucleoFreqs[3]; /*ATA->ATT*/ rateMatrix[15][12] = nucleoFreqs[0]; /*ATT->ATA*/
        rateMatrix[20][21] = nucleoFreqs[1]; /*CCA->CCC*/ rateMatrix[21][20] = nucleoFreqs[0]; /*CCC->CCA*/
        rateMatrix[20][23] = nucleoFreqs[3]; /*CCA->CCT*/ rateMatrix[23][20] = nucleoFreqs[0]; /*CCT->CCA*/
        rateMatrix[21][22] = nucleoFreqs[2]; /*CCC->CCG*/ rateMatrix[22][21] = nucleoFreqs[1]; /*CCG->CCC*/
        rateMatrix[22][23] = nucleoFreqs[3]; /*CCG->CCT*/ rateMatrix[23][22] = nucleoFreqs[2]; /*CCT->CCG*/
        rateMatrix[24][25] = nucleoFreqs[1]; /*CGA->CGC*/ rateMatrix[25][24] = nucleoFreqs[0]; /*CGC->CGA*/
        rateMatrix[24][27] = nucleoFreqs[3]; /*CGA->CGT*/ rateMatrix[27][24] = nucleoFreqs[0]; /*CGT->CGA*/
        rateMatrix[25][26] = nucleoFreqs[2]; /*CGC->CGG*/ rateMatrix[26][25] = nucleoFreqs[1]; /*CGG->CGC*/
        rateMatrix[26][27] = nucleoFreqs[3]; /*CGG->CGT*/ rateMatrix[27][26] = nucleoFreqs[2]; /*CGT->CGG*/
        rateMatrix[28][29] = nucleoFreqs[1]; /*CTA->CTC*/ rateMatrix[29][28] = nucleoFreqs[0]; /*CTC->CTA*/
        rateMatrix[28][31] = nucleoFreqs[3]; /*CTA->CTT*/ rateMatrix[31][28] = nucleoFreqs[0]; /*CTT->CTA*/
        rateMatrix[29][30] = nucleoFreqs[2]; /*CTC->CTG*/ rateMatrix[30][29] = nucleoFreqs[1]; /*CTG->CTC*/
        rateMatrix[30][31] = nucleoFreqs[3]; /*CTG->CTT*/ rateMatrix[31][30] = nucleoFreqs[2]; /*CTT->CTG*/
        rateMatrix[36][37] = nucleoFreqs[1]; /*GCA->GCC*/ rateMatrix[37][36] = nucleoFreqs[0]; /*GCC->GCA*/
        rateMatrix[36][39] = nucleoFreqs[3]; /*GCA->GCT*/ rateMatrix[39][36] = nucleoFreqs[0]; /*GCT->GCA*/
        rateMatrix[37][38] = nucleoFreqs[2]; /*GCC->GCG*/ rateMatrix[38][37] = nucleoFreqs[1]; /*GCG->GCC*/
        rateMatrix[38][39] = nucleoFreqs[3]; /*GCG->GCT*/ rateMatrix[39][38] = nucleoFreqs[2]; /*GCT->GCG*/
        rateMatrix[40][41] = nucleoFreqs[1]; /*GGA->GGC*/ rateMatrix[41][40] = nucleoFreqs[0]; /*GGC->GGA*/
        rateMatrix[40][43] = nucleoFreqs[3]; /*GGA->GGT*/ rateMatrix[43][40] = nucleoFreqs[0]; /*GGT->GGA*/
        rateMatrix[41][42] = nucleoFreqs[2]; /*GGC->GGG*/ rateMatrix[42][41] = nucleoFreqs[1]; /*GGG->GGC*/
        rateMatrix[42][43] = nucleoFreqs[3]; /*GGG->GGT*/ rateMatrix[43][42] = nucleoFreqs[2]; /*GGT->GGG*/
        rateMatrix[44][45] = nucleoFreqs[1]; /*GTA->GTC*/ rateMatrix[45][44] = nucleoFreqs[0]; /*GTC->GTA*/
        rateMatrix[44][47] = nucleoFreqs[3]; /*GTA->GTT*/ rateMatrix[47][44] = nucleoFreqs[0]; /*GTT->GTA*/
        rateMatrix[45][46] = nucleoFreqs[2]; /*GTC->GTG*/ rateMatrix[46][45] = nucleoFreqs[1]; /*GTG->GTC*/
        rateMatrix[46][47] = nucleoFreqs[3]; /*GTG->GTT*/ rateMatrix[47][46] = nucleoFreqs[2]; /*GTT->GTG*/
        rateMatrix[50][51] = nucleoFreqs[1]; /*TCA->TCC*/ rateMatrix[51][50] = nucleoFreqs[0]; /*TCC->TCA*/
        rateMatrix[50][53] = nucleoFreqs[3]; /*TCA->TCT*/ rateMatrix[53][50] = nucleoFreqs[0]; /*TCT->TCA*/
        rateMatrix[51][52] = nucleoFreqs[2]; /*TCC->TCG*/ rateMatrix[52][51] = nucleoFreqs[1]; /*TCG->TCC*/
        rateMatrix[52][53] = nucleoFreqs[3]; /*TCG->TCT*/ rateMatrix[53][52] = nucleoFreqs[2]; /*TCT->TCG*/
    	//syn, transition
        rateMatrix[0][2] = nucleoFreqs[2] * kappa; /*AAA->AAG*/ rateMatrix[2][0] = nucleoFreqs[0] * kappa; /*AAG->AAA*/
        rateMatrix[1][3] = nucleoFreqs[3] * kappa; /*AAC->AAT*/ rateMatrix[3][1] = nucleoFreqs[1] * kappa; /*AAT->AAC*/
        rateMatrix[4][6] = nucleoFreqs[2] * kappa; /*ACA->ACG*/ rateMatrix[6][4] = nucleoFreqs[0] * kappa; /*ACG->ACA*/
        rateMatrix[5][7] = nucleoFreqs[3] * kappa; /*ACC->ACT*/ rateMatrix[7][5] = nucleoFreqs[1] * kappa; /*ACT->ACC*/
        rateMatrix[8][10] = nucleoFreqs[2] * kappa; /*AGA->AGG*/ rateMatrix[10][8] = nucleoFreqs[0] * kappa; /*AGG->AGA*/
        rateMatrix[9][11] = nucleoFreqs[3] * kappa; /*AGC->AGT*/ rateMatrix[11][9] = nucleoFreqs[1] * kappa; /*AGT->AGC*/
        rateMatrix[13][15] = nucleoFreqs[3] * kappa; /*ATC->ATT*/ rateMatrix[15][13] = nucleoFreqs[1] * kappa; /*ATT->ATC*/
        rateMatrix[16][18] = nucleoFreqs[2] * kappa; /*CAA->CAG*/ rateMatrix[18][16] = nucleoFreqs[0] * kappa; /*CAG->CAA*/
        rateMatrix[17][19] = nucleoFreqs[3] * kappa; /*CAC->CAT*/ rateMatrix[19][17] = nucleoFreqs[1] * kappa; /*CAT->CAC*/
        rateMatrix[20][22] = nucleoFreqs[2] * kappa; /*CCA->CCG*/ rateMatrix[22][20] = nucleoFreqs[0] * kappa; /*CCG->CCA*/
        rateMatrix[21][23] = nucleoFreqs[3] * kappa; /*CCC->CCT*/ rateMatrix[23][21] = nucleoFreqs[1] * kappa; /*CCT->CCC*/
        rateMatrix[24][26] = nucleoFreqs[2] * kappa; /*CGA->CGG*/ rateMatrix[26][24] = nucleoFreqs[0] * kappa; /*CGG->CGA*/
        rateMatrix[25][27] = nucleoFreqs[3] * kappa; /*CGC->CGT*/ rateMatrix[27][25] = nucleoFreqs[1] * kappa; /*CGT->CGC*/
        rateMatrix[28][30] = nucleoFreqs[2] * kappa; /*CTA->CTG*/ rateMatrix[30][28] = nucleoFreqs[0] * kappa; /*CTG->CTA*/
        rateMatrix[28][57] = nucleoFreqs[3] * kappa; /*CTA->TTA*/ rateMatrix[57][28] = nucleoFreqs[1] * kappa; /*TTA->CTA*/
        rateMatrix[29][31] = nucleoFreqs[3] * kappa; /*CTC->CTT*/ rateMatrix[31][29] = nucleoFreqs[1] * kappa; /*CTT->CTC*/
        rateMatrix[30][59] = nucleoFreqs[3] * kappa; /*CTG->TTG*/ rateMatrix[59][30] = nucleoFreqs[1] * kappa; /*TTG->CTG*/
        rateMatrix[32][34] = nucleoFreqs[2] * kappa; /*GAA->GAG*/ rateMatrix[34][32] = nucleoFreqs[0] * kappa; /*GAG->GAA*/
        rateMatrix[33][35] = nucleoFreqs[3] * kappa; /*GAC->GAT*/ rateMatrix[35][33] = nucleoFreqs[1] * kappa; /*GAT->GAC*/
        rateMatrix[36][38] = nucleoFreqs[2] * kappa; /*GCA->GCG*/ rateMatrix[38][36] = nucleoFreqs[0] * kappa; /*GCG->GCA*/
        rateMatrix[37][39] = nucleoFreqs[3] * kappa; /*GCC->GCT*/ rateMatrix[39][37] = nucleoFreqs[1] * kappa; /*GCT->GCC*/
        rateMatrix[40][42] = nucleoFreqs[2] * kappa; /*GGA->GGG*/ rateMatrix[42][40] = nucleoFreqs[0] * kappa; /*GGG->GGA*/
        rateMatrix[41][43] = nucleoFreqs[3] * kappa; /*GGC->GGT*/ rateMatrix[43][41] = nucleoFreqs[1] * kappa; /*GGT->GGC*/
        rateMatrix[44][46] = nucleoFreqs[2] * kappa; /*GTA->GTG*/ rateMatrix[46][44] = nucleoFreqs[0] * kappa; /*GTG->GTA*/
        rateMatrix[45][47] = nucleoFreqs[3] * kappa; /*GTC->GTT*/ rateMatrix[47][45] = nucleoFreqs[1] * kappa; /*GTT->GTC*/
        rateMatrix[48][49] = nucleoFreqs[3] * kappa; /*TAC->TAT*/ rateMatrix[49][48] = nucleoFreqs[1] * kappa; /*TAT->TAC*/
        rateMatrix[50][52] = nucleoFreqs[2] * kappa; /*TCA->TCG*/ rateMatrix[52][50] = nucleoFreqs[0] * kappa; /*TCG->TCA*/
        rateMatrix[51][53] = nucleoFreqs[3] * kappa; /*TCC->TCT*/ rateMatrix[53][51] = nucleoFreqs[1] * kappa; /*TCT->TCC*/
        rateMatrix[54][56] = nucleoFreqs[3] * kappa; /*TGC->TGT*/ rateMatrix[56][54] = nucleoFreqs[1] * kappa; /*TGT->TGC*/
        rateMatrix[57][59] = nucleoFreqs[2] * kappa; /*TTA->TTG*/ rateMatrix[59][57] = nucleoFreqs[0] * kappa; /*TTG->TTA*/
        rateMatrix[58][60] = nucleoFreqs[3] * kappa; /*TTC->TTT*/ rateMatrix[60][58] = nucleoFreqs[1] * kappa; /*TTT->TTC*/
    	//non-syn, transversion
        rateMatrix[0][1] = nucleoFreqs[1] * omega; /*AAA->AAC*/ rateMatrix[1][0] = nucleoFreqs[0] * omega; /*AAC->AAA*/
        rateMatrix[0][3] = nucleoFreqs[3] * omega; /*AAA->AAT*/ rateMatrix[3][0] = nucleoFreqs[0] * omega; /*AAT->AAA*/
        rateMatrix[0][4] = nucleoFreqs[1] * omega; /*AAA->ACA*/ rateMatrix[4][0] = nucleoFreqs[0] * omega; /*ACA->AAA*/
        rateMatrix[0][12] = nucleoFreqs[3] * omega; /*AAA->ATA*/ rateMatrix[12][0] = nucleoFreqs[0] * omega; /*ATA->AAA*/
        rateMatrix[0][16] = nucleoFreqs[1] * omega; /*AAA->CAA*/ rateMatrix[16][0] = nucleoFreqs[0] * omega; /*CAA->AAA*/
        rateMatrix[1][2] = nucleoFreqs[2] * omega; /*AAC->AAG*/ rateMatrix[2][1] = nucleoFreqs[1] * omega; /*AAG->AAC*/
        rateMatrix[1][5] = nucleoFreqs[1] * omega; /*AAC->ACC*/ rateMatrix[5][1] = nucleoFreqs[0] * omega; /*ACC->AAC*/
        rateMatrix[1][13] = nucleoFreqs[3] * omega; /*AAC->ATC*/ rateMatrix[13][1] = nucleoFreqs[0] * omega; /*ATC->AAC*/
        rateMatrix[1][17] = nucleoFreqs[1] * omega; /*AAC->CAC*/ rateMatrix[17][1] = nucleoFreqs[0] * omega; /*CAC->AAC*/
        rateMatrix[1][48] = nucleoFreqs[3] * omega; /*AAC->TAC*/ rateMatrix[48][1] = nucleoFreqs[0] * omega; /*TAC->AAC*/
        rateMatrix[2][3] = nucleoFreqs[3] * omega; /*AAG->AAT*/ rateMatrix[3][2] = nucleoFreqs[2] * omega; /*AAT->AAG*/
        rateMatrix[2][6] = nucleoFreqs[1] * omega; /*AAG->ACG*/ rateMatrix[6][2] = nucleoFreqs[0] * omega; /*ACG->AAG*/
        rateMatrix[2][14] = nucleoFreqs[3] * omega; /*AAG->ATG*/ rateMatrix[14][2] = nucleoFreqs[0] * omega; /*ATG->AAG*/
        rateMatrix[2][18] = nucleoFreqs[1] * omega; /*AAG->CAG*/ rateMatrix[18][2] = nucleoFreqs[0] * omega; /*CAG->AAG*/
        rateMatrix[3][7] = nucleoFreqs[1] * omega; /*AAT->ACT*/ rateMatrix[7][3] = nucleoFreqs[0] * omega; /*ACT->AAT*/
        rateMatrix[3][15] = nucleoFreqs[3] * omega; /*AAT->ATT*/ rateMatrix[15][3] = nucleoFreqs[0] * omega; /*ATT->AAT*/
        rateMatrix[3][19] = nucleoFreqs[1] * omega; /*AAT->CAT*/ rateMatrix[19][3] = nucleoFreqs[0] * omega; /*CAT->AAT*/
        rateMatrix[3][49] = nucleoFreqs[3] * omega; /*AAT->TAT*/ rateMatrix[49][3] = nucleoFreqs[0] * omega; /*TAT->AAT*/
        rateMatrix[4][8] = nucleoFreqs[2] * omega; /*ACA->AGA*/ rateMatrix[8][4] = nucleoFreqs[1] * omega; /*AGA->ACA*/
        rateMatrix[4][20] = nucleoFreqs[1] * omega; /*ACA->CCA*/ rateMatrix[20][4] = nucleoFreqs[0] * omega; /*CCA->ACA*/
        rateMatrix[4][50] = nucleoFreqs[3] * omega; /*ACA->TCA*/ rateMatrix[50][4] = nucleoFreqs[0] * omega; /*TCA->ACA*/
        rateMatrix[5][9] = nucleoFreqs[2] * omega; /*ACC->AGC*/ rateMatrix[9][5] = nucleoFreqs[1] * omega; /*AGC->ACC*/
        rateMatrix[5][21] = nucleoFreqs[1] * omega; /*ACC->CCC*/ rateMatrix[21][5] = nucleoFreqs[0] * omega; /*CCC->ACC*/
        rateMatrix[5][51] = nucleoFreqs[3] * omega; /*ACC->TCC*/ rateMatrix[51][5] = nucleoFreqs[0] * omega; /*TCC->ACC*/
        rateMatrix[6][10] = nucleoFreqs[2] * omega; /*ACG->AGG*/ rateMatrix[10][6] = nucleoFreqs[1] * omega; /*AGG->ACG*/
        rateMatrix[6][22] = nucleoFreqs[1] * omega; /*ACG->CCG*/ rateMatrix[22][6] = nucleoFreqs[0] * omega; /*CCG->ACG*/
        rateMatrix[6][52] = nucleoFreqs[3] * omega; /*ACG->TCG*/ rateMatrix[52][6] = nucleoFreqs[0] * omega; /*TCG->ACG*/
        rateMatrix[7][11] = nucleoFreqs[2] * omega; /*ACT->AGT*/ rateMatrix[11][7] = nucleoFreqs[1] * omega; /*AGT->ACT*/
        rateMatrix[7][23] = nucleoFreqs[1] * omega; /*ACT->CCT*/ rateMatrix[23][7] = nucleoFreqs[0] * omega; /*CCT->ACT*/
        rateMatrix[7][53] = nucleoFreqs[3] * omega; /*ACT->TCT*/ rateMatrix[53][7] = nucleoFreqs[0] * omega; /*TCT->ACT*/
        rateMatrix[8][9] = nucleoFreqs[1] * omega; /*AGA->AGC*/ rateMatrix[9][8] = nucleoFreqs[0] * omega; /*AGC->AGA*/
        rateMatrix[8][11] = nucleoFreqs[3] * omega; /*AGA->AGT*/ rateMatrix[11][8] = nucleoFreqs[0] * omega; /*AGT->AGA*/
        rateMatrix[8][12] = nucleoFreqs[3] * omega; /*AGA->ATA*/ rateMatrix[12][8] = nucleoFreqs[2] * omega; /*ATA->AGA*/
        rateMatrix[9][10] = nucleoFreqs[2] * omega; /*AGC->AGG*/ rateMatrix[10][9] = nucleoFreqs[1] * omega; /*AGG->AGC*/
        rateMatrix[9][13] = nucleoFreqs[3] * omega; /*AGC->ATC*/ rateMatrix[13][9] = nucleoFreqs[2] * omega; /*ATC->AGC*/
        rateMatrix[9][25] = nucleoFreqs[1] * omega; /*AGC->CGC*/ rateMatrix[25][9] = nucleoFreqs[0] * omega; /*CGC->AGC*/
        rateMatrix[9][54] = nucleoFreqs[3] * omega; /*AGC->TGC*/ rateMatrix[54][9] = nucleoFreqs[0] * omega; /*TGC->AGC*/
        rateMatrix[10][11] = nucleoFreqs[3] * omega; /*AGG->AGT*/ rateMatrix[11][10] = nucleoFreqs[2] * omega; /*AGT->AGG*/
        rateMatrix[10][14] = nucleoFreqs[3] * omega; /*AGG->ATG*/ rateMatrix[14][10] = nucleoFreqs[2] * omega; /*ATG->AGG*/
        rateMatrix[10][55] = nucleoFreqs[3] * omega; /*AGG->TGG*/ rateMatrix[55][10] = nucleoFreqs[0] * omega; /*TGG->AGG*/
        rateMatrix[11][15] = nucleoFreqs[3] * omega; /*AGT->ATT*/ rateMatrix[15][11] = nucleoFreqs[2] * omega; /*ATT->AGT*/
        rateMatrix[11][27] = nucleoFreqs[1] * omega; /*AGT->CGT*/ rateMatrix[27][11] = nucleoFreqs[0] * omega; /*CGT->AGT*/
        rateMatrix[11][56] = nucleoFreqs[3] * omega; /*AGT->TGT*/ rateMatrix[56][11] = nucleoFreqs[0] * omega; /*TGT->AGT*/
        rateMatrix[12][28] = nucleoFreqs[1] * omega; /*ATA->CTA*/ rateMatrix[28][12] = nucleoFreqs[0] * omega; /*CTA->ATA*/
        rateMatrix[12][57] = nucleoFreqs[3] * omega; /*ATA->TTA*/ rateMatrix[57][12] = nucleoFreqs[0] * omega; /*TTA->ATA*/
        rateMatrix[13][14] = nucleoFreqs[2] * omega; /*ATC->ATG*/ rateMatrix[14][13] = nucleoFreqs[1] * omega; /*ATG->ATC*/
        rateMatrix[13][29] = nucleoFreqs[1] * omega; /*ATC->CTC*/ rateMatrix[29][13] = nucleoFreqs[0] * omega; /*CTC->ATC*/
        rateMatrix[13][58] = nucleoFreqs[3] * omega; /*ATC->TTC*/ rateMatrix[58][13] = nucleoFreqs[0] * omega; /*TTC->ATC*/
        rateMatrix[14][15] = nucleoFreqs[3] * omega; /*ATG->ATT*/ rateMatrix[15][14] = nucleoFreqs[2] * omega; /*ATT->ATG*/
        rateMatrix[14][30] = nucleoFreqs[1] * omega; /*ATG->CTG*/ rateMatrix[30][14] = nucleoFreqs[0] * omega; /*CTG->ATG*/
        rateMatrix[14][59] = nucleoFreqs[3] * omega; /*ATG->TTG*/ rateMatrix[59][14] = nucleoFreqs[0] * omega; /*TTG->ATG*/
        rateMatrix[15][31] = nucleoFreqs[1] * omega; /*ATT->CTT*/ rateMatrix[31][15] = nucleoFreqs[0] * omega; /*CTT->ATT*/
        rateMatrix[15][60] = nucleoFreqs[3] * omega; /*ATT->TTT*/ rateMatrix[60][15] = nucleoFreqs[0] * omega; /*TTT->ATT*/
        rateMatrix[16][17] = nucleoFreqs[1] * omega; /*CAA->CAC*/ rateMatrix[17][16] = nucleoFreqs[0] * omega; /*CAC->CAA*/
        rateMatrix[16][19] = nucleoFreqs[3] * omega; /*CAA->CAT*/ rateMatrix[19][16] = nucleoFreqs[0] * omega; /*CAT->CAA*/
        rateMatrix[16][20] = nucleoFreqs[1] * omega; /*CAA->CCA*/ rateMatrix[20][16] = nucleoFreqs[0] * omega; /*CCA->CAA*/
        rateMatrix[16][28] = nucleoFreqs[3] * omega; /*CAA->CTA*/ rateMatrix[28][16] = nucleoFreqs[0] * omega; /*CTA->CAA*/
        rateMatrix[16][32] = nucleoFreqs[2] * omega; /*CAA->GAA*/ rateMatrix[32][16] = nucleoFreqs[1] * omega; /*GAA->CAA*/
        rateMatrix[17][18] = nucleoFreqs[2] * omega; /*CAC->CAG*/ rateMatrix[18][17] = nucleoFreqs[1] * omega; /*CAG->CAC*/
        rateMatrix[17][21] = nucleoFreqs[1] * omega; /*CAC->CCC*/ rateMatrix[21][17] = nucleoFreqs[0] * omega; /*CCC->CAC*/
        rateMatrix[17][29] = nucleoFreqs[3] * omega; /*CAC->CTC*/ rateMatrix[29][17] = nucleoFreqs[0] * omega; /*CTC->CAC*/
        rateMatrix[17][33] = nucleoFreqs[2] * omega; /*CAC->GAC*/ rateMatrix[33][17] = nucleoFreqs[1] * omega; /*GAC->CAC*/
        rateMatrix[18][19] = nucleoFreqs[3] * omega; /*CAG->CAT*/ rateMatrix[19][18] = nucleoFreqs[2] * omega; /*CAT->CAG*/
        rateMatrix[18][22] = nucleoFreqs[1] * omega; /*CAG->CCG*/ rateMatrix[22][18] = nucleoFreqs[0] * omega; /*CCG->CAG*/
        rateMatrix[18][30] = nucleoFreqs[3] * omega; /*CAG->CTG*/ rateMatrix[30][18] = nucleoFreqs[0] * omega; /*CTG->CAG*/
        rateMatrix[18][34] = nucleoFreqs[2] * omega; /*CAG->GAG*/ rateMatrix[34][18] = nucleoFreqs[1] * omega; /*GAG->CAG*/
        rateMatrix[19][23] = nucleoFreqs[1] * omega; /*CAT->CCT*/ rateMatrix[23][19] = nucleoFreqs[0] * omega; /*CCT->CAT*/
        rateMatrix[19][31] = nucleoFreqs[3] * omega; /*CAT->CTT*/ rateMatrix[31][19] = nucleoFreqs[0] * omega; /*CTT->CAT*/
        rateMatrix[19][35] = nucleoFreqs[2] * omega; /*CAT->GAT*/ rateMatrix[35][19] = nucleoFreqs[1] * omega; /*GAT->CAT*/
        rateMatrix[20][24] = nucleoFreqs[2] * omega; /*CCA->CGA*/ rateMatrix[24][20] = nucleoFreqs[1] * omega; /*CGA->CCA*/
        rateMatrix[20][36] = nucleoFreqs[2] * omega; /*CCA->GCA*/ rateMatrix[36][20] = nucleoFreqs[1] * omega; /*GCA->CCA*/
        rateMatrix[21][25] = nucleoFreqs[2] * omega; /*CCC->CGC*/ rateMatrix[25][21] = nucleoFreqs[1] * omega; /*CGC->CCC*/
        rateMatrix[21][37] = nucleoFreqs[2] * omega; /*CCC->GCC*/ rateMatrix[37][21] = nucleoFreqs[1] * omega; /*GCC->CCC*/
        rateMatrix[22][26] = nucleoFreqs[2] * omega; /*CCG->CGG*/ rateMatrix[26][22] = nucleoFreqs[1] * omega; /*CGG->CCG*/
        rateMatrix[22][38] = nucleoFreqs[2] * omega; /*CCG->GCG*/ rateMatrix[38][22] = nucleoFreqs[1] * omega; /*GCG->CCG*/
        rateMatrix[23][27] = nucleoFreqs[2] * omega; /*CCT->CGT*/ rateMatrix[27][23] = nucleoFreqs[1] * omega; /*CGT->CCT*/
        rateMatrix[23][39] = nucleoFreqs[2] * omega; /*CCT->GCT*/ rateMatrix[39][23] = nucleoFreqs[1] * omega; /*GCT->CCT*/
        rateMatrix[24][28] = nucleoFreqs[3] * omega; /*CGA->CTA*/ rateMatrix[28][24] = nucleoFreqs[2] * omega; /*CTA->CGA*/
        rateMatrix[24][40] = nucleoFreqs[2] * omega; /*CGA->GGA*/ rateMatrix[40][24] = nucleoFreqs[1] * omega; /*GGA->CGA*/
        rateMatrix[25][29] = nucleoFreqs[3] * omega; /*CGC->CTC*/ rateMatrix[29][25] = nucleoFreqs[2] * omega; /*CTC->CGC*/
        rateMatrix[25][41] = nucleoFreqs[2] * omega; /*CGC->GGC*/ rateMatrix[41][25] = nucleoFreqs[1] * omega; /*GGC->CGC*/
        rateMatrix[26][30] = nucleoFreqs[3] * omega; /*CGG->CTG*/ rateMatrix[30][26] = nucleoFreqs[2] * omega; /*CTG->CGG*/
        rateMatrix[26][42] = nucleoFreqs[2] * omega; /*CGG->GGG*/ rateMatrix[42][26] = nucleoFreqs[1] * omega; /*GGG->CGG*/
        rateMatrix[27][31] = nucleoFreqs[3] * omega; /*CGT->CTT*/ rateMatrix[31][27] = nucleoFreqs[2] * omega; /*CTT->CGT*/
        rateMatrix[27][43] = nucleoFreqs[2] * omega; /*CGT->GGT*/ rateMatrix[43][27] = nucleoFreqs[1] * omega; /*GGT->CGT*/
        rateMatrix[28][44] = nucleoFreqs[2] * omega; /*CTA->GTA*/ rateMatrix[44][28] = nucleoFreqs[1] * omega; /*GTA->CTA*/
        rateMatrix[29][45] = nucleoFreqs[2] * omega; /*CTC->GTC*/ rateMatrix[45][29] = nucleoFreqs[1] * omega; /*GTC->CTC*/
        rateMatrix[30][46] = nucleoFreqs[2] * omega; /*CTG->GTG*/ rateMatrix[46][30] = nucleoFreqs[1] * omega; /*GTG->CTG*/
        rateMatrix[31][47] = nucleoFreqs[2] * omega; /*CTT->GTT*/ rateMatrix[47][31] = nucleoFreqs[1] * omega; /*GTT->CTT*/
        rateMatrix[32][33] = nucleoFreqs[1] * omega; /*GAA->GAC*/ rateMatrix[33][32] = nucleoFreqs[0] * omega; /*GAC->GAA*/
        rateMatrix[32][35] = nucleoFreqs[3] * omega; /*GAA->GAT*/ rateMatrix[35][32] = nucleoFreqs[0] * omega; /*GAT->GAA*/
        rateMatrix[32][36] = nucleoFreqs[1] * omega; /*GAA->GCA*/ rateMatrix[36][32] = nucleoFreqs[0] * omega; /*GCA->GAA*/
        rateMatrix[32][44] = nucleoFreqs[3] * omega; /*GAA->GTA*/ rateMatrix[44][32] = nucleoFreqs[0] * omega; /*GTA->GAA*/
        rateMatrix[33][34] = nucleoFreqs[2] * omega; /*GAC->GAG*/ rateMatrix[34][33] = nucleoFreqs[1] * omega; /*GAG->GAC*/
        rateMatrix[33][37] = nucleoFreqs[1] * omega; /*GAC->GCC*/ rateMatrix[37][33] = nucleoFreqs[0] * omega; /*GCC->GAC*/
        rateMatrix[33][45] = nucleoFreqs[3] * omega; /*GAC->GTC*/ rateMatrix[45][33] = nucleoFreqs[0] * omega; /*GTC->GAC*/
        rateMatrix[33][48] = nucleoFreqs[3] * omega; /*GAC->TAC*/ rateMatrix[48][33] = nucleoFreqs[2] * omega; /*TAC->GAC*/
        rateMatrix[34][35] = nucleoFreqs[3] * omega; /*GAG->GAT*/ rateMatrix[35][34] = nucleoFreqs[2] * omega; /*GAT->GAG*/
        rateMatrix[34][38] = nucleoFreqs[1] * omega; /*GAG->GCG*/ rateMatrix[38][34] = nucleoFreqs[0] * omega; /*GCG->GAG*/
        rateMatrix[34][46] = nucleoFreqs[3] * omega; /*GAG->GTG*/ rateMatrix[46][34] = nucleoFreqs[0] * omega; /*GTG->GAG*/
        rateMatrix[35][39] = nucleoFreqs[1] * omega; /*GAT->GCT*/ rateMatrix[39][35] = nucleoFreqs[0] * omega; /*GCT->GAT*/
        rateMatrix[35][47] = nucleoFreqs[3] * omega; /*GAT->GTT*/ rateMatrix[47][35] = nucleoFreqs[0] * omega; /*GTT->GAT*/
        rateMatrix[35][49] = nucleoFreqs[3] * omega; /*GAT->TAT*/ rateMatrix[49][35] = nucleoFreqs[2] * omega; /*TAT->GAT*/
        rateMatrix[36][40] = nucleoFreqs[2] * omega; /*GCA->GGA*/ rateMatrix[40][36] = nucleoFreqs[1] * omega; /*GGA->GCA*/
        rateMatrix[36][50] = nucleoFreqs[3] * omega; /*GCA->TCA*/ rateMatrix[50][36] = nucleoFreqs[2] * omega; /*TCA->GCA*/
        rateMatrix[37][41] = nucleoFreqs[2] * omega; /*GCC->GGC*/ rateMatrix[41][37] = nucleoFreqs[1] * omega; /*GGC->GCC*/
        rateMatrix[37][51] = nucleoFreqs[3] * omega; /*GCC->TCC*/ rateMatrix[51][37] = nucleoFreqs[2] * omega; /*TCC->GCC*/
        rateMatrix[38][42] = nucleoFreqs[2] * omega; /*GCG->GGG*/ rateMatrix[42][38] = nucleoFreqs[1] * omega; /*GGG->GCG*/
        rateMatrix[38][52] = nucleoFreqs[3] * omega; /*GCG->TCG*/ rateMatrix[52][38] = nucleoFreqs[2] * omega; /*TCG->GCG*/
        rateMatrix[39][43] = nucleoFreqs[2] * omega; /*GCT->GGT*/ rateMatrix[43][39] = nucleoFreqs[1] * omega; /*GGT->GCT*/
        rateMatrix[39][53] = nucleoFreqs[3] * omega; /*GCT->TCT*/ rateMatrix[53][39] = nucleoFreqs[2] * omega; /*TCT->GCT*/
        rateMatrix[40][44] = nucleoFreqs[3] * omega; /*GGA->GTA*/ rateMatrix[44][40] = nucleoFreqs[2] * omega; /*GTA->GGA*/
        rateMatrix[41][45] = nucleoFreqs[3] * omega; /*GGC->GTC*/ rateMatrix[45][41] = nucleoFreqs[2] * omega; /*GTC->GGC*/
        rateMatrix[41][54] = nucleoFreqs[3] * omega; /*GGC->TGC*/ rateMatrix[54][41] = nucleoFreqs[2] * omega; /*TGC->GGC*/
        rateMatrix[42][46] = nucleoFreqs[3] * omega; /*GGG->GTG*/ rateMatrix[46][42] = nucleoFreqs[2] * omega; /*GTG->GGG*/
        rateMatrix[42][55] = nucleoFreqs[3] * omega; /*GGG->TGG*/ rateMatrix[55][42] = nucleoFreqs[2] * omega; /*TGG->GGG*/
        rateMatrix[43][47] = nucleoFreqs[3] * omega; /*GGT->GTT*/ rateMatrix[47][43] = nucleoFreqs[2] * omega; /*GTT->GGT*/
        rateMatrix[43][56] = nucleoFreqs[3] * omega; /*GGT->TGT*/ rateMatrix[56][43] = nucleoFreqs[2] * omega; /*TGT->GGT*/
        rateMatrix[44][57] = nucleoFreqs[3] * omega; /*GTA->TTA*/ rateMatrix[57][44] = nucleoFreqs[2] * omega; /*TTA->GTA*/
        rateMatrix[45][58] = nucleoFreqs[3] * omega; /*GTC->TTC*/ rateMatrix[58][45] = nucleoFreqs[2] * omega; /*TTC->GTC*/
        rateMatrix[46][59] = nucleoFreqs[3] * omega; /*GTG->TTG*/ rateMatrix[59][46] = nucleoFreqs[2] * omega; /*TTG->GTG*/
        rateMatrix[47][60] = nucleoFreqs[3] * omega; /*GTT->TTT*/ rateMatrix[60][47] = nucleoFreqs[2] * omega; /*TTT->GTT*/
        rateMatrix[48][51] = nucleoFreqs[1] * omega; /*TAC->TCC*/ rateMatrix[51][48] = nucleoFreqs[0] * omega; /*TCC->TAC*/
        rateMatrix[48][58] = nucleoFreqs[3] * omega; /*TAC->TTC*/ rateMatrix[58][48] = nucleoFreqs[0] * omega; /*TTC->TAC*/
        rateMatrix[49][53] = nucleoFreqs[1] * omega; /*TAT->TCT*/ rateMatrix[53][49] = nucleoFreqs[0] * omega; /*TCT->TAT*/
        rateMatrix[49][60] = nucleoFreqs[3] * omega; /*TAT->TTT*/ rateMatrix[60][49] = nucleoFreqs[0] * omega; /*TTT->TAT*/
        rateMatrix[51][54] = nucleoFreqs[2] * omega; /*TCC->TGC*/ rateMatrix[54][51] = nucleoFreqs[1] * omega; /*TGC->TCC*/
        rateMatrix[52][55] = nucleoFreqs[2] * omega; /*TCG->TGG*/ rateMatrix[55][52] = nucleoFreqs[1] * omega; /*TGG->TCG*/
        rateMatrix[53][56] = nucleoFreqs[2] * omega; /*TCT->TGT*/ rateMatrix[56][53] = nucleoFreqs[1] * omega; /*TGT->TCT*/
        rateMatrix[54][55] = nucleoFreqs[2] * omega; /*TGC->TGG*/ rateMatrix[55][54] = nucleoFreqs[1] * omega; /*TGG->TGC*/
        rateMatrix[54][58] = nucleoFreqs[3] * omega; /*TGC->TTC*/ rateMatrix[58][54] = nucleoFreqs[2] * omega; /*TTC->TGC*/
        rateMatrix[55][56] = nucleoFreqs[3] * omega; /*TGG->TGT*/ rateMatrix[56][55] = nucleoFreqs[2] * omega; /*TGT->TGG*/
        rateMatrix[55][59] = nucleoFreqs[3] * omega; /*TGG->TTG*/ rateMatrix[59][55] = nucleoFreqs[2] * omega; /*TTG->TGG*/
        rateMatrix[56][60] = nucleoFreqs[3] * omega; /*TGT->TTT*/ rateMatrix[60][56] = nucleoFreqs[2] * omega; /*TTT->TGT*/
        rateMatrix[57][58] = nucleoFreqs[1] * omega; /*TTA->TTC*/ rateMatrix[58][57] = nucleoFreqs[0] * omega; /*TTC->TTA*/
        rateMatrix[57][60] = nucleoFreqs[3] * omega; /*TTA->TTT*/ rateMatrix[60][57] = nucleoFreqs[0] * omega; /*TTT->TTA*/
        rateMatrix[58][59] = nucleoFreqs[2] * omega; /*TTC->TTG*/ rateMatrix[59][58] = nucleoFreqs[1] * omega; /*TTG->TTC*/
        rateMatrix[59][60] = nucleoFreqs[3] * omega; /*TTG->TTT*/ rateMatrix[60][59] = nucleoFreqs[2] * omega; /*TTT->TTG*/
    	//non-syn, transition
        rateMatrix[0][8] = nucleoFreqs[2] * kappa * omega; /*AAA->AGA*/ rateMatrix[8][0] = nucleoFreqs[0] * kappa * omega; /*AGA->AAA*/
        rateMatrix[0][32] = nucleoFreqs[2] * kappa * omega; /*AAA->GAA*/ rateMatrix[32][0] = nucleoFreqs[0] * kappa * omega; /*GAA->AAA*/
        rateMatrix[1][9] = nucleoFreqs[2] * kappa * omega; /*AAC->AGC*/ rateMatrix[9][1] = nucleoFreqs[0] * kappa * omega; /*AGC->AAC*/
        rateMatrix[1][33] = nucleoFreqs[2] * kappa * omega; /*AAC->GAC*/ rateMatrix[33][1] = nucleoFreqs[0] * kappa * omega; /*GAC->AAC*/
        rateMatrix[2][10] = nucleoFreqs[2] * kappa * omega; /*AAG->AGG*/ rateMatrix[10][2] = nucleoFreqs[0] * kappa * omega; /*AGG->AAG*/
        rateMatrix[2][34] = nucleoFreqs[2] * kappa * omega; /*AAG->GAG*/ rateMatrix[34][2] = nucleoFreqs[0] * kappa * omega; /*GAG->AAG*/
        rateMatrix[3][11] = nucleoFreqs[2] * kappa * omega; /*AAT->AGT*/ rateMatrix[11][3] = nucleoFreqs[0] * kappa * omega; /*AGT->AAT*/
        rateMatrix[3][35] = nucleoFreqs[2] * kappa * omega; /*AAT->GAT*/ rateMatrix[35][3] = nucleoFreqs[0] * kappa * omega; /*GAT->AAT*/
        rateMatrix[4][12] = nucleoFreqs[3] * kappa * omega; /*ACA->ATA*/ rateMatrix[12][4] = nucleoFreqs[1] * kappa * omega; /*ATA->ACA*/
        rateMatrix[4][36] = nucleoFreqs[2] * kappa * omega; /*ACA->GCA*/ rateMatrix[36][4] = nucleoFreqs[0] * kappa * omega; /*GCA->ACA*/
        rateMatrix[5][13] = nucleoFreqs[3] * kappa * omega; /*ACC->ATC*/ rateMatrix[13][5] = nucleoFreqs[1] * kappa * omega; /*ATC->ACC*/
        rateMatrix[5][37] = nucleoFreqs[2] * kappa * omega; /*ACC->GCC*/ rateMatrix[37][5] = nucleoFreqs[0] * kappa * omega; /*GCC->ACC*/
        rateMatrix[6][14] = nucleoFreqs[3] * kappa * omega; /*ACG->ATG*/ rateMatrix[14][6] = nucleoFreqs[1] * kappa * omega; /*ATG->ACG*/
        rateMatrix[6][38] = nucleoFreqs[2] * kappa * omega; /*ACG->GCG*/ rateMatrix[38][6] = nucleoFreqs[0] * kappa * omega; /*GCG->ACG*/
        rateMatrix[7][15] = nucleoFreqs[3] * kappa * omega; /*ACT->ATT*/ rateMatrix[15][7] = nucleoFreqs[1] * kappa * omega; /*ATT->ACT*/
        rateMatrix[7][39] = nucleoFreqs[2] * kappa * omega; /*ACT->GCT*/ rateMatrix[39][7] = nucleoFreqs[0] * kappa * omega; /*GCT->ACT*/
        rateMatrix[8][40] = nucleoFreqs[2] * kappa * omega; /*AGA->GGA*/ rateMatrix[40][8] = nucleoFreqs[0] * kappa * omega; /*GGA->AGA*/
        rateMatrix[9][41] = nucleoFreqs[2] * kappa * omega; /*AGC->GGC*/ rateMatrix[41][9] = nucleoFreqs[0] * kappa * omega; /*GGC->AGC*/
        rateMatrix[10][42] = nucleoFreqs[2] * kappa * omega; /*AGG->GGG*/ rateMatrix[42][10] = nucleoFreqs[0] * kappa * omega; /*GGG->AGG*/
        rateMatrix[11][43] = nucleoFreqs[2] * kappa * omega; /*AGT->GGT*/ rateMatrix[43][11] = nucleoFreqs[0] * kappa * omega; /*GGT->AGT*/
        rateMatrix[12][14] = nucleoFreqs[2] * kappa * omega; /*ATA->ATG*/ rateMatrix[14][12] = nucleoFreqs[0] * kappa * omega; /*ATG->ATA*/
        rateMatrix[12][44] = nucleoFreqs[2] * kappa * omega; /*ATA->GTA*/ rateMatrix[44][12] = nucleoFreqs[0] * kappa * omega; /*GTA->ATA*/
        rateMatrix[13][45] = nucleoFreqs[2] * kappa * omega; /*ATC->GTC*/ rateMatrix[45][13] = nucleoFreqs[0] * kappa * omega; /*GTC->ATC*/
        rateMatrix[14][46] = nucleoFreqs[2] * kappa * omega; /*ATG->GTG*/ rateMatrix[46][14] = nucleoFreqs[0] * kappa * omega; /*GTG->ATG*/
        rateMatrix[15][47] = nucleoFreqs[2] * kappa * omega; /*ATT->GTT*/ rateMatrix[47][15] = nucleoFreqs[0] * kappa * omega; /*GTT->ATT*/
        rateMatrix[16][24] = nucleoFreqs[2] * kappa * omega; /*CAA->CGA*/ rateMatrix[24][16] = nucleoFreqs[0] * kappa * omega; /*CGA->CAA*/
        rateMatrix[17][25] = nucleoFreqs[2] * kappa * omega; /*CAC->CGC*/ rateMatrix[25][17] = nucleoFreqs[0] * kappa * omega; /*CGC->CAC*/
        rateMatrix[17][48] = nucleoFreqs[3] * kappa * omega; /*CAC->TAC*/ rateMatrix[48][17] = nucleoFreqs[1] * kappa * omega; /*TAC->CAC*/
        rateMatrix[18][26] = nucleoFreqs[2] * kappa * omega; /*CAG->CGG*/ rateMatrix[26][18] = nucleoFreqs[0] * kappa * omega; /*CGG->CAG*/
        rateMatrix[19][27] = nucleoFreqs[2] * kappa * omega; /*CAT->CGT*/ rateMatrix[27][19] = nucleoFreqs[0] * kappa * omega; /*CGT->CAT*/
        rateMatrix[19][49] = nucleoFreqs[3] * kappa * omega; /*CAT->TAT*/ rateMatrix[49][19] = nucleoFreqs[1] * kappa * omega; /*TAT->CAT*/
        rateMatrix[20][28] = nucleoFreqs[3] * kappa * omega; /*CCA->CTA*/ rateMatrix[28][20] = nucleoFreqs[1] * kappa * omega; /*CTA->CCA*/
        rateMatrix[20][50] = nucleoFreqs[3] * kappa * omega; /*CCA->TCA*/ rateMatrix[50][20] = nucleoFreqs[1] * kappa * omega; /*TCA->CCA*/
        rateMatrix[21][29] = nucleoFreqs[3] * kappa * omega; /*CCC->CTC*/ rateMatrix[29][21] = nucleoFreqs[1] * kappa * omega; /*CTC->CCC*/
        rateMatrix[21][51] = nucleoFreqs[3] * kappa * omega; /*CCC->TCC*/ rateMatrix[51][21] = nucleoFreqs[1] * kappa * omega; /*TCC->CCC*/
        rateMatrix[22][30] = nucleoFreqs[3] * kappa * omega; /*CCG->CTG*/ rateMatrix[30][22] = nucleoFreqs[1] * kappa * omega; /*CTG->CCG*/
        rateMatrix[22][52] = nucleoFreqs[3] * kappa * omega; /*CCG->TCG*/ rateMatrix[52][22] = nucleoFreqs[1] * kappa * omega; /*TCG->CCG*/
        rateMatrix[23][31] = nucleoFreqs[3] * kappa * omega; /*CCT->CTT*/ rateMatrix[31][23] = nucleoFreqs[1] * kappa * omega; /*CTT->CCT*/
        rateMatrix[23][53] = nucleoFreqs[3] * kappa * omega; /*CCT->TCT*/ rateMatrix[53][23] = nucleoFreqs[1] * kappa * omega; /*TCT->CCT*/
        rateMatrix[25][54] = nucleoFreqs[3] * kappa * omega; /*CGC->TGC*/ rateMatrix[54][25] = nucleoFreqs[1] * kappa * omega; /*TGC->CGC*/
        rateMatrix[26][55] = nucleoFreqs[3] * kappa * omega; /*CGG->TGG*/ rateMatrix[55][26] = nucleoFreqs[1] * kappa * omega; /*TGG->CGG*/
        rateMatrix[27][56] = nucleoFreqs[3] * kappa * omega; /*CGT->TGT*/ rateMatrix[56][27] = nucleoFreqs[1] * kappa * omega; /*TGT->CGT*/
        rateMatrix[29][58] = nucleoFreqs[3] * kappa * omega; /*CTC->TTC*/ rateMatrix[58][29] = nucleoFreqs[1] * kappa * omega; /*TTC->CTC*/
        rateMatrix[31][60] = nucleoFreqs[3] * kappa * omega; /*CTT->TTT*/ rateMatrix[60][31] = nucleoFreqs[1] * kappa * omega; /*TTT->CTT*/
        rateMatrix[32][40] = nucleoFreqs[2] * kappa * omega; /*GAA->GGA*/ rateMatrix[40][32] = nucleoFreqs[0] * kappa * omega; /*GGA->GAA*/
        rateMatrix[33][41] = nucleoFreqs[2] * kappa * omega; /*GAC->GGC*/ rateMatrix[41][33] = nucleoFreqs[0] * kappa * omega; /*GGC->GAC*/
        rateMatrix[34][42] = nucleoFreqs[2] * kappa * omega; /*GAG->GGG*/ rateMatrix[42][34] = nucleoFreqs[0] * kappa * omega; /*GGG->GAG*/
        rateMatrix[35][43] = nucleoFreqs[2] * kappa * omega; /*GAT->GGT*/ rateMatrix[43][35] = nucleoFreqs[0] * kappa * omega; /*GGT->GAT*/
        rateMatrix[36][44] = nucleoFreqs[3] * kappa * omega; /*GCA->GTA*/ rateMatrix[44][36] = nucleoFreqs[1] * kappa * omega; /*GTA->GCA*/
        rateMatrix[37][45] = nucleoFreqs[3] * kappa * omega; /*GCC->GTC*/ rateMatrix[45][37] = nucleoFreqs[1] * kappa * omega; /*GTC->GCC*/
        rateMatrix[38][46] = nucleoFreqs[3] * kappa * omega; /*GCG->GTG*/ rateMatrix[46][38] = nucleoFreqs[1] * kappa * omega; /*GTG->GCG*/
        rateMatrix[39][47] = nucleoFreqs[3] * kappa * omega; /*GCT->GTT*/ rateMatrix[47][39] = nucleoFreqs[1] * kappa * omega; /*GTT->GCT*/
        rateMatrix[48][54] = nucleoFreqs[2] * kappa * omega; /*TAC->TGC*/ rateMatrix[54][48] = nucleoFreqs[0] * kappa * omega; /*TGC->TAC*/
        rateMatrix[49][56] = nucleoFreqs[2] * kappa * omega; /*TAT->TGT*/ rateMatrix[56][49] = nucleoFreqs[0] * kappa * omega; /*TGT->TAT*/
        rateMatrix[50][57] = nucleoFreqs[3] * kappa * omega; /*TCA->TTA*/ rateMatrix[57][50] = nucleoFreqs[1] * kappa * omega; /*TTA->TCA*/
        rateMatrix[51][58] = nucleoFreqs[3] * kappa * omega; /*TCC->TTC*/ rateMatrix[58][51] = nucleoFreqs[1] * kappa * omega; /*TTC->TCC*/
        rateMatrix[52][59] = nucleoFreqs[3] * kappa * omega; /*TCG->TTG*/ rateMatrix[59][52] = nucleoFreqs[1] * kappa * omega; /*TTG->TCG*/
        rateMatrix[53][60] = nucleoFreqs[3] * kappa * omega; /*TCT->TTT*/ rateMatrix[60][53] = nucleoFreqs[1] * kappa * omega; /*TTT->TCT*/
    
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
        double[] fFreqs = getFrequencies();
        
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * fFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
        
    }// setupRateMatrix
	
    @Override
    public double[] getFrequencies() {
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	double [] codonFreqs = new double[nrOfStates];
    	
    	//total probability without stop codons: TAA, TAG, TGA
    	double totalProb = 1.0 - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[0] - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[2] - nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[0];

    	codonFreqs[0] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // AAA
    	codonFreqs[1] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // AAC
    	codonFreqs[2] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // AAG
    	codonFreqs[3] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // AAT
    	codonFreqs[4] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // ACA
    	codonFreqs[5] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // ACC
    	codonFreqs[6] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // ACG
    	codonFreqs[7] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // ACT
    	codonFreqs[8] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // AGA
    	codonFreqs[9] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // AGC
    	codonFreqs[10] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // AGG
    	codonFreqs[11] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // AGT
    	codonFreqs[12] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // ATA
    	codonFreqs[13] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // ATC
    	codonFreqs[14] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // ATG
    	codonFreqs[15] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // ATT
    	codonFreqs[16] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // CAA
    	codonFreqs[17] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // CAC
    	codonFreqs[18] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // CAG
    	codonFreqs[19] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // CAT
    	codonFreqs[20] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // CCA
    	codonFreqs[21] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // CCC
    	codonFreqs[22] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // CCG
    	codonFreqs[23] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // CCT
    	codonFreqs[24] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // CGA
    	codonFreqs[25] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // CGC
    	codonFreqs[26] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // CGG
    	codonFreqs[27] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // CGT
    	codonFreqs[28] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // CTA
    	codonFreqs[29] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // CTC
    	codonFreqs[30] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // CTG
    	codonFreqs[31] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // CTT
    	codonFreqs[32] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // GAA
    	codonFreqs[33] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // GAC
    	codonFreqs[34] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // GAG
    	codonFreqs[35] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // GAT
    	codonFreqs[36] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // GCA
    	codonFreqs[37] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // GCC
    	codonFreqs[38] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // GCG
    	codonFreqs[39] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // GCT
    	codonFreqs[40] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // GGA
    	codonFreqs[41] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // GGC
    	codonFreqs[42] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // GGG
    	codonFreqs[43] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // GGT
    	codonFreqs[44] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // GTA
    	codonFreqs[45] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // GTC
    	codonFreqs[46] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // GTG
    	codonFreqs[47] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // GTT
    	codonFreqs[48] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // TAC
    	codonFreqs[49] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // TAT
    	codonFreqs[50] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // TCA
    	codonFreqs[51] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // TCC
    	codonFreqs[52] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // TCG
    	codonFreqs[53] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // TCT
    	codonFreqs[54] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // TGC
    	codonFreqs[55] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // TGG
    	codonFreqs[56] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // TGT
    	codonFreqs[57] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // TTA
    	codonFreqs[58] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // TTC
    	codonFreqs[59] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // TTG
    	codonFreqs[60] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // TTT
    	
    	return codonFreqs;
    }
    
    /**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        storedUpdateMatrix = updateMatrix;
        storedEigenDecomposition = eigenDecomposition.copy();
//        System.arraycopy(relativeRates, 0, storedRelativeRates, 0, relativeRates.length);

        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {

        updateMatrix = storedUpdateMatrix;

        // To restore all this stuff just swap the pointers...
//        double[] tmp1 = storedRelativeRates;
//        storedRelativeRates = relativeRates;
//        relativeRates = tmp1;

        EigenDecomposition tmp = storedEigenDecomposition;
        storedEigenDecomposition = eigenDecomposition;
        eigenDecomposition = tmp;
        super.restore();

    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrix = true;
        return true;
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
                setupRateMatrix();
                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) throws Exception {
        if (dataType instanceof Codon) {
            return true;
        }
        throw new Exception("Can only handle codon data");
    }

    /***************************************************************************************/
    //for testing purpose, need to be disabled after testing
    /***************************************************************************************/    
    /**
     * access to (copy of) rate matrix *
     */
    public double[][] getRateMatrix() {
        return rateMatrix.clone();
    }
    
    public void prepareMatricesForTest(){
        setupRateMatrix();
    }
}
