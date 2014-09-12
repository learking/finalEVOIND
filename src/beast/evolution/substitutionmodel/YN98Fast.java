package beast.evolution.substitutionmodel;

import java.util.Arrays;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.UserDataType;
import beast.evolution.tree.Node;

public class YN98Fast extends SubstitutionModel.Base{

    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "omega parameter in YN98 model", Validate.REQUIRED);
    //frequencies for pi_A, pi_C, pi_G and pi_T
    public Input<Frequencies> nucleoFreqInput =
            new Input<Frequencies>("nucleoFrequencies", "substitution model equilibrium state frequencies", Validate.REQUIRED);
    
    double[][] rateMatrix;
    DoubleMatrix doubleRateMatrix;
    
    Frequencies nucleoFrequencies;
    double[][] symmMatrix;
    double[] diagMatrix;
    
    protected boolean updateMatrix = true;
    private boolean storedUpdateMatrix = true;
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public YN98Fast() {
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
        
        symmMatrix = new double[nrOfStates][nrOfStates];
        diagMatrix = new double[nrOfStates];
        rateMatrix = new double[nrOfStates][nrOfStates];
        doubleRateMatrix = new DoubleMatrix();
    }
    
    /**
     * sets up rate matrix *
     */
    void setupSymmMatrix(){
    	
    	final double k = kappaInput.get().getValue();
    	final double omega = omegaInput.get().getValue();
    	
    	//0:A 1:C 2:G 3:T
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	//construct symmetric codon substitution matrix (61X61):
    	//0: diagonal positions, codons differ at more than one position
    	symmMatrix[0][1] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->AAC: IS NonSynonymous;
    	symmMatrix[0][2] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * k ; 	//AAA->AAG: IS Transition;
    	symmMatrix[0][3] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->AAT: IS NonSynonymous;
    	symmMatrix[0][4] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->ACA: IS NonSynonymous;
    	symmMatrix[0][8] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//AAA->AGA: IS Transition;IS NonSynonymous;
    	symmMatrix[0][12] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->ATA: IS NonSynonymous;
    	symmMatrix[0][16] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAA->CAA: IS NonSynonymous;
    	symmMatrix[0][32] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//AAA->GAA: IS Transition;IS NonSynonymous;
    	symmMatrix[1][2] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAC->AAG: IS NonSynonymous;
    	symmMatrix[1][3] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * k ; 	//AAC->AAT: IS Transition;
    	symmMatrix[1][5] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->ACC: IS NonSynonymous;
    	symmMatrix[1][9] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//AAC->AGC: IS Transition;IS NonSynonymous;
    	symmMatrix[1][13] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->ATC: IS NonSynonymous;
    	symmMatrix[1][17] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->CAC: IS NonSynonymous;
    	symmMatrix[1][33] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//AAC->GAC: IS Transition;IS NonSynonymous;
    	symmMatrix[1][48] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AAC->TAC: IS NonSynonymous;
    	symmMatrix[2][3] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AAG->AAT: IS NonSynonymous;
    	symmMatrix[2][6] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->ACG: IS NonSynonymous;
    	symmMatrix[2][10] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//AAG->AGG: IS Transition;IS NonSynonymous;
    	symmMatrix[2][14] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->ATG: IS NonSynonymous;
    	symmMatrix[2][18] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AAG->CAG: IS NonSynonymous;
    	symmMatrix[2][34] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//AAG->GAG: IS Transition;IS NonSynonymous;
    	symmMatrix[3][7] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->ACT: IS NonSynonymous;
    	symmMatrix[3][11] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//AAT->AGT: IS Transition;IS NonSynonymous;
    	symmMatrix[3][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->ATT: IS NonSynonymous;
    	symmMatrix[3][19] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->CAT: IS NonSynonymous;
    	symmMatrix[3][35] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//AAT->GAT: IS Transition;IS NonSynonymous;
    	symmMatrix[3][49] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AAT->TAT: IS NonSynonymous;
    	symmMatrix[4][5] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACA->ACC: 
    	symmMatrix[4][6] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k ; 	//ACA->ACG: IS Transition;
    	symmMatrix[4][7] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACA->ACT: 
    	symmMatrix[4][8] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//ACA->AGA: IS NonSynonymous;
    	symmMatrix[4][12] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * k * omega ; 	//ACA->ATA: IS Transition;IS NonSynonymous;
    	symmMatrix[4][20] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//ACA->CCA: IS NonSynonymous;
    	symmMatrix[4][36] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//ACA->GCA: IS Transition;IS NonSynonymous;
    	symmMatrix[4][50] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//ACA->TCA: IS NonSynonymous;
    	symmMatrix[5][6] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACC->ACG: 
    	symmMatrix[5][7] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k ; 	//ACC->ACT: IS Transition;
    	symmMatrix[5][9] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//ACC->AGC: IS NonSynonymous;
    	symmMatrix[5][13] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//ACC->ATC: IS Transition;IS NonSynonymous;
    	symmMatrix[5][21] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//ACC->CCC: IS NonSynonymous;
    	symmMatrix[5][37] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//ACC->GCC: IS Transition;IS NonSynonymous;
    	symmMatrix[5][51] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//ACC->TCC: IS NonSynonymous;
    	symmMatrix[6][7] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) ; 	//ACG->ACT: 
    	symmMatrix[6][10] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//ACG->AGG: IS NonSynonymous;
    	symmMatrix[6][14] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * k * omega ; 	//ACG->ATG: IS Transition;IS NonSynonymous;
    	symmMatrix[6][22] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//ACG->CCG: IS NonSynonymous;
    	symmMatrix[6][38] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//ACG->GCG: IS Transition;IS NonSynonymous;
    	symmMatrix[6][52] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//ACG->TCG: IS NonSynonymous;
    	symmMatrix[7][11] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ACT->AGT: IS NonSynonymous;
    	symmMatrix[7][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//ACT->ATT: IS Transition;IS NonSynonymous;
    	symmMatrix[7][23] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//ACT->CCT: IS NonSynonymous;
    	symmMatrix[7][39] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//ACT->GCT: IS Transition;IS NonSynonymous;
    	symmMatrix[7][53] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//ACT->TCT: IS NonSynonymous;
    	symmMatrix[8][9] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGA->AGC: IS NonSynonymous;
    	symmMatrix[8][10] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * k ; 	//AGA->AGG: IS Transition;
    	symmMatrix[8][11] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGA->AGT: IS NonSynonymous;
    	symmMatrix[8][12] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//AGA->ATA: IS NonSynonymous;
    	symmMatrix[8][24] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) ; 	//AGA->CGA: 
    	symmMatrix[8][40] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//AGA->GGA: IS Transition;IS NonSynonymous;
    	symmMatrix[9][10] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGC->AGG: IS NonSynonymous;
    	symmMatrix[9][11] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * k ; 	//AGC->AGT: IS Transition;
    	symmMatrix[9][13] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//AGC->ATC: IS NonSynonymous;
    	symmMatrix[9][25] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//AGC->CGC: IS NonSynonymous;
    	symmMatrix[9][41] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//AGC->GGC: IS Transition;IS NonSynonymous;
    	symmMatrix[9][54] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//AGC->TGC: IS NonSynonymous;
    	symmMatrix[10][11] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGG->AGT: IS NonSynonymous;
    	symmMatrix[10][14] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//AGG->ATG: IS NonSynonymous;
    	symmMatrix[10][26] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//AGG->CGG: 
    	symmMatrix[10][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//AGG->GGG: IS Transition;IS NonSynonymous;
    	symmMatrix[10][55] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//AGG->TGG: IS NonSynonymous;
    	symmMatrix[11][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//AGT->ATT: IS NonSynonymous;
    	symmMatrix[11][27] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//AGT->CGT: IS NonSynonymous;
    	symmMatrix[11][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//AGT->GGT: IS Transition;IS NonSynonymous;
    	symmMatrix[11][56] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//AGT->TGT: IS NonSynonymous;
    	symmMatrix[12][13] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) ; 	//ATA->ATC: 
    	symmMatrix[12][14] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//ATA->ATG: IS Transition;IS NonSynonymous;
    	symmMatrix[12][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) ; 	//ATA->ATT: 
    	symmMatrix[12][28] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//ATA->CTA: IS NonSynonymous;
    	symmMatrix[12][44] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * k * omega ; 	//ATA->GTA: IS Transition;IS NonSynonymous;
    	symmMatrix[12][57] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//ATA->TTA: IS NonSynonymous;
    	symmMatrix[13][14] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ATC->ATG: IS NonSynonymous;
    	symmMatrix[13][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k ; 	//ATC->ATT: IS Transition;
    	symmMatrix[13][29] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//ATC->CTC: IS NonSynonymous;
    	symmMatrix[13][45] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//ATC->GTC: IS Transition;IS NonSynonymous;
    	symmMatrix[13][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//ATC->TTC: IS NonSynonymous;
    	symmMatrix[14][15] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//ATG->ATT: IS NonSynonymous;
    	symmMatrix[14][30] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//ATG->CTG: IS NonSynonymous;
    	symmMatrix[14][46] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * k * omega ; 	//ATG->GTG: IS Transition;IS NonSynonymous;
    	symmMatrix[14][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//ATG->TTG: IS NonSynonymous;
    	symmMatrix[15][31] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//ATT->CTT: IS NonSynonymous;
    	symmMatrix[15][47] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//ATT->GTT: IS Transition;IS NonSynonymous;
    	symmMatrix[15][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//ATT->TTT: IS NonSynonymous;
    	symmMatrix[16][17] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CAC: IS NonSynonymous;
    	symmMatrix[16][18] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k ; 	//CAA->CAG: IS Transition;
    	symmMatrix[16][19] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CAT: IS NonSynonymous;
    	symmMatrix[16][20] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CCA: IS NonSynonymous;
    	symmMatrix[16][24] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CAA->CGA: IS Transition;IS NonSynonymous;
    	symmMatrix[16][28] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAA->CTA: IS NonSynonymous;
    	symmMatrix[16][32] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[0]) * omega ; 	//CAA->GAA: IS NonSynonymous;
    	symmMatrix[17][18] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAC->CAG: IS NonSynonymous;
    	symmMatrix[17][19] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k ; 	//CAC->CAT: IS Transition;
    	symmMatrix[17][21] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CAC->CCC: IS NonSynonymous;
    	symmMatrix[17][25] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CAC->CGC: IS Transition;IS NonSynonymous;
    	symmMatrix[17][29] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CAC->CTC: IS NonSynonymous;
    	symmMatrix[17][33] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//CAC->GAC: IS NonSynonymous;
    	symmMatrix[17][48] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * k * omega ; 	//CAC->TAC: IS Transition;IS NonSynonymous;
    	symmMatrix[18][19] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CAG->CAT: IS NonSynonymous;
    	symmMatrix[18][22] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CAG->CCG: IS NonSynonymous;
    	symmMatrix[18][26] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CAG->CGG: IS Transition;IS NonSynonymous;
    	symmMatrix[18][30] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CAG->CTG: IS NonSynonymous;
    	symmMatrix[18][34] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[2]) * omega ; 	//CAG->GAG: IS NonSynonymous;
    	symmMatrix[19][23] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CAT->CCT: IS NonSynonymous;
    	symmMatrix[19][27] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CAT->CGT: IS Transition;IS NonSynonymous;
    	symmMatrix[19][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CAT->CTT: IS NonSynonymous;
    	symmMatrix[19][35] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//CAT->GAT: IS NonSynonymous;
    	symmMatrix[19][49] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * k * omega ; 	//CAT->TAT: IS Transition;IS NonSynonymous;
    	symmMatrix[20][21] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCA->CCC: 
    	symmMatrix[20][22] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k ; 	//CCA->CCG: IS Transition;
    	symmMatrix[20][23] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCA->CCT: 
    	symmMatrix[20][24] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CCA->CGA: IS NonSynonymous;
    	symmMatrix[20][28] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CCA->CTA: IS Transition;IS NonSynonymous;
    	symmMatrix[20][36] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CCA->GCA: IS NonSynonymous;
    	symmMatrix[20][50] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * k * omega ; 	//CCA->TCA: IS Transition;IS NonSynonymous;
    	symmMatrix[21][22] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCC->CCG: 
    	symmMatrix[21][23] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k ; 	//CCC->CCT: IS Transition;
    	symmMatrix[21][25] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CCC->CGC: IS NonSynonymous;
    	symmMatrix[21][29] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CCC->CTC: IS Transition;IS NonSynonymous;
    	symmMatrix[21][37] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CCC->GCC: IS NonSynonymous;
    	symmMatrix[21][51] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * k * omega ; 	//CCC->TCC: IS Transition;IS NonSynonymous;
    	symmMatrix[22][23] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) ; 	//CCG->CCT: 
    	symmMatrix[22][26] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CCG->CGG: IS NonSynonymous;
    	symmMatrix[22][30] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CCG->CTG: IS Transition;IS NonSynonymous;
    	symmMatrix[22][38] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CCG->GCG: IS NonSynonymous;
    	symmMatrix[22][52] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k * omega ; 	//CCG->TCG: IS Transition;IS NonSynonymous;
    	symmMatrix[23][27] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CCT->CGT: IS NonSynonymous;
    	symmMatrix[23][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CCT->CTT: IS Transition;IS NonSynonymous;
    	symmMatrix[23][39] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CCT->GCT: IS NonSynonymous;
    	symmMatrix[23][53] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k * omega ; 	//CCT->TCT: IS Transition;IS NonSynonymous;
    	symmMatrix[24][25] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGA->CGC: 
    	symmMatrix[24][26] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k ; 	//CGA->CGG: IS Transition;
    	symmMatrix[24][27] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGA->CGT: 
    	symmMatrix[24][28] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//CGA->CTA: IS NonSynonymous;
    	symmMatrix[24][40] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//CGA->GGA: IS NonSynonymous;
    	symmMatrix[25][26] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGC->CGG: 
    	symmMatrix[25][27] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * k ; 	//CGC->CGT: IS Transition;
    	symmMatrix[25][29] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//CGC->CTC: IS NonSynonymous;
    	symmMatrix[25][41] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//CGC->GGC: IS NonSynonymous;
    	symmMatrix[25][54] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//CGC->TGC: IS Transition;IS NonSynonymous;
    	symmMatrix[26][27] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) ; 	//CGG->CGT: 
    	symmMatrix[26][30] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//CGG->CTG: IS NonSynonymous;
    	symmMatrix[26][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//CGG->GGG: IS NonSynonymous;
    	symmMatrix[26][55] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//CGG->TGG: IS Transition;IS NonSynonymous;
    	symmMatrix[27][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//CGT->CTT: IS NonSynonymous;
    	symmMatrix[27][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//CGT->GGT: IS NonSynonymous;
    	symmMatrix[27][56] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//CGT->TGT: IS Transition;IS NonSynonymous;
    	symmMatrix[28][29] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTA->CTC: 
    	symmMatrix[28][30] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k ; 	//CTA->CTG: IS Transition;
    	symmMatrix[28][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTA->CTT: 
    	symmMatrix[28][44] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//CTA->GTA: IS NonSynonymous;
    	symmMatrix[28][57] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * k ; 	//CTA->TTA: IS Transition;
    	symmMatrix[29][30] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTC->CTG: 
    	symmMatrix[29][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * k ; 	//CTC->CTT: IS Transition;
    	symmMatrix[29][45] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//CTC->GTC: IS NonSynonymous;
    	symmMatrix[29][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//CTC->TTC: IS Transition;IS NonSynonymous;
    	symmMatrix[30][31] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) ; 	//CTG->CTT: 
    	symmMatrix[30][46] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//CTG->GTG: IS NonSynonymous;
    	symmMatrix[30][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * k ; 	//CTG->TTG: IS Transition;
    	symmMatrix[31][47] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//CTT->GTT: IS NonSynonymous;
    	symmMatrix[31][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//CTT->TTT: IS Transition;IS NonSynonymous;
    	symmMatrix[32][33] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GAC: IS NonSynonymous;
    	symmMatrix[32][34] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * k ; 	//GAA->GAG: IS Transition;
    	symmMatrix[32][35] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GAT: IS NonSynonymous;
    	symmMatrix[32][36] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GCA: IS NonSynonymous;
    	symmMatrix[32][40] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//GAA->GGA: IS Transition;IS NonSynonymous;
    	symmMatrix[32][44] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAA->GTA: IS NonSynonymous;
    	symmMatrix[33][34] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAC->GAG: IS NonSynonymous;
    	symmMatrix[33][35] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * k ; 	//GAC->GAT: IS Transition;
    	symmMatrix[33][37] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GAC->GCC: IS NonSynonymous;
    	symmMatrix[33][41] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//GAC->GGC: IS Transition;IS NonSynonymous;
    	symmMatrix[33][45] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GAC->GTC: IS NonSynonymous;
    	symmMatrix[33][48] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[1]) * omega ; 	//GAC->TAC: IS NonSynonymous;
    	symmMatrix[34][35] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GAG->GAT: IS NonSynonymous;
    	symmMatrix[34][38] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GAG->GCG: IS NonSynonymous;
    	symmMatrix[34][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//GAG->GGG: IS Transition;IS NonSynonymous;
    	symmMatrix[34][46] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GAG->GTG: IS NonSynonymous;
    	symmMatrix[35][39] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GAT->GCT: IS NonSynonymous;
    	symmMatrix[35][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//GAT->GGT: IS Transition;IS NonSynonymous;
    	symmMatrix[35][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GAT->GTT: IS NonSynonymous;
    	symmMatrix[35][49] = 1.0 / (nucleoFreqs[0] * nucleoFreqs[3]) * omega ; 	//GAT->TAT: IS NonSynonymous;
    	symmMatrix[36][37] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCA->GCC: 
    	symmMatrix[36][38] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k ; 	//GCA->GCG: IS Transition;
    	symmMatrix[36][39] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCA->GCT: 
    	symmMatrix[36][40] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GCA->GGA: IS NonSynonymous;
    	symmMatrix[36][44] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * k * omega ; 	//GCA->GTA: IS Transition;IS NonSynonymous;
    	symmMatrix[36][50] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[0]) * omega ; 	//GCA->TCA: IS NonSynonymous;
    	symmMatrix[37][38] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCC->GCG: 
    	symmMatrix[37][39] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k ; 	//GCC->GCT: IS Transition;
    	symmMatrix[37][41] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GCC->GGC: IS NonSynonymous;
    	symmMatrix[37][45] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * k * omega ; 	//GCC->GTC: IS Transition;IS NonSynonymous;
    	symmMatrix[37][51] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[1]) * omega ; 	//GCC->TCC: IS NonSynonymous;
    	symmMatrix[38][39] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) ; 	//GCG->GCT: 
    	symmMatrix[38][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GCG->GGG: IS NonSynonymous;
    	symmMatrix[38][46] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k * omega ; 	//GCG->GTG: IS Transition;IS NonSynonymous;
    	symmMatrix[38][52] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[2]) * omega ; 	//GCG->TCG: IS NonSynonymous;
    	symmMatrix[39][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GCT->GGT: IS NonSynonymous;
    	symmMatrix[39][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k * omega ; 	//GCT->GTT: IS Transition;IS NonSynonymous;
    	symmMatrix[39][53] = 1.0 / (nucleoFreqs[1] * nucleoFreqs[3]) * omega ; 	//GCT->TCT: IS NonSynonymous;
    	symmMatrix[40][41] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGA->GGC: 
    	symmMatrix[40][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k ; 	//GGA->GGG: IS Transition;
    	symmMatrix[40][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGA->GGT: 
    	symmMatrix[40][44] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[0]) * omega ; 	//GGA->GTA: IS NonSynonymous;
    	symmMatrix[41][42] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGC->GGG: 
    	symmMatrix[41][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * k ; 	//GGC->GGT: IS Transition;
    	symmMatrix[41][45] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GGC->GTC: IS NonSynonymous;
    	symmMatrix[41][54] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[1]) * omega ; 	//GGC->TGC: IS NonSynonymous;
    	symmMatrix[42][43] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) ; 	//GGG->GGT: 
    	symmMatrix[42][46] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GGG->GTG: IS NonSynonymous;
    	symmMatrix[42][55] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[2]) * omega ; 	//GGG->TGG: IS NonSynonymous;
    	symmMatrix[43][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GGT->GTT: IS NonSynonymous;
    	symmMatrix[43][56] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * omega ; 	//GGT->TGT: IS NonSynonymous;
    	symmMatrix[44][45] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTA->GTC: 
    	symmMatrix[44][46] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k ; 	//GTA->GTG: IS Transition;
    	symmMatrix[44][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTA->GTT: 
    	symmMatrix[44][57] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * omega ; 	//GTA->TTA: IS NonSynonymous;
    	symmMatrix[45][46] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTC->GTG: 
    	symmMatrix[45][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) * k ; 	//GTC->GTT: IS Transition;
    	symmMatrix[45][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//GTC->TTC: IS NonSynonymous;
    	symmMatrix[46][47] = 1.0 / (nucleoFreqs[2] * nucleoFreqs[3]) ; 	//GTG->GTT: 
    	symmMatrix[46][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//GTG->TTG: IS NonSynonymous;
    	symmMatrix[47][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//GTT->TTT: IS NonSynonymous;
    	symmMatrix[48][49] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * k ; 	//TAC->TAT: IS Transition;
    	symmMatrix[48][51] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TAC->TCC: IS NonSynonymous;
    	symmMatrix[48][54] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//TAC->TGC: IS Transition;IS NonSynonymous;
    	symmMatrix[48][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TAC->TTC: IS NonSynonymous;
    	symmMatrix[49][53] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TAT->TCT: IS NonSynonymous;
    	symmMatrix[49][56] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//TAT->TGT: IS Transition;IS NonSynonymous;
    	symmMatrix[49][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TAT->TTT: IS NonSynonymous;
    	symmMatrix[50][51] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCA->TCC: 
    	symmMatrix[50][52] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k ; 	//TCA->TCG: IS Transition;
    	symmMatrix[50][53] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCA->TCT: 
    	symmMatrix[50][57] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[0]) * k * omega ; 	//TCA->TTA: IS Transition;IS NonSynonymous;
    	symmMatrix[51][52] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCC->TCG: 
    	symmMatrix[51][53] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k ; 	//TCC->TCT: IS Transition;
    	symmMatrix[51][54] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TCC->TGC: IS NonSynonymous;
    	symmMatrix[51][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * k * omega ; 	//TCC->TTC: IS Transition;IS NonSynonymous;
    	symmMatrix[52][53] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) ; 	//TCG->TCT: 
    	symmMatrix[52][55] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TCG->TGG: IS NonSynonymous;
    	symmMatrix[52][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * k * omega ; 	//TCG->TTG: IS Transition;IS NonSynonymous;
    	symmMatrix[53][56] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TCT->TGT: IS NonSynonymous;
    	symmMatrix[53][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k * omega ; 	//TCT->TTT: IS Transition;IS NonSynonymous;
    	symmMatrix[54][55] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGC->TGG: IS NonSynonymous;
    	symmMatrix[54][56] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * k ; 	//TGC->TGT: IS Transition;
    	symmMatrix[54][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[1]) * omega ; 	//TGC->TTC: IS NonSynonymous;
    	symmMatrix[55][56] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGG->TGT: IS NonSynonymous;
    	symmMatrix[55][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[2]) * omega ; 	//TGG->TTG: IS NonSynonymous;
    	symmMatrix[56][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TGT->TTT: IS NonSynonymous;
    	symmMatrix[57][58] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTA->TTC: IS NonSynonymous;
    	symmMatrix[57][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k ; 	//TTA->TTG: IS Transition;
    	symmMatrix[57][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTA->TTT: IS NonSynonymous;
    	symmMatrix[58][59] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTC->TTG: IS NonSynonymous;
    	symmMatrix[58][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * k ; 	//TTC->TTT: IS Transition;
    	symmMatrix[59][60] = 1.0 / (nucleoFreqs[3] * nucleoFreqs[3]) * omega ; 	//TTG->TTT: IS NonSynonymous;


    	//let lower triangle equal to upper
        for (int i = 0; i < nrOfStates; i++) {
        	for (int j = 0; j < nrOfStates; j++) {
    			if (i < j){
    				symmMatrix[j][i] = symmMatrix[i][j];
    			}
        	}
        }
    }

    void setupDiagMatrix(){
    	//0:A 1:C 2:G 3:T
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	
    	double totalProb = 1.0 - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[0] - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[2] - nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[0];
    	
    	//example:
    	diagMatrix[0] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // AAA
    	diagMatrix[1] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // AAC
    	diagMatrix[2] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // AAG
    	diagMatrix[3] = nucleoFreqs[0] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // AAT
    	diagMatrix[4] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // ACA
    	diagMatrix[5] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // ACC
    	diagMatrix[6] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // ACG
    	diagMatrix[7] = nucleoFreqs[0] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // ACT
    	diagMatrix[8] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // AGA
    	diagMatrix[9] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // AGC
    	diagMatrix[10] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // AGG
    	diagMatrix[11] = nucleoFreqs[0] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // AGT
    	diagMatrix[12] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // ATA
    	diagMatrix[13] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // ATC
    	diagMatrix[14] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // ATG
    	diagMatrix[15] = nucleoFreqs[0] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // ATT
    	diagMatrix[16] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // CAA
    	diagMatrix[17] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // CAC
    	diagMatrix[18] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // CAG
    	diagMatrix[19] = nucleoFreqs[1] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // CAT
    	diagMatrix[20] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // CCA
    	diagMatrix[21] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // CCC
    	diagMatrix[22] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // CCG
    	diagMatrix[23] = nucleoFreqs[1] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // CCT
    	diagMatrix[24] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // CGA
    	diagMatrix[25] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // CGC
    	diagMatrix[26] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // CGG
    	diagMatrix[27] = nucleoFreqs[1] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // CGT
    	diagMatrix[28] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // CTA
    	diagMatrix[29] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // CTC
    	diagMatrix[30] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // CTG
    	diagMatrix[31] = nucleoFreqs[1] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // CTT
    	diagMatrix[32] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[0] / totalProb; // GAA
    	diagMatrix[33] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // GAC
    	diagMatrix[34] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[2] / totalProb; // GAG
    	diagMatrix[35] = nucleoFreqs[2] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // GAT
    	diagMatrix[36] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // GCA
    	diagMatrix[37] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // GCC
    	diagMatrix[38] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // GCG
    	diagMatrix[39] = nucleoFreqs[2] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // GCT
    	diagMatrix[40] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[0] / totalProb; // GGA
    	diagMatrix[41] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // GGC
    	diagMatrix[42] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // GGG
    	diagMatrix[43] = nucleoFreqs[2] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // GGT
    	diagMatrix[44] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // GTA
    	diagMatrix[45] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // GTC
    	diagMatrix[46] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // GTG
    	diagMatrix[47] = nucleoFreqs[2] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // GTT
    	diagMatrix[48] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[1] / totalProb; // TAC
    	diagMatrix[49] = nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[3] / totalProb; // TAT
    	diagMatrix[50] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[0] / totalProb; // TCA
    	diagMatrix[51] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[1] / totalProb; // TCC
    	diagMatrix[52] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[2] / totalProb; // TCG
    	diagMatrix[53] = nucleoFreqs[3] * nucleoFreqs[1] * nucleoFreqs[3] / totalProb; // TCT
    	diagMatrix[54] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[1] / totalProb; // TGC
    	diagMatrix[55] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[2] / totalProb; // TGG
    	diagMatrix[56] = nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[3] / totalProb; // TGT
    	diagMatrix[57] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[0] / totalProb; // TTA
    	diagMatrix[58] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[1] / totalProb; // TTC
    	diagMatrix[59] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[2] / totalProb; // TTG
    	diagMatrix[60] = nucleoFreqs[3] * nucleoFreqs[3] * nucleoFreqs[3] / totalProb; // TTT
    }
    
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
    
    @Override
    public double[] getFrequencies() {
    	return diagMatrix;
    	/*
    	//this function should return stationary (should use up-to-date frequencies)
    	// in YN98 case: pi_j1 * pi_j2 * pi_j3 / (1 - P(stop codon 1) - P(stop codon 2) - P(stop codon 3))
    	double [] codonFreqs = new double[nrOfStates];
    	//can we compute codonFreqs from diagMatrix?
    	
    	//0:A 1:C 2:G 3:T
    	double[] nucleoFreqs = nucleoFrequencies.getFreqs();
    	
    	//total probability without stop codons
    	double totalProb = 1.0 - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[0] - nucleoFreqs[3] * nucleoFreqs[0] * nucleoFreqs[2] - nucleoFreqs[3] * nucleoFreqs[2] * nucleoFreqs[0];

    	//example:
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
        */
    }
    
	@Override
	public void getTransitionProbabilities(Node node, double fStartTime,
			double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;
        //System.out.println("distance:" + distance);

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupRateMatrix();
                doubleRateMatrix = new DoubleMatrix(rateMatrix);
                //System.out.println("doubleM:");
                //doubleRateMatrix.print();
                //System.out.println("rateM:" + Arrays.deepToString(rateMatrix));
                updateMatrix = false;
            }
        }

        //get exponentiated matrix and assign it to matrix
        //first: determine a reliable (and fast) method to do expm
        double[][] tmpMatrix = MatrixFunctions.expm(doubleRateMatrix.mul(distance)).toArray2();
        //System.out.println("tmpMatrix:" + Arrays.deepToString(tmpMatrix));	
        
        //copy values        
        int i,j;
        int u = 0;
        for (i = 0; i < nrOfStates; i++) {
            for (j = 0; j < nrOfStates; j++) {
                matrix[u] = tmpMatrix[i][j];
                u++;
            }
        }
  
        //System.out.println("probabilities:" + Arrays.toString(matrix));		
	}

	//should be useless method
	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

    /**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        storedUpdateMatrix = updateMatrix;
        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {
        updateMatrix = storedUpdateMatrix;
        super.restore();

    }
    
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrix = true;
        return true;
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
    public double[][] getRateMatrix() {
        return rateMatrix.clone();
    }
    
    public double[][] getSymmMatrix() {
        return symmMatrix.clone();
    }  
    
    public double[] getDiagMatrix() {
        return diagMatrix.clone();
    }
    public void prepareMatricesForTest(){
        setupRateMatrix();
    }
}
