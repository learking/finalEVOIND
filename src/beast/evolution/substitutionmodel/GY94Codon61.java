package beast.evolution.substitutionmodel;

import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class GY94Codon61  extends GeneralSubstitutionModel {
    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "omega parameter in YN98 model", Validate.REQUIRED);

    public Input<RealParameter> codonProbInput = new Input<RealParameter>("codonProb", "probabilities for each codon in OneStruct model", Validate.REQUIRED);    
    double[] codonProb;
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public GY94Codon61() {
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
        
        updateMatrix = true;
        nrOfStates = 61;

        eigenSystem = createEigenSystem();
        
        rateMatrix = new double[nrOfStates][nrOfStates];
    }
    
    @Override
    protected void setupRelativeRates() {}
    
    @Override
    protected void setupRateMatrix() {
    	final double k = kappaInput.get().getValue();
    	final double omega = omegaInput.get().getValue();
    	
    	double[] fFreqs = getFrequencies();
    	// fill in rateMatrix directly
    	rateMatrix[0][1] = fFreqs[1] * omega ; 	//AAA->AAC: IS NonSynonymous;
    	rateMatrix[0][2] = fFreqs[2] * k ; 	//AAA->AAG: IS Transition;
    	rateMatrix[0][3] = fFreqs[3] * omega ; 	//AAA->AAT: IS NonSynonymous;
    	rateMatrix[0][4] = fFreqs[4] * omega ; 	//AAA->ACA: IS NonSynonymous;
    	rateMatrix[0][8] = fFreqs[8] * k * omega ; 	//AAA->AGA: IS Transition;IS NonSynonymous;
    	rateMatrix[0][12] = fFreqs[12] * omega ; 	//AAA->ATA: IS NonSynonymous;
    	rateMatrix[0][16] = fFreqs[16] * omega ; 	//AAA->CAA: IS NonSynonymous;
    	rateMatrix[0][32] = fFreqs[32] * k * omega ; 	//AAA->GAA: IS Transition;IS NonSynonymous;
    	rateMatrix[1][0] = fFreqs[0] * omega ; 	//AAC->AAA: IS NonSynonymous;
    	rateMatrix[1][2] = fFreqs[2] * omega ; 	//AAC->AAG: IS NonSynonymous;
    	rateMatrix[1][3] = fFreqs[3] * k ; 	//AAC->AAT: IS Transition;
    	rateMatrix[1][5] = fFreqs[5] * omega ; 	//AAC->ACC: IS NonSynonymous;
    	rateMatrix[1][9] = fFreqs[9] * k * omega ; 	//AAC->AGC: IS Transition;IS NonSynonymous;
    	rateMatrix[1][13] = fFreqs[13] * omega ; 	//AAC->ATC: IS NonSynonymous;
    	rateMatrix[1][17] = fFreqs[17] * omega ; 	//AAC->CAC: IS NonSynonymous;
    	rateMatrix[1][33] = fFreqs[33] * k * omega ; 	//AAC->GAC: IS Transition;IS NonSynonymous;
    	rateMatrix[1][48] = fFreqs[48] * omega ; 	//AAC->TAC: IS NonSynonymous;
    	rateMatrix[2][0] = fFreqs[0] * k ; 	//AAG->AAA: IS Transition;
    	rateMatrix[2][1] = fFreqs[1] * omega ; 	//AAG->AAC: IS NonSynonymous;
    	rateMatrix[2][3] = fFreqs[3] * omega ; 	//AAG->AAT: IS NonSynonymous;
    	rateMatrix[2][6] = fFreqs[6] * omega ; 	//AAG->ACG: IS NonSynonymous;
    	rateMatrix[2][10] = fFreqs[10] * k * omega ; 	//AAG->AGG: IS Transition;IS NonSynonymous;
    	rateMatrix[2][14] = fFreqs[14] * omega ; 	//AAG->ATG: IS NonSynonymous;
    	rateMatrix[2][18] = fFreqs[18] * omega ; 	//AAG->CAG: IS NonSynonymous;
    	rateMatrix[2][34] = fFreqs[34] * k * omega ; 	//AAG->GAG: IS Transition;IS NonSynonymous;
    	rateMatrix[3][0] = fFreqs[0] * omega ; 	//AAT->AAA: IS NonSynonymous;
    	rateMatrix[3][1] = fFreqs[1] * k ; 	//AAT->AAC: IS Transition;
    	rateMatrix[3][2] = fFreqs[2] * omega ; 	//AAT->AAG: IS NonSynonymous;
    	rateMatrix[3][7] = fFreqs[7] * omega ; 	//AAT->ACT: IS NonSynonymous;
    	rateMatrix[3][11] = fFreqs[11] * k * omega ; 	//AAT->AGT: IS Transition;IS NonSynonymous;
    	rateMatrix[3][15] = fFreqs[15] * omega ; 	//AAT->ATT: IS NonSynonymous;
    	rateMatrix[3][19] = fFreqs[19] * omega ; 	//AAT->CAT: IS NonSynonymous;
    	rateMatrix[3][35] = fFreqs[35] * k * omega ; 	//AAT->GAT: IS Transition;IS NonSynonymous;
    	rateMatrix[3][49] = fFreqs[49] * omega ; 	//AAT->TAT: IS NonSynonymous;
    	rateMatrix[4][0] = fFreqs[0] * omega ; 	//ACA->AAA: IS NonSynonymous;
    	rateMatrix[4][5] = fFreqs[5] ; 	//ACA->ACC: 
    	rateMatrix[4][6] = fFreqs[6] * k ; 	//ACA->ACG: IS Transition;
    	rateMatrix[4][7] = fFreqs[7] ; 	//ACA->ACT: 
    	rateMatrix[4][8] = fFreqs[8] * omega ; 	//ACA->AGA: IS NonSynonymous;
    	rateMatrix[4][12] = fFreqs[12] * k * omega ; 	//ACA->ATA: IS Transition;IS NonSynonymous;
    	rateMatrix[4][20] = fFreqs[20] * omega ; 	//ACA->CCA: IS NonSynonymous;
    	rateMatrix[4][36] = fFreqs[36] * k * omega ; 	//ACA->GCA: IS Transition;IS NonSynonymous;
    	rateMatrix[4][50] = fFreqs[50] * omega ; 	//ACA->TCA: IS NonSynonymous;
    	rateMatrix[5][1] = fFreqs[1] * omega ; 	//ACC->AAC: IS NonSynonymous;
    	rateMatrix[5][4] = fFreqs[4] ; 	//ACC->ACA: 
    	rateMatrix[5][6] = fFreqs[6] ; 	//ACC->ACG: 
    	rateMatrix[5][7] = fFreqs[7] * k ; 	//ACC->ACT: IS Transition;
    	rateMatrix[5][9] = fFreqs[9] * omega ; 	//ACC->AGC: IS NonSynonymous;
    	rateMatrix[5][13] = fFreqs[13] * k * omega ; 	//ACC->ATC: IS Transition;IS NonSynonymous;
    	rateMatrix[5][21] = fFreqs[21] * omega ; 	//ACC->CCC: IS NonSynonymous;
    	rateMatrix[5][37] = fFreqs[37] * k * omega ; 	//ACC->GCC: IS Transition;IS NonSynonymous;
    	rateMatrix[5][51] = fFreqs[51] * omega ; 	//ACC->TCC: IS NonSynonymous;
    	rateMatrix[6][2] = fFreqs[2] * omega ; 	//ACG->AAG: IS NonSynonymous;
    	rateMatrix[6][4] = fFreqs[4] * k ; 	//ACG->ACA: IS Transition;
    	rateMatrix[6][5] = fFreqs[5] ; 	//ACG->ACC: 
    	rateMatrix[6][7] = fFreqs[7] ; 	//ACG->ACT: 
    	rateMatrix[6][10] = fFreqs[10] * omega ; 	//ACG->AGG: IS NonSynonymous;
    	rateMatrix[6][14] = fFreqs[14] * k * omega ; 	//ACG->ATG: IS Transition;IS NonSynonymous;
    	rateMatrix[6][22] = fFreqs[22] * omega ; 	//ACG->CCG: IS NonSynonymous;
    	rateMatrix[6][38] = fFreqs[38] * k * omega ; 	//ACG->GCG: IS Transition;IS NonSynonymous;
    	rateMatrix[6][52] = fFreqs[52] * omega ; 	//ACG->TCG: IS NonSynonymous;
    	rateMatrix[7][3] = fFreqs[3] * omega ; 	//ACT->AAT: IS NonSynonymous;
    	rateMatrix[7][4] = fFreqs[4] ; 	//ACT->ACA: 
    	rateMatrix[7][5] = fFreqs[5] * k ; 	//ACT->ACC: IS Transition;
    	rateMatrix[7][6] = fFreqs[6] ; 	//ACT->ACG: 
    	rateMatrix[7][11] = fFreqs[11] * omega ; 	//ACT->AGT: IS NonSynonymous;
    	rateMatrix[7][15] = fFreqs[15] * k * omega ; 	//ACT->ATT: IS Transition;IS NonSynonymous;
    	rateMatrix[7][23] = fFreqs[23] * omega ; 	//ACT->CCT: IS NonSynonymous;
    	rateMatrix[7][39] = fFreqs[39] * k * omega ; 	//ACT->GCT: IS Transition;IS NonSynonymous;
    	rateMatrix[7][53] = fFreqs[53] * omega ; 	//ACT->TCT: IS NonSynonymous;
    	rateMatrix[8][0] = fFreqs[0] * k * omega ; 	//AGA->AAA: IS Transition;IS NonSynonymous;
    	rateMatrix[8][4] = fFreqs[4] * omega ; 	//AGA->ACA: IS NonSynonymous;
    	rateMatrix[8][9] = fFreqs[9] * omega ; 	//AGA->AGC: IS NonSynonymous;
    	rateMatrix[8][10] = fFreqs[10] * k ; 	//AGA->AGG: IS Transition;
    	rateMatrix[8][11] = fFreqs[11] * omega ; 	//AGA->AGT: IS NonSynonymous;
    	rateMatrix[8][12] = fFreqs[12] * omega ; 	//AGA->ATA: IS NonSynonymous;
    	rateMatrix[8][24] = fFreqs[24] ; 	//AGA->CGA: 
    	rateMatrix[8][40] = fFreqs[40] * k * omega ; 	//AGA->GGA: IS Transition;IS NonSynonymous;
    	rateMatrix[9][1] = fFreqs[1] * k * omega ; 	//AGC->AAC: IS Transition;IS NonSynonymous;
    	rateMatrix[9][5] = fFreqs[5] * omega ; 	//AGC->ACC: IS NonSynonymous;
    	rateMatrix[9][8] = fFreqs[8] * omega ; 	//AGC->AGA: IS NonSynonymous;
    	rateMatrix[9][10] = fFreqs[10] * omega ; 	//AGC->AGG: IS NonSynonymous;
    	rateMatrix[9][11] = fFreqs[11] * k ; 	//AGC->AGT: IS Transition;
    	rateMatrix[9][13] = fFreqs[13] * omega ; 	//AGC->ATC: IS NonSynonymous;
    	rateMatrix[9][25] = fFreqs[25] * omega ; 	//AGC->CGC: IS NonSynonymous;
    	rateMatrix[9][41] = fFreqs[41] * k * omega ; 	//AGC->GGC: IS Transition;IS NonSynonymous;
    	rateMatrix[9][54] = fFreqs[54] * omega ; 	//AGC->TGC: IS NonSynonymous;
    	rateMatrix[10][2] = fFreqs[2] * k * omega ; 	//AGG->AAG: IS Transition;IS NonSynonymous;
    	rateMatrix[10][6] = fFreqs[6] * omega ; 	//AGG->ACG: IS NonSynonymous;
    	rateMatrix[10][8] = fFreqs[8] * k ; 	//AGG->AGA: IS Transition;
    	rateMatrix[10][9] = fFreqs[9] * omega ; 	//AGG->AGC: IS NonSynonymous;
    	rateMatrix[10][11] = fFreqs[11] * omega ; 	//AGG->AGT: IS NonSynonymous;
    	rateMatrix[10][14] = fFreqs[14] * omega ; 	//AGG->ATG: IS NonSynonymous;
    	rateMatrix[10][26] = fFreqs[26] ; 	//AGG->CGG: 
    	rateMatrix[10][42] = fFreqs[42] * k * omega ; 	//AGG->GGG: IS Transition;IS NonSynonymous;
    	rateMatrix[10][55] = fFreqs[55] * omega ; 	//AGG->TGG: IS NonSynonymous;
    	rateMatrix[11][3] = fFreqs[3] * k * omega ; 	//AGT->AAT: IS Transition;IS NonSynonymous;
    	rateMatrix[11][7] = fFreqs[7] * omega ; 	//AGT->ACT: IS NonSynonymous;
    	rateMatrix[11][8] = fFreqs[8] * omega ; 	//AGT->AGA: IS NonSynonymous;
    	rateMatrix[11][9] = fFreqs[9] * k ; 	//AGT->AGC: IS Transition;
    	rateMatrix[11][10] = fFreqs[10] * omega ; 	//AGT->AGG: IS NonSynonymous;
    	rateMatrix[11][15] = fFreqs[15] * omega ; 	//AGT->ATT: IS NonSynonymous;
    	rateMatrix[11][27] = fFreqs[27] * omega ; 	//AGT->CGT: IS NonSynonymous;
    	rateMatrix[11][43] = fFreqs[43] * k * omega ; 	//AGT->GGT: IS Transition;IS NonSynonymous;
    	rateMatrix[11][56] = fFreqs[56] * omega ; 	//AGT->TGT: IS NonSynonymous;
    	rateMatrix[12][0] = fFreqs[0] * omega ; 	//ATA->AAA: IS NonSynonymous;
    	rateMatrix[12][4] = fFreqs[4] * k * omega ; 	//ATA->ACA: IS Transition;IS NonSynonymous;
    	rateMatrix[12][8] = fFreqs[8] * omega ; 	//ATA->AGA: IS NonSynonymous;
    	rateMatrix[12][13] = fFreqs[13] ; 	//ATA->ATC: 
    	rateMatrix[12][14] = fFreqs[14] * k * omega ; 	//ATA->ATG: IS Transition;IS NonSynonymous;
    	rateMatrix[12][15] = fFreqs[15] ; 	//ATA->ATT: 
    	rateMatrix[12][28] = fFreqs[28] * omega ; 	//ATA->CTA: IS NonSynonymous;
    	rateMatrix[12][44] = fFreqs[44] * k * omega ; 	//ATA->GTA: IS Transition;IS NonSynonymous;
    	rateMatrix[12][57] = fFreqs[57] * omega ; 	//ATA->TTA: IS NonSynonymous;
    	rateMatrix[13][1] = fFreqs[1] * omega ; 	//ATC->AAC: IS NonSynonymous;
    	rateMatrix[13][5] = fFreqs[5] * k * omega ; 	//ATC->ACC: IS Transition;IS NonSynonymous;
    	rateMatrix[13][9] = fFreqs[9] * omega ; 	//ATC->AGC: IS NonSynonymous;
    	rateMatrix[13][12] = fFreqs[12] ; 	//ATC->ATA: 
    	rateMatrix[13][14] = fFreqs[14] * omega ; 	//ATC->ATG: IS NonSynonymous;
    	rateMatrix[13][15] = fFreqs[15] * k ; 	//ATC->ATT: IS Transition;
    	rateMatrix[13][29] = fFreqs[29] * omega ; 	//ATC->CTC: IS NonSynonymous;
    	rateMatrix[13][45] = fFreqs[45] * k * omega ; 	//ATC->GTC: IS Transition;IS NonSynonymous;
    	rateMatrix[13][58] = fFreqs[58] * omega ; 	//ATC->TTC: IS NonSynonymous;
    	rateMatrix[14][2] = fFreqs[2] * omega ; 	//ATG->AAG: IS NonSynonymous;
    	rateMatrix[14][6] = fFreqs[6] * k * omega ; 	//ATG->ACG: IS Transition;IS NonSynonymous;
    	rateMatrix[14][10] = fFreqs[10] * omega ; 	//ATG->AGG: IS NonSynonymous;
    	rateMatrix[14][12] = fFreqs[12] * k * omega ; 	//ATG->ATA: IS Transition;IS NonSynonymous;
    	rateMatrix[14][13] = fFreqs[13] * omega ; 	//ATG->ATC: IS NonSynonymous;
    	rateMatrix[14][15] = fFreqs[15] * omega ; 	//ATG->ATT: IS NonSynonymous;
    	rateMatrix[14][30] = fFreqs[30] * omega ; 	//ATG->CTG: IS NonSynonymous;
    	rateMatrix[14][46] = fFreqs[46] * k * omega ; 	//ATG->GTG: IS Transition;IS NonSynonymous;
    	rateMatrix[14][59] = fFreqs[59] * omega ; 	//ATG->TTG: IS NonSynonymous;
    	rateMatrix[15][3] = fFreqs[3] * omega ; 	//ATT->AAT: IS NonSynonymous;
    	rateMatrix[15][7] = fFreqs[7] * k * omega ; 	//ATT->ACT: IS Transition;IS NonSynonymous;
    	rateMatrix[15][11] = fFreqs[11] * omega ; 	//ATT->AGT: IS NonSynonymous;
    	rateMatrix[15][12] = fFreqs[12] ; 	//ATT->ATA: 
    	rateMatrix[15][13] = fFreqs[13] * k ; 	//ATT->ATC: IS Transition;
    	rateMatrix[15][14] = fFreqs[14] * omega ; 	//ATT->ATG: IS NonSynonymous;
    	rateMatrix[15][31] = fFreqs[31] * omega ; 	//ATT->CTT: IS NonSynonymous;
    	rateMatrix[15][47] = fFreqs[47] * k * omega ; 	//ATT->GTT: IS Transition;IS NonSynonymous;
    	rateMatrix[15][60] = fFreqs[60] * omega ; 	//ATT->TTT: IS NonSynonymous;
    	rateMatrix[16][0] = fFreqs[0] * omega ; 	//CAA->AAA: IS NonSynonymous;
    	rateMatrix[16][17] = fFreqs[17] * omega ; 	//CAA->CAC: IS NonSynonymous;
    	rateMatrix[16][18] = fFreqs[18] * k ; 	//CAA->CAG: IS Transition;
    	rateMatrix[16][19] = fFreqs[19] * omega ; 	//CAA->CAT: IS NonSynonymous;
    	rateMatrix[16][20] = fFreqs[20] * omega ; 	//CAA->CCA: IS NonSynonymous;
    	rateMatrix[16][24] = fFreqs[24] * k * omega ; 	//CAA->CGA: IS Transition;IS NonSynonymous;
    	rateMatrix[16][28] = fFreqs[28] * omega ; 	//CAA->CTA: IS NonSynonymous;
    	rateMatrix[16][32] = fFreqs[32] * omega ; 	//CAA->GAA: IS NonSynonymous;
    	rateMatrix[17][1] = fFreqs[1] * omega ; 	//CAC->AAC: IS NonSynonymous;
    	rateMatrix[17][16] = fFreqs[16] * omega ; 	//CAC->CAA: IS NonSynonymous;
    	rateMatrix[17][18] = fFreqs[18] * omega ; 	//CAC->CAG: IS NonSynonymous;
    	rateMatrix[17][19] = fFreqs[19] * k ; 	//CAC->CAT: IS Transition;
    	rateMatrix[17][21] = fFreqs[21] * omega ; 	//CAC->CCC: IS NonSynonymous;
    	rateMatrix[17][25] = fFreqs[25] * k * omega ; 	//CAC->CGC: IS Transition;IS NonSynonymous;
    	rateMatrix[17][29] = fFreqs[29] * omega ; 	//CAC->CTC: IS NonSynonymous;
    	rateMatrix[17][33] = fFreqs[33] * omega ; 	//CAC->GAC: IS NonSynonymous;
    	rateMatrix[17][48] = fFreqs[48] * k * omega ; 	//CAC->TAC: IS Transition;IS NonSynonymous;
    	rateMatrix[18][2] = fFreqs[2] * omega ; 	//CAG->AAG: IS NonSynonymous;
    	rateMatrix[18][16] = fFreqs[16] * k ; 	//CAG->CAA: IS Transition;
    	rateMatrix[18][17] = fFreqs[17] * omega ; 	//CAG->CAC: IS NonSynonymous;
    	rateMatrix[18][19] = fFreqs[19] * omega ; 	//CAG->CAT: IS NonSynonymous;
    	rateMatrix[18][22] = fFreqs[22] * omega ; 	//CAG->CCG: IS NonSynonymous;
    	rateMatrix[18][26] = fFreqs[26] * k * omega ; 	//CAG->CGG: IS Transition;IS NonSynonymous;
    	rateMatrix[18][30] = fFreqs[30] * omega ; 	//CAG->CTG: IS NonSynonymous;
    	rateMatrix[18][34] = fFreqs[34] * omega ; 	//CAG->GAG: IS NonSynonymous;
    	rateMatrix[19][3] = fFreqs[3] * omega ; 	//CAT->AAT: IS NonSynonymous;
    	rateMatrix[19][16] = fFreqs[16] * omega ; 	//CAT->CAA: IS NonSynonymous;
    	rateMatrix[19][17] = fFreqs[17] * k ; 	//CAT->CAC: IS Transition;
    	rateMatrix[19][18] = fFreqs[18] * omega ; 	//CAT->CAG: IS NonSynonymous;
    	rateMatrix[19][23] = fFreqs[23] * omega ; 	//CAT->CCT: IS NonSynonymous;
    	rateMatrix[19][27] = fFreqs[27] * k * omega ; 	//CAT->CGT: IS Transition;IS NonSynonymous;
    	rateMatrix[19][31] = fFreqs[31] * omega ; 	//CAT->CTT: IS NonSynonymous;
    	rateMatrix[19][35] = fFreqs[35] * omega ; 	//CAT->GAT: IS NonSynonymous;
    	rateMatrix[19][49] = fFreqs[49] * k * omega ; 	//CAT->TAT: IS Transition;IS NonSynonymous;
    	rateMatrix[20][4] = fFreqs[4] * omega ; 	//CCA->ACA: IS NonSynonymous;
    	rateMatrix[20][16] = fFreqs[16] * omega ; 	//CCA->CAA: IS NonSynonymous;
    	rateMatrix[20][21] = fFreqs[21] ; 	//CCA->CCC: 
    	rateMatrix[20][22] = fFreqs[22] * k ; 	//CCA->CCG: IS Transition;
    	rateMatrix[20][23] = fFreqs[23] ; 	//CCA->CCT: 
    	rateMatrix[20][24] = fFreqs[24] * omega ; 	//CCA->CGA: IS NonSynonymous;
    	rateMatrix[20][28] = fFreqs[28] * k * omega ; 	//CCA->CTA: IS Transition;IS NonSynonymous;
    	rateMatrix[20][36] = fFreqs[36] * omega ; 	//CCA->GCA: IS NonSynonymous;
    	rateMatrix[20][50] = fFreqs[50] * k * omega ; 	//CCA->TCA: IS Transition;IS NonSynonymous;
    	rateMatrix[21][5] = fFreqs[5] * omega ; 	//CCC->ACC: IS NonSynonymous;
    	rateMatrix[21][17] = fFreqs[17] * omega ; 	//CCC->CAC: IS NonSynonymous;
    	rateMatrix[21][20] = fFreqs[20] ; 	//CCC->CCA: 
    	rateMatrix[21][22] = fFreqs[22] ; 	//CCC->CCG: 
    	rateMatrix[21][23] = fFreqs[23] * k ; 	//CCC->CCT: IS Transition;
    	rateMatrix[21][25] = fFreqs[25] * omega ; 	//CCC->CGC: IS NonSynonymous;
    	rateMatrix[21][29] = fFreqs[29] * k * omega ; 	//CCC->CTC: IS Transition;IS NonSynonymous;
    	rateMatrix[21][37] = fFreqs[37] * omega ; 	//CCC->GCC: IS NonSynonymous;
    	rateMatrix[21][51] = fFreqs[51] * k * omega ; 	//CCC->TCC: IS Transition;IS NonSynonymous;
    	rateMatrix[22][6] = fFreqs[6] * omega ; 	//CCG->ACG: IS NonSynonymous;
    	rateMatrix[22][18] = fFreqs[18] * omega ; 	//CCG->CAG: IS NonSynonymous;
    	rateMatrix[22][20] = fFreqs[20] * k ; 	//CCG->CCA: IS Transition;
    	rateMatrix[22][21] = fFreqs[21] ; 	//CCG->CCC: 
    	rateMatrix[22][23] = fFreqs[23] ; 	//CCG->CCT: 
    	rateMatrix[22][26] = fFreqs[26] * omega ; 	//CCG->CGG: IS NonSynonymous;
    	rateMatrix[22][30] = fFreqs[30] * k * omega ; 	//CCG->CTG: IS Transition;IS NonSynonymous;
    	rateMatrix[22][38] = fFreqs[38] * omega ; 	//CCG->GCG: IS NonSynonymous;
    	rateMatrix[22][52] = fFreqs[52] * k * omega ; 	//CCG->TCG: IS Transition;IS NonSynonymous;
    	rateMatrix[23][7] = fFreqs[7] * omega ; 	//CCT->ACT: IS NonSynonymous;
    	rateMatrix[23][19] = fFreqs[19] * omega ; 	//CCT->CAT: IS NonSynonymous;
    	rateMatrix[23][20] = fFreqs[20] ; 	//CCT->CCA: 
    	rateMatrix[23][21] = fFreqs[21] * k ; 	//CCT->CCC: IS Transition;
    	rateMatrix[23][22] = fFreqs[22] ; 	//CCT->CCG: 
    	rateMatrix[23][27] = fFreqs[27] * omega ; 	//CCT->CGT: IS NonSynonymous;
    	rateMatrix[23][31] = fFreqs[31] * k * omega ; 	//CCT->CTT: IS Transition;IS NonSynonymous;
    	rateMatrix[23][39] = fFreqs[39] * omega ; 	//CCT->GCT: IS NonSynonymous;
    	rateMatrix[23][53] = fFreqs[53] * k * omega ; 	//CCT->TCT: IS Transition;IS NonSynonymous;
    	rateMatrix[24][8] = fFreqs[8] ; 	//CGA->AGA: 
    	rateMatrix[24][16] = fFreqs[16] * k * omega ; 	//CGA->CAA: IS Transition;IS NonSynonymous;
    	rateMatrix[24][20] = fFreqs[20] * omega ; 	//CGA->CCA: IS NonSynonymous;
    	rateMatrix[24][25] = fFreqs[25] ; 	//CGA->CGC: 
    	rateMatrix[24][26] = fFreqs[26] * k ; 	//CGA->CGG: IS Transition;
    	rateMatrix[24][27] = fFreqs[27] ; 	//CGA->CGT: 
    	rateMatrix[24][28] = fFreqs[28] * omega ; 	//CGA->CTA: IS NonSynonymous;
    	rateMatrix[24][40] = fFreqs[40] * omega ; 	//CGA->GGA: IS NonSynonymous;
    	rateMatrix[25][9] = fFreqs[9] * omega ; 	//CGC->AGC: IS NonSynonymous;
    	rateMatrix[25][17] = fFreqs[17] * k * omega ; 	//CGC->CAC: IS Transition;IS NonSynonymous;
    	rateMatrix[25][21] = fFreqs[21] * omega ; 	//CGC->CCC: IS NonSynonymous;
    	rateMatrix[25][24] = fFreqs[24] ; 	//CGC->CGA: 
    	rateMatrix[25][26] = fFreqs[26] ; 	//CGC->CGG: 
    	rateMatrix[25][27] = fFreqs[27] * k ; 	//CGC->CGT: IS Transition;
    	rateMatrix[25][29] = fFreqs[29] * omega ; 	//CGC->CTC: IS NonSynonymous;
    	rateMatrix[25][41] = fFreqs[41] * omega ; 	//CGC->GGC: IS NonSynonymous;
    	rateMatrix[25][54] = fFreqs[54] * k * omega ; 	//CGC->TGC: IS Transition;IS NonSynonymous;
    	rateMatrix[26][10] = fFreqs[10] ; 	//CGG->AGG: 
    	rateMatrix[26][18] = fFreqs[18] * k * omega ; 	//CGG->CAG: IS Transition;IS NonSynonymous;
    	rateMatrix[26][22] = fFreqs[22] * omega ; 	//CGG->CCG: IS NonSynonymous;
    	rateMatrix[26][24] = fFreqs[24] * k ; 	//CGG->CGA: IS Transition;
    	rateMatrix[26][25] = fFreqs[25] ; 	//CGG->CGC: 
    	rateMatrix[26][27] = fFreqs[27] ; 	//CGG->CGT: 
    	rateMatrix[26][30] = fFreqs[30] * omega ; 	//CGG->CTG: IS NonSynonymous;
    	rateMatrix[26][42] = fFreqs[42] * omega ; 	//CGG->GGG: IS NonSynonymous;
    	rateMatrix[26][55] = fFreqs[55] * k * omega ; 	//CGG->TGG: IS Transition;IS NonSynonymous;
    	rateMatrix[27][11] = fFreqs[11] * omega ; 	//CGT->AGT: IS NonSynonymous;
    	rateMatrix[27][19] = fFreqs[19] * k * omega ; 	//CGT->CAT: IS Transition;IS NonSynonymous;
    	rateMatrix[27][23] = fFreqs[23] * omega ; 	//CGT->CCT: IS NonSynonymous;
    	rateMatrix[27][24] = fFreqs[24] ; 	//CGT->CGA: 
    	rateMatrix[27][25] = fFreqs[25] * k ; 	//CGT->CGC: IS Transition;
    	rateMatrix[27][26] = fFreqs[26] ; 	//CGT->CGG: 
    	rateMatrix[27][31] = fFreqs[31] * omega ; 	//CGT->CTT: IS NonSynonymous;
    	rateMatrix[27][43] = fFreqs[43] * omega ; 	//CGT->GGT: IS NonSynonymous;
    	rateMatrix[27][56] = fFreqs[56] * k * omega ; 	//CGT->TGT: IS Transition;IS NonSynonymous;
    	rateMatrix[28][12] = fFreqs[12] * omega ; 	//CTA->ATA: IS NonSynonymous;
    	rateMatrix[28][16] = fFreqs[16] * omega ; 	//CTA->CAA: IS NonSynonymous;
    	rateMatrix[28][20] = fFreqs[20] * k * omega ; 	//CTA->CCA: IS Transition;IS NonSynonymous;
    	rateMatrix[28][24] = fFreqs[24] * omega ; 	//CTA->CGA: IS NonSynonymous;
    	rateMatrix[28][29] = fFreqs[29] ; 	//CTA->CTC: 
    	rateMatrix[28][30] = fFreqs[30] * k ; 	//CTA->CTG: IS Transition;
    	rateMatrix[28][31] = fFreqs[31] ; 	//CTA->CTT: 
    	rateMatrix[28][44] = fFreqs[44] * omega ; 	//CTA->GTA: IS NonSynonymous;
    	rateMatrix[28][57] = fFreqs[57] * k ; 	//CTA->TTA: IS Transition;
    	rateMatrix[29][13] = fFreqs[13] * omega ; 	//CTC->ATC: IS NonSynonymous;
    	rateMatrix[29][17] = fFreqs[17] * omega ; 	//CTC->CAC: IS NonSynonymous;
    	rateMatrix[29][21] = fFreqs[21] * k * omega ; 	//CTC->CCC: IS Transition;IS NonSynonymous;
    	rateMatrix[29][25] = fFreqs[25] * omega ; 	//CTC->CGC: IS NonSynonymous;
    	rateMatrix[29][28] = fFreqs[28] ; 	//CTC->CTA: 
    	rateMatrix[29][30] = fFreqs[30] ; 	//CTC->CTG: 
    	rateMatrix[29][31] = fFreqs[31] * k ; 	//CTC->CTT: IS Transition;
    	rateMatrix[29][45] = fFreqs[45] * omega ; 	//CTC->GTC: IS NonSynonymous;
    	rateMatrix[29][58] = fFreqs[58] * k * omega ; 	//CTC->TTC: IS Transition;IS NonSynonymous;
    	rateMatrix[30][14] = fFreqs[14] * omega ; 	//CTG->ATG: IS NonSynonymous;
    	rateMatrix[30][18] = fFreqs[18] * omega ; 	//CTG->CAG: IS NonSynonymous;
    	rateMatrix[30][22] = fFreqs[22] * k * omega ; 	//CTG->CCG: IS Transition;IS NonSynonymous;
    	rateMatrix[30][26] = fFreqs[26] * omega ; 	//CTG->CGG: IS NonSynonymous;
    	rateMatrix[30][28] = fFreqs[28] * k ; 	//CTG->CTA: IS Transition;
    	rateMatrix[30][29] = fFreqs[29] ; 	//CTG->CTC: 
    	rateMatrix[30][31] = fFreqs[31] ; 	//CTG->CTT: 
    	rateMatrix[30][46] = fFreqs[46] * omega ; 	//CTG->GTG: IS NonSynonymous;
    	rateMatrix[30][59] = fFreqs[59] * k ; 	//CTG->TTG: IS Transition;
    	rateMatrix[31][15] = fFreqs[15] * omega ; 	//CTT->ATT: IS NonSynonymous;
    	rateMatrix[31][19] = fFreqs[19] * omega ; 	//CTT->CAT: IS NonSynonymous;
    	rateMatrix[31][23] = fFreqs[23] * k * omega ; 	//CTT->CCT: IS Transition;IS NonSynonymous;
    	rateMatrix[31][27] = fFreqs[27] * omega ; 	//CTT->CGT: IS NonSynonymous;
    	rateMatrix[31][28] = fFreqs[28] ; 	//CTT->CTA: 
    	rateMatrix[31][29] = fFreqs[29] * k ; 	//CTT->CTC: IS Transition;
    	rateMatrix[31][30] = fFreqs[30] ; 	//CTT->CTG: 
    	rateMatrix[31][47] = fFreqs[47] * omega ; 	//CTT->GTT: IS NonSynonymous;
    	rateMatrix[31][60] = fFreqs[60] * k * omega ; 	//CTT->TTT: IS Transition;IS NonSynonymous;
    	rateMatrix[32][0] = fFreqs[0] * k * omega ; 	//GAA->AAA: IS Transition;IS NonSynonymous;
    	rateMatrix[32][16] = fFreqs[16] * omega ; 	//GAA->CAA: IS NonSynonymous;
    	rateMatrix[32][33] = fFreqs[33] * omega ; 	//GAA->GAC: IS NonSynonymous;
    	rateMatrix[32][34] = fFreqs[34] * k ; 	//GAA->GAG: IS Transition;
    	rateMatrix[32][35] = fFreqs[35] * omega ; 	//GAA->GAT: IS NonSynonymous;
    	rateMatrix[32][36] = fFreqs[36] * omega ; 	//GAA->GCA: IS NonSynonymous;
    	rateMatrix[32][40] = fFreqs[40] * k * omega ; 	//GAA->GGA: IS Transition;IS NonSynonymous;
    	rateMatrix[32][44] = fFreqs[44] * omega ; 	//GAA->GTA: IS NonSynonymous;
    	rateMatrix[33][1] = fFreqs[1] * k * omega ; 	//GAC->AAC: IS Transition;IS NonSynonymous;
    	rateMatrix[33][17] = fFreqs[17] * omega ; 	//GAC->CAC: IS NonSynonymous;
    	rateMatrix[33][32] = fFreqs[32] * omega ; 	//GAC->GAA: IS NonSynonymous;
    	rateMatrix[33][34] = fFreqs[34] * omega ; 	//GAC->GAG: IS NonSynonymous;
    	rateMatrix[33][35] = fFreqs[35] * k ; 	//GAC->GAT: IS Transition;
    	rateMatrix[33][37] = fFreqs[37] * omega ; 	//GAC->GCC: IS NonSynonymous;
    	rateMatrix[33][41] = fFreqs[41] * k * omega ; 	//GAC->GGC: IS Transition;IS NonSynonymous;
    	rateMatrix[33][45] = fFreqs[45] * omega ; 	//GAC->GTC: IS NonSynonymous;
    	rateMatrix[33][48] = fFreqs[48] * omega ; 	//GAC->TAC: IS NonSynonymous;
    	rateMatrix[34][2] = fFreqs[2] * k * omega ; 	//GAG->AAG: IS Transition;IS NonSynonymous;
    	rateMatrix[34][18] = fFreqs[18] * omega ; 	//GAG->CAG: IS NonSynonymous;
    	rateMatrix[34][32] = fFreqs[32] * k ; 	//GAG->GAA: IS Transition;
    	rateMatrix[34][33] = fFreqs[33] * omega ; 	//GAG->GAC: IS NonSynonymous;
    	rateMatrix[34][35] = fFreqs[35] * omega ; 	//GAG->GAT: IS NonSynonymous;
    	rateMatrix[34][38] = fFreqs[38] * omega ; 	//GAG->GCG: IS NonSynonymous;
    	rateMatrix[34][42] = fFreqs[42] * k * omega ; 	//GAG->GGG: IS Transition;IS NonSynonymous;
    	rateMatrix[34][46] = fFreqs[46] * omega ; 	//GAG->GTG: IS NonSynonymous;
    	rateMatrix[35][3] = fFreqs[3] * k * omega ; 	//GAT->AAT: IS Transition;IS NonSynonymous;
    	rateMatrix[35][19] = fFreqs[19] * omega ; 	//GAT->CAT: IS NonSynonymous;
    	rateMatrix[35][32] = fFreqs[32] * omega ; 	//GAT->GAA: IS NonSynonymous;
    	rateMatrix[35][33] = fFreqs[33] * k ; 	//GAT->GAC: IS Transition;
    	rateMatrix[35][34] = fFreqs[34] * omega ; 	//GAT->GAG: IS NonSynonymous;
    	rateMatrix[35][39] = fFreqs[39] * omega ; 	//GAT->GCT: IS NonSynonymous;
    	rateMatrix[35][43] = fFreqs[43] * k * omega ; 	//GAT->GGT: IS Transition;IS NonSynonymous;
    	rateMatrix[35][47] = fFreqs[47] * omega ; 	//GAT->GTT: IS NonSynonymous;
    	rateMatrix[35][49] = fFreqs[49] * omega ; 	//GAT->TAT: IS NonSynonymous;
    	rateMatrix[36][4] = fFreqs[4] * k * omega ; 	//GCA->ACA: IS Transition;IS NonSynonymous;
    	rateMatrix[36][20] = fFreqs[20] * omega ; 	//GCA->CCA: IS NonSynonymous;
    	rateMatrix[36][32] = fFreqs[32] * omega ; 	//GCA->GAA: IS NonSynonymous;
    	rateMatrix[36][37] = fFreqs[37] ; 	//GCA->GCC: 
    	rateMatrix[36][38] = fFreqs[38] * k ; 	//GCA->GCG: IS Transition;
    	rateMatrix[36][39] = fFreqs[39] ; 	//GCA->GCT: 
    	rateMatrix[36][40] = fFreqs[40] * omega ; 	//GCA->GGA: IS NonSynonymous;
    	rateMatrix[36][44] = fFreqs[44] * k * omega ; 	//GCA->GTA: IS Transition;IS NonSynonymous;
    	rateMatrix[36][50] = fFreqs[50] * omega ; 	//GCA->TCA: IS NonSynonymous;
    	rateMatrix[37][5] = fFreqs[5] * k * omega ; 	//GCC->ACC: IS Transition;IS NonSynonymous;
    	rateMatrix[37][21] = fFreqs[21] * omega ; 	//GCC->CCC: IS NonSynonymous;
    	rateMatrix[37][33] = fFreqs[33] * omega ; 	//GCC->GAC: IS NonSynonymous;
    	rateMatrix[37][36] = fFreqs[36] ; 	//GCC->GCA: 
    	rateMatrix[37][38] = fFreqs[38] ; 	//GCC->GCG: 
    	rateMatrix[37][39] = fFreqs[39] * k ; 	//GCC->GCT: IS Transition;
    	rateMatrix[37][41] = fFreqs[41] * omega ; 	//GCC->GGC: IS NonSynonymous;
    	rateMatrix[37][45] = fFreqs[45] * k * omega ; 	//GCC->GTC: IS Transition;IS NonSynonymous;
    	rateMatrix[37][51] = fFreqs[51] * omega ; 	//GCC->TCC: IS NonSynonymous;
    	rateMatrix[38][6] = fFreqs[6] * k * omega ; 	//GCG->ACG: IS Transition;IS NonSynonymous;
    	rateMatrix[38][22] = fFreqs[22] * omega ; 	//GCG->CCG: IS NonSynonymous;
    	rateMatrix[38][34] = fFreqs[34] * omega ; 	//GCG->GAG: IS NonSynonymous;
    	rateMatrix[38][36] = fFreqs[36] * k ; 	//GCG->GCA: IS Transition;
    	rateMatrix[38][37] = fFreqs[37] ; 	//GCG->GCC: 
    	rateMatrix[38][39] = fFreqs[39] ; 	//GCG->GCT: 
    	rateMatrix[38][42] = fFreqs[42] * omega ; 	//GCG->GGG: IS NonSynonymous;
    	rateMatrix[38][46] = fFreqs[46] * k * omega ; 	//GCG->GTG: IS Transition;IS NonSynonymous;
    	rateMatrix[38][52] = fFreqs[52] * omega ; 	//GCG->TCG: IS NonSynonymous;
    	rateMatrix[39][7] = fFreqs[7] * k * omega ; 	//GCT->ACT: IS Transition;IS NonSynonymous;
    	rateMatrix[39][23] = fFreqs[23] * omega ; 	//GCT->CCT: IS NonSynonymous;
    	rateMatrix[39][35] = fFreqs[35] * omega ; 	//GCT->GAT: IS NonSynonymous;
    	rateMatrix[39][36] = fFreqs[36] ; 	//GCT->GCA: 
    	rateMatrix[39][37] = fFreqs[37] * k ; 	//GCT->GCC: IS Transition;
    	rateMatrix[39][38] = fFreqs[38] ; 	//GCT->GCG: 
    	rateMatrix[39][43] = fFreqs[43] * omega ; 	//GCT->GGT: IS NonSynonymous;
    	rateMatrix[39][47] = fFreqs[47] * k * omega ; 	//GCT->GTT: IS Transition;IS NonSynonymous;
    	rateMatrix[39][53] = fFreqs[53] * omega ; 	//GCT->TCT: IS NonSynonymous;
    	rateMatrix[40][8] = fFreqs[8] * k * omega ; 	//GGA->AGA: IS Transition;IS NonSynonymous;
    	rateMatrix[40][24] = fFreqs[24] * omega ; 	//GGA->CGA: IS NonSynonymous;
    	rateMatrix[40][32] = fFreqs[32] * k * omega ; 	//GGA->GAA: IS Transition;IS NonSynonymous;
    	rateMatrix[40][36] = fFreqs[36] * omega ; 	//GGA->GCA: IS NonSynonymous;
    	rateMatrix[40][41] = fFreqs[41] ; 	//GGA->GGC: 
    	rateMatrix[40][42] = fFreqs[42] * k ; 	//GGA->GGG: IS Transition;
    	rateMatrix[40][43] = fFreqs[43] ; 	//GGA->GGT: 
    	rateMatrix[40][44] = fFreqs[44] * omega ; 	//GGA->GTA: IS NonSynonymous;
    	rateMatrix[41][9] = fFreqs[9] * k * omega ; 	//GGC->AGC: IS Transition;IS NonSynonymous;
    	rateMatrix[41][25] = fFreqs[25] * omega ; 	//GGC->CGC: IS NonSynonymous;
    	rateMatrix[41][33] = fFreqs[33] * k * omega ; 	//GGC->GAC: IS Transition;IS NonSynonymous;
    	rateMatrix[41][37] = fFreqs[37] * omega ; 	//GGC->GCC: IS NonSynonymous;
    	rateMatrix[41][40] = fFreqs[40] ; 	//GGC->GGA: 
    	rateMatrix[41][42] = fFreqs[42] ; 	//GGC->GGG: 
    	rateMatrix[41][43] = fFreqs[43] * k ; 	//GGC->GGT: IS Transition;
    	rateMatrix[41][45] = fFreqs[45] * omega ; 	//GGC->GTC: IS NonSynonymous;
    	rateMatrix[41][54] = fFreqs[54] * omega ; 	//GGC->TGC: IS NonSynonymous;
    	rateMatrix[42][10] = fFreqs[10] * k * omega ; 	//GGG->AGG: IS Transition;IS NonSynonymous;
    	rateMatrix[42][26] = fFreqs[26] * omega ; 	//GGG->CGG: IS NonSynonymous;
    	rateMatrix[42][34] = fFreqs[34] * k * omega ; 	//GGG->GAG: IS Transition;IS NonSynonymous;
    	rateMatrix[42][38] = fFreqs[38] * omega ; 	//GGG->GCG: IS NonSynonymous;
    	rateMatrix[42][40] = fFreqs[40] * k ; 	//GGG->GGA: IS Transition;
    	rateMatrix[42][41] = fFreqs[41] ; 	//GGG->GGC: 
    	rateMatrix[42][43] = fFreqs[43] ; 	//GGG->GGT: 
    	rateMatrix[42][46] = fFreqs[46] * omega ; 	//GGG->GTG: IS NonSynonymous;
    	rateMatrix[42][55] = fFreqs[55] * omega ; 	//GGG->TGG: IS NonSynonymous;
    	rateMatrix[43][11] = fFreqs[11] * k * omega ; 	//GGT->AGT: IS Transition;IS NonSynonymous;
    	rateMatrix[43][27] = fFreqs[27] * omega ; 	//GGT->CGT: IS NonSynonymous;
    	rateMatrix[43][35] = fFreqs[35] * k * omega ; 	//GGT->GAT: IS Transition;IS NonSynonymous;
    	rateMatrix[43][39] = fFreqs[39] * omega ; 	//GGT->GCT: IS NonSynonymous;
    	rateMatrix[43][40] = fFreqs[40] ; 	//GGT->GGA: 
    	rateMatrix[43][41] = fFreqs[41] * k ; 	//GGT->GGC: IS Transition;
    	rateMatrix[43][42] = fFreqs[42] ; 	//GGT->GGG: 
    	rateMatrix[43][47] = fFreqs[47] * omega ; 	//GGT->GTT: IS NonSynonymous;
    	rateMatrix[43][56] = fFreqs[56] * omega ; 	//GGT->TGT: IS NonSynonymous;
    	rateMatrix[44][12] = fFreqs[12] * k * omega ; 	//GTA->ATA: IS Transition;IS NonSynonymous;
    	rateMatrix[44][28] = fFreqs[28] * omega ; 	//GTA->CTA: IS NonSynonymous;
    	rateMatrix[44][32] = fFreqs[32] * omega ; 	//GTA->GAA: IS NonSynonymous;
    	rateMatrix[44][36] = fFreqs[36] * k * omega ; 	//GTA->GCA: IS Transition;IS NonSynonymous;
    	rateMatrix[44][40] = fFreqs[40] * omega ; 	//GTA->GGA: IS NonSynonymous;
    	rateMatrix[44][45] = fFreqs[45] ; 	//GTA->GTC: 
    	rateMatrix[44][46] = fFreqs[46] * k ; 	//GTA->GTG: IS Transition;
    	rateMatrix[44][47] = fFreqs[47] ; 	//GTA->GTT: 
    	rateMatrix[44][57] = fFreqs[57] * omega ; 	//GTA->TTA: IS NonSynonymous;
    	rateMatrix[45][13] = fFreqs[13] * k * omega ; 	//GTC->ATC: IS Transition;IS NonSynonymous;
    	rateMatrix[45][29] = fFreqs[29] * omega ; 	//GTC->CTC: IS NonSynonymous;
    	rateMatrix[45][33] = fFreqs[33] * omega ; 	//GTC->GAC: IS NonSynonymous;
    	rateMatrix[45][37] = fFreqs[37] * k * omega ; 	//GTC->GCC: IS Transition;IS NonSynonymous;
    	rateMatrix[45][41] = fFreqs[41] * omega ; 	//GTC->GGC: IS NonSynonymous;
    	rateMatrix[45][44] = fFreqs[44] ; 	//GTC->GTA: 
    	rateMatrix[45][46] = fFreqs[46] ; 	//GTC->GTG: 
    	rateMatrix[45][47] = fFreqs[47] * k ; 	//GTC->GTT: IS Transition;
    	rateMatrix[45][58] = fFreqs[58] * omega ; 	//GTC->TTC: IS NonSynonymous;
    	rateMatrix[46][14] = fFreqs[14] * k * omega ; 	//GTG->ATG: IS Transition;IS NonSynonymous;
    	rateMatrix[46][30] = fFreqs[30] * omega ; 	//GTG->CTG: IS NonSynonymous;
    	rateMatrix[46][34] = fFreqs[34] * omega ; 	//GTG->GAG: IS NonSynonymous;
    	rateMatrix[46][38] = fFreqs[38] * k * omega ; 	//GTG->GCG: IS Transition;IS NonSynonymous;
    	rateMatrix[46][42] = fFreqs[42] * omega ; 	//GTG->GGG: IS NonSynonymous;
    	rateMatrix[46][44] = fFreqs[44] * k ; 	//GTG->GTA: IS Transition;
    	rateMatrix[46][45] = fFreqs[45] ; 	//GTG->GTC: 
    	rateMatrix[46][47] = fFreqs[47] ; 	//GTG->GTT: 
    	rateMatrix[46][59] = fFreqs[59] * omega ; 	//GTG->TTG: IS NonSynonymous;
    	rateMatrix[47][15] = fFreqs[15] * k * omega ; 	//GTT->ATT: IS Transition;IS NonSynonymous;
    	rateMatrix[47][31] = fFreqs[31] * omega ; 	//GTT->CTT: IS NonSynonymous;
    	rateMatrix[47][35] = fFreqs[35] * omega ; 	//GTT->GAT: IS NonSynonymous;
    	rateMatrix[47][39] = fFreqs[39] * k * omega ; 	//GTT->GCT: IS Transition;IS NonSynonymous;
    	rateMatrix[47][43] = fFreqs[43] * omega ; 	//GTT->GGT: IS NonSynonymous;
    	rateMatrix[47][44] = fFreqs[44] ; 	//GTT->GTA: 
    	rateMatrix[47][45] = fFreqs[45] * k ; 	//GTT->GTC: IS Transition;
    	rateMatrix[47][46] = fFreqs[46] ; 	//GTT->GTG: 
    	rateMatrix[47][60] = fFreqs[60] * omega ; 	//GTT->TTT: IS NonSynonymous;
    	rateMatrix[48][1] = fFreqs[1] * omega ; 	//TAC->AAC: IS NonSynonymous;
    	rateMatrix[48][17] = fFreqs[17] * k * omega ; 	//TAC->CAC: IS Transition;IS NonSynonymous;
    	rateMatrix[48][33] = fFreqs[33] * omega ; 	//TAC->GAC: IS NonSynonymous;
    	rateMatrix[48][49] = fFreqs[49] * k ; 	//TAC->TAT: IS Transition;
    	rateMatrix[48][51] = fFreqs[51] * omega ; 	//TAC->TCC: IS NonSynonymous;
    	rateMatrix[48][54] = fFreqs[54] * k * omega ; 	//TAC->TGC: IS Transition;IS NonSynonymous;
    	rateMatrix[48][58] = fFreqs[58] * omega ; 	//TAC->TTC: IS NonSynonymous;
    	rateMatrix[49][3] = fFreqs[3] * omega ; 	//TAT->AAT: IS NonSynonymous;
    	rateMatrix[49][19] = fFreqs[19] * k * omega ; 	//TAT->CAT: IS Transition;IS NonSynonymous;
    	rateMatrix[49][35] = fFreqs[35] * omega ; 	//TAT->GAT: IS NonSynonymous;
    	rateMatrix[49][48] = fFreqs[48] * k ; 	//TAT->TAC: IS Transition;
    	rateMatrix[49][53] = fFreqs[53] * omega ; 	//TAT->TCT: IS NonSynonymous;
    	rateMatrix[49][56] = fFreqs[56] * k * omega ; 	//TAT->TGT: IS Transition;IS NonSynonymous;
    	rateMatrix[49][60] = fFreqs[60] * omega ; 	//TAT->TTT: IS NonSynonymous;
    	rateMatrix[50][4] = fFreqs[4] * omega ; 	//TCA->ACA: IS NonSynonymous;
    	rateMatrix[50][20] = fFreqs[20] * k * omega ; 	//TCA->CCA: IS Transition;IS NonSynonymous;
    	rateMatrix[50][36] = fFreqs[36] * omega ; 	//TCA->GCA: IS NonSynonymous;
    	rateMatrix[50][51] = fFreqs[51] ; 	//TCA->TCC: 
    	rateMatrix[50][52] = fFreqs[52] * k ; 	//TCA->TCG: IS Transition;
    	rateMatrix[50][53] = fFreqs[53] ; 	//TCA->TCT: 
    	rateMatrix[50][57] = fFreqs[57] * k * omega ; 	//TCA->TTA: IS Transition;IS NonSynonymous;
    	rateMatrix[51][5] = fFreqs[5] * omega ; 	//TCC->ACC: IS NonSynonymous;
    	rateMatrix[51][21] = fFreqs[21] * k * omega ; 	//TCC->CCC: IS Transition;IS NonSynonymous;
    	rateMatrix[51][37] = fFreqs[37] * omega ; 	//TCC->GCC: IS NonSynonymous;
    	rateMatrix[51][48] = fFreqs[48] * omega ; 	//TCC->TAC: IS NonSynonymous;
    	rateMatrix[51][50] = fFreqs[50] ; 	//TCC->TCA: 
    	rateMatrix[51][52] = fFreqs[52] ; 	//TCC->TCG: 
    	rateMatrix[51][53] = fFreqs[53] * k ; 	//TCC->TCT: IS Transition;
    	rateMatrix[51][54] = fFreqs[54] * omega ; 	//TCC->TGC: IS NonSynonymous;
    	rateMatrix[51][58] = fFreqs[58] * k * omega ; 	//TCC->TTC: IS Transition;IS NonSynonymous;
    	rateMatrix[52][6] = fFreqs[6] * omega ; 	//TCG->ACG: IS NonSynonymous;
    	rateMatrix[52][22] = fFreqs[22] * k * omega ; 	//TCG->CCG: IS Transition;IS NonSynonymous;
    	rateMatrix[52][38] = fFreqs[38] * omega ; 	//TCG->GCG: IS NonSynonymous;
    	rateMatrix[52][50] = fFreqs[50] * k ; 	//TCG->TCA: IS Transition;
    	rateMatrix[52][51] = fFreqs[51] ; 	//TCG->TCC: 
    	rateMatrix[52][53] = fFreqs[53] ; 	//TCG->TCT: 
    	rateMatrix[52][55] = fFreqs[55] * omega ; 	//TCG->TGG: IS NonSynonymous;
    	rateMatrix[52][59] = fFreqs[59] * k * omega ; 	//TCG->TTG: IS Transition;IS NonSynonymous;
    	rateMatrix[53][7] = fFreqs[7] * omega ; 	//TCT->ACT: IS NonSynonymous;
    	rateMatrix[53][23] = fFreqs[23] * k * omega ; 	//TCT->CCT: IS Transition;IS NonSynonymous;
    	rateMatrix[53][39] = fFreqs[39] * omega ; 	//TCT->GCT: IS NonSynonymous;
    	rateMatrix[53][49] = fFreqs[49] * omega ; 	//TCT->TAT: IS NonSynonymous;
    	rateMatrix[53][50] = fFreqs[50] ; 	//TCT->TCA: 
    	rateMatrix[53][51] = fFreqs[51] * k ; 	//TCT->TCC: IS Transition;
    	rateMatrix[53][52] = fFreqs[52] ; 	//TCT->TCG: 
    	rateMatrix[53][56] = fFreqs[56] * omega ; 	//TCT->TGT: IS NonSynonymous;
    	rateMatrix[53][60] = fFreqs[60] * k * omega ; 	//TCT->TTT: IS Transition;IS NonSynonymous;
    	rateMatrix[54][9] = fFreqs[9] * omega ; 	//TGC->AGC: IS NonSynonymous;
    	rateMatrix[54][25] = fFreqs[25] * k * omega ; 	//TGC->CGC: IS Transition;IS NonSynonymous;
    	rateMatrix[54][41] = fFreqs[41] * omega ; 	//TGC->GGC: IS NonSynonymous;
    	rateMatrix[54][48] = fFreqs[48] * k * omega ; 	//TGC->TAC: IS Transition;IS NonSynonymous;
    	rateMatrix[54][51] = fFreqs[51] * omega ; 	//TGC->TCC: IS NonSynonymous;
    	rateMatrix[54][55] = fFreqs[55] * omega ; 	//TGC->TGG: IS NonSynonymous;
    	rateMatrix[54][56] = fFreqs[56] * k ; 	//TGC->TGT: IS Transition;
    	rateMatrix[54][58] = fFreqs[58] * omega ; 	//TGC->TTC: IS NonSynonymous;
    	rateMatrix[55][10] = fFreqs[10] * omega ; 	//TGG->AGG: IS NonSynonymous;
    	rateMatrix[55][26] = fFreqs[26] * k * omega ; 	//TGG->CGG: IS Transition;IS NonSynonymous;
    	rateMatrix[55][42] = fFreqs[42] * omega ; 	//TGG->GGG: IS NonSynonymous;
    	rateMatrix[55][52] = fFreqs[52] * omega ; 	//TGG->TCG: IS NonSynonymous;
    	rateMatrix[55][54] = fFreqs[54] * omega ; 	//TGG->TGC: IS NonSynonymous;
    	rateMatrix[55][56] = fFreqs[56] * omega ; 	//TGG->TGT: IS NonSynonymous;
    	rateMatrix[55][59] = fFreqs[59] * omega ; 	//TGG->TTG: IS NonSynonymous;
    	rateMatrix[56][11] = fFreqs[11] * omega ; 	//TGT->AGT: IS NonSynonymous;
    	rateMatrix[56][27] = fFreqs[27] * k * omega ; 	//TGT->CGT: IS Transition;IS NonSynonymous;
    	rateMatrix[56][43] = fFreqs[43] * omega ; 	//TGT->GGT: IS NonSynonymous;
    	rateMatrix[56][49] = fFreqs[49] * k * omega ; 	//TGT->TAT: IS Transition;IS NonSynonymous;
    	rateMatrix[56][53] = fFreqs[53] * omega ; 	//TGT->TCT: IS NonSynonymous;
    	rateMatrix[56][54] = fFreqs[54] * k ; 	//TGT->TGC: IS Transition;
    	rateMatrix[56][55] = fFreqs[55] * omega ; 	//TGT->TGG: IS NonSynonymous;
    	rateMatrix[56][60] = fFreqs[60] * omega ; 	//TGT->TTT: IS NonSynonymous;
    	rateMatrix[57][12] = fFreqs[12] * omega ; 	//TTA->ATA: IS NonSynonymous;
    	rateMatrix[57][28] = fFreqs[28] * k ; 	//TTA->CTA: IS Transition;
    	rateMatrix[57][44] = fFreqs[44] * omega ; 	//TTA->GTA: IS NonSynonymous;
    	rateMatrix[57][50] = fFreqs[50] * k * omega ; 	//TTA->TCA: IS Transition;IS NonSynonymous;
    	rateMatrix[57][58] = fFreqs[58] * omega ; 	//TTA->TTC: IS NonSynonymous;
    	rateMatrix[57][59] = fFreqs[59] * k ; 	//TTA->TTG: IS Transition;
    	rateMatrix[57][60] = fFreqs[60] * omega ; 	//TTA->TTT: IS NonSynonymous;
    	rateMatrix[58][13] = fFreqs[13] * omega ; 	//TTC->ATC: IS NonSynonymous;
    	rateMatrix[58][29] = fFreqs[29] * k * omega ; 	//TTC->CTC: IS Transition;IS NonSynonymous;
    	rateMatrix[58][45] = fFreqs[45] * omega ; 	//TTC->GTC: IS NonSynonymous;
    	rateMatrix[58][48] = fFreqs[48] * omega ; 	//TTC->TAC: IS NonSynonymous;
    	rateMatrix[58][51] = fFreqs[51] * k * omega ; 	//TTC->TCC: IS Transition;IS NonSynonymous;
    	rateMatrix[58][54] = fFreqs[54] * omega ; 	//TTC->TGC: IS NonSynonymous;
    	rateMatrix[58][57] = fFreqs[57] * omega ; 	//TTC->TTA: IS NonSynonymous;
    	rateMatrix[58][59] = fFreqs[59] * omega ; 	//TTC->TTG: IS NonSynonymous;
    	rateMatrix[58][60] = fFreqs[60] * k ; 	//TTC->TTT: IS Transition;
    	rateMatrix[59][14] = fFreqs[14] * omega ; 	//TTG->ATG: IS NonSynonymous;
    	rateMatrix[59][30] = fFreqs[30] * k ; 	//TTG->CTG: IS Transition;
    	rateMatrix[59][46] = fFreqs[46] * omega ; 	//TTG->GTG: IS NonSynonymous;
    	rateMatrix[59][52] = fFreqs[52] * k * omega ; 	//TTG->TCG: IS Transition;IS NonSynonymous;
    	rateMatrix[59][55] = fFreqs[55] * omega ; 	//TTG->TGG: IS NonSynonymous;
    	rateMatrix[59][57] = fFreqs[57] * k ; 	//TTG->TTA: IS Transition;
    	rateMatrix[59][58] = fFreqs[58] * omega ; 	//TTG->TTC: IS NonSynonymous;
    	rateMatrix[59][60] = fFreqs[60] * omega ; 	//TTG->TTT: IS NonSynonymous;
    	rateMatrix[60][15] = fFreqs[15] * omega ; 	//TTT->ATT: IS NonSynonymous;
    	rateMatrix[60][31] = fFreqs[31] * k * omega ; 	//TTT->CTT: IS Transition;IS NonSynonymous;
    	rateMatrix[60][47] = fFreqs[47] * omega ; 	//TTT->GTT: IS NonSynonymous;
    	rateMatrix[60][49] = fFreqs[49] * omega ; 	//TTT->TAT: IS NonSynonymous;
    	rateMatrix[60][53] = fFreqs[53] * k * omega ; 	//TTT->TCT: IS Transition;IS NonSynonymous;
    	rateMatrix[60][56] = fFreqs[56] * omega ; 	//TTT->TGT: IS NonSynonymous;
    	rateMatrix[60][57] = fFreqs[57] * omega ; 	//TTT->TTA: IS NonSynonymous;
    	rateMatrix[60][58] = fFreqs[58] * k ; 	//TTT->TTC: IS Transition;
    	rateMatrix[60][59] = fFreqs[59] * omega ; 	//TTT->TTG: IS NonSynonymous;
    	
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }
        
        System.out.println(Arrays.deepToString(rateMatrix));
        
        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * fFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
    }
    
    double[] getCodonProb(){
    	Double[] codonProb = codonProbInput.get().getValues();
    	double[] codonProbResult = new double[codonProb.length];
    	for (int i = 0; i < codonProb.length; i++) {
    		codonProbResult[i] = codonProb[i];
    	}
    	return codonProbResult;
    }
    
    @Override
    public double[] getFrequencies() {
    	double[] codonP = getCodonProb();
    	return codonP;
    }
    
}