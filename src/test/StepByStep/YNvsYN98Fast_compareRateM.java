package test.StepByStep;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN;
import beast.evolution.substitutionmodel.YN98Fast;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class YNvsYN98Fast_compareRateM extends TestCase {

    protected TreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new TreeLikelihood();
    }
    
    static public Alignment getCodonAlignment() throws Exception {
        Sequence human = new Sequence("human", "ATGACGGAATATAAGCTGGTGGTGGTGGGCGCCGGCGGTGTGGGCAAGAGTGCGCTGACCATCCAGCTGATCCAGAACCATTTTGTGGACGAATACGACCCCACTATAGAGGATTCCTACCGGAAGCAGGTGGTCATTGATGGGGAGACGTGCCTGTTGGACATCCTGGATACCGCCGGCCAGGAGGAGTACAGCGCCATGCGGGACCAGTACATGCGCACCGGGGAGGGCTTCCTGTGTGTGTTTGCCATCAACAACACCAAGTCTTTTGAGGACATCCACCAGTACAGGGAGCAGATCAAACGGGTGAAGGACTCGGATGACGTGCCCATGGTGCTGGTGGGGAACAAGTGTGACCTGGCTGCACGCACTGTGGAATCTCGGCAGGCTCAGGACCTCGCCCGAAGCTACGGCATCCCCTACATCGAGACCTCGGCCAAGACCCGGCAGGGAGTGGAGGATGCCTTCTACACGTTGGTGCGTGAGATCCGGCAGCAC");
        Sequence fly = new Sequence("fly", "ATGACGGAATACAAATTGGTTGTTGTTGGTGCGGGAGGCGTTGGCAAATCGGCGTTGACCATCCAACTAATTCAGAATCATTTTGTTGACGAATACGATCCCACAATCGAGGACTCGTACCGAAAGCAAGTGGTCATTGATGGAGAAACCTGCCTTCTGGACATCTTGGATACCGCTGGACAGGAGGAGTACTCGGCTATGCGGGATCAGTATATGCGCACGGGCGAGGGCTTCCTGTTAGTCTTTGCCGTAAATAGTGCAAAATCCTTTGAAGACATCGGCACATACCGCGAGCAGATCAAACGAGTCAAGGATGCCGAGGAGGTGCCAATGGTGCTAGTGGGCAATAAGTGTGACTTGACCACGTGGAACGTTAAAAACGAACAGGCAAGAGAGGTGGCCAAACAATACGGCATTCCATACATTGAGACATCAGCCAAGACGCGCATGGGCGTTGATGATGCATTTTACACACTCGTGCGCGAGATCCGAAAGGAC");
        Sequence chicken = new Sequence("chicken", "ATGACTGAGTATAAGCTTGTTGTCGTTGGAGCTGGTGGTGTGGGCAAGAGCGCCTTGACAATACAGCTCATTCAGAACCACTTTGTGGATGAGTATGACCCTACCATAGAGGATTCCTACAGAAAGCAAGTAGTAATTGATGGGGAAACCTGTCTCTTGGATATTCTTGATACAGCAGGTCAAGAAGAATATAGTGCAATGAGGGACCAATATATGAGAACAGGAGAAGGCTTTCTGTGTGTTTTTGCTATAAACAATACAAAATCTTTTGAAGATATTCACCATTATAGGGAACAAATAAAGAGAGTTAAAGACTCTGAAGATGTCCCAATGGTGCTAGTAGGAAACAAATGTGATTTGCCTTCCAGAACAGTAGATACAAAACAAGCTCAGGATTTAGCAAGAAGTTATGGAATTCCTTTTATTGAAACATCAGCAAAGACAAGACAGGGTGTTGATGATGCCTTCTATACATTAGTTCGAGAAATCAGAAAACAC");
        Sequence liza = new Sequence("liza", "ATGACGGAATATAAGCTGGTTGTGGTAGGAGCTGGAGGTGTTGGCAAGAGCGCACTTACTATTCAGCTCATCCAGAATCACTTTGTGGACGAATATGACCCCACAATTGAGGACTCCTACAGAAAGCAGGTAGTTATTGACGGAGAGACGTGTCTCTTGGACATCCTGGACACTGCAGGTCAAGAGGAGTACAGCGCCATGAGAGATCAGTACATGAGGACAGGGGAGGGCTTTCTCTGTGTCTTTGCCATCAACAACACCAAGTCCTTCGAGGACATTCACCACTATAGAGAACAGATTAAGCGGGTGAAGGACTCTGAGGACGTCCCCATGGTGTTGGTGGGGAACAAGTGTGACCTCCCGTCCCGGACAGTGGACACCAAGCAGGCTCAGGACTTAGCACGCAGCTACGGCATTCCCTTTATTGAGACCTCAGCCAAAACCAGACAGGGCGTTGATGATGCCTTTTACACGTTAGTGCGAGAAATCCGCAAGCAT");
        Sequence tribolium = new Sequence("tribolium", "ATGACTGAATACAAACTAGTAGTAGTTGGAGCAGGTGGTGTCGGCAAATCAGCTTTGACCATACAATTAATCCAAAATCACTTCGTCGACGAATACGACCCTACCATTGAAGACTCCTATCGAAAACAAGTAGTCATCGATGGGGAAACGTGTTTACTGGATATTTTGGATACGGCAGGACAGGAAGAATACAGTGCCATGCGAGACCAGTACATGAGGACAGGGGAAGGTTTCCTTTTGGTTTTCGCCGTTAATTCAGCTAAAAGTTTCGAAGACATTGGAACATACAGGGAACAAATTAAAAGGGTTAAAGATGCCGAAGTCGTACCAATGGTACTCGTAGGAAACAAATGCGACCTCACTTCGTGGGCTGTAGACATGAACCAAGCCAGAGAGGTGGCGCGGCAGTACGGGATCCCGTTCGTGGAGACGTCGGCGAAGACCAGGATGGGTGTGGACGAGGCATTTTACACGTTAGTTAGAGAAATACGTAAGGAC");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", fly, "sequence", chicken, "sequence", liza, "sequence", tribolium,
                "userDataType", codon);

        return data;
    }
    
    static public Tree getTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((human:1.486397618641998,(fly:0.7972225185644256,tribolium:0.7972225185644256):0.6891751000775724):0.9323295960287876,(chicken:1.181379831183816,liza:1.181379831183816):1.2373473834869697):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
	@Test
	public void testRateMatrix_allEqual() throws Exception {			
        RealParameter f = new RealParameter(new Double[]{0.2140365728337246,0.29317049897596853,0.10983236706946545,0.38296056112084154});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        yn98fast.prepareMatricesForTest();
        double[][] yn98fast_rateM = yn98fast.getRateMatrix();
        
        YN yn = new YN();
        yn.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        yn.prepareMatricesForTest();
        double[][] yn_rateM = yn.getRateMatrix();
        
        //compare rate matrix
        for (int i = 0; i < 61; i++) {
        	for (int j = 0; j < 61; j++) {	
                assertEquals(yn98fast_rateM[i][j], yn_rateM[i][j], 1e-10);
        	}
        }
        
        //compare transition probabilities
        double distance = 1.1;
        
        double[] mat_yn98fast = new double[61 * 61];
        double[] mat_yn = new double[61 * 61];
        
        yn98fast.getTransitionProbabilities(null, distance, 0, 1, mat_yn98fast);
        yn.getTransitionProbabilities(null, distance, 0, 1, mat_yn);
        
        for (int i = 0; i < mat_yn98fast.length; i++) {
                assertEquals(mat_yn98fast[i], mat_yn[i], 1e-10);
        }
        
        
        Alignment data = getCodonAlignment();
        Tree tree = getTree(data);
        
        SiteModel siteModel_yn98fast = new SiteModel();
        siteModel_yn98fast.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_yn98fast = newTreeLikelihood();
        likelihood_yn98fast.initByName("data", data, "tree", tree, "siteModel", siteModel_yn98fast);
        double fLogP_yn98fast = likelihood_yn98fast.calculateLogP();
        
        //System.out.println("yn98fast LogP:" + fLogP_yn98fast);
        
        SiteModel siteModel_yn = new SiteModel();
        siteModel_yn.initByName("substModel", yn);      
        TreeLikelihood likelihood_yn = newTreeLikelihood();
        likelihood_yn.initByName("data", data, "tree", tree, "siteModel", siteModel_yn);
        double fLogP_yn = likelihood_yn.calculateLogP();
        
        assertEquals(fLogP_yn98fast, fLogP_yn, 1e-10);
        
	}
}
