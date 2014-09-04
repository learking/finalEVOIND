package test.beast.evolution.likelihood;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

/**
 * This test aims to check whether YN98 is compatible with the rest of BEAST2 program (before spending time to 
 * compile my plugin, construct .xml input file and test using real-world sequence alignment)
 * *
 */

public class YN98_1_Test extends TestCase {
	
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
    
    static public Tree getTree_sample0(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "(((human:1.6667276805951001,liza:1.6667276805951001):1.7137605617323786,(chicken:1.6780683322532925,tribolium:1.6780683322532925):1.7024199100741861):93.03370812766852,fly:96.41419636999599):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
    @Test
    public void testYN98_sample0() throws Exception {
    	
        Alignment data = getCodonAlignment();
        Tree tree = getTree_sample0(data);

        RealParameter f = new RealParameter(new Double[]{0.21232809709968758,0.292956974771638,0.11086021329597975,0.38385471483269473});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98 yn98 = new YN98();
        yn98.initByName("kappa", "2.3046949838646555", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", yn98);

        TreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);

        double fLogP = likelihood.calculateLogP();
        assertEquals(fLogP, -838.9757864452677, BEASTTestCase.PRECISION);
    }
    
}