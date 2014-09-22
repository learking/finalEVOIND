package test.StepByStep;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.likelihood.TreeLikelihoodSimplified;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98Fast;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class PAML2SeqLikelihood extends TestCase {
	
    static public Alignment getData_2seq() throws Exception {
        Sequence seqA = new Sequence("human", "ATGACGGAATATAAGCTGGTGGTGGTGGGCGCCGGCGGTGTGGGCAAGAGTGCGCTGACCATCCAGCTGATCCAGAACCATTTTGTGGACGAATACGACCCCACTATAGAGGATTCCTACCGGAAGCAGGTGGTCATTGATGGGGAGACGTGCCTGTTGGACATCCTGGATACCGCCGGCCAGGAGGAGTACAGCGCCATGCGGGACCAGTACATGCGCACCGGGGAGGGCTTCCTGTGTGTGTTTGCCATCAACAACACCAAGTCTTTTGAGGACATCCACCAGTACAGGGAGCAGATCAAACGGGTGAAGGACTCGGATGACGTGCCCATGGTGCTGGTGGGGAACAAGTGTGACCTGGCTGCACGCACTGTGGAATCTCGGCAGGCTCAGGACCTCGCCCGAAGCTACGGCATCCCCTACATCGAGACCTCGGCCAAGACCCGGCAGGGAGTGGAGGATGCCTTCTACACGTTGGTGCGTGAGATCCGGCAGCAC");
        Sequence seqB = new Sequence("tribolium", "ATGACTGAATACAAACTAGTAGTAGTTGGAGCAGGTGGTGTCGGCAAATCAGCTTTGACCATACAATTAATCCAAAATCACTTCGTCGACGAATACGACCCTACCATTGAAGACTCCTATCGAAAACAAGTAGTCATCGATGGGGAAACGTGTTTACTGGATATTTTGGATACGGCAGGACAGGAAGAATACAGTGCCATGCGAGACCAGTACATGAGGACAGGGGAAGGTTTCCTTTTGGTTTTCGCCGTTAATTCAGCTAAAAGTTTCGAAGACATTGGAACATACAGGGAACAAATTAAAAGGGTTAAAGATGCCGAAGTCGTACCAATGGTACTCGTAGGAAACAAATGCGACCTCACTTCGTGGGCTGTAGACATGAACCAAGCCAGAGAGGTGGCGCGGCAGTACGGGATCCCGTTCGTGGAGACGTCGGCGAAGACCAGGATGGGTGTGGACGAGGCATTTTACACGTTAGTTAGAGAAATACGTAAGGAC");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", seqA, "sequence", seqB, "userDataType", codon);

        return data;
    }
    
    static public Tree getTree_2seq(Alignment data) throws Exception {
    	TreeParser tree = new TreeParser();
    	tree.initByName("taxa", data,
            "newick", "(human:35.311075, tribolium:35.362197);",
            "IsLabelledNewick", true, "adjustTipHeights", false);
    	return tree;
    }
    
    static public Tree getTree_2seq_long(Alignment data) throws Exception {
    	TreeParser tree = new TreeParser();
    	tree.initByName("taxa", data,
            "newick", "(human:826.2991454953574, tribolium:4645.405319290922);",
            "IsLabelledNewick", true, "adjustTipHeights", false);
    	return tree;
    }
    
	public void testTreeLikelihood_2seq() throws Exception {
        Alignment data_2seq = getData_2seq();
        Tree tree = getTree_2seq(data_2seq);
        
        RealParameter f = new RealParameter(new Double[]{0.25,0.25,0.25,0.25});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "0.94325", "omega", "0.0012", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_2seq = new SiteModel();
        siteModel_2seq.initByName("substModel", yn98fast);      
        TreeLikelihoodSimplified likelihood_2seq = new TreeLikelihoodSimplified();
        likelihood_2seq.initByName("data", data_2seq, "tree", tree, "siteModel", siteModel_2seq);
        double fLogP_2seq = likelihood_2seq.calculateLogP();
        
        /*
        SiteModel siteModel_3sites_original = new SiteModel();
        siteModel_3sites_original.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_3sites_original = new TreeLikelihood();
        likelihood_3sites_original.initByName("data", data_3sites, "tree", tree, "siteModel", siteModel_3sites);
        double fLogP_3sites_original = likelihood_3sites.calculateLogP();
        */
        
        //assertEquals(fLogP_3sites, fLogP_3sites_original, 1e-10);
        System.out.println("likelihood:" + fLogP_2seq);
	}
    
	public void testTreeLikelihood_2seq_longTree() throws Exception {
        Alignment data_2seq = getData_2seq();
        Tree longtree = getTree_2seq_long(data_2seq);
        
        RealParameter f = new RealParameter(new Double[]{0.25,0.25,0.25,0.25});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "1.2431800243209412", "omega", "0.000018564105316624207", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_2seq = new SiteModel();
        siteModel_2seq.initByName("substModel", yn98fast);      
        TreeLikelihoodSimplified likelihood_2seq = new TreeLikelihoodSimplified();
        likelihood_2seq.initByName("data", data_2seq, "tree", longtree, "siteModel", siteModel_2seq);
        double fLogP_2seq = likelihood_2seq.calculateLogP();

        System.out.println("likelihood:" + fLogP_2seq);
	}
}
