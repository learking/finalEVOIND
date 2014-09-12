package test.StepByStep;

import org.junit.Test;

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

public class TestTreeLikelihoodYN98Fast extends TestCase {

    static public Alignment getData_1site() throws Exception {
        Sequence seqA = new Sequence("A", "AAA");
        Sequence seqB = new Sequence("B", "TTT");
        Sequence seqC = new Sequence("C", "AAA");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", seqA, "sequence", seqB, "sequence", seqC, "userDataType", codon);

        return data;
    }
    
    static public Alignment getData_2sitesDup() throws Exception {
        Sequence seqA = new Sequence("A", "AAAAAA");
        Sequence seqB = new Sequence("B", "TTTTTT");
        Sequence seqC = new Sequence("C", "AAAAAA");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", seqA, "sequence", seqB, "sequence", seqC, "userDataType", codon);

        return data;
    }
    
    static public Alignment getData_3sites() throws Exception {
        Sequence seqA = new Sequence("A", "AAAAATCCC");
        Sequence seqB = new Sequence("B", "TTTTTTCCC");
        Sequence seqC = new Sequence("C", "ATAAAACCC");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", seqA, "sequence", seqB, "sequence", seqC, "userDataType", codon);

        return data;
    }
    
        static public Tree getTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((A:1.0, C:1.0):1.0, B:2.0);",
                "IsLabelledNewick", true);
        return tree;
    }
    
    static public Tree getTree_unRooted(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((A:1.0, C:1.0):0.5, B:2.5);","adjustTipHeights", false,
                "IsLabelledNewick", true);
        return tree;
    }
   
	public void testTreeLikelihood_3sitesPattern() throws Exception {
        Alignment data_3sites = getData_3sites();
        Tree tree = getTree_unRooted(data_3sites);
        
        RealParameter f = new RealParameter(new Double[]{0.21,0.29,0.12,0.38});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_3sites = new SiteModel();
        siteModel_3sites.initByName("substModel", yn98fast);      
        TreeLikelihoodSimplified likelihood_3sites = new TreeLikelihoodSimplified();
        likelihood_3sites.initByName("data", data_3sites, "tree", tree, "siteModel", siteModel_3sites);
        double fLogP_3sites = likelihood_3sites.calculateLogP();
        
        SiteModel siteModel_3sites_original = new SiteModel();
        siteModel_3sites_original.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_3sites_original = new TreeLikelihood();
        likelihood_3sites_original.initByName("data", data_3sites, "tree", tree, "siteModel", siteModel_3sites);
        double fLogP_3sites_original = likelihood_3sites.calculateLogP();
        
        assertEquals(fLogP_3sites, fLogP_3sites_original, 1e-10);
        //System.out.println("One site dup likelihood:" + fLogP_2sitesDup);
	}
    
    /*
	@Test
	public void testTreeLikelihood_1siteVS2sitesDup() throws Exception {
        Alignment data_1site = getData_1site();
        Tree tree_1site = getTree(data_1site);
        
        RealParameter f = new RealParameter(new Double[]{0.21,0.29,0.12,0.38});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_1site = new SiteModel();
        siteModel_1site.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_1site = new TreeLikelihood();
        likelihood_1site.initByName("data", data_1site, "tree", tree_1site, "siteModel", siteModel_1site);
        double fLogP_1site = likelihood_1site.calculateLogP();
        
        Alignment data_2sitesDup = getData_2sitesDup();
        Tree tree_2sitesDup = getTree(data_2sitesDup);
        SiteModel siteModel_2sitesDup = new SiteModel();
        siteModel_2sitesDup.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_2sitesDup = new TreeLikelihood();
        likelihood_2sitesDup.initByName("data", data_2sitesDup, "tree", tree_2sitesDup, "siteModel", siteModel_2sitesDup);
        double fLogP_2sitesDup = likelihood_2sitesDup.calculateLogP();
        
        assertEquals(fLogP_1site * 2.0 , fLogP_2sitesDup, 1e-10);
        //System.out.println("One site dup likelihood:" + fLogP_2sitesDup);
	}
	
	@Test
	public void testTreeLikelihood_1siteDiffTrees() throws Exception {
        Alignment data_1site = getData_1site();
        Tree tree_1site = getTree(data_1site);
        
        RealParameter f = new RealParameter(new Double[]{0.21,0.29,0.12,0.38});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_1site = new SiteModel();
        siteModel_1site.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_1site = new TreeLikelihood();
        likelihood_1site.initByName("data", data_1site, "tree", tree_1site, "siteModel", siteModel_1site);
        double fLogP_1site = likelihood_1site.calculateLogP();
        
        Tree tree_1siteUnRooted = getTree_unRooted(data_1site);
        SiteModel siteModel_1siteUnRooted = new SiteModel();
        siteModel_1siteUnRooted.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_1siteUnRooted = new TreeLikelihood();
        likelihood_1siteUnRooted.initByName("data", data_1site, "tree", tree_1siteUnRooted, "siteModel", siteModel_1siteUnRooted);
        double fLogP_1siteUnRooted = likelihood_1siteUnRooted.calculateLogP();
        
        assertEquals(fLogP_1site , fLogP_1siteUnRooted, 1e-10);
        //System.out.println("One site dup likelihood:" + fLogP_2sitesDup);
	}
	*/
	
    /*
	@Test
	public void testTreeLikelihood_originalVSsimple() throws Exception {
        Alignment data = getData_1site();
        Tree tree_unRooted = getTree_unRooted(data);
        
        RealParameter f = new RealParameter(new Double[]{0.21,0.29,0.12,0.38});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_original = new SiteModel();
        siteModel_original.initByName("substModel", yn98fast);      
        TreeLikelihood likelihood_original = new TreeLikelihood();
        likelihood_original.initByName("data", data, "tree", tree_unRooted, "siteModel", siteModel_original);
        double fLogP_original = likelihood_original.calculateLogP();
        
        SiteModel siteModel_simple = new SiteModel();
        siteModel_simple.initByName("substModel", yn98fast);      
        TreeLikelihoodSimplified likelihood_simple = new TreeLikelihoodSimplified();
        likelihood_simple.initByName("data", data, "tree", tree_unRooted, "siteModel", siteModel_simple);
        double fLogP_simple = likelihood_simple.calculateLogP();
        
        assertEquals(fLogP_original , fLogP_simple, 1e-10);
        assertEquals(likelihood_simple.getTreeNodeCount(), 5);
        
	}
	*/
    
    
}
