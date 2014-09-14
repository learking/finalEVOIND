package test.StepByStep;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihoodSimplified;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98Fast;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class AlexLikelihoodVerify_YN98Fast extends TestCase {

    static public Alignment getData_1site() throws Exception {
        Sequence seqA = new Sequence("A", "AAATGTAAA");
        Sequence seqB = new Sequence("B", "TTTGGGTTT");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", seqA, "sequence", seqB, "userDataType", codon);

        return data;
    }
    
    static public Tree getTree_2branch_unRooted(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "(A:1.0,B:2.5);","adjustTipHeights", false,
                "IsLabelledNewick", true);
        return tree;
    }
    
	@Test
	public void testTreeLikelihood_1site() throws Exception {
        Alignment data_1site = getData_1site();
        Tree tree_1site = getTree_2branch_unRooted(data_1site);
        
        RealParameter f = new RealParameter(new Double[]{0.21,0.29,0.12,0.38});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        YN98Fast yn98fast = new YN98Fast();
        yn98fast.initByName("kappa", "2.1", "omega", "0.05", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel siteModel_1site = new SiteModel();
        siteModel_1site.initByName("substModel", yn98fast);      
        TreeLikelihoodSimplified likelihood_1site = new TreeLikelihoodSimplified();
        likelihood_1site.initByName("data", data_1site, "tree", tree_1site, "siteModel", siteModel_1site);
        double fLogP_1site = likelihood_1site.calculateLogP();
        
        //print out rateM (for entire tree) and root dist
        PrintStream ps = System.out; //backup
        
        System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98Fast_rateM_Sep13.txt")));
        System.out.println(Arrays.deepToString(yn98fast.getRateMatrix()));   
        System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/YN98Fast_rootDist_Sep13.txt")));
        System.out.println(Arrays.toString(yn98fast.getFrequencies()));
        
        System.setOut(ps);
        System.out.println("Likelihood:" + fLogP_1site);
        //assertEquals(fLogP_1site * 2.0 , fLogP_2sitesDup, 1e-10);
	}
    
}
