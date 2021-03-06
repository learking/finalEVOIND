package test.beast.evolution.likelihood;

import java.io.FileOutputStream;
import java.io.PrintStream;

import junit.framework.TestCase;

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

public class YN98Fast_compareLikelihood  extends TestCase {
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
    
    static public Tree getTree_1(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((human:1.486397618641998,(fly:0.7972225185644256,tribolium:0.7972225185644256):0.6891751000775724):0.9323295960287876,(chicken:1.181379831183816,liza:1.181379831183816):1.2373473834869697):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
    @Test
    public void testYN98Fast_1() throws Exception {
    	
        Alignment data = getCodonAlignment();
        Tree tree = getTree_1(data);

        RealParameter f = new RealParameter(new Double[]{0.2140365728337246,0.29317049897596853,0.10983236706946545,0.38296056112084154});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98Fast yn98 = new YN98Fast();
        yn98.initByName("kappa", "2.368504291751431", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", yn98);
        
        //on mac
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/debugging_YN98Fast.txt")));
        //System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/debugging_YN98.txt")));
        TreeLikelihoodSimplified likelihoodSimple = new TreeLikelihoodSimplified();
        likelihoodSimple.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double fLogP = likelihoodSimple.calculateLogP();
        System.out.println("TreeLikelihood is:" + fLogP);
    }

    static public Tree getTree_2(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((human:14.86397618641998,(fly:7.972225185644256,tribolium:7.972225185644256):6.891751000775724):9.323295960287876,(chicken:11.81379831183816,liza:11.81379831183816):12.373473834869697):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
    @Test
    public void testYN98Fast_2() throws Exception {
    	
        Alignment data = getCodonAlignment();
        Tree tree = getTree_2(data);

        RealParameter f = new RealParameter(new Double[]{0.2140365728337246,0.29317049897596853,0.10983236706946545,0.38296056112084154});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98Fast yn98 = new YN98Fast();
        yn98.initByName("kappa", "2.368504291751431", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", yn98);
        
        //on mac
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/debugging_YN98Fast.txt")));
        //System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/debugging_YN98.txt")));
        TreeLikelihoodSimplified likelihoodSimple = new TreeLikelihoodSimplified();
        likelihoodSimple.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double fLogP = likelihoodSimple.calculateLogP();
        System.out.println("TreeLikelihood is:" + fLogP);
    }
    
    static public Tree getTree_3(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((((human:0.2875822683162259,liza:0.2875822683162259):0.11712789441759525,chicken:0.40471016273382115):0.1462622073342687,fly:0.5509723700680899):2.9378739494747114E12,tribolium:2.937873949475262E12):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
    @Test
    public void testYN98Fast_3() throws Exception {
    	
        Alignment data = getCodonAlignment();
        Tree tree = getTree_3(data);

        RealParameter f = new RealParameter(new Double[]{0.25,0.25,0.25,0.25});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98Fast yn98 = new YN98Fast();
        yn98.initByName("kappa", "1.0", "omega", "1.0", "nucleoFrequencies", nucleoFrequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", yn98);
        
        //on mac
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/debugging_YN98Fast.txt")));
        //System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/debugging_YN98.txt")));
        TreeLikelihoodSimplified likelihoodSimple = new TreeLikelihoodSimplified();
        likelihoodSimple.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double fLogP = likelihoodSimple.calculateLogP();
        System.out.println("TreeLikelihood is:" + fLogP);
    }
}
