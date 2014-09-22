package test.beast.evolution.substmodel;

import java.io.FileOutputStream;
import java.io.PrintStream;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.likelihood.TreeLikelihoodSimplified;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.OneStruct;
import beast.evolution.substitutionmodel.OneStructFast;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class OneStructVsOneStructFast_rateMTest extends TestCase {
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
    
    static public Tree getTree_sample500000(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((human:1.486397618641998,(fly:0.7972225185644256,tribolium:0.7972225185644256):0.6891751000775724):0.9323295960287876,(chicken:1.181379831183816,liza:1.181379831183816):1.2373473834869697):0.0;",
                "IsLabelledNewick", true);
        return tree;
    }
    
    @Test
    public void testOneStructEachSiteLikelihood() throws Exception {
        
        Alignment data = getCodonAlignment();
        Tree tree = getTree_sample500000(data);
        
        //use FilteredAlignment to create a sub data for each site
        FilteredAlignment site1 = new FilteredAlignment();
        site1.initByName("data", data, "filter", "1");
        
        RealParameter f = new RealParameter(new Double[]{0.2, 0.1, 0.3, 0.4});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        
        RealParameter codonProb1 = new RealParameter(new Double[]{
        		0.031350262,0.016678685,0.024086503,0.013927975,0.025667266,0.019513173,0.001036515,
        		0.019661785,0.004234123,0.021709737,0.023116887,0.008893605,0.023613895,0.016321916,
        		0.021248437,0.034586834,0.013469791,0.013974397,0.019040362,0.023127154,0.012276834,
        		0.016818336,0.013005605,0.016724293,0.031102112,0.025201710,0.002783288,0.028634640,
        		0.011845951,0.035306332,0.007706424,0.016484928,0.001266775,0.015199627,0.024995155,
        		0.027195620,0.022472747,0.008187028,0.029989578,0.003142386,0.027191536,0.010992749,
        		0.010295569,0.005584570,0.029299030,0.006562422,0.004537246,0.009528992,0.002336504,
        		0.002458360,0.004808568,0.014593801,0.006875164,0.009996637,0.019877895,0.012218437,
        		0.026721705,0.015533992,0.020767555,0.021984981,0.012235612});
        
        OneStruct site1_oneStruct = new OneStruct();
        
        site1_oneStruct.initByName("kappa", "1.5", "omega", "0.9", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb1);
        
        SiteModel site1_siteModel = new SiteModel();
        site1_siteModel.initByName("substModel", site1_oneStruct);
        
        TreeLikelihoodSimplified site1_likelihood = new TreeLikelihoodSimplified();
        site1_likelihood.initByName("data", site1, "tree", tree, "siteModel", site1_siteModel);
        
        //System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/debugging_OneStruct.txt")));
        //on mac
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/debugging_OneStruct.txt")));        
        site1_likelihood.calculateLogP();
    }
    
    @Test
    public void testOneStructFastEachSiteLikelihood() throws Exception {
        
        Alignment data = getCodonAlignment();
        Tree tree = getTree_sample500000(data);
        
        //use FilteredAlignment to create a sub data for each site
        FilteredAlignment site1 = new FilteredAlignment();
        site1.initByName("data", data, "filter", "1");
        
        RealParameter f = new RealParameter(new Double[]{0.2, 0.1, 0.3, 0.4});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        
        RealParameter codonProb1 = new RealParameter(new Double[]{
        		0.031350262,0.016678685,0.024086503,0.013927975,0.025667266,0.019513173,0.001036515,
        		0.019661785,0.004234123,0.021709737,0.023116887,0.008893605,0.023613895,0.016321916,
        		0.021248437,0.034586834,0.013469791,0.013974397,0.019040362,0.023127154,0.012276834,
        		0.016818336,0.013005605,0.016724293,0.031102112,0.025201710,0.002783288,0.028634640,
        		0.011845951,0.035306332,0.007706424,0.016484928,0.001266775,0.015199627,0.024995155,
        		0.027195620,0.022472747,0.008187028,0.029989578,0.003142386,0.027191536,0.010992749,
        		0.010295569,0.005584570,0.029299030,0.006562422,0.004537246,0.009528992,0.002336504,
        		0.002458360,0.004808568,0.014593801,0.006875164,0.009996637,0.019877895,0.012218437,
        		0.026721705,0.015533992,0.020767555,0.021984981,0.012235612});
        
        OneStructFast site1_oneStruct = new OneStructFast();
        
        site1_oneStruct.initByName("kappa", "1.5", "omega", "0.9", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb1);
        
        SiteModel site1_siteModel = new SiteModel();
        site1_siteModel.initByName("substModel", site1_oneStruct);
        
        TreeLikelihoodSimplified site1_likelihood = new TreeLikelihoodSimplified();
        site1_likelihood.initByName("data", site1, "tree", tree, "siteModel", site1_siteModel);
        
        //System.setOut(new PrintStream(new FileOutputStream("/home/kuangyu/Desktop/debugging_OneStruct.txt")));
        //on mac
        //System.setOut(new PrintStream(new FileOutputStream("/Users/kwang2/Desktop/debugging_OneStructFast.txt")));        
        site1_likelihood.calculateLogP();
    }
}