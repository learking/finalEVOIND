package test.beast.evolution.tree;

import java.io.File;

import javax.xml.parsers.DocumentBuilderFactory;

import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.tree.BranchTree;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class BranchTreeBasicTest extends TestCase {

	Alignment data;
	Tree tree;
	
    static public Alignment getSimpleAlignment() throws Exception {
        Sequence Human = new Sequence("Human", "CCA");
        Sequence Chimpanzee = new Sequence("Chimpanzee", "CCC");
        Sequence Gibbon = new Sequence("Gibbon", "CTA");
        Sequence Gorilla = new Sequence("Gorilla", "CCA");

        Alignment data = new Alignment();
        data.initByName("sequence", Human, "sequence", Chimpanzee, "sequence", Gorilla, "sequence", Gibbon,
                "dataType", "nucleotide"
        );
        return data;
    }
	
    static public Tree getSimpleTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "(((Human: 0.1, Chimpanzee: 0.05): 0.3, Gibbon:0.5):0.2 , Gorilla: 0.7)", "IsLabelledNewick", true, "adjustTipHeights", false);
        return tree;
    }
	
    @Override
    protected void setUp() throws Exception {
    	super.setUp();
        data = getSimpleAlignment();
        tree = getSimpleTree(data);
    }
    
    @Test
    public void testToString() throws Exception{
    	String treeStr = tree.toString();
    	tree.log(0, System.out);
    	
    	//System.out.println(treeStr);
    	
		BranchTree branchTree = new BranchTree();
		branchTree.initByName("initial", tree);
		
		String branchTreeStr = branchTree.toString();
		branchTree.log(0, System.out);
		
		assertEquals(treeStr, branchTreeStr);
    }
    
    @Test
    public void testFromXML() throws Exception{
    	//prepare parser
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        
        //under linux
        //String stateFileName = "/home/kuangyu/Dropbox/research/meeting_memo/Jan_2014/pathTree.xml.state";
        //under mac
        //String stateFileName = "/Users/kwang2/Desktop/Dropbox/research/meeting_memo/Jan_2014/pathTree.xml.state";
        
        //real state file
        String stateFileName = "/Users/kwang2/Desktop/Dropbox/research/graduationProcess/papers/evoproteinAnalysis/utilityScripts/branchTreeState/branchTree.xml.state";
        
        //String stateFileName = "/home/kuangyu/workspace/beast2/realData_twoStruct_gapped.xml.state";
        Document doc = factory.newDocumentBuilder().parse(new File(stateFileName));
        doc.normalize();
        //get node that belongs to PathTree
        final NodeList nodes = doc.getElementsByTagName("*");
        final Node topNode = nodes.item(0);
        final NodeList children = topNode.getChildNodes();  
        final Node child = children.item(1);

        //verify its ID:
        //final String sID = child.getAttributes().getNamedItem("id").getNodeValue();
        //System.out.println(sID);

        //show its text content:
        BranchTree branchTree = new BranchTree();
        branchTree.fromXML(child);

		System.out.print(branchTree.toString());
    }
    
}
