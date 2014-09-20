package beast.evolution.tree;

import java.io.PrintStream;

import beast.core.StateNodeInitialiser;
import beast.core.Input.Validate;
import beast.util.TreeParser;

public class BranchTree extends Tree {
	
	protected double[] branchLengths;
	protected double[] storedBranchLengths;
	
    //Always need an initial tree (From TreeParser) to start with
    public BranchTree() {
    	m_initial.setRule(Validate.REQUIRED);
    }
	
    @Override
    public void initAndValidate() throws Exception {
        if (m_initial.get() != null  && !(this instanceof StateNodeInitialiser)) {
        	final Tree other = m_initial.get();
        	root = other.root.copy();
        	nodeCount = other.nodeCount;
        	internalNodeCount = other.internalNodeCount;
        	leafNodeCount = other.leafNodeCount;
        	
        	//init "tree as array" representation
            if (nodeCount >= 0) {
                initArrays();
            }
            
            //init branchLengths
            branchLengths = new double[nodeCount];
            storedBranchLengths = new double[nodeCount];

			for(int i=0; i<nodeCount; i++) {
				branchLengths[i] = m_nodes[i].getLength();
				storedBranchLengths[i] = m_nodes[i].getLength();
			}
        	
        }else{
        	throw new RuntimeException("Always need an initial tree (From TreeParser) to start with");        	
        }
    }
    
    /**
     * StateNode implementation
     */
    @Override
    public String toString() {
    	//part I: tree
    	return toShortNewick(root, true);
    }
    
    /**
     * @return beast.tree in Newick format, with length and meta data
     *         information. Unlike toNewick(), here Nodes are numbered, instead of
     *         using the node labels.
     *         If there are internal nodes with non-null IDs then their numbers are also printed.
     *         Also, all internal nodes are labelled if bPrintInternalNodeNumbers
     *         is set true. This is useful for example when storing a State to file
     *         so that it can be restored.
     */
    public String toShortNewick(Node node, final boolean bPrintInternalNodeNumbers) {
        final StringBuilder buf = new StringBuilder();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toShortNewick(node.getLeft(), bPrintInternalNodeNumbers));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toShortNewick(node.getRight(), bPrintInternalNodeNumbers));
            }
            buf.append(")");
            if (getID() != null) {
                buf.append(node.getNr());
            } else if (bPrintInternalNodeNumbers) {
                buf.append(node.getNr());
            }

        } else {
            buf.append(node.getNr());
        }
        //buf.append(node.getNewickMetaData());
        buf.append(":").append(branchLengths[node.getNr()]);
        return buf.toString();
    }
    
    public void init(PrintStream out) throws Exception {
    }

    public void log(int nSample, PrintStream out) {
        Tree tree = (Tree) getCurrent();
        out.print("tree STATE_" + nSample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String sNewick = toSortedNewick(tree.getRoot(), dummy);
        out.print(sNewick);
        out.print(";");
    }

    public String toSortedNewick(Node node, int[] iMaxNodeInClade) {
        StringBuilder buf = new StringBuilder();
        if (node.getLeft() != null) {
            buf.append("(");
            String sChild1 = toSortedNewick(node.getLeft(), iMaxNodeInClade);
            int iChild1 = iMaxNodeInClade[0];
            if (node.getRight() != null) {
                String sChild2 = toSortedNewick(node.getRight(), iMaxNodeInClade);
                int iChild2 = iMaxNodeInClade[0];
                if (iChild1 > iChild2) {
                    buf.append(sChild2);
                    buf.append(",");
                    buf.append(sChild1);
                } else {
                    buf.append(sChild1);
                    buf.append(",");
                    buf.append(sChild2);
                    iMaxNodeInClade[0] = iChild1;
                }
            } else {
                buf.append(sChild1);
            }
            buf.append(")");
            if (getID() != null) {
                buf.append(node.getNr()+1);
            }
        } else {
            iMaxNodeInClade[0] = node.getNr();
            buf.append(node.getNr() + 1);
        }
        
        buf.append(":").append(branchLengths[node.getNr()]);
        return buf.toString();
    }
    
    /**
     * @see beast.core.Loggable *
     */
    public void close(PrintStream out) {
        out.print("End;");
    }
    
    /**
     * reconstruct PathTree from XML fragment in the form of a DOM node *
     */
    
    @Override
    public void fromXML(final org.w3c.dom.Node node) {
        final String sNewick = node.getTextContent();
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            setRoot(parser.parseNewick(sNewick));
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        initArrays();
        
        //init branchLengths
        branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

		for(int i=0; i<nodeCount; i++) {
			branchLengths[i] = m_nodes[i].getLength();
			storedBranchLengths[i] = m_nodes[i].getLength();
		}
    }
    
    /**
     * StateNode implementation 
     * @throws Exception *
     */
    @Override
    protected void store() {
    	super.store();
    	
    	// store branche lengths
    	double[] tmp = new double[nodeCount];
		for(int i=0; i<nodeCount; i++) {
			tmp[i] = branchLengths[i];
		}
		storedBranchLengths = tmp;
    }
    
    @Override
    public void restore() {
    	super.restore();
    	// restore branche lengths
    	// store branche lengths
    	double[] tmp = new double[nodeCount];
		for(int i=0; i<nodeCount; i++) {
			tmp[i] = storedBranchLengths[i];
		}
		branchLengths = tmp;
    }
    
    /*
     * getters
     * */
    public double[] getCurrentBranchLengths(){
    	return branchLengths;
    }
    
    /******
     * for testing purpose
     * *****/
    //tmp (for testing)
    public void fakeStore(){
    	store();
    }
    
}
