package beast.evolution.operators;

import beast.core.Input;
import beast.evolution.tree.BranchTree;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;


public class CChangeOperator extends TreeOperator {
	
    public Input<Double> sizeInput = new Input<Double>("size", "size of the slide, default 1.0", 1.0);

    // shadows size
    double fSize;

    @Override
    public void initAndValidate() {
        fSize = sizeInput.get();
    }
    
	@Override
	public double proposal() {
        final BranchTree branchTree = (BranchTree) treeInput.get(this);

        double logq = 0;

        Node i;

        // 1. choose a random node avoiding root
        do {
            i = branchTree.getNode(Randomizer.nextInt(branchTree.getNodeCount()));
        } while (i.isRoot());
        
        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldBranchLength = branchTree.getBranchLengthI(i.getNr());
        final double newBranchLength = Math.abs(oldBranchLength + delta);
        

        // 3. set new branch length
        branchTree.setBranchLengthI(i.getNr(), newBranchLength);
        //always zero
		return logq;
	}

    private double getDelta() {
    	    //return Randomizer.nextGaussian() * fSize;
            return Randomizer.nextDouble() * fSize;
    }
    
}
