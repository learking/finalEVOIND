package beast.evolution.operators;

import beast.core.Input;
import beast.evolution.tree.BranchTree;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

public class BranchScaler extends TreeOperator {
    public final Input<Double> scaleFactorInput = new Input<Double>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);
    public Input<Double> scaleUpperLimit = new Input<Double>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    public Input<Double> scaleLowerLimit = new Input<Double>("lower", "Lower limit of scale factor", 1e-8);
    
    /**
     * shadows input *
     */
    private double m_fScaleFactor;

    //private double upper, lower;
    
    @Override
    public void initAndValidate() throws Exception {
        m_fScaleFactor = scaleFactorInput.get();
        //upper = scaleUpperLimit.get();
        //lower = scaleLowerLimit.get();
    }
    
    protected double getScaler() {
        return (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
    }
    
	@Override
	public double proposal() {
        final BranchTree branchTree = (BranchTree) treeInput.get(this);
        Node i;

        // 1. choose a random node avoiding root
        do {
            i = branchTree.getNode(Randomizer.nextInt(branchTree.getNodeCount()));
        } while (i.isRoot());
        
        final double scale = getScaler();
        final double oldBranchLength = branchTree.getBranchLengthI(i.getNr());
        final double newBranchLength = oldBranchLength*scale;
        
        // 3. set new branch length
        branchTree.setBranchLengthI(i.getNr(), newBranchLength);
        
		return -Math.log(scale);
	}

}
