package beast.evolution.tree;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;

@Description("Logger to report length (sum of all branch lengths) of a tree")
public class TreeLengthLogger  extends CalculationNode implements Loggable, Function {
    public Input<BranchTree> treeInput = new Input<BranchTree>("tree", "tree to report entireLength for.", Validate.REQUIRED);
    
    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) throws Exception {
        final Tree tree = treeInput.get();
        if (getID() == null || getID().matches("\\s*")) {
            out.print(tree.getID() + ".length\t");
        } else {
            out.print(getID() + "\t");
        }
    }

	@Override
	public void log(int nSample, PrintStream out) {
        out.print(treeInput.get().getEntireTreeLength() + "\t");		
	}

	@Override
	public void close(PrintStream out) {
        // nothing to do		
	}
	
	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue() {
		return treeInput.get().getEntireTreeLength();
	}

	@Override
	public double getArrayValue(int iDim) {
		return treeInput.get().getEntireTreeLength();
	}    
    
}
