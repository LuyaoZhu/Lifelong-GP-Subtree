package yimei.jss.algorithm.lifelongGPSubtree;

import ec.*;
import ec.gp.*;
import yimei.jss.gp.data.DoubleData;

public class ConstNode extends GPNode {

    public double value;

    public ConstNode(double value) {
        this.value = value;
    }

    @Override
    public String toString() {
        return String.format("%.6f", value);
    }

    @Override
    public int expectedChildren() {
        return 0;
    }

    @Override
    public void eval(final EvolutionState state,
                     final int thread,
                     final GPData input,
                     final ADFStack stack,
                     final GPIndividual individual,
                     final Problem problem) {
        ((DoubleData) input).value = value;
    }
}
