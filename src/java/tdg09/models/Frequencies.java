package tdg09.models;

import org.apache.commons.lang.ArrayUtils;

/**
 * A parameter class to hold an array of amino acid frequencies.
 *
 * @author Asif Tamuri
 * @version 1.1
 */
public class Frequencies extends Parameter {
    public Frequencies(double[] frequencies, boolean optimise) {
        setValue(frequencies);
        setOptimiseValue(optimise);
    }

    public double[] getValue() {
        return (double[]) super.getValue();
    }

    public String toString() {
        return "models.Frequencies = " + ArrayUtils.toString(getValue());
    }
}
