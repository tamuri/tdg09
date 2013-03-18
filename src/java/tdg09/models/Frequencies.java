package tdg09.models;

import org.apache.commons.lang.ArrayUtils;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
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
