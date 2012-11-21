package models;

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
public class ParametersForMinimisation {
    private double[] parameters;
    private double[] stepSize;
    private double[] lowerBounds;
    private double[] upperBounds;

    public ParametersForMinimisation(double[] parameters, double[] stepSize, double[] lowerBounds, double[] upperBounds) {
        this.parameters = parameters;
        this.stepSize = stepSize;
        this.lowerBounds = lowerBounds;
        this.upperBounds = upperBounds;
    }

    public double[] getParameters() {
        return parameters;
    }

    public double[] getStepSize() {
        return stepSize;
    }

    public double getLowerBounds(int i) {
        return lowerBounds[i];
    }

    public double getUpperBounds(int i) {
        return upperBounds[i];
    }

    @Override
    public String toString() {
        return "models.ParametersForMinimisation{" +
                "parameters=" + (parameters == null ? null : ArrayUtils.toString(parameters)) +
                ", stepSize=" + (stepSize == null ? null : ArrayUtils.toString(stepSize)) +
                ", lowerBounds=" + (lowerBounds == null ? null : ArrayUtils.toString(lowerBounds)) +
                ", upperBounds=" + (upperBounds == null ? null : ArrayUtils.toString(upperBounds)) +
                '}';
    }
}
