package tdg09.models;

import org.apache.commons.lang.ArrayUtils;

/**
 * Groups together model parameters for the optimiser
 *
 * @author Asif Tamuri
 * @version 1.1
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
