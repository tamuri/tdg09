package tdg09.models;

/**
 * A parameter class to hold the rate scaling factor
 *
 * @author Asif Tamuri
 * @version 1.1
 */
public class Rate extends Parameter {
    public Rate(Double rate, boolean optimise) {
        setValue(rate);
        setOptimiseValue(optimise);
    }

    public Double getValue() {
        return (Double) super.getValue();
    }

    public String toString() {
        return "models.Rate = " + getValue();
    }
}
