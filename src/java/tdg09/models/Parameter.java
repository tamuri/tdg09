package tdg09.models;

/**
 * The interface that all model parameters should extend
 *
 * @author Asif Tamuri
 * @version 1.1
 */
public abstract class Parameter {
    private Object value;
    private boolean optimiseValue;

    public Object getValue() {
        return this.value;
    }

    public void setValue(Object value) {
        this.value = value;
    }

    public boolean isOptimiseValue() {
        return optimiseValue;
    }

    public void setOptimiseValue(boolean optimiseValue) {
        this.optimiseValue = optimiseValue;
    }
}
