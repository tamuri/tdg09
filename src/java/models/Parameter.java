package models;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
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
