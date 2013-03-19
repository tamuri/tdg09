package tdg09;

import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 19/03/2013 10:56
 */
public class Model {
    public final String name;
    public final int parameters;
    public final double rate;
    public final double lnL;
    public final List<Map<Character,Double>> frequencies;
    public final boolean converged;

    public Model(String name, int parameters, double rate, double lnL, List<Map<Character, Double>> frequencies, boolean converged) {
        this.name = name;
        this.parameters = parameters;
        this.rate = rate;
        this.lnL = lnL;
        this.frequencies = frequencies;
        this.converged = converged;
    }
}
