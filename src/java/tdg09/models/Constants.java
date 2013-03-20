package tdg09.models;

/**
 * Program defaults
 *
 * @author Asif Tamuri
 * @author Richard Goldstein
 * @version 1.1
 */

public interface Constants {
    static char[] aaNames = "ARNDCQEGHILKMFPSTWYV".toCharArray();
    static double DEFAULT_TOL = 5.0e-6; // tolerance used in determining the convergence of the minimization
    static boolean DELETE_SINGLES = true; // whether to remove residues that occur only once
    static double DEFAULT_ETA = 1.0; // default rate
    static double ETA_STEPSIZE = 0.5;
    static double FREQ_STEPSIZE = 1.0;
    static double PARAM_UPPER_BOUND_ETA = 100.0;
    static double PARAM_LOWER_BOUND_ETA = 1.0E-6;
    static double PARAM_UPPER_BOUND_FREQ = 10.0;
    static double PARAM_LOWER_BOUND_FREQ = -500.0;
    static String DEFAULT_MODEL_NAME = "DEFAULT_MODEL_NAME";
}

