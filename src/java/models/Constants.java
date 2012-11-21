package models;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @author Richard Goldstein <rgoldst@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
 */

public interface Constants {
    static char[] aaNames = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
    static double DEFAULT_TOL = 1.0E-8; // tolerance used in determining the convergence of the minimization
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

