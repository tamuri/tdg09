package tdg09;

import flanagan.math.Minimisation;
import org.apache.commons.lang.ArrayUtils;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg09.models.*;

import java.util.concurrent.Callable;

/**
* Author: Asif Tamuri (atamuri@ebi.ac.uk)
* Date: 18/03/2013 21:54
*/
public class SiteAnalyser implements Callable<Void> {
    final Alignment alignment;
    final Tree tree;
    final int site;
    final static String LOCK = "LOCK";

    SiteAnalyser(Alignment alignment, Tree tree, int site) {
        this.alignment = alignment;
        this.tree = tree;
        this.site = site;
    }

        public Void call() throws Exception {
            PhyloAlignment phyloAlignment = new PhyloAlignment(alignment, tree, site, Constants.DELETE_SINGLES);

            if (phyloAlignment.getNumberOfActiveAminoAcids() == 0) {
                System.out.println("No active residues");
                return null;
            }

            Minimisation min = new Minimisation();

            // HOMOGENEOUS (model B)
            LikelihoodCalculator modelWAGssF = new LikelihoodCalculator(phyloAlignment);
            // Rate nuB = new Rate(nuA.getValue(), true);
            Rate nuB = new Rate(Constants.DEFAULT_ETA, true);
            Frequencies piB = new Frequencies(phyloAlignment.getEquilibriumFrequencies(), true);
            modelWAGssF.setParameters(nuB, piB);
            modelWAGssF.setDefaultModel(new CladeModel(nuB, piB, true));
            runMinimisation(modelWAGssF, min);
            printResult("WAG+ssF", min.getParamValues().length, -min.getMinimum(), nuB, piB);


            // HOST-SPECIFIC FREQUENCIES, SHARED RATE SCALING FACTOR (model C)
            LikelihoodCalculator modelWAGlssF = new LikelihoodCalculator(phyloAlignment);
            Rate nuC = new Rate(nuB.getValue(), true);
            Frequencies piC1 = new Frequencies(piB.getValue(), true);
            Frequencies piC2 = new Frequencies(piB.getValue(), true);
            modelWAGlssF.setParameters(nuC, piC1, piC2);
            modelWAGlssF.addCladeModel("Av", new CladeModel(nuC, piC1, true));
            modelWAGlssF.addCladeModel("Hu", new CladeModel(nuC, piC2, true));
            runMinimisation(modelWAGlssF, min);
            printResult("WAG+lssF", min.getParamValues().length, -min.getMinimum(), nuC, piC1, piC2);
            System.out.printf("\nAR\t%s\n", ArrayUtils.toString(modelWAGlssF.getAncestralRecon()));
            return null;
        }


    private void runMinimisation(LikelihoodCalculator model, Minimisation min) {
        ParametersForMinimisation pfm = model.getParametersForMinimisation();
        addConstraints(min, pfm);
        min.nelderMead(model, pfm.getParameters(), pfm.getStepSize(), Constants.DEFAULT_TOL);
        model.updateParameters(min.getParamValues());
    }

    private void printResult(String modelName, int parameterCount, double logLikelihood, Parameter... parameters) {
        printModelTitle(modelName);
        printParameters(parameterCount, parameters);
        printLogLikelihood(logLikelihood);
    }

    private void printLogLikelihood(double minimum) {
        System.out.printf("Log-likelihood = %s\n", minimum);
    }

    private void printModelTitle(String title) {
        System.out.printf("\n---%s---\n", title);
    }

    private void printParameters(int parameterCount, Parameter... params) {
        System.out.printf("Parameters     = %s\n", parameterCount);
        for (Parameter p : params) {
            if (p.getClass() == Rate.class && p.isOptimiseValue()) {
                System.out.printf("Rate           = %s\n", p.getValue());
            } else if (p.getClass() == Frequencies.class && p.isOptimiseValue()) {
                System.out.printf("Frequency      = ");
                Frequencies f = (Frequencies) p;
                for (int i = 0; i < f.getValue().length; i++) {
                    if (f.getValue()[i] > 0) {
                        System.out.printf("%s:%s ", Constants.aaNames[i], f.getValue()[i]);
                    }
                }
                System.out.println();
            }
        }
    }

    private void addConstraints(Minimisation min, ParametersForMinimisation pfm) {
        min.removeConstraints();
        for (int i = 0; i < pfm.getParameters().length; i++) {
            min.addConstraint(i, -1, pfm.getLowerBounds(i));
            min.addConstraint(i, 1, pfm.getUpperBounds(i));
        }
    }
}
