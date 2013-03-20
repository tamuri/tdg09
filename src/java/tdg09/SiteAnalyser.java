package tdg09;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import flanagan.math.Minimisation;
import pal.alignment.Alignment;
import pal.tree.Tree;
import tdg09.models.*;

import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

/**
 * Given an alignment, tree and location, runs the homogeneous and non-homogeneous models. Is Callable and run in an
 * ExecutionService by the Analyse class.
 *
 * @author Asif Tamuri
 * @version 1.1
 */
public class SiteAnalyser implements Callable<Result> {
    final Alignment alignment;
    final Tree tree;
    final int site;
    final List<String> groups;

    SiteAnalyser(Alignment alignment, Tree tree, int site, List<String> groups) {
        this.alignment = alignment;
        this.tree = tree;
        this.site = site;
        this.groups = groups;
    }

    public Result call() throws Exception {

        PhyloAlignment phyloAlignment = new PhyloAlignment(alignment, tree, site, Constants.DELETE_SINGLES);

        if (phyloAlignment.getNumberOfActiveAminoAcids() == 0) {
            return new Result(site, phyloAlignment.getResidues(), null);
        }

        List<Result.Model> models = Lists.newArrayList();

        Minimisation min = new Minimisation();
        min.suppressNoConvergenceMessage();

        // HOMOGENEOUS (model B)
        LikelihoodCalculator modelWAGssF = new LikelihoodCalculator(phyloAlignment);
        Rate nuB = new Rate(Constants.DEFAULT_ETA, true);
        Frequencies piB = new Frequencies(phyloAlignment.getEquilibriumFrequencies(), true);
        modelWAGssF.setParameters(nuB, piB);
        modelWAGssF.setDefaultModel(new CladeModel(nuB, piB, true));
        runMinimisation(modelWAGssF, min);
        models.add(getModel("WAG+ssF", min.getConvStatus(), min.getParamValues().length, -min.getMinimum(), Lists.newArrayList(nuB, piB)));

        // HOST-SPECIFIC FREQUENCIES, SHARED RATE SCALING FACTOR (model C)
        LikelihoodCalculator modelWAGlssF = new LikelihoodCalculator(phyloAlignment);
        Rate nuC = new Rate(nuB.getValue(), true);
        List<Parameter> groupParameters = Lists.newArrayList();
        groupParameters.add(nuC);
        for (String group : groups) {
            Frequencies f = new Frequencies(piB.getValue(), true);
            groupParameters.add(f);
            modelWAGlssF.addCladeModel(group, new CladeModel(nuC, f, true));
        }
        modelWAGlssF.setParameters(groupParameters.toArray(new Parameter[groupParameters.size()]));
        runMinimisation(modelWAGlssF, min);
        models.add(getModel("WAG+lssF", min.getConvStatus(), min.getParamValues().length, -min.getMinimum(), groupParameters));

        return new Result(site, phyloAlignment.getResidues(), models);
    }

    private Result.Model getModel(String name, boolean converged, int n, double logLikelihood, List<Parameter> params) {

        double rate = Double.NaN;
        List<Map<Character, Double>> frequencies = Lists.newArrayList();

        for (Parameter p : params) {
            if (p.getClass() == Rate.class && p.isOptimiseValue()) {
                rate = ((Rate) p).getValue();
            } else if (p.getClass() == Frequencies.class && p.isOptimiseValue()) {
                Frequencies f = (Frequencies) p;
                Map<Character, Double> thisF = Maps.newLinkedHashMap();
                for (int i = 0; i < f.getValue().length; i++) {
                    if (f.getValue()[i] > 0) {
                        thisF.put(Constants.aaNames[i], f.getValue()[i]);
                    }
                }
                frequencies.add(thisF);
            }
        }


        return new Result.Model(name, n, rate, logLikelihood, frequencies, converged);
    }

    private void runMinimisation(LikelihoodCalculator model, Minimisation min) {
        ParametersForMinimisation pfm = model.getParametersForMinimisation();
        addConstraints(min, pfm);
        min.nelderMead(model, pfm.getParameters(), pfm.getStepSize(), Constants.DEFAULT_TOL);
        model.updateParameters(min.getParamValues());
    }

    private void addConstraints(Minimisation min, ParametersForMinimisation pfm) {
        min.removeConstraints();
        for (int i = 0; i < pfm.getParameters().length; i++) {
            min.addConstraint(i, -1, pfm.getLowerBounds(i));
            min.addConstraint(i, 1, pfm.getUpperBounds(i));
        }
    }


}
