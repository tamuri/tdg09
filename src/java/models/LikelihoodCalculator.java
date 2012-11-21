package models;

import org.apache.commons.lang.ArrayUtils;

import java.util.*;

import flanagan.math.MinimisationFunction;
import pal.tree.Node;

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
public class LikelihoodCalculator implements MinimisationFunction{
    private Map<String,CladeModel> cladeModels = new HashMap<String, CladeModel>();
    private List<String> orderedCladeKeys = new ArrayList<String>();
    private PhyloAlignment phyloAlignment;
    private String rootCladeModelKey = null;
    private List<Parameter> parameters;
    private double[] ancestralRecon = new double[20];
    private double[][] probMatrix = new double[20][20];
    private double[][] probMatrix0 = new double[20][20];
    private double[][] probMatrix1 = new double[20][20];


    public LikelihoodCalculator(PhyloAlignment phyloAlignment) {
        this.phyloAlignment = phyloAlignment;
    }

    public void setDefaultModel(CladeModel cladeModel) {
        cladeModels.put(Constants.DEFAULT_MODEL_NAME, cladeModel);
        rootCladeModelKey = Constants.DEFAULT_MODEL_NAME;
    }

    public double[] getAncestralRecon() {
        return ancestralRecon;
    }

    public ParametersForMinimisation getParametersForMinimisation() {
        double[] paramsForMinimisation = ArrayUtils.EMPTY_DOUBLE_ARRAY;
        double[] stepsize = ArrayUtils.EMPTY_DOUBLE_ARRAY;
        double[] lowerBounds = ArrayUtils.EMPTY_DOUBLE_ARRAY;
        double[] upperBounds = ArrayUtils.EMPTY_DOUBLE_ARRAY;

        for (Parameter p : parameters) {
            if (p.getClass() == Rate.class && p.isOptimiseValue()) {
                paramsForMinimisation = ArrayUtils.add(paramsForMinimisation, ((Rate) p).getValue());
                stepsize = ArrayUtils.add(stepsize, Constants.ETA_STEPSIZE);
                lowerBounds = ArrayUtils.add(lowerBounds, Constants.PARAM_LOWER_BOUND_ETA);
                upperBounds = ArrayUtils.add(upperBounds, Constants.PARAM_UPPER_BOUND_ETA);
            } else if (p.getClass() == Frequencies.class && p.isOptimiseValue()) {
                double[] freq = aminoFrequenciesToParameters(((Frequencies) p).getValue());
                paramsForMinimisation = ArrayUtils.addAll(paramsForMinimisation, freq);
                for (int i = 0; i < freq.length; i++) {
                    stepsize = ArrayUtils.add(stepsize, Constants.FREQ_STEPSIZE);
                    lowerBounds = ArrayUtils.add(lowerBounds, Constants.PARAM_LOWER_BOUND_FREQ);
                    upperBounds = ArrayUtils.add(upperBounds, Constants.PARAM_UPPER_BOUND_FREQ);
                }
            }
        }
        ParametersForMinimisation pfm = new ParametersForMinimisation(paramsForMinimisation, stepsize, lowerBounds, upperBounds);

        return pfm;
    }

    public double function(double[] parameters) {
        updateParameters(parameters);
        double l = calculate();
		return -1.0 * l; // Flanagan's optimization *minimises*
    }

    public void updateParameters(double[] params) {
        double[] tmp = ArrayUtils.clone(params);

        for (Parameter p : parameters) {
            if (p.getClass() == Rate.class && p.isOptimiseValue()) {
                p.setValue(tmp[0]);
                tmp = ArrayUtils.remove(tmp, 0);
            } else if (p.getClass() == Frequencies.class && p.isOptimiseValue()) {
                p.setValue(parametersToAminoFrequencies(ArrayUtils.subarray(tmp, 0, phyloAlignment.getNumberOfActiveAminoAcids())));
                tmp = ArrayUtils.subarray(tmp, phyloAlignment.getNumberOfActiveAminoAcids(), tmp.length);
            }
        }

        // We've updated the parameters. Notify every clade model to create new models.
        for (String key : cladeModels.keySet()) {
            cladeModels.get(key).updateAminoAcidModel();
        }
    }

    private double calculate() {
        Node root = phyloAlignment.getTree().getRoot();
		double[] conditionals = downTree(root);
		double sum = 0.0;
		for (int i = 0; i < 20; i++) {
            sum += conditionals[i] * cladeModels.get(rootCladeModelKey).getFrequencies().getValue()[i];

		}
        for (int i = 0; i < 20; i++) {
            ancestralRecon[i] = (conditionals[i] * cladeModels.get(rootCladeModelKey).getFrequencies().getValue()[i]) / sum;

        }
		return Math.log(sum);
    }   

    private double[] downTree(Node parent) {
		double[] conditionals = new double[20];
        String parentName, childName;

        parentName = parent.getIdentifier().getName();
		if (parent.isLeaf()) {

			int iAA = phyloAlignment.getAminoAcidForSequence(parentName);

			if ((iAA >= 0) && (iAA < 20)) {
				conditionals[iAA] = 1.0;
			} else {
				for (int iRes = 0; iRes < 20; iRes++) {
					conditionals[iRes] = 1.0;
				}
			}
		} else {
            // this is an internal node - initialise the conditionals array for this node
			for (int i = 0; i < 20; i++) {
				conditionals[i] = 1.0;
			}

            // for each child from this node
			for (int i = 0; i < parent.getChildCount(); i++) {
				Node child = parent.getChild(i);

                // recurse down the tree to get the conditionals for this child
                double[] lowerConditional = downTree(child);

                // if we have a single clade model
                if (cladeModels.size() == 1) {
                    // homogeneous model
                    cladeModels.get(Constants.DEFAULT_MODEL_NAME).getProbabilityMatrix(probMatrix, child.getBranchLength());
                    updateIntrahostConditionals(lowerConditional, conditionals, probMatrix);
                } else {
                    // clade-specific model
                    childName = child.getIdentifier().getName();

                    if (parentName.length() == 0 // the root of the tree is a parent without a label
                            || childName.substring(0,2).equals(parentName.substring(0,2))) {
                        // we're on a branch in a common clade
                        cladeModels.get(childName.substring(0,2)).getProbabilityMatrix(probMatrix, child.getBranchLength());
                        updateIntrahostConditionals(lowerConditional, conditionals, probMatrix);
                    } else {
                        // we're on a branch switching clades
                        cladeModels.get(parentName.substring(0,2)).getProbabilityMatrix(probMatrix0, child.getBranchLength() * 0.5);
                        cladeModels.get(childName.substring(0,2)).getProbabilityMatrix(probMatrix1, child.getBranchLength() * 0.5);
                        updateInterhostConditionals(lowerConditional, conditionals, probMatrix0, probMatrix1);
                    }
                }

			}

		}
		return conditionals;
	}

    private void updateInterhostConditionals(double[] lowerConditional, double[] conditionals, double[][] probMatrix0, double[][] probMatrix1) {
        double[][] probMatrix = new double[20][20];
        for (int i = 0; i < 20; i++) {
            double branchProb = 0.0;
            for (int j = 0; j < 20; j++) {
                for (int k = 0; k < 20; k++) {
                    probMatrix[i][j] += probMatrix0[i][k] * probMatrix1[k][j];
                }
                branchProb += lowerConditional[j] * probMatrix[i][j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private void updateIntrahostConditionals(double[] lowerConditional, double[] conditionals, double[][] probMatrix) {
        for (int i = 0; i < 20; i++) {
            double branchProb = 0.0;
            for (int j = 0; j < 20; j++) {
                branchProb += lowerConditional[j] * probMatrix[i][j];
            }
            conditionals[i] *= branchProb;
        }
    }


    private double[] parametersToAminoFrequencies(double[] parameters) {
        double[] freq = new double[20];

		double sum = 0.0;
        int paramPos = 0;

		for (int i = 0; i < 20; i++) {
			if (phyloAlignment.isAminoAcidActive(i)) {
				freq[i] = Math.exp(parameters[paramPos]);
				paramPos++;
			} else if (i == phyloAlignment.getMostCommonRes()) {
				freq[i] = 1.0;
			} else {
				freq[i] = 0.0;
			}
			sum += freq[i];
		}

		for (int i = 0; i < 20; i++)
			freq[i] /= sum;

        return freq;
    }

    private double[] aminoFrequenciesToParameters(double[] freqs) {
        double params[] = ArrayUtils.EMPTY_DOUBLE_ARRAY;
        double x = 1 / freqs[phyloAlignment.getMostCommonRes()];
        for (int i = 0; i < 20; i++) {
            if (phyloAlignment.isAminoAcidActive(i)) {
                params = ArrayUtils.add(params, Math.log(freqs[i] * x));
            }
        }
        return params;
    }

    public void addCladeModel(String labelPrefix, CladeModel cladeModel) {
        if (cladeModels.size() == 0) {
            rootCladeModelKey = labelPrefix;
        }
        orderedCladeKeys.add(labelPrefix);
        cladeModels.put(labelPrefix, cladeModel);
    }

    public void setParameters(Parameter... parameters) {
        this.parameters = Arrays.asList(parameters);
    }

}
