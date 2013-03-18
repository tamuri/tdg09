package tdg09.models;

import pal.substmodel.RateMatrix;
import pal.substmodel.WAG;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
 */
public class CladeModel {
    private Rate rate;
    private Frequencies frequencies;
    RateMatrix substitutionModel;
    int size = 0;
    int[] aminoAcidPos;
    double[][] Pt;
    boolean collapse;

    public CladeModel(Rate rate, Frequencies frequencies, boolean collapse) {
        this.rate = rate;
        this.frequencies = frequencies;
        this.collapse = collapse;

        if (this.collapse) {
            double[] f = frequencies.getValue();
            int collapsed_f_size = 0;
            for (int i = 0; i < 20; i++) {
                if (f[i] > 0) {
                    collapsed_f_size++;
                }
            }
            size = collapsed_f_size;
            Pt = new double[size][size];

            aminoAcidPos = new int[size];
            int new_f_pos = 0;
            for (int i = 0; i < 20; i++) {
                if (frequencies.getValue()[i] > 0) {
                    aminoAcidPos[new_f_pos] = i;
                    new_f_pos++;
                }
            }
        }
    }

    public void updateAminoAcidModel() {
        if (this.collapse) {
            double[] new_f = new double[size];
            int new_f_pos = 0;
            for (int i = 0; i < 20; i++) {
                if (frequencies.getValue()[i] > 0) {
                    new_f[new_f_pos] = frequencies.getValue()[i];
                    aminoAcidPos[new_f_pos] = i;
                    new_f_pos++;
                }
            }
            substitutionModel = new CollapsedAminoAcidModel(new_f, aminoAcidPos);
        } else {
            substitutionModel = new WAG(frequencies.getValue());
        }
    }

    public void getProbabilityMatrix(double[][] probMatrix, double branchLength) {
        substitutionModel.setDistance(rate.getValue() * branchLength);

        // If the requested probMatrix has the same dimension as our clade model's rate matrix (i.e. when size = 20)
        if (!this.collapse || size == probMatrix.length) {
            // Return the transition probabilities directly from the model
            substitutionModel.getTransitionProbabilities(probMatrix);
        } else {
            // Otherwise the model will return a rate matrix of a reduced size (e.g. for 4 amino acids only)
            substitutionModel.getTransitionProbabilities(Pt);
            // Fill the full 20x20 matrix with entries for the amino acids of interest
            for (int i = 0; i < Pt.length; i++) {
                for (int j = 0; j < Pt[0].length; j++) {
                    probMatrix[aminoAcidPos[i]][aminoAcidPos[j]] = Pt[i][j];
                }
            }
        }
    }

    public Frequencies getFrequencies() {
        return frequencies;
    }
}
