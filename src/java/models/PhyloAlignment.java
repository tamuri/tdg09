package models;

import pal.alignment.Alignment;
import pal.tree.SimpleTree;
import pal.tree.Tree;
import pal.datatype.MolecularDataType;
import pal.datatype.DataTypeTool;

import java.util.*;

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
public class PhyloAlignment {

    private SimpleTree tree = null;
    private Map<String, Integer> seqData = new HashMap<String, Integer>();
    private int[] aaCount = new int[21];
    private int mostCommonRes;
    private boolean[] activeAminoAcids = new boolean[20];
    private int activeAminoAcidCount = 0;

    public PhyloAlignment(Alignment alignment, Tree tree, int site, boolean deleteSingles) {
        this.tree = (SimpleTree) tree;
        loadAlignmentAtSite(alignment, site, deleteSingles);
    }

    public void loadAlignmentAtSite(Alignment alignment, int site, boolean deleteSingles) {

        // Read in each sequence and store in seqData. Count occurrence of each residue
        MolecularDataType molDataType = DataTypeTool.getUniverisalAminoAcids();

        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            int iSeq = molDataType.getState(alignment.getData(i, site - 1));
            if ((iSeq > 19) || (iSeq < 0)) {
                iSeq = 20;
            }
            aaCount[iSeq]++;
            seqData.put(alignment.getIdentifier(i).getName(), iSeq);
        }

        // Find the most common residue
        int maxCount = aaCount[0];
        for (int i = 1; i < 20; i++) {
            if (aaCount[i] > maxCount) {
                maxCount = aaCount[i];
                mostCommonRes = i;
            }
        }

        if (deleteSingles) {
            Vector<Integer> deleteThese = new Vector<Integer>();
            for (int i = 0; i < 20; i++) {
                if (aaCount[i] == 1) {
                    deleteThese.add(i);
                }
            }

            for (String nodeName : seqData.keySet()) {
                if (deleteThese.contains(seqData.get(nodeName))) {
                    System.out.printf("Deleting single residue from %s (%s)\n", nodeName, molDataType.getChar(seqData.get(nodeName)));
                    aaCount[seqData.get(nodeName)]--;
                    seqData.put(nodeName, 20);
                }
            }
        }

        // set (and output) the number of active amino acids in this alignment
        System.out.printf("Residues: [");
        for (int i = 0; i < 20; i++) {
            activeAminoAcids[i] = false;
            if ((aaCount[i] > 0) && (i != mostCommonRes)) {
                activeAminoAcids[i] = true;
                activeAminoAcidCount++;
            }
            if (aaCount[i] > 0) {
                System.out.printf("%s:%s, ", Constants.aaNames[i], aaCount[i]);
            }
        }

        System.out.println("]");

        System.out.println("Most common residue is " + mostCommonRes + " " + Constants.aaNames[mostCommonRes]);
        System.out.println("Active residues : " + activeAminoAcidCount);

    }

    protected int getNumberOfActiveAminoAcids() {
        return activeAminoAcidCount;
    }

    protected int getAminoAcidCount(int i) {
        return aaCount[i];
    }

    protected boolean isAminoAcidActive(int pos) {
        return activeAminoAcids[pos];
    }

    protected int getMostCommonRes() {
        return mostCommonRes;
    }

    public Tree getTree() {
        return tree;
    }

    public int getAminoAcidForSequence(String name) {
        return seqData.get(name);
    }

    public double[] getEquilibriumFrequencies() {
        double[] freq = new double[20];
        double sum = 0;
        for (int i = 0; i < 20; i++) {
            freq[i] = getAminoAcidCount(i);
            sum += freq[i];
        }

        for (int i = 0; i < 20; i++) {
            freq[i] = freq[i] /sum ;
        }

        return freq;
    }
}

