package tdg09.models;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import pal.alignment.Alignment;
import pal.datatype.DataTypeTool;
import pal.datatype.MolecularDataType;
import pal.tree.SimpleTree;
import pal.tree.Tree;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Prepares the site pattern for analysis by the LikelihoodCalculator
 *
 * @author Asif Tamuri
 * @author Richard Goldstein
 * @version 1.1
 */
public class PhyloAlignment {

    private SimpleTree tree = null;
    private Map<String, Integer> seqData = new HashMap<String, Integer>();
    private int[] aaCount = new int[21];
    private int mostCommonRes;
    private boolean[] activeAminoAcids = new boolean[20];
    private int activeAminoAcidCount = 0;
    private Map<Character, Integer> residues = Maps.newLinkedHashMap();

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
            List<Integer> deleteThese = Lists.newArrayList();
            for (int i = 0; i < 20; i++) {
                if (aaCount[i] == 1) {
                    deleteThese.add(i);
                }
            }

            for (String nodeName : seqData.keySet()) {
                if (deleteThese.contains(seqData.get(nodeName))) {
                    aaCount[seqData.get(nodeName)]--;
                    seqData.put(nodeName, 20);
                }
            }
        }

        for (int i = 0; i < 20; i++) {
            activeAminoAcids[i] = false;
            if ((aaCount[i] > 0) && (i != mostCommonRes)) {
                activeAminoAcids[i] = true;
                activeAminoAcidCount++;
            }
            if (aaCount[i] > 0) {
                residues.put(Constants.aaNames[i], aaCount[i]);
            }
        }


    }

    public Map<Character, Integer> getResidues() {
        return residues;
    }

    public int getNumberOfActiveAminoAcids() {
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

