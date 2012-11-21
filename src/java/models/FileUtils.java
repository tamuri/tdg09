package models;

import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.DataTypeTool;
import pal.tree.ReadTree;
import pal.tree.Tree;

import java.io.BufferedReader;
import java.io.FileReader;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
 */
public class FileUtils {
    public static Tree readTree(String filePath) {
        Tree tree;
        try {
            tree = new ReadTree(filePath);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return tree;
    }

    public static Alignment readAlignment(String filePath) {
        Alignment a;
        try {
            a = AlignmentReaders.readPhylipClustalAlignment(new BufferedReader(new FileReader(filePath)), DataTypeTool.getUniverisalAminoAcids());
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return a;
    }
}
