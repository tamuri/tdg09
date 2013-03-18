package tdg09.utils;

import com.google.common.collect.Sets;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.DataTypeTool;
import pal.tree.ReadTree;
import pal.tree.Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Set;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 18/03/2013 16:48
 */
public class CoreUtils {
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

    public static boolean isTreeAndAlignmentValid(Tree tree, Alignment alignment) {
        Set<String> nodes = Sets.newHashSet();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            nodes.add(tree.getExternalNode(i).getIdentifier().getName());
        }

        Set<String> seqs = Sets.newHashSet();
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            seqs.add(alignment.getIdentifier(i).getName());
        }

        return !(seqs.size() != nodes.size() || Sets.symmetricDifference(seqs, nodes).size() > 0);
    }

}
