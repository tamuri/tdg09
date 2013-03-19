package tdg09;

import com.google.common.collect.Sets;
import pal.alignment.Alignment;
import pal.alignment.AlignmentReaders;
import pal.datatype.DataTypeTool;
import pal.tree.ReadTree;
import pal.tree.Tree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.Comparator;
import java.util.Set;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 18/03/2013 16:48
 */
public class Utils {
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

    /**
     * Returns a generics Comparator for comparing a double field in class T.
     * Usage: Collections.sort(listOfMyClass, Utils.doubleComparator("myDoubleField", MyClass.class));
     */
    public static <T> Comparator<T> doubleComparator(final String name, final Class<T> clazz) {
        return new Comparator<T>() {
            public int compare(T o1, T o2) {
                try {
                    Field f = clazz.getField(name);

                    double d1 = f.getDouble(o1);
                    double d2 = f.getDouble(o2);

                    if (d1 < d2) {
                        return -1;
                    } else if (d1 == d2) {
                        return 0;
                    } else {
                        return 1;
                    }

                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }
        };
    }

}
