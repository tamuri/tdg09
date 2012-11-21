package trees;

import pal.tree.*;

import java.io.*;
import java.util.*;

/**
 * @author Asif Tamuri <atamuri@nimr.mrc.ac.uk>
 * @version 1.0
 *
 * Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009)
 * Identifying Changes in Selective Constraints: Host Shifts in Influenza.
 * PLoS Comput Biol 5(11): e1000564.
 * doi:10.1371/journal.pcbi.1000564
 */
public class TreeNodeLabeler {
        public static void main(String[] args) {
        TreeNodeLabeler m = new TreeNodeLabeler();
        m.run(args[0]);
    }

    private void run(String fileIn) {
        Tree t;
        try {
            t = TreeTool.readTree(new BufferedReader(new FileReader(fileIn)));
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException();
        }

        SimpleTree st = new SimpleTree(t);
        List<Node> unknownNodes = new ArrayList<Node>();
        List<Node> knownNodes = new ArrayList<Node>();

        // to begin with all internal nodes are unnamed
        for (int i=0; i< st.getInternalNodeCount(); i++) {
            unknownNodes.add(st.getInternalNode(i));
        }

        int preUnknownNodeCount = unknownNodes.size();
        int postUnknownNodeCount = -1;

        while (unknownNodes.size() > 0 && preUnknownNodeCount != postUnknownNodeCount) {
            System.out.printf("%s != %s\n", preUnknownNodeCount, postUnknownNodeCount);
            preUnknownNodeCount = unknownNodes.size();
            // loop through each of the unknown nodes
            for (Node n : unknownNodes) {
                // get a distinct list of all hosts for the node
                Set<String> s = new HashSet<String>();

                // loop over every child node
                for (int j=0; j<n.getChildCount(); j++) {
                    Node c = n.getChild(j);
                    String childNodeName = c.getIdentifier().getName();
                    if (c.isLeaf()) {
                        s.add(childNodeName.substring(0,2));
                    } else {
                        // one of the children is an internal node with a name
                        if (childNodeName != null && childNodeName.length() > 0) {
                            s.add(childNodeName.substring(0,2));
                        }
                    }
                }

                // list all the possible names for this node
                System.out.printf("%s = ", n.getNumber());
                printStringSet(s);
                System.out.printf("\n");

                if (s.size() == 1) {
                    n.getIdentifier().setName(s.iterator().next());
                    knownNodes.add(n);
                }
            }

            // remove all known nodes
            for (Node n : knownNodes) {
                unknownNodes.remove(n);
            }

            postUnknownNodeCount = unknownNodes.size();
        }

        // now the only remaining unknown nodes are those that neighbour two different clades (and the root node)
        for (Node n : unknownNodes) {
            // the name of this hostshift node = parent node name (remember, could be the root node = leave unlabelled)
            if (n.getParent() != null) {
                String parentNodeName = n.getParent().getIdentifier().getName();
                if (parentNodeName != null && parentNodeName.length() > 0) {
                    n.getIdentifier().setName(parentNodeName.substring(0,2) + "_HS");
                }
            } else {
                n.getIdentifier().setName("ROOT");
            }
        }
     
        FileWriter fw;
        try {
            fw = new FileWriter(fileIn + ".out");
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException();
        }

        PrintWriter pw = new PrintWriter(fw);
        TreeUtils.printNH(st, pw);
        pw.close();
    }

    private void printStringSet(Set<String> set) {
        for (String s : set) {
            System.out.printf("%s ", s);
        }
    }
}
