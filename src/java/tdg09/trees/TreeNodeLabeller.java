package tdg09.trees;

import pal.tree.*;
import tdg09.Utils;

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
public class TreeNodeLabeller {
    public static void main(String[] args) throws Exception {
        TreeNodeLabeller tnl = new TreeNodeLabeller();

        Tree t = Utils.readTree(args[0]);
        Tree out = tnl.label(t);

        FileWriter fw  = new FileWriter(args[0] + ".out");
        PrintWriter pw = new PrintWriter(fw);
        TreeUtils.printNH(out, pw);
        fw.close();
        pw.close();
    }

    public SimpleTree label(Tree t) {

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
            // System.out.printf("Branch %s..%s different.\n", preUnknownNodeCount, postUnknownNodeCount);
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
                if (s.size() > 1) {
                    System.out.printf("# Node %s from %s resolved\n", n.getNumber(), s.toString());
                }

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
                    n.getIdentifier().setName(parentNodeName.substring(0,2) + "_GS");
                }
            } else {
                // n.getIdentifier().setName("ROOT");
            }
        }


        return st;


    }

    private void printStringSet(Set<String> set) {
        for (String s : set) {
            System.out.printf("[%s] ", s);
        }
    }
}
