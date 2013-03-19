package tdg09;

import com.beust.jcommander.JCommander;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.yaml.snakeyaml.Yaml;
import pal.alignment.Alignment;
import pal.statistics.LikelihoodRatioTest;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import tdg09.trees.TreeNodeLabeler;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.sql.Timestamp;
import java.util.*;
import java.util.concurrent.*;

/**
 * The main entry point for tdg09 analysis.
 *
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 18/03/2013 16:31
 */
public class Analyse {
    Options options = new Options();

    public static void main(String[] args) {
        Analyse a = new Analyse();
        a.run(args);
    }

    private void run(String[] args) {
        // Parse the command-line options
        JCommander jc = new JCommander(options);
        try {
            jc.parse(args);
        } catch (Exception e) {
            jc.setProgramName("java -cp tdg09.jar tdg09.Analyse");
            jc.usage();
            System.exit(1);
        }

        // TODO: print output to console and also a file with the 'name' prefix, ala RAxML
        System.out.printf("StartTime: %s\n", new Timestamp(System.currentTimeMillis()));
        System.out.printf("WorkingDirectory: %s\n", System.getProperty("user.dir"));
        System.out.printf("Options: %s\n\n", Joiner.on(" ").join(Lists.newArrayList(args)));


        // Load the tree and alignment
        System.out.printf("TreeFile: %s\n", new File(options.treePath).getAbsolutePath());
        Tree tree = Utils.readTree(options.treePath);

        System.out.printf("AlignmentFile: %s\n\n", new File(options.alignmentPath).getAbsolutePath());
        Alignment alignment = Utils.readAlignment(options.alignmentPath);


        validate(tree, alignment);

        // Label the tree's internal nodes if they haven't already
        boolean labelled = true;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            Node node = tree.getInternalNode(i);
            if (!node.isRoot() && node.getIdentifier().getName().length() == 0) {
                labelled = false;
                break;
            }
        }

        if (!labelled) {
            System.out.println("# The internal nodes of the tree are not labelled. Labelling...");
            TreeNodeLabeler labeler = new TreeNodeLabeler();
            tree = labeler.label(tree);
        }

        // Output the labelling of internal nodes
        printGroupChanges(tree.getRoot());
        System.out.println();

        StringWriter sw = new StringWriter();
        TreeUtils.printNH(tree, new PrintWriter(sw));
        System.out.printf("LabelledTree:\n\t%s\n", sw.toString().replaceAll("\n", ""));
        System.out.println();

        // We're now ready to run
        ExecutorService pool = Executors.newFixedThreadPool(options.threads);
        List<Future<Result>> futures = Lists.newArrayList();

        for (int i = 1; i <= alignment.getSiteCount(); i++) {
        //for (int i = 1; i <= 10; i++) {
            futures.add(pool.submit(new SiteAnalyser(alignment, tree, i, options.groups)));
        }

        List<Result> results = Lists.newArrayList();

        try {
            for (Future<Result> future : futures) {
                results.add(future.get());
                future.get().toYaml();
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        pool.shutdown();

        System.out.print("\n---\n\n");
        doLRT(results);

        System.out.print("\n---\n\n");
        printSummary(results);

        System.out.println();
        System.out.printf("EndTime: %s\n", new Timestamp(System.currentTimeMillis()));

        System.out.println("...");
    }

    private void printSummary(List<Result> results) {
        System.out.println("# Other information");
        System.out.println();

        List<Integer> sites = Lists.newArrayList();
        for (Result r : results) if (r.models == null) sites.add(r.site);
        System.out.println("ConservedPositions:");
        System.out.printf ("    Count: %s\n", sites.size());
        System.out.printf ("    Sites: %s\n", new Yaml().dump(sites));
    }

    private void doLRT(List<Result> results) {
        System.out.println("# Likelihood ratio test:");
        System.out.println();

        List<Result> polymorphicSites = Lists.newArrayList();

        for (Result r : results) {
            if (r.models != null) {
                polymorphicSites.add(r);

                r.lrt = LikelihoodRatioTest.getSignificance(
                        r.models.get(1).lnL - r.models.get(0).lnL, // delta lnL
                        r.models.get(1).parameters - r.models.get(0).parameters); // degrees of freedom
            }
        }

        // order by likelihood ratio test p-value
        Collections.sort(polymorphicSites, Utils.getDoubleComparator("lrt", results.get(0)));

        // calculate false discovery rate (naive method)
        for (int i = 0; i < polymorphicSites.size(); i++) {
            int rank = i + 1;
            Result r = polymorphicSites.get(i);
            r.fdr = r.lrt * polymorphicSites.size() / rank;
        }

        // order by false discovery rate
        Collections.sort(polymorphicSites, Utils.getDoubleComparator("fdr", results.get(0)));

        System.out.println("# Site, delta lnL,  dof, P-value,   FDR");
        for (Result r : polymorphicSites) {
            System.out.printf("- %1$4d, %2$3.7f,  %3$s,   %4$.7f, %5$.7f%n",
                    r.site,
                    r.models.get(1).lnL - r.models.get(0).lnL,
                    r.models.get(1).parameters - r.models.get(0).parameters,
                    r.lrt,
                    r.fdr);
        }
    }

    private void validate(Tree tree, Alignment alignment) {
        // 1. Check the tree and alignment are in agreement
        if (!Utils.isTreeAndAlignmentValid(tree, alignment)) {
            System.out.println("ERROR: The tree and alignment do not have the same taxa.");
            System.exit(1);
        } else {
            System.out.printf("Alignment:\n\tSequenceCount: %s\n\tSiteCount: %s\n\n", alignment.getSequenceCount(), alignment.getSiteCount());
        }

        // 2. Check that all taxa have an assigned group and all groups are used
        List<String> taxa = Lists.newArrayList();
        Set<String> usedGroups = Sets.newHashSet();

        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            taxa.add(alignment.getIdentifier(i).getName());
        }

        for (String group : options.groups) {
            Iterator<String> it = taxa.iterator();
            while (it.hasNext()) {
                if (it.next().startsWith(group)) {
                    it.remove();
                    usedGroups.add(group);
                }
            }
        }

        if (taxa.size() > 0) {
            System.out.printf("ERROR: %s taxa are not allocated to a group { %s }\n", taxa.size(), taxa.toString());
            System.exit(1);
        }

        Set<String> unusedGroups = Sets.difference(Sets.newHashSet(options.groups), usedGroups);
        if (unusedGroups.size() > 0) {
            System.out.printf("# WARNING: group(s) %s unused and will be ignored.\n", unusedGroups.toString());
            for (String s : unusedGroups) {
                options.groups.remove(s);
            }
        }

        System.out.printf("Groups: [%s]\n\n", Joiner.on(", ").join(options.groups));

    }

    private void printGroupChanges(Node node) {
        String nodeGroup;
        if (node.isRoot()) {
            System.out.printf("# Assuming that root of tree is in group [%s]\n", options.groups.get(0));
            nodeGroup = options.groups.get(0);
        } else {
            nodeGroup = getGroup(node);
        }

        for (int i = 0; i < node.getChildCount(); i++) {
            Node child = node.getChild(i);
            String childGroup = getGroup(child);
            if (!nodeGroup.equals(childGroup)) {
                System.out.printf("# Switching from group [%s] to [%s] at branch %s..%s\n", nodeGroup, childGroup, node.getNumber(), child.getNumber());
            }
            printGroupChanges(child);
        }

    }

    private String getGroup(Node node) {
        for (String g : options.groups) {
            if (node.getIdentifier().getName().startsWith(g)) return g;
        }
        return null;
    }


}


