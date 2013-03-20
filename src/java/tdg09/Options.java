package tdg09;

import com.beust.jcommander.Parameter;

import java.util.List;

/**
 * A container for command-line options, used by JCommander.
 * @author Asif Tamuri
 * @version 1.1
 */
public class Options {
    @Parameter(required = true, names = "-alignment", description = "Alignment in PHYLIP format")
    public String alignmentPath;

    @Parameter(required = true, names = "-tree", description = "Tree in NEWICK format")
    public String treePath;

    @Parameter(required = false, names = "-threads", description = "Number of threads to use")
    public int threads = 1;

    @Parameter(required = true, names = "-groups", variableArity = true, description = "Group labels to partition tree e.g. Av Hu")
    public List<String> groups;
}
