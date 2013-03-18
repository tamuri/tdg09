package tdg09;

import com.beust.jcommander.Parameter;

import java.util.List;

/**
 * Command-line options
 *
 * User: atamuri
 * Date: 18/03/2013 16:32
 */
public class Options {
    @Parameter(required = true, names = "-alignment")
    public String alignmentPath;

    @Parameter(required = true, names = "-tree")
    public String treePath;

    @Parameter(required = false, names = "-threads")
    public int threads = 1;

    @Parameter(required = true, names = "-groups", variableArity = true)
    public List<String> groups;

    // TODO: have a 'name' for each run, like RAxML
}
