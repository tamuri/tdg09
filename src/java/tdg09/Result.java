package tdg09;

import com.google.common.collect.Lists;
import org.yaml.snakeyaml.introspector.BeanAccess;
import org.yaml.snakeyaml.introspector.Property;
import org.yaml.snakeyaml.introspector.PropertyUtils;
import org.yaml.snakeyaml.representer.Representer;
import tdg09.models.Constants;

import java.beans.IntrospectionException;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Holds all the relevant results from the analysis of a site
 *
 * @author Asif Tamuri
 * @version 1.1
 */
public class Result {
    public final int site;
    public final Map<Character, Integer> residues;
    public final List<Model> models;
    public double lrt;
    public double fdr;

    @Override
    public String toString() {
        return "Result{" +
                "site=" + site +
                ", residues=" + residues +
                ", models=" + models +
                ", lrt=" + lrt +
                ", fdr=" + fdr +
                '}';
    }

    public Result(int site, Map<Character, Integer> residues, List<Model> models) {
        this.site = site;
        this.residues = residues;
        this.models = models;
    }

    public static class Model {
        public final String name;
        public final int parameters;
        public final double rate;
        public final double lnL;
        public final List<Map<Character,Double>> frequencies;
        public final boolean converged;

        public Model(String name, int parameters, double rate, double lnL, List<Map<Character, Double>> frequencies, boolean converged) {
            this.name = name;
            this.parameters = parameters;
            this.rate = rate;
            this.lnL = lnL;
            this.frequencies = frequencies;
            this.converged = converged;
        }

        @Override
        public String toString() {
            return "Model{" +
                    "name='" + name + '\'' +
                    ", parameters=" + parameters +
                    ", rate=" + rate +
                    ", lnL=" + lnL +
                    ", frequencies=" + frequencies +
                    ", converged=" + converged +
                    '}';
        }

        public List<List<Double>> getFlatFrequencies() {
            List<List<Double>> f = Lists.newArrayList();

            for (Map<Character, Double> mf : this.frequencies) {
                List<Double> d = Lists.newArrayList();
                for (int i = 0; i < Constants.aaNames.length; i++) {
                    if (mf.containsKey(Constants.aaNames[i])) {
                        d.add(mf.get(Constants.aaNames[i]));
                    } else {
                        d.add(0.0);
                    }
                }
                f.add(d);
            }

            return f;
        }
    }

    public static Representer getYamlRepresenter() {

        Representer r = new Representer();
        r.setPropertyUtils(new PropertyUtils() {
            @Override
            protected Set<Property> createPropertySet(Class<? extends Object> type, BeanAccess bAccess) throws IntrospectionException {
                Set<Property> result = new LinkedHashSet<Property>();

                if (type == Result.class) {
                    result.add(super.getProperty(type, "site"));
                    result.add(super.getProperty(type, "residues"));
                    result.add(super.getProperty(type, "models"));
                } else if (type == Result.Model.class) {
                    result.add(super.getProperty(type, "name"));
                    result.add(super.getProperty(type, "converged"));
                    result.add(super.getProperty(type, "parameters"));
                    result.add(super.getProperty(type, "rate"));
                    result.add(super.getProperty(type, "frequencies"));
                    result.add(super.getProperty(type, "lnL"));
                } else {
                    result.addAll(super.getProperties(type, bAccess));
                }

                return result;
            }
        });

        return r;
    }
}

