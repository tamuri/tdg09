package tdg09;

import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.introspector.BeanAccess;
import org.yaml.snakeyaml.introspector.Property;
import org.yaml.snakeyaml.introspector.PropertyUtils;
import org.yaml.snakeyaml.representer.Representer;

import java.beans.IntrospectionException;
import java.util.*;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 19/03/2013 09:43
 */
public class Result {
    public final int site;
    public final Map<Character, Integer> residues;
    public final List<Model> models;
    public double lrt;
    public double fdr;

    public Result(int site, Map<Character, Integer> residues, List<Model> models) {
        this.site = site;
        this.residues = residues;
        this.models = models;
    }

    public void toYaml() {
        class PropertyOrder extends PropertyUtils {
            @Override
            protected Set<Property> createPropertySet(Class<? extends Object> type, BeanAccess bAccess) throws IntrospectionException {
                Set<Property> result = new LinkedHashSet<Property>();

                if (type == Result.class) {
                    result.add(super.getProperty(type, "site"));
                    result.add(super.getProperty(type, "residues"));
                    result.add(super.getProperty(type, "models"));
                } else if (type == Model.class) {
                    result.add(super.getProperty(type, "name"));
                    result.add(super.getProperty(type, "converged"));
                    result.add(super.getProperty(type, "parameters"));
                    result.add(super.getProperty(type, "rate"));
                    result.add(super.getProperty(type, "frequencies"));
                    result.add(super.getProperty(type, "lnL"));
                }

                return result;
            }
        }

        Representer representer = new Representer();
        representer.setPropertyUtils(new PropertyOrder());

        Yaml yaml = new Yaml(representer);
        System.out.println(yaml.dump(this));
    }
}

