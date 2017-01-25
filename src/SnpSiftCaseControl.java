package edu.carnegiescience.dpb.evanslab;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.Map;
import java.util.HashMap;

/**
 * Read the VCF output from SnpSift caseControl. Format is:
 *
 * CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Samp1    Samp2    ...    SampN
 */
public class SnpSiftCaseControl {

    public static void main(String[] args) {

        if (args.length!=2) {
            System.out.println("Usage: SnpSiftCaseControl <caseControl-annoted VCF file> <max CC_REC value>");
            System.exit(1);
        }
        
        try {
            
            File file = new File(args[0]);
            if (!file.exists()) {
                System.err.println("Error: "+args[0]+" does not exist.");
                System.exit(1);
            }

            Double maxCCRec = Double.parseDouble(args[1]);

            // header

            System.out.println("Chr\tPos\tRef\tAlt\tCC_TREND\tCC_GENO\tCC_ALL\tCC_DOM\tCC_REC\tType1\tName1\tType2\tName2\tType3\tName3\tType4\tName4\tType5\tName5");
            
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line;
            while ((line=reader.readLine())!=null) {
                if (!line.startsWith("#")) {

                    // 0      1       2       3       4       5       6       7       8       9        10     ...  N+8    
                    // CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Samp1    Samp2  ...  SampN
                    String[] parts = line.split("\t");
                    String chrom = parts[0];
                    int pos = Integer.parseInt(parts[1]);
                    String id = parts[2];
                    char ref = parts[3].charAt(0);
                    char alt = parts[4].charAt(0);
                    double qual = Double.parseDouble(parts[5]);
                    String filter = parts[6];
                    String info = parts[7];
                    String format = parts[8];
                    // nullify dot values
                    if (id.equals(".")) id = null;
                    if (filter.equals(".")) filter = null;
                    // load the samples
                    String[] samples = new String[parts.length-9];
                    for (int i=0; i<samples.length; i++) {
                        samples[i] = parts[i+9];
                        if (samples[i].equals(".")) samples[i] = null;
                    }

                    // parse info string into a map
                    Map<String,String> infoMap = new HashMap<String,String>();
                    String[] items = info.split(";");
                    for (int i=0; i<items.length; i++) {
                        String[] pieces = items[i].split("=");
                        infoMap.put(pieces[0], pieces[1]);
                    }

                    // extract some annotation
                    // 0|1                |2       |3                            |4                            |5                |6                            |||9           ||||||
                    // G|intergenic_region|MODIFIER|Zm00001d034986-Zm00001d034987|Zm00001d034986-Zm00001d034987|intergenic_region|Zm00001d034986-Zm00001d034987|||n.1011743T>G||||||
                    Map<String,String> annMap = new HashMap<String,String>();
                    String[] annParts = infoMap.get("ANN").split(",");
                    for (int i=0; i<annParts.length; i++) {
                        String[] pieces = annParts[i].split("\\|");
                        if (pieces.length>3) annMap.put(pieces[1],pieces[3]); // should limit to gene names
                    }

                    // extract p values
                    double ccTrend = Double.parseDouble(infoMap.get("CC_TREND"));
                    double ccGeno = Double.parseDouble(infoMap.get("CC_GENO"));
                    double ccAll = Double.parseDouble(infoMap.get("CC_ALL"));
                    double ccDom = Double.parseDouble(infoMap.get("CC_DOM"));
                    double ccRec = Double.parseDouble(infoMap.get("CC_REC"));

                    // only print out loci with small enough p values
                    if (!Double.isNaN(ccRec) && ccRec<maxCCRec) {
                        System.out.print(chrom+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+ccTrend+"\t"+ccGeno+"\t"+ccAll+"\t"+ccDom+"\t"+ccRec);
                        for (String annFeature : annMap.keySet()) {
                            System.out.print("\t"+annFeature+"\t"+annMap.get(annFeature));
                        }
                        System.out.println("");
                    }
                    
                }

            }

        } catch (Exception ex) {

            System.err.println(ex.toString());
            System.exit(1);

        }

    }


}
