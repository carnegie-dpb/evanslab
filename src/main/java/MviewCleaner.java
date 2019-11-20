import java.io.*;

/**
 * Clean out empty lines from Mview or ClustalW text output.
 */
/*
  CLUSTAL O(1.2.4) multiple sequence alignment
  
012                   22                                                         81
  |                   |                                                          |
  Zm00014a002653      CCAATCAGACCCAATAGCATTCAGTGGGCCTTAACTTCCTGCAGTTCCTGGGCTGTGGAA        60
  97991               ------------------------------------------------------------        0
  844                 ------------------------------------------------------------        0
 */
public class MviewCleaner {

    static String txtDelimiter = "[";
    static String htmlDelimiter = "<SPAN style=\"color:#666666\">";
    
    public static void main(String[] args) throws IOException {
        if (args.length==0) {
            System.err.println("Usage: MviewCleaner mview.txt|mview.html|clustalo-*.clustal_num");
            System.exit(0);
        }
        boolean isTxt = args[0].endsWith("txt");
        boolean isHtml = args[0].endsWith("html");
        boolean isClustal = args[0].endsWith("clustal_num");
        if (!isTxt && !isHtml && !isClustal) {
            System.err.println("Input file must end in .txt or .html or .clustal_num.");
            System.exit(1);
        }
        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        int start = -1;
        String line = null;
        while ((line=reader.readLine())!=null) {
            boolean skipit = false;
            if (isTxt && line.indexOf(txtDelimiter)>0) start = line.indexOf(txtDelimiter);
            if (isHtml && line.indexOf(htmlDelimiter)>0) start = line.indexOf(htmlDelimiter);
            if (isClustal) start = 20;
            if (start>0) {
                String content = line.replace(htmlDelimiter,"").replace("</SPAN>","");
                if (content.length()>start) {
                    int jEnd = content.length();
                    if (isClustal) jEnd = Math.min(jEnd,80);
                    skipit = content.charAt(start)=='-';
                    for (int j=start; j<jEnd; j++) {
                        skipit = skipit && (content.charAt(j)=='-' || content.charAt(j)==' ');
                    }
                }
            }
            if (!skipit) {
                System.out.println(line);
            }
        }
    }
}
