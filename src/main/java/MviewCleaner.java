import java.io.*;

/**
 * Clean out empty lines from MVIEW text output.
 */
public class MviewCleaner {

    static String txtDelimiter = "[";
    static String htmlDelimiter = "<SPAN style=\"color:#666666\">";
    
    public static void main(String[] args) throws IOException {
        if (args.length==0) {
            System.err.println("Usage: MviewCleaner mview.txt|mview.html");
            System.exit(0);
        }
        boolean isTxt = args[0].endsWith("txt");
        boolean isHtml = args[0].endsWith("html");
        if (!isTxt && !isHtml) {
            System.err.println("Input file must end in .txt or .html.");
            System.exit(1);
        }
        BufferedReader reader = new BufferedReader(new FileReader(args[0]));
        int start = -1;
        String line = null;
        while ((line=reader.readLine())!=null) {
            boolean skipit = false;
            if (isTxt && line.indexOf(txtDelimiter)>0) start = line.indexOf(txtDelimiter);
            if (isHtml && line.indexOf(htmlDelimiter)>0) start = line.indexOf(htmlDelimiter);
            if (start>0) {
                String content = line.replace(htmlDelimiter,"").replace("</SPAN>","");
                if (content.length()>start) {
                    skipit = content.charAt(start)=='-';
                    for (int j=start; j<content.length(); j++) {
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
