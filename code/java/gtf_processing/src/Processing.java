import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import java.util.List;
import java.util.ArrayList;
import java.util.Scanner;

/***
 * This is a class to process the Mouse GTF file.
 * The processing is 2 step 
 * 1. Change the chromosome column to chrChromsomenumber for eg. '1' to 'chr1'
 * 2. Split the file of features chromsomewise. 
 * 
 * 
 * The resulting files are then used for Genewise sync file creation using the Popoolation2 tools.
 * @author nehiljain
 *
 */


public class Processing {
  
  public static void main(String... aArgs) throws IOException{

	Processing text = new Processing();
	
    //treat as a large file - use some buffering
    text.readLargerTextFile(FILE_NAME);   
  }
  /***
   * The FILE_NAME is the input filename
   * OUTPUT_Path is th output directory where all the different split files are written
   * Encoding is standard
   * rows is the arraylist which stores all the lines as one String element which is cleared after the change of chromosome
   * row is each row itself which stores the line read in and the modification
   * chrosomeNumber stores the previous chromosomeNumber
   */
  final static String FILE_NAME = "/Users/nehiljain/Documents/workspace/gtf_processing/Mus_musculus.GRCm38.75.gtf";
  final static String OUTPUT_Path = "/Users/nehiljain/Documents/workspace/gtf_processing/data_files/";
  final static Charset ENCODING = StandardCharsets.UTF_8;
  static List<String> rows = new ArrayList<String>();
  static StringBuffer row;
  static StringBuffer chromosomeNumber = new StringBuffer("");
  //For larger files
  
  void readLargerTextFile(String aFileName ) throws IOException {
	  //creates a oath object for Input File
	  Path inPath = Paths.get(aFileName);
	  //start a scanner for INPUT FILE
    try (Scanner scanner =  new Scanner(inPath, ENCODING.name())){
    	//Scan each line
    		while (scanner.hasNextLine()){
    			//process each line in some way
    			row = new StringBuffer(scanner.nextLine());
    			//checks for starting and change of chromosome
    			if (chromosomeNumber.equals("") ||
    					!chromosomeNumber.toString().equals(row.substring(0, 
    												row.indexOf("\t")))) {
    				//Write the arraylist to a file
    				writeLargerTextFile(OUTPUT_Path + chromosomeNumber.toString() + "_Mod_Mus.gtf", rows);
    				//EMoty the array list to save space 
    				rows.clear();
    				//store the new chrosome number
    				chromosomeNumber = new StringBuffer(row.substring(0, 
							row.indexOf("\t")));
    				
    				log("Changed!! -" + chromosomeNumber);
    			
    			}
    			//do the modification of each row
    			row.replace(0, row.indexOf("\t"), "chr" + chromosomeNumber);
    			rows.add(row.toString());
    	}      
    }
  }
  
  void writeLargerTextFile(String aFileName, List<String> aLines) throws IOException {
	    Path path = Paths.get(aFileName);
	    try (BufferedWriter writer = Files.newBufferedWriter(path, ENCODING)){
	      for(String line : aLines){
	        writer.write(line);
	        writer.newLine();
	      }
	    }
	  }


  private static void log(Object aMsg){
    System.out.println(String.valueOf(aMsg));
  }
  
} 