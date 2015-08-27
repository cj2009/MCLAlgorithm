package mclalgorithm;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import static mclalgorithm.MCLAlgorithm.exp;
import static mclalgorithm.MCLAlgorithm.inflate;
import static mclalgorithm.MCLAlgorithm.normalize;
import static mclalgorithm.MCLAlgorithm.readFile;
import org.jblas.DoubleMatrix;

/*
 * This is a program that performs the MCL algorithm in Java, while also comparing
 * the matrix from each iteration to the matrix produced by our Hadoop program 
 * in each iteration. This helps to verify that both algorithms are producing 
 * exactly the same output after each iteration.
 */

/**
 *
 * @author chris_joseph
 */
public class MCLTesting {
    
    public static void main(String args[]) throws Exception {
        int e = 2; // expansion parameter; responsible for the strengthening/weaking of currents
        int r = 2; // inflation parameter; controls the granularity of the clusters
        

        // This is the initial array:
        // The input file must follow certain a format. See the comments associated
        // with the 'readFile' method.
        DoubleMatrix matrix = readFile("edges.txt");
        System.out.println("Successfully read edge file from disk!");
        
        
        // Add self edges. This is part of the MCL algorithm.
        for (int i = 0; i < matrix.rows; i++) {
            matrix.put(i, i, 1.0);
        }
        
        
        // Initialization step: normalize the matrix.
        normalize(matrix);
        System.out.println("Markov matrix created!");

        
        // Now we alternate between expanding and inflating this matrix, until
        // a steady state is obtained:
        // According to Macropol (UCSB), most graphs need anywhere from 10
        // to 100 steps to converge (based on empirical results).
        for (int i = 0; i < 30; i++) {
            System.out.println("CYCLE " + i + ":");
            matrix = exp(matrix, e);
            System.out.println("\t\t->Expansion complete");
            inflate(matrix, r);
            System.out.println("\t\t->Inflation Complete");
            
            String filePath = "testing/out-MR3-run0-cycle" + i + "/part-r-00000";
            DoubleMatrix matrix2 = readFile2(filePath, matrix.rows);
            System.out.println("\t\t->Successfully read Hadoop output file from disk!");
            
            
            
            // Compare the results with the results of our Hadoop algorithm.
            // We're going to write the results of these comparisons to a file.
            File file = new File("testing_results/results-" + i + ".txt");
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);


            bw.write("These are the results from cycle " + i + " of the expansion/inflation process.\n\n");
            bw.write("Java Matrix: \n");
            writeMatrix(matrix, bw);
            bw.write("\n\nHadoop Matrix: \n");
            writeMatrix(matrix2, bw);

            boolean b = matrixEquals(matrix, matrix2);
            if (b) {
//                bw.write("\nThese matrices are the same.");
                System.out.println("\t\t->Matrices are the same!\n");
            }
                
            else {
//                bw.write("These matrices are NOT THE SAME!");
                System.out.println("\t\t->Matrices are NOT THE SAME!\n");
            }
                

            bw.close();
            fw.close();
        }
        
        

    }
    
    
    public static DoubleMatrix readFile2(String fileName, int n) throws IOException {
        ArrayList<Integer[]> edges = new ArrayList();
        
        Path path = Paths.get(fileName);
        BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8);
        String line = null;
        String parts[];
        String entries[];
        String s[];
        int v1, v2;

        DoubleMatrix matrix = new DoubleMatrix(n, n);
        
        while ((line = br.readLine()) != null) {
            // Store each edge into the arraylist.
            
            parts = line.split("\t");
            entries = parts[1].split(",");
            v1 = Integer.parseInt(parts[0]);
            for (String entry : entries) {
                s = entry.split(":");
                v2 = Integer.parseInt(s[0]);
                matrix.put(v2, v1, Double.parseDouble(s[1]));
            }
            

        }
        br.close();
        
        // After we finish reading the file, return the matrix:
        return matrix;
    }
    
    public static void writeMatrix(DoubleMatrix matrix, BufferedWriter bw) throws IOException {
        DecimalFormat df = new DecimalFormat("#0.00000");
        for (int i = 0; i < matrix.rows; i++) {
            String s = "";
            
            for (int j = 0; j < matrix.columns - 1; j++) {
                s += df.format(matrix.get(i, j)) + "\t";
            }
            s += df.format(matrix.get(i, matrix.columns - 1));
            
            
            bw.write(s + "\n");
        }
    }
    
    
    /**
     * Tests if two matrices are the same, within 5 decimal places of precision.
     */
    public static boolean matrixEquals(DoubleMatrix m1, DoubleMatrix m2) {
        boolean returnValue = true;
        for (int i = 0; i < m1.rows; i++) {
            for (int j = 0; j < m1.columns; j++) {
                double diff = Math.abs(m1.get(i, j) - m2.get(i, j));
                if (diff > 0.000001) {
                    System.out.printf("\t\t  NO! correct[%d][%d] = %.5f, m2[%d][%d] = %.5f\n", i, j, m1.get(i, j), i, j, m2.get(i, j));
                    returnValue = false;
                }
            }
        }
        return returnValue;
    }
}
