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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Set;
import org.jblas.DoubleMatrix;

/**
 * This is a class that performs the MCL algorithm on a graph. It works only on 
 * undirected, unweighted graphs.
 * Based on the original MCL Algorithm by Stijn van Dongen.
 * Also based on the presentation by Kathy Macropol (UCSB):
 * (www.cs.ucsb.edu/~xyan/classes/CS595D-2009winter/MCL_Presentation2.pdf)
 * 
 * @author chris_joseph
 */
public class MCLAlgorithm {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        int cycles = 20;
        int e = 2; // expansion parameter; responsible for the strengthening/weaking of currents
        double r = 2.0; // inflation parameter; controls the granularity of the clusters
        // See the PDF mentioned above for detailed terminology
        

        // This is the initial array:
        // The input file must follow certain a format. See the comments associated
        // with the 'readFile' method below.
        DoubleMatrix matrix = readFile("edges.txt");
        
        System.out.println("Graph successfully loaded from file.");
        DoubleMatrix matrix2 = matrix.dup(); // keep a backup, which we'll need later
        
        
        // Add self edges. This is part of the MCL algorithm.
        for (int i = 0; i < matrix.columns; i++) {
            matrix.put(i, i, 1.0);
        }
        
        
        // Initialization step: normalize the matrix.
        normalize(matrix);
        System.out.println("Initial normalization complete.");

        
        // Now we alternate between expanding and inflating this matrix, until
        // a steady state is obtained:
        // According to Macropol (UCSB), most graphs need anywhere from 10
        // to 100 steps to converge (based on empirical results).
        for (int i = 0; i < cycles; i++) {
            System.out.println("Cycle " + i + ":");
            matrix = exp(matrix, e);
            System.out.println("\t-->Expansion complete");
            inflate(matrix, r);
            System.out.println("\t-->Inflation complete");
        }
        System.out.println("MCL complete. \n");
        
        
        // Now, interpret the results!
        // This is the list of clusters. Each list contains all the vertices
        // that belong to one particular cluster.
        ArrayList<ArrayList<Integer>> clusters = identifyClusters(matrix);
        
        
        // Print all the clusters:
        for (ArrayList<Integer> cluster : clusters) {
            System.out.print("Cluster: ");
            for (Integer i : cluster)
                System.out.print(i + ", ");
            System.out.println("\b\b");
        }
        System.out.println("Total of " + clusters.size() + " clusters identified.");
        
        
        // At this point, we have collected all the clusters into a list.
        // Now output any edges from the original edges that are incident to
        // any 2 vertices from the same cluster:
        ArrayList<String> clusterEdges = identifyClusterEdges(clusters, matrix2);        
        
        
        // Finally, write all these edges to a file:
        File file = new File("cluster_edges.txt");

        
        // If file doesnt exist, then create it
        if (!file.exists()) {
                file.createNewFile();
        }
        

        FileWriter fw = new FileWriter(file.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fw);
        for (String edge : clusterEdges)
            bw.write(edge);
        bw.close();
        
    }
    
    
    /**
     * Normalizes the columns of the matrix, as explained in the PDF.
     * @param matrix 
     */
    public static void normalize(DoubleMatrix matrix) {
        // Normalize the columns:
        for (int j = 0; j < matrix.columns; j++) {
            // First, determine the sum of the column:
            double sum = 0.0;
            for (int i = 0; i < matrix.rows; i++)
                sum += matrix.get(i, j);

            // Now, divide each value by the sum to normalize it:
            for (int i = 0; i < matrix.rows; i++) {
                double avg = matrix.get(i, j) / sum;
                matrix.put(i, j, avg);
            }
        }
    }
    
    /**
     * Exponentiates each column to a given power r and then normalizes them.
     * See the PDF for more info.
     * @param matrix
     * @param r 
     */
    public static void inflate(DoubleMatrix matrix, double r) {
        // Outer loop iterates thru the columns:
        for (int j = 0; j < matrix.columns; j++) {
            
            // Go thru the elements of the current column, raising them to
            // the r-th power.
            for (int i = 0; i < matrix.rows; i++) {
                double d = Math.pow(matrix.get(i, j), r);
                matrix.put(i, j, d);
            }
        }
        
        // Last step is to normalize all the columns:
        normalize(matrix);
    }
    
    
    /**
     * THIS METHOD IS NO LONGER BEING USED BECAUSE OF THE EXTERNAL MATRIX 
     * MULTIPLICATION LIBRARY.
     * 
     * Simple method for cell-by-cell multiplication of 2 identically sized 
     * matrices. This method is not very robust, because it assumes that both
     * matrices are square matrices of the same size. Will crash on any other
     * types of inputs.
     * @param a
     * @param b
     * @return 
     */
    public static double[][] multiply(double a[][], double b[][]) {
        double c[][] = new double[a.length][a.length];
        
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a.length; j++) {
                
                double sum = 0.0;
                for (int k = 0; k < a.length; k++)
                    sum += a[i][k] * b[k][j];
                c[i][j] = sum;
                
            }
        }
        
        return c;
    }
    
    /**
     * Exponentiates a matrix to the p-th power. This is equivalent to multiplying
     * the matrix by itself p times. This method uses a simple form of recursion
     * to perform the exponentiation (divide and conquer).
     * @param matrix
     * @param p
     * @return
     * @throws Exception 
     */
    public static DoubleMatrix exp(DoubleMatrix matrix, int p) throws Exception {
        if (p < 1)
            throw new Exception("Non-positive exponents not supported.");
        if (p == 2) {
            DoubleMatrix result = matrix.mmul(matrix);
            matrix = null;
            return result;
        }
        
        if (p % 2 == 0) {
            DoubleMatrix m = exp(matrix, p/2);
            matrix = null;
            return m.mmul(m);
        }
        else {
            DoubleMatrix m = exp(matrix, (p-1));
            DoubleMatrix result = m.mmul(matrix);
            matrix = null;
            return result;
        }
    }
    
    /**
     * Opens a text file containing a list of edges, and then creates a matrix
     * out of that. Edges will be treated as undirected, and duplicate edges will
     * have no effect on the output of this function. Input format: one edge per
     * line; each edge must be of the form a\tb, where a and b are the numbers
     * (integers) of 2 vertices and \t is the tab character.
     * @param fileName
     * @return
     */
    public static DoubleMatrix readFile(String fileName) throws IOException {
        ArrayList<Integer[]> edges = new ArrayList();
        int largestNode = 0; // largest node encountered so far.
        
        Path path = Paths.get(fileName);
        BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8);
        String line = null;
        String vertices[];
        int v1, v2;
        
        while ((line = br.readLine()) != null) {
            // Store each edge into the arraylist. Also keep track of what is the 
            // largest node value we encounter, because that's how big our matrix
            // will be.
            if (line.contains(" "))
                vertices = line.split(" ");
            else
                vertices = line.split("\t");
            
            v1 = Integer.parseInt(vertices[0]);
            v2 = Integer.parseInt(vertices[1]);
            Integer vertices2[] = {v1, v2};
            edges.add(vertices2);
            
            if (v1 > largestNode)
                largestNode = v1;
            if (v2 > largestNode)
                largestNode = v2;
        }
        br.close();
        
        // After we finish reading the file, it's time to create the matrix and
        // initialize the edges:
        DoubleMatrix matrix = new DoubleMatrix(largestNode + 1, largestNode + 1);
        for (Integer[] i : edges) {
            v1 = i[0];
            v2 = i[1];
            matrix.put(v1, v2, 1.0);
            matrix.put(v2, v1, 1.0);
        }
        return matrix;
    }
    
    
    /**
     * Given the final matrix that has reached a steady state, this method will
     * interpret it to identify the clusters. Returns an ArrayList of ArrayLists.
     * Nodes that belong to the same cluster will be contained in the same list.
     * See the PDF for how this is performed.
     * @param matrix
     * @return 
     */
    public static ArrayList<ArrayList<Integer>> identifyClusters(DoubleMatrix matrix) {
        
        ArrayList<ArrayList<Integer>> result = new ArrayList<ArrayList<Integer>>();
        
        
        for (int i = 0; i < matrix.rows; i++) {
            
            // We use a hashtable to gather the cluster from each row:
            HashMap<String, ArrayList<Integer>> hashtable = new HashMap<String, ArrayList<Integer>>();
            
            for (int j = 0; j < matrix.columns; j++) {
                if (Math.abs(matrix.get(i, j)) < 0.00001)
                    continue;
                
                String d = "" + round(matrix.get(i, j), 5);
                if (hashtable.containsKey(d)) {
                    // update the list in the hashtable:
                    ArrayList<Integer> list = hashtable.get(d);
                    if (!list.contains(j))
                        list.add(j);
                }
                else {
                    // create a new list in the hashtable:
                    ArrayList<Integer> list = new ArrayList<Integer>();
                    list.add(j);
                    hashtable.put(d, list);
                }
            }
            
            // After we're done traversing this row, take the clusters identified 
            // so far and add them to the result to be returned:
            Set<String> keySet = hashtable.keySet();
            for (String k : keySet) {
                ArrayList<Integer> list = hashtable.get(k);
                
                // Output all the numbers in this list.
                // They comprise a single cluster.

                Collections.sort(list);

                ArrayList<Integer> currentCluster = new ArrayList<Integer>();
                for (int x = 0; x < list.size(); x++)
                    currentCluster.add(list.get(x));
                result.add(currentCluster);
            }
        }
        
        // It's possible that the same cluster was identified more than once 
        // (because it might have appeared on more than 1 row of the final 
        // matrix). To remove duplicates, we can sort the list of clusters and
        // then then look for any duplicates (which will be adjacent to each 
        // other in the sorted list).
        sortClusters(result);
        removeDuplicateClusters(result);
        
        return result;
    }
    
    
    /**
    A simple function that rounds a decimal to the specified number of places.
    Adapted from:
    http://stackoverflow.com/questions/2808535/round-a-double-to-2-decimal-places
    */
    public static double round(double value, int places) {
        long factor = (long) Math.pow(10, places);
        value = value * factor;
        long tmp = Math.round(value);
        return (double) tmp / factor;
    }
    
    
    /**
     * Given the list of clusters and the original incidence matrix, this method 
     * returns an ArrayList that contains all the edges that should be outputted.
     * These edges are the ones that are part of the cluster.
     * @param clusters
     * @param incidenceMatrix
     * @return 
     */
    public static ArrayList<String> identifyClusterEdges(ArrayList<ArrayList<Integer>> clusters, DoubleMatrix incidenceMatrix) {
        ArrayList<String> clusterEdges = new ArrayList();
        
        for (ArrayList<Integer> cluster : clusters) {
            for (int i = 0; i < cluster.size() - 1; i++) {
                for (int j = i + 1; j < cluster.size(); j++) {
                    int v1 = cluster.get(i);
                    int v2 = cluster.get(j);
                    
                    if (Math.abs(incidenceMatrix.get(v1, v2) - 1.0) < 0.0000001) {
                        clusterEdges.add(v1 + "\t" + v2 + "\n");
                    }
                }
            }
        }
        return clusterEdges;
    }
    
    /**
     * Sorts the list of clusters in order to easily allow us to detect and 
     * remove duplicates.
     * @param clusters 
     */
    public static void sortClusters(ArrayList<ArrayList<Integer>> clusters) {
        Collections.sort(clusters, new Comparator<ArrayList<Integer>>() {

            @Override
            public int compare(ArrayList<Integer> o1, ArrayList<Integer> o2) {
                if (o1.size() < o2.size())
                    return -1;
                if (o2.size() < o1.size())
                    return 1;
                else {
                    // both clusters have same size. check every element manually:
                    for (int i = 0; i < o1.size(); i++) {
                        if (o1.get(i) < o2.get(i))
                            return -1;
                        if (o2.get(i) < o1.get(i))
                            return 1;
                    }
                    return 0; // in the end, if no mismatches were found, 
                    // then we know these two clusters are the same.
                }
            }
            
        });
    }
    
    /**
     * Given a list of SORTED clusters, this method checks if any two adjacent 
     * clusters are equivalent. If so, then they are duplicates, and the 
     * cluster with the latter index is removed.
     * 
     * Note: the numbers within each cluster itself must be sorted for this 
     * method to work.
     * @param clusters 
     */
    public static void removeDuplicateClusters(ArrayList<ArrayList<Integer>> clusters) {
        int i = 0;
        while (i < clusters.size() - 1) {
            ArrayList<Integer> c1 = clusters.get(i);
            ArrayList<Integer> c2 = clusters.get(i+1);
            
            if (areClustersEqual(c1, c2)) {
                clusters.remove(c2);
            }
            else {
                i++;
            }
            
        }
    }
    
    /**
     * Returns true if the 2 given clusters are equal to each other. The numbers 
     * within each cluster itself must be sorted for this method to work.
     * @param c1
     * @param c2
     * @return 
     */
    public static boolean areClustersEqual(ArrayList<Integer> c1, ArrayList<Integer> c2) {
        if (c1.size() == c2.size()) {
            for (int i = 0; i < c1.size(); i++) {
                if (c1.get(i) != c2.get(i))
                    return false;
            }
            return true;
        }
        else
            return false;
    }
}