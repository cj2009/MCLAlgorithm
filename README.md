# MCLAlgorithm

This NetBeans project contains my implementation of the popular Markov Clustering (MCL) algorithm in Java.

MCLAlgorithm:
 * This is a class that performs the MCL algorithm on a graph.
 * It works only on undirected, unweighted graphs.
 * Based on the original MCL Algorithm by Stijn van Dongen.
 * Also based on the presentation by Kathy Macropol (UCSB): (www.cs.ucsb.edu/~xyan/classes/CS595D-2009winter/MCL_Presentation2.pdf)
 
MCLTesting:
 * This is a program that performs the MCL algorithm in Java, while also comparing the matrix from each iteration to the matrix produced by our Hadoop program in each iteration. This helps to verify that both algorithms are producing exactly the same output after each iteration.
