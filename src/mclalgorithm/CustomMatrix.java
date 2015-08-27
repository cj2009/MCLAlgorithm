package mclalgorithm;

import org.jblas.DoubleMatrix;

/**
 * This is a custom class that acts as a simple 2-dimensional matrix for double 
 * values. This is created as a template that should be implemented with any 
 * matrix library of one's choice. This makes the project more modular and allows 
 * one to easily swap in and out any matrix libraries for testing, without 
 * needing to change the MCAlgorithm class.
 * 
 * This class will require the following essential methods:
 *  + put(int row, int column, double value)
 *  + get(int row, int value)
 *  + zero()
 *  + square()
 *  + multiply()
 *  + clone()
 *  + rows()
 *  + cols()
 * 
 * The current implementation is using the JBLAS library (www.jblas.org).
 * @author chris_joseph
 */
public class CustomMatrix {
    private DoubleMatrix m;
    
    public CustomMatrix(int rows, int columns) {
        m = DoubleMatrix.zeros(rows, columns);
    }
    
    private CustomMatrix() {
        
    }
    
    public double get(int row, int col) {
        return m.get(row, col);
    }
    
    public void put(int row, int col, double value) {
        m.put(row, col, value);
    }
    
    public void zero() {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                m.put(i, j, 0.0);
            }
        }
    }
    
    public CustomMatrix square() {
        CustomMatrix result = new CustomMatrix();        
        DoubleMatrix dm = m.mmul(m);
        result.m = dm;
        return result;
    }
    
    public CustomMatrix clone() {
        CustomMatrix result = new CustomMatrix();
        DoubleMatrix dm = m.dup();
        result.m = dm;
        return result;
    }
    
    public int rows() {
        return m.rows;
    }
    
    public int cols() {
        return m.columns;
    }
    
    public CustomMatrix multiply(CustomMatrix matrix) {
        CustomMatrix result = new CustomMatrix();
        
        DoubleMatrix dm = m.mmul(matrix.m);
        result.m = dm;
        return result;
    }
    
    public void print() {
        for (int i = 0; i < m.rows; i++) {
            for (int j = 0; j < m.columns; j++) {
                System.out.printf("%.3f\t", m.get(i, j));
            }
            System.out.printf("\b\n");
        }
    }
}
