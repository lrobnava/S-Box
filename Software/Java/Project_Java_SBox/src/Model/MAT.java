/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Model;

/**
 *
 * @author skyli
 */
public class MAT {

    private int n;
    private int m;
    private int vector[][];
    
    public void createMAT(int n, int m){
        this.n = n;        
        this.m = m;
        this.vector = new int[n][m];
    }
    
    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public int getM() {
        return m;
    }

    public void setM(int m) {
        this.m = m;
    }

    public int[][] getVector() {
        return vector;
    }

    public void setVector(int[][] vector) {
        this.vector = vector;
    }
    
}
