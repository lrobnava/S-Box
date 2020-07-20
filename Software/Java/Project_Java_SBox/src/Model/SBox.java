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
public class SBox {

    private int[] sbox;
    private int nonlinearity;
    private int hammingweight;
    private int deltauniformity;
    private double robustness;    
    
    public int getNonlinearity() {
        return nonlinearity;
    }

    public void setNonlinearity(int nonlinearity) {
        this.nonlinearity = nonlinearity;
    }

    public int getHammingweight() {
        return hammingweight;
    }

    public void setHammingweight(int hammingweight) {
        this.hammingweight = hammingweight;
    }

    public int getDeltauniformity() {
        return deltauniformity;
    }

    public void setDeltauniformity(int deltauniformity) {
        this.deltauniformity = deltauniformity;
    }

    public double getRobustness() {
        return robustness;
    }

    public void setRobustness(double robustness) {
        this.robustness = robustness;
    }
        
    public void setSbox(int[] sbox){
        this.sbox = sbox;
    }
    
    public int[] getSbox(){
        return this.sbox;
    }
    
}
