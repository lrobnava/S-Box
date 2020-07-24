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
public class Calculator {
    
    public void calculateRepresentationSBox(SBox sbox){
        
        int n = sbox.getN();
        int m = sbox.getM();
        int input[] = sbox.getSbox();

        Representation representation = new Representation(n, m);
        representation.setInput(input);
        representation.calculateTruthTable();
        representation.calculateLinearCombinations();
        representation.calculateWalshTransform();
        representation.calculateAutocorrelationFast();
        sbox.setRepresentation(representation);
        calculateHammingWeights(sbox);
    }
    
    public int calculateNonLinearity (SBox sbox){
        
        int  i, j, rows, columns, value = 0, temp = 0, min = Integer.MAX_VALUE;
        
        rows = 1 << sbox.getM();
        columns = 1 << sbox.getN() - 1;
        
        Representation represention = sbox.getRepresentation();
        
        int mshift = rows;
        int walshtransformation[][] = represention.getWalshtransform();
        
        int components[] = new int[columns];
        
        for (int k = 0; k < columns; k++) {
            components[k] = 0;
        }
        
        for (j = 0; j < columns; j++){
                if (walshtransformation[0][j] < 0)
                    temp = walshtransformation[0][j] * (-1);
                
                for (i = 1; i < rows; i++){
                    if (Math.abs(walshtransformation[i][j]) > temp)
                            temp = Math.abs(walshtransformation[i][j]);
                }
                
                components[j] = (mshift - temp) >> 1;
                
                if (components[j] < min){
                    min = components[j];
                }
        }
        sbox.setNonlinearity(min);
        return min;
    }  
    
    public void calculateHammingWeights(SBox sbox){
        int m = sbox.getM();
        int mshift = 1 << m;
        int hamming[] = new int[mshift];
        for(int i = 0; i < mshift; i++){
		hamming[i] = calculateHammingWeight(i);            
        }
        sbox.setHammingweight(hamming);
    }
    
    public int calculateHammingWeight(int x){
        int res;
        for (res = 0; x > 0; x = x >> 1)
            res = res + (x & 0x01);
        return res;
    }
    
    public int calculateCorrelationImmunity(SBox sbox){
        int  i, j, columns, order, rows, min = Integer.MAX_VALUE;
        int m =  sbox.getM();
        int n = sbox.getN();
        columns = 1 << n - 1;
        rows = 1 << m;
        int components[] = new int[columns];
        
        Representation representation = sbox.getRepresentation();
        int walshtransform[][];
        walshtransform = representation.getWalshtransform();
        
        int hamming[] = sbox.getHammingweight();
        
        for (j = 0; j < columns; j++){
            order = 1;
            for (i = 1; i < rows; i++){
                if (order == hamming[i] && walshtransform[i][j] != 0){
                        components[j] = order - 1;
                        break;
                }
                if (i == (rows - 1) && order <= m){
                        i = 1;
                        order++;
                }
                components[j] = order - 2;
            }
            if (components[j] < min)
                min = components[j];
        }
        return min;
    }    
    
    public int calculateAbsoluteIndicator(SBox sbox) {
        int i, j, rows, columns;
        int max = 0, temp = 0, temp2 = 0;

        int m = sbox.getM();
        int n = sbox.getN();
        rows = 1 << m;
        columns = 1 << n - 1;
        
        Representation representation = sbox.getRepresentation();
        int autocorrelationfast[][] = representation.getAutocorrelationfast();
        
        int components[] = new int[columns];

        for (j = 0; j < columns; j++)
        {
                temp = Math.abs(autocorrelationfast[1][j]); //disregard first value since it is 2^n
                for (i = 2; i < rows; i++)	
                {
                        temp2 = Math.abs(autocorrelationfast[i][j]);
                        if (temp2 > temp)
                                temp = temp2;
                }
                components[j] = temp;
                if (temp > max)
                        max = temp;
        }
        return max;
    }
    
    public int calculateSumOfSquareIndicator (SBox sbox) //the smaller, the bettter
    {
        int i, j, rows, columns;
        int max = 0, sum = 0;

        int m = sbox.getM();
        int n = sbox.getN();
        
        rows = 1 << m;
        columns = 1 << n - 1;
        
        Representation representation = sbox.getRepresentation();
        int autocorrelationfast[][] = representation.getAutocorrelationfast();
        
        int components[] = new int[columns];
        
        for (j = 0; j < columns; j++)
        {
            for (i = 0; i < rows; i++)	
            {
                sum += autocorrelationfast[i][j] * autocorrelationfast[i][j];
            }
            components[j] = sum;
            if (sum > max)
                max = sum;
            sum = 0;
        }
        return max;
    }

    public int calculateAlgebraicDegree (SBox sbox)
    {
        int  tmp, weight, deg;
        int i, j, rows, columns, max = 0;

        Representation representation = sbox.getRepresentation();
        int algebraicnormalform[][] = representation.getAlgebraicnormalform();
        
        int m = sbox.getM();
        int n = sbox.getN();
        
        rows = 1 << m;
        columns = 1 << n - 1;
        
        int components[] = new int[columns];

        for (j = 0; j < columns; j++)
        {
            if (algebraicnormalform[rows - 1][j] != 0)
                deg = m;
            else
            {
                for (deg = 0, i = 1; i < (rows - 1); ++i)
                    if (algebraicnormalform[i][j] != 0) 
                    {
                        for (weight = 0, tmp = i; tmp > 0; tmp >>= 1)
                            weight = weight + tmp%2;
                        if (weight > deg)
                            deg = weight;
                    }
            }
            components[j] = deg;
            if (components[j] > max)
                max = components[j];
        }
        return max;
    }
}

