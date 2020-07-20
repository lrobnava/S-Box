/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Model;

public class Representation {
    
    private int m;
    private int n;
    private int mshift;
    private int nshift;
    private int truthtable [][];
    private final int input[];
    private final int linearcombinations[][];
    private final int walshtransform[][];
    private final int autocorrelationfast[][];
    private final int algebraicnormalform[][];
    
    public Representation(int n, int m){
        this.n = n;
        this.m = m;
        this.truthtable = new int[1 << m][1 << n];
        this.input = new int[1 << m];
        this.linearcombinations = new int[1 << m][1 << n - 1];
        this.walshtransform = new int[1 << m][1 << n - 1];
        this.autocorrelationfast = new int[1 << m][1 << n - 1];
        this.algebraicnormalform = new int[1 << m][1 << n - 1];
    }
    
    public int[][] getTruthTable(){
        return this.truthtable;
    }
    
    public void setTruthTable(int truthtable[][]){
        this.truthtable = truthtable;
    }
    
    public void calculateTruthTable(int m_shift, int n){
    int i, j;
    for (i = 0; i < m_shift; i++)
        for (j = 0; j < n; j++)
            this.truthtable[i][j] = this.input[i] >> (n - j - 1) & 0x01;
    }
    
    public void linear_combinations ()
    {
        int m_shift = 1 << this.m;
        int n_shift = 1 << this.n;
	int i, j, shift;
        
        HelperFunctions hf = new HelperFunctions();
        
	if (n == 1) //for Boolean function there is no linear combinations
            for (i = 0; i < m_shift; i++)
                this.linearcombinations[i][0] = this.truthtable[i][0]; 
	else
	{
            for (i = 1; i < n_shift; i++){ // start from 1, 0 does not have any interesting case, all linear combinations
                for (j = 0; j < n; j++){
                    shift = i >> j & 0x01;
                    if (1==shift) {
                        hf.xor_elements(this.linearcombinations, 
                                this.truthtable, 
                                this.linearcombinations, 
                                m, 
                                n_shift, 
                                m_shift);
                    }
                }
            }
	}
    }
    
    public int walsh_transform()
    {
        int m, halfm, t1, t2, a, b, r, j;
	int i, z, rows, columns;
	rows = this.mshift;
	columns = this.nshift - 1;

	for (z = 0; z < columns; z++)
	{
            for (i = 0; i < rows; ++i)
                this.walshtransform[i][z] = (this.linearcombinations[i][z] == 0) ? 1 : -1;

            for (i = 1; i <= this.m; ++i) 
            {
                m  = (1 << i);
                halfm = m/2;
                for (r = 0; r < (int) rows; r += m) 
                {
                    t1 = r;
                    t2 = r + halfm;
                    for (j = 0; j < halfm; ++j, ++t1, ++t2) 
                    {
                        a = this.walshtransform[t1][z];
                        b = this.walshtransform[t2][z];
                        this.walshtransform[t1][z] = a + b;
                        this.walshtransform[t2][z] = a - b;
                    }
                }
            }
	}
	return 1;
    }
    
}
