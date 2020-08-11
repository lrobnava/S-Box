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
    private int input[];
    private int linearcombinations[][];
    private int walshtransform[][];
    private int autocorrelationfast[][];
    private int algebraicnormalform[][];
    
    public Representation(int n, int m){
        this.n = n;
        this.m = m;
        this.nshift = 1 << n;
        this.mshift = 1 << m;
        this.truthtable = new int[this.mshift][this.n];
        this.input = new int[this.mshift];
        this.linearcombinations = new int[this.mshift][this.nshift - 1];
        this.walshtransform = new int[this.mshift][this.nshift - 1];
        this.autocorrelationfast = new int[this.mshift][this.nshift - 1];
        this.algebraicnormalform = new int[this.mshift][this.nshift - 1];
    }
    
    public void calculateTruthTable(){
    int i, j;
    for (i = 0; i < this.mshift; i++)
        for (j = 0; j < this.n; j++)
            this.truthtable[i][j] = this.input[i] >> (this.n - j - 1) & 0x01;
    }
    
    public void calculateLinearCombinations ()
    {
	int i, j, shift;      
        HelperFunctions hf = new HelperFunctions();
        
	if (this.n == 1) //for Boolean function there is no linear combinations
            for (i = 0; i < this.mshift; i++)
                this.linearcombinations[i][0] = this.truthtable[i][0]; 
	else
	{
            for (i = 1; i < this.nshift; i++){ // start from 1, 0 does not have any interesting case, all linear combinations
                for (j = 0; j < this.n; j++){
                    shift = (i >> j) & 0x01;
                    if (0 < shift) {
                        hf.xor_elements( this.linearcombinations, 
                        this.truthtable, 
                        this.linearcombinations, 
                        i - 1, 
                        j, 
                        this.mshift );
                    }
                }
            }
	}
    }
    
    public int calculateWalshTransform(){
        int m, halfm, t1, t2, a, b, r, j;
	int i, z, rows, columns;
	rows = this.mshift;
	columns = this.nshift - 1;

	for (z = 0; z < columns; z++){
            for (i = 0; i < rows; ++i)
                this.walshtransform[i][z] = (this.linearcombinations[i][z] == 0) ? 1 : -1;

            for (i = 1; i <= this.m; ++i) {
                m  = (1 << i);
                halfm = m/2;
                for (r = 0; r < (int) rows; r += m) {
                    t1 = r;
                    t2 = r + halfm;
                    for (j = 0; j < halfm; ++j, ++t1, ++t2) {
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
    
    public int calculateAutocorrelationFast (){
        int m, halfm, t1, t2, r, a, b, j;
        int i, z, rows, columns;
        rows = this.mshift;
        columns = this.nshift - 1;

        for (z = 0; z < columns; z++){
            for (i = 0; i < rows; i++)
                this.autocorrelationfast[i][z] = -1* this.walshtransform[i][z]*this.walshtransform[i][z];

            for (i = 1; i <= this.m; ++i) {
                m  = (1 << i);
                halfm = m/2;
                for (r = 0; r < (int) rows; r += m) {
                    t1 = r;
                    t2 = r + halfm;
                    for (j = 0; j < halfm; ++j, ++t1, ++t2) {
                        a = this.autocorrelationfast[t1][z];
                        b = this.autocorrelationfast[t2][z];
                        this.autocorrelationfast[t1][z] = a + b;
                        this.autocorrelationfast[t2][z] = a - b;
                    }
                }
            }
            for (i = 0; i < rows; i++)
                this.autocorrelationfast[i][z] /= (1 << this.m)*(-1);
        }
        return 1;
    }    
 
    public int calculateAlgebraicNormalForm (){
        int tmp = 0, res = 0;
        int i, j, k, z, rows, columns;

        rows = 1 << (this.m - 1);
        columns = this.nshift -1;

        int t[] = new int[rows]; 
        int u[] = new int[rows];

        for (z = 0; z < columns; z++){
            for (i = 0; i < (rows << 1); ++i)
                this.algebraicnormalform[i][z] = (int) this.linearcombinations[i][z];
            for (i = 0; i < rows; ++i)
                u[i] = t[i] = 0;

            for (i = 0; i < this.m; ++i){
                for (j = 0; j < rows; ++j){
                    t[j] = this.algebraicnormalform[2*j][z];
                    u[j] = (this.algebraicnormalform[2*j][z] == this.algebraicnormalform[2*j+1][z]) ? 0 : 1;
                }
                for (k = 0; k < rows; ++k){
                    this.algebraicnormalform[k][z] = t[k];
                    this.algebraicnormalform[rows + k][z] = u[k];
                }
            }
        }
        return 1;
    }
    
    public int[][] getTruthTable(){
        return this.truthtable;
    }
    
    public void setTruthTable(int truthtable[][]){
        this.truthtable = truthtable;
    }    

    public int getM() {
        return m;
    }

    public void setM(int m) {
        this.m = m;
    }

    public int getN() {
        return n;
    }

    public void setN(int n) {
        this.n = n;
    }

    public int getMshift() {
        return mshift;
    }

    public void setMshift(int mshift) {
        this.mshift = mshift;
    }

    public int getNshift() {
        return nshift;
    }

    public void setNshift(int nshift) {
        this.nshift = nshift;
    }

    public int[] getInput() {
        return input;
    }

    public void setInput(int[] input) {
        this.input = input;
    }

    public int[][] getLinearcombinations() {
        return linearcombinations;
    }

    public void setLinearcombinations(int[][] linearcombinations) {
        this.linearcombinations = linearcombinations;
    }

    public int[][] getWalshtransform() {
        return walshtransform;
    }

    public void setWalshtransform(int[][] walshtransform) {
        this.walshtransform = walshtransform;
    }

    public int[][] getAutocorrelationfast() {
        return autocorrelationfast;
    }

    public void setAutocorrelationfast(int[][] autocorrelationfast) {
        this.autocorrelationfast = autocorrelationfast;
    }

    public int[][] getAlgebraicnormalform() {
        return algebraicnormalform;
    }

    public void setAlgebraicnormalform(int[][] algebraicnormalform) {
        this.algebraicnormalform = algebraicnormalform;
    }

}
