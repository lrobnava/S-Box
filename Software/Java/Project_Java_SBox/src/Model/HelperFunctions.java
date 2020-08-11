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
public class HelperFunctions {
    
    public void xor_elements (int _matrix1[][], int truth_table[][], int _matrix2[][], int column_ll, int column_tt, int m_shift){
        int i;
        for (i = 0; i < m_shift; i++)
            _matrix2[i][column_ll] = _matrix1[i][column_ll] ^ truth_table[i][column_tt];
    }

    //by Frederic Laffite, from boolfun
    public int choose (int n, int k){
        
        int i, num = 1, den = 1;

        if (k < 0 || k > n) 
            return 0;
        for (i = 0; i < k; ++i) {
            num *= n--;
            den *= (k-i); 
        }
        
        return (num/den);
    }    
    
    public int[] sortIncreasingDeg (int v[], int hamming[])
    {
        int i,j,tmp,len;
        len = v.length;
        for (i = 0; i < len-1; ++i)
            for (j = i+1; j < len; ++j)
                //if(hamming_weight(v[j]) < hamming_weight(v[i])) 
                if(hamming[v[j]] < hamming[v[i]]) 
                {
                    tmp = v[j];
                    v[j] = v[i];
                    v[i] = tmp;
                }
        return v;
    }    
    
    public MAT getMatrix (int tt[][], int m, MAT mat, int monomials[], int Nm, int ai, int b)
    {
        int Ns, i, j, len;
        len = 1 << m;
        int[] support = getSupport(tt, m, b);
        Ns = support.length;
        int[][] v = null;
        if (Ns == 0 || Ns == len)
            mat = null;
        else 
            {
            if(Nm > Ns) 
                mat.createMAT(Nm, Nm);
            else
                mat.createMAT(Nm, Ns);
            v = mat.getVector();
            for (i = 0; i < Nm; ++i)
                for (j = 0; j < Ns; ++j)
                    v[i][j] = preceq(monomials[i], support[j]);
            mat.setVector(v);
        }
        
        return mat;
    }    
    
    public int[] getSupport (int[][] tt, int n, int b)
    {
        int i, k, len, N;

        len = 1 << n;
        for (N = 0, i = 0; i < len; ++i)
            if(tt[i][0] != b)
                N = N + 1;
        int[] res = new int[N];
        for (k = 0, i = 0; i < len; ++i)
            if (tt[i][0] != b)
                res[k++] = i;
        return res;
    }    
    
    public int preceq (int a, int b)
    {
        int res = 1;

        while ((a > 0 || b > 0) && (res ==1)) 
        {
            if ((a & 1) > (b & 1)) 
                res = 0;
            a >>= 1; 
            b >>= 1; 
        }
        return res;
    }
    
    public int solveMatrix (MAT m, int[] hamming, int[] monomials, int b)
    {
        int _n = m.getN();
        int _m = m.getM();
        int _v[][] = m.getVector();
        int i, j, l, res, processed_lines, zero_lines;
        int deg[] = new int[_n];
        
        for (res = 0, i = 0; i < _n; ++i) 
        {
            //deg[i] = hamming_weight(monomials[i]);
            deg[i] = hamming[monomials[i]];
            if (deg[i] > res) 
                res = deg[i];
        }
        processed_lines = zero_lines = 0;
        for (i = 0; i < _n; ++i) 
        {
            for (j = 0; j < _m && _v[i][j] == 0; ++j);
            if (j == _m) 
            {
                ++zero_lines;
                if (deg[i] < res && deg[i] != 0)
                    res = deg[i];
            } 
            else 
            {
                ++processed_lines;
                if (i != j && i < _m && j < _m)
                    m = swap_columns(m, i, j);
                for (l = i+1; l < _n && i < _m; ++l) 
                {
                    if (i < _m && _v[l][i] != 0) 
                    {
                        m = addLine(m, l, i);
                        deg[l] = (deg[i] > deg[l]) ? deg[i] : deg[l];
                    }
                }
            }
        }

        return res;
    }

    public MAT swap_columns (MAT mat, int a, int b)
    {
        int i;
        int _n = mat.getN();
        int _v[][] = mat.getVector();
        int tmp[] = new int[_n];
        for (i = 0; i < _n; ++i) 
            tmp[i] = _v[i][a];
        for (i = 0; i < _n; ++i) 
            _v[i][a] = _v[i][b];
        for (i = 0; i < _n; ++i) 
            _v[i][b] = tmp[i];
        mat.setVector(_v);
        return mat;
    }

    public MAT addLine (MAT mat, int dst, int src)
    {
        int j;
        int _m = mat.getM();
        int[][] _v = mat.getVector();
        for (j = 0; j < _m; ++j)
            _v[dst][j] = (_v[dst][j] + _v[src][j]) & 1;
        mat.setVector(_v);
        return mat;
    }
    
    public int innerProduct (int a, int b, int M)
    {
        int i;
        int res = 0;
        for (i = 0; i < M; i++)
        {
            res ^=  (((a >> i) & 0x01) & ((b >> i) & 0x01));
        }
        return res;
    }    
    
}
