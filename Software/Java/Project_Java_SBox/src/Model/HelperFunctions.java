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
    
    public void xor_elements (int _matrix1[][], int truth_table[][], int _matrix2[][], int column_ll, int column_tt, int m_shift)
    {
        int i;
        for (i = 0; i < m_shift; i++)
            _matrix2[i][column_ll] = _matrix1[i][column_ll] ^ truth_table[i][column_tt];
    }

}
