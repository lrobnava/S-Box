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
    
    public void calculateRepresentationSBox(SBox sbox) {
        int n = sbox.getN();
        int m = sbox.getM();
        int input[] = sbox.getSbox();

        Representation representation = new Representation(n, m);
        representation.setInput(input);
        representation.calculateTruthTable();
        representation.calculateLinearCombinations();
        representation.calculateWalshTransform();
        representation.calculateAlgebraicNormalForm();
        representation.calculateAutocorrelationFast();
        sbox.setRepresentation(representation);
        calculateHammingWeights(sbox);
    }
    
    public int calculateNonLinearity (SBox sbox) {
        int  i, j; 
        int rows, columns, mshift; 
        int temp = 0, min = Integer.MAX_VALUE;
        int walshtransformation[][], components[];
        
        Representation represention = sbox.getRepresentation();
        rows = 1 << sbox.getM();
        columns = 1 << sbox.getN() - 1;
        
        mshift = rows;
        walshtransformation = represention.getWalshtransform();
        components = new int[columns];
        
        for (int k = 0; k < columns; k++) 
            components[k] = 0;
        
        for (j = 0; j < columns; j++) {

            if (walshtransformation[0][j] < 0)
                temp = walshtransformation[0][j] * (-1);

            for (i = 1; i < rows; i++) {
                if (Math.abs(walshtransformation[i][j]) > temp)
                    temp = Math.abs(walshtransformation[i][j]);
            }

            components[j] = (mshift - temp) >> 1;

            if (components[j] < min)
                min = components[j];

        }

        sbox.setNonlinearity(min);
        return min;
    }  
    
    public void calculateHammingWeights(SBox sbox) {
        int m = sbox.getM();
        int mshift = 1 << m;
        int hamming[] = new int[mshift];
        
        for(int i = 0; i < mshift; i++)
            hamming[i] = calculateHammingWeight(i);            

        sbox.setHammingweight(hamming);
    }
    
    public int calculateHammingWeight(int x) {
        int res;

        for (res = 0; x > 0; x = x >> 1)
            res = res + (x & 0x01);

        return res;
    }
    
    public int calculateCorrelationImmunity(SBox sbox) {
        int  i, j, columns, order, rows;
        int min = Integer.MAX_VALUE, m =  sbox.getM(), n = sbox.getN();
        int hamming[], components[];
        int walshtransform[][];
        Representation representation = sbox.getRepresentation();
        
        columns = 1 << n - 1;
        rows = 1 << m;
        
        components = new int[columns];
        walshtransform = representation.getWalshtransform();
        hamming = sbox.getHammingweight();
        
        for (j = 0; j < columns; j++) {
            order = 1;
            for (i = 1; i < rows; i++) {
                if (order == hamming[i] && walshtransform[i][j] != 0){
                    components[j] = order - 1;
                    break;
                }
                if (i == (rows - 1) && order <= m) {
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
        int max = 0, temp = 0, temp2 = 0, m = sbox.getM(), n = sbox.getN();
        int components[];
        int autocorrelationfast[][];
        Representation representation = sbox.getRepresentation();
        
        rows = 1 << m;
        columns = 1 << n - 1;
        
        autocorrelationfast = representation.getAutocorrelationfast();
        components = new int[columns];

        for (j = 0; j < columns; j++) {
            
            temp = Math.abs(autocorrelationfast[1][j]); //disregard first value since it is 2^n
            
            for (i = 2; i < rows; i++) {
                
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
    
    public int calculateSumOfSquareIndicator (SBox sbox) { //the smaller, the bettter
        int i, j, rows, columns;
        int max = 0, sum = 0, m = sbox.getM(), n = sbox.getN();
        int components[];
        int autocorrelationfast[][];
        Representation representation = sbox.getRepresentation();
        
        rows = 1 << m;
        columns = 1 << n - 1;
        
        autocorrelationfast = representation.getAutocorrelationfast();
        components = new int[columns];
        
        for (j = 0; j < columns; j++) {
            
            for (i = 0; i < rows; i++)                 
                sum += autocorrelationfast[i][j] * autocorrelationfast[i][j];

            components[j] = sum;

            if (sum > max)
                max = sum;

            sum = 0;

        }

        return max;
    }

    public int calculateAlgebraicDegree (SBox sbox) {
        int i, j, rows, columns, tmp, weight, deg; 
        int max = 0, m = sbox.getM(), n = sbox.getN();
        int components[];
        int algebraicnormalform[][];
        Representation representation = sbox.getRepresentation();
        algebraicnormalform = representation.getAlgebraicnormalform();
        
        rows = 1 << m;
        columns = 1 << n - 1;
        
        components = new int[columns];

        for (j = 0; j < columns; j++) {
            
            if (algebraicnormalform[rows - 1][j] != 0)
                deg = m;
            else {
                
                for (deg = 0, i = 1; i < (rows - 1); ++i) {
                    
                    if (algebraicnormalform[i][j] != 0) {
                        
                        for (weight = 0, tmp = i; tmp > 0; tmp >>= 1)
                            weight = weight + tmp%2;
                        
                        if (weight > deg)
                            deg = weight;
                    
                    }
                
                }
            
            }
            
            components[j] = deg;
            
            if (components[j] > max)
                max = components[j];
        
        }
        
        return max;
    }
    
    public int calculateAlgebraicImmunity(SBox sbox) {//u8 **ll, uint *components)
        int a, b, deg, i, columns, rows, Nm; 
        int res = 0, min = Integer.MAX_VALUE, n = sbox.getN(), m = sbox.getM();
        int hamming[], components[], monomials[];
        int algebraicnormalform[][];
        HelperFunctions helperfunction = new HelperFunctions();
        Representation representation = sbox.getRepresentation();
        
        hamming = sbox.getHammingweight();
        algebraicnormalform = representation.getAlgebraicnormalform();
        
        MAT m0 = new MAT();
        MAT m1 = new MAT();

        rows = 1 << m;
        columns = (1 << n) - 1;

        components = new int[columns];

        deg = (m >> 1) + (m & 1);
        monomials = calculateMonomials(m, deg, hamming);
        Nm = monomials.length;
        
        monomials = helperfunction.sortIncreasingDeg(monomials, hamming);

        for (i = 0; i < columns; i++) {
            
            m0 = helperfunction.getMatrix(algebraicnormalform, m, m0, monomials, Nm, deg, 0);
            
            if (m0 == null)
                res = 0;
            
            else {
                
                m1 = helperfunction.getMatrix(algebraicnormalform, m, m1, monomials, Nm, deg, 1);
                a = helperfunction.solveMatrix(m0, hamming, monomials, 0);
                b = helperfunction.solveMatrix(m1, hamming, monomials, 1);
                res = (a < b) ? a : b;
            
            }
            
            components[i] = res;
            
            if (res < min)
                min = res;
        
        }

        return min;
    }    
    
    public int[] calculateMonomials (int n, int d, int hamming[]) {
        int i, k, N;//Tenemos que devolver N mayuscula
        int res[]; //Tenemos que devolver res
        HelperFunctions helperfunctions = new HelperFunctions();
        
        for (N = 0, k = 0; k <= d; ++k)
            N = N + helperfunctions.choose(n, k);
        
        res = new int[N];
        
        for (k = 0, i = 0; i<(1<<n); ++i)
            //if (hamming_weight(i) <= d)
            if (hamming[i] <= d)
               res[k++] = i;
        
        return res;
    }    
    
    public int calculatePropagationCharacteristics(SBox sbox) {//the higher the better{
        int i, j, z, order = 1, count = 0;
        int help[];
        int mshift = 1 << sbox.getM();
        int hamming[] = sbox.getHammingweight();
        int inputarray[] = sbox.getSbox();
        int N = sbox.getN();
        
        help = new int[mshift];

        do {
            for (i = 1; i < mshift; i++) {
                //if (order == hamming_weight(i))
                if (order == hamming[i]) {
                    
                    calculateDerivative(inputarray, help, i, mshift);
                    //could be done with tt, wt, balance
                    for (j = 0; j < mshift; j++) {
                        
                        for (z = 0; z < mshift; z++) {
                            
                            if (help[z] == j)
                                count++;
                        
                        }
                        
                        if (count != 1)
                            return order - 1;
                        
                        count = 0;
                    }
                }
            }
            order++;
        } while (order <= N);
        return order-1;
    }    
    
    public void calculateDerivative (int input_array[], int output_array[], int shift, int columns) {
        
        int i;
        
        for (i = 0; i < columns; i++)
            output_array[i] = input_array[i] ^ input_array[i ^ shift];
    
    }
    
    public int calculateNumFixedPoints (SBox sbox) {
        
        int inputarray[] = sbox.getSbox();
        int i, res = 0, rows = 1 << sbox.getM();

        for (i = 0; i < rows; i++) {
            
            if ((inputarray[i] ^ i) == 0) //if (output)
                res++; //System.out.println("Fixed point is %x on position %d.\n", inputarray[i], i);

        }
        return res;
    }
    
    public int calculateNumOppositeFixedPoints (SBox sbox) {
        
        int i, res = 0, rows = 1 << sbox.getM();
        int inputarray[] = sbox.getSbox();
        
        for (i = 0; i < rows; i++) {
            
            if (inputarray[i] == (~i & (rows - 1))) //if (output)
                res++; //printf("Opposite fixed point is %x on position %d.\n", input_array[i], i);

        }

        return res;

    }    
    
    
    public float calculateComputeKappaCPA(int inputBits, int outputBits, int keyBits, SBox sbox) { 

        int inputSize = 1 << inputBits;
        int outputSize = 1 << outputBits;
        int keySize = 1 << keyBits;
        int coefficientSize = (keySize * (keySize-1)) >> 1; // combinatorics on key picks i and j
        int confusionCounter = 0; //# of times we observe confusion
        int outi,outj; // Sbox output	
        int coefficientCounter = 0; //how many coefficients have we computed
        int keyi, keyj, input, i = 0, j = 0, k = 0, coefficientSum = 0, temp = 0, fcounter = 0, in = 0;
        float mean = 0, var = 0;

        int inputarray[] = sbox.getSbox();
        int hamming[] = sbox.getHammingweight();
        //float reducedCoefficients[] = new float[coefficientSize];
        //float frequency[] = new float[coefficientSize];
        float confusionCharacteristic[] = new float[coefficientSize]; 

        for (keyi = 0; keyi < keySize; keyi++) {

            for (keyj = keyi + 1; keyj < keySize; keyj++) {
                
                for (input = 0; input < inputSize ; input++) {

                    outi = inputarray[keyi ^ input];
                    outj = inputarray[keyj ^ input];

                    temp = hamming[outi] - hamming[outj]; 
                    temp *=  temp;
                    coefficientSum += temp;

                }	
                //input set is over, lets compute confusion coefficient for (keyi, keyj)
                confusionCharacteristic[coefficientCounter]=(float)coefficientSum / (float)inputSize;
                coefficientCounter ++;
                confusionCounter = 0; //dpa
                coefficientSum = 0; //cpa
            
            }
        }

        fcounter = 0;
        for (i = 0; i < coefficientSize; i++)
        {
                if(confusionCharacteristic[10] == confusionCharacteristic[i])
                {
                        fcounter++;
                }
                mean = mean + confusionCharacteristic[i];
        }

        mean = mean / (float)coefficientSize;

        for (i = 0; i < coefficientSize; i++)
                var = (float) (var + (float) Math.pow((confusionCharacteristic[i] - (float) mean), 2));

        var = var / (float)coefficientSize;

        return var;
    }    
    
    public float calculateSNR_DPA (SBox sbox)
    {
        HelperFunctions helperfunctions = new HelperFunctions();
        int tt[][] = sbox.getRepresentation().getTruthTable();
        int mshift = 1 << sbox.getM();
        int M = sbox.getM();
        int N = sbox.getN();
        float sum1 = 0, sum2 = 0, sum3 = 0, help = 0;
        int k, i, x, temp;

        temp = (1 << (2 * M)) * N;

        for (k = 0; k < mshift; k++)
        {
            sum2 = 0;
            for (i = 0; i < N; i++)
            {
                sum3 = 0;
                for (x = 0; x < mshift; x++)
                {
                    sum3 += ((1 - 2 * helperfunctions.innerProduct(x, k, M)) * (1 - 2 * tt[x][i]));
                }
                sum2 += sum3;
            }
            sum2 = sum2 * sum2 * sum2 * sum2;
            sum1 += sum2;
        }
        sum1 = (float) Math.pow (sum1, (float) - 0.5);
        return sum1 * temp;
    }    
    
    public int calculateDeltaUniformity (SBox sbox)
    {
        
        int i, j, delta = 0, R = 0;
        int mshift = 1 << sbox.getM();
        int nshift = 1 << sbox.getN();
        int DDT[][] = new int[mshift][nshift];
        int inputarray[] = sbox.getSbox();

        for (i = 0; i < mshift; i++)
        {
            for (j = 0 ; j < nshift; j++)
            {
                DDT [i ^ j] [inputarray[i] ^ inputarray[j]]++;			
            }
        }

        for (i = 0; i < mshift; i++)
        {
            for (j = 0 ; j < nshift; j++)
            {
                if (DDT [i][j] > delta && (i != 0 && j != 0))
                    delta = DDT [i][j];	
            }
        }

        for (i = 1; i < nshift; i++)
        {
            if (DDT [i][0] != 0)
                R++;
        }

        float robustness = (1 - (R / (float)mshift)) * (float)(1 - (delta / (float)mshift));

        return delta;
    }
    
    public float calculateRobustness (SBox sbox)
    {
        
        int i, j, delta = 0, R = 0;
        int mshift = 1 << sbox.getM();
        int nshift = 1 << sbox.getN();
        int DDT[][] = new int[mshift][nshift];
        int inputarray[] = sbox.getSbox();

        for (i = 0; i < mshift; i++)
        {
            for (j = 0 ; j < nshift; j++)
            {
                DDT [i ^ j] [inputarray[i] ^ inputarray[j]]++;			
            }
        }

        for (i = 0; i < mshift; i++)
        {
            for (j = 0 ; j < nshift; j++)
            {
                if (DDT [i][j] > delta && (i != 0 && j != 0))
                    delta = DDT [i][j];	
            }
        }

        for (i = 1; i < nshift; i++)
        {
            if (DDT [i][0] != 0)
                R++;
        }

        float robustness = (1 - (R / (float)mshift)) * (float)(1 - (delta / (float)mshift));
        
        return robustness;
    }    
    
    
    public float calculateTransparency (SBox sbox)
    {
        int N = sbox.getN();
        int M = sbox.getM();
        int nshift = 1 << N;
        int mshift = 1 << M;
        int hamming[] = sbox.getHammingweight();
        int inputarray[] = sbox.getSbox();
        int b, a, v, x, tmp = 0;
        float res = 0;
        float sigma1, sigma2, sigma3 = 0, C;
        HelperFunctions helperfunctions = new HelperFunctions();

        if (N > 1)
            return calculateTransparencyFast (sbox);

        C = (float) ((1 << (2 * M)) - mshift);

        for (b = 0; b < nshift; b++)
        {
            sigma1 = 0;
            for (a = 1; a < mshift; a++)
            {
                    sigma2 = 0;
                    for (v = 1; v < nshift; v = v << 1) //can go from 1 since we need Hamming weight to be 1
                    {
                            sigma3 = 0;
                            for (x = 0; x < mshift; x++)
                            {
                                    sigma3 += (float) (1 - 2 * helperfunctions.innerProduct(v, inputarray[x]^inputarray[x^a], M));
                            }
                            sigma2 += (float) (1 - 2 * helperfunctions.innerProduct(v, b, M)) * sigma3;
                    }

                    if (sigma2 < 0)
                            sigma2 = sigma2 * (-1);

                    sigma1 += sigma2;
            }
            //tmp = 2 * hamming_weight(b);
            tmp = 2 * hamming[b];
            if (res < Math.abs((int)N - tmp) - sigma1/C)
                res = Math.abs((int)N - tmp) - sigma1/C;
        }
        return res;
    }    
    
    public float calculateTransparencyFast(SBox sbox)
    {
            int tmp = 0;
            int b, a, x;
            int N = sbox.getN();
            int M = sbox.getM();
            int mshift = 1 << M;
            int nshift = 1 << N;
            int hamming[] = sbox.getHammingweight();
            int inputarray[] = sbox.getSbox();
            float res = 0, z = 0, temp = 0, sigma1 = 0, sigma2 = 0, treshold = 0;
            float C = (float) ((1 << (2 * N)) - nshift);
            float K = (float) (N * nshift);

            for (b = 0; b < nshift; b++)
            {
                    treshold = (float) (N - 2 * hamming[b]);
                    if (treshold < 0)
                            treshold = treshold * (-1);
                    treshold = (treshold - res) * C;
                    if (treshold >= 0)
                    {
                            sigma2 = 0;
                            for (a = 1; a < nshift; a++)
                            {
                                    sigma1 = 0;
                                    for (x = 0; x < nshift; x++)
                                    {
                                            //sigma1 = sigma1 + hamming_weight(b ^ (evaluate_box(tt,x,a, N)));
                                            //sigma1 = sigma1 + hamming[b ^ (evaluate_box(tt, x, a))];
                                            sigma1 = sigma1 + hamming[b ^ (inputarray[x]^inputarray[x^a])];
                                    }
                                    z = K - 2 * sigma1;
                                    if (sigma1 > treshold)
                                            break;
                                    temp = treshold - (sigma1/C);

                                    if (z < 0)
                                            z = z * (-1);
                                    sigma2 += z;
                            }
                            //tmp = 2 * hamming_weight(b);
                            tmp = 2 * hamming[b];
                            if (res < Math.abs((int)N - tmp) - sigma2/C)
                                    res = Math.abs((int)N - tmp) - sigma2/C;
                    }
            }
            return res;
    }    
    
    
}

