/***************************************************
 *  This is a hybrid algorithm of Strassen and Coppersmith-Winograd 
 *  This applies Winograd threshold to reduce recursion overwhelming at small inputs
 *  Could be reduced only to 4 multiplications but i dont know why it gets a wrong output 
 *  Thus i gave up on it and accepted the 4s time at hard case to solve the Paskets in O(N)
 *  
 *  This code belongs to Yousef Elkammar and only available for AA&D staff to view.
 *
 **************************************************/

using System;
using System.Threading.Tasks;

namespace Problem
{
    public static class MatrixMultiplication
    {
        public const int WINOGRAD_THRESHOLD = 128;

        static public int[,] MatrixMultiply(int[,] mat1, int[,] mat2, int size)
        {
            // Apply thresholding - faster than Strassen in small matrices
            if (size <= WINOGRAD_THRESHOLD) return StandardMatrixMultiply(mat1, mat2, size);
            
            /*                                    DIVIDE           */
            int h = size / 2;

            //Strassen's parameters.
            //dividing into submatrices 
            int[,] A11 = GetSubmatrix(mat1, 0, 0, h);
            int[,] A12 = GetSubmatrix(mat1, 0, h, h);
            int[,] A21 = GetSubmatrix(mat1, h, 0, h);
            int[,] A22 = GetSubmatrix(mat1, h, h, h);

            int[,] B11 = GetSubmatrix(mat2, 0, 0, h);
            int[,] B12 = GetSubmatrix(mat2, 0, h, h);
            int[,] B21 = GetSubmatrix(mat2, h, 0, h);
            int[,] B22 = GetSubmatrix(mat2, h, h, h);

            // Calculate intermediate matrices
            int[,] S1 = Subtract(B12, B22, h);
            int[,] S2 = Add(A11, A12, h);
            int[,] S3 = Add(A21, A22, h);
            int[,] S4 = Subtract(B21, B11, h);
            int[,] S5 = Add(A11, A22, h);
            int[,] S6 = Add(B11, B22, h);
            int[,] S7 = Subtract(A12, A22, h);
            int[,] S8 = Add(B21, B22, h);
            int[,] S9 = Subtract(A11, A21, h);
            int[,] S10 = Add(B11, B12, h);

            /*                                       CONQUER                     */

            // Getting submatrix using parallel processing
            var tasks = new[] {
                Task.Run(() => MatrixMultiply(A11, S1, h)),
                Task.Run(() => MatrixMultiply(S2, B22, h)),
                Task.Run(() => MatrixMultiply(S3, B11, h)),
                Task.Run(() => MatrixMultiply(A22, S4, h)),
                Task.Run(() => MatrixMultiply(S5, S6, h)),
                Task.Run(() => MatrixMultiply(S7, S8, h)),
                Task.Run(() => MatrixMultiply(S9, S10, h))
            };

            Task.WhenAll(tasks);

            // Get the results of the tasks
            int[,] product1 = tasks[0].Result;
            int[,] product2 = tasks[1].Result;
            int[,] product3 = tasks[2].Result;
            int[,] product4 = tasks[3].Result;
            int[,] product5 = tasks[4].Result;
            int[,] product6 = tasks[5].Result;
            int[,] product7 = tasks[6].Result;

            /*                                     COMBINE                     */
            // Generate output matrix
            int[,] comb11 = Add(Subtract(Add(product5, product4, h), product2, h), product6, h);
            int[,] comb12 = Add(product1, product2, h);
            int[,] comb13 = Add(product3, product4, h);
            int[,] comb14 = Subtract(Subtract(Add(product5, product1, h), product3, h), product7, h);

            int[,] comb = new int[size, size];

            //loop to combine the outout
            Parallel.For(0, h, i =>
            {
                for (int j = 0; j < h; j++)
                {
                    comb[i, j] = comb11[i, j];
                    comb[i, j + h] = comb12[i, j];
                    comb[i + h, j] = comb13[i, j];
                    comb[i + h, j + h] = comb14[i, j];
                }
            });
            return comb;
        }

        // Normal O(N^3), Faster in small matrices based on Winograd hypothesis.
        static unsafe private int[,] StandardMatrixMultiply(int[,] M1, int[,] M2, int N)
        {
            int[,] reslt = new int[N, N];

            fixed (int* pM1 = M1, pM2 = M2, pC = reslt)
            {
                int* pRowM1 = pM1, pRowM2 = pM2, pRowC = pC;

                for (int i = 0; i < N; i++, pRowM1 += N, pRowC += N)
                {
                    for (int j = 0; j < N; j++, pRowM2++)
                    {
                        int sum = 0;
                        int k = 0;

                        for (; k < N - 1; k += 2)
                        {
                            sum += pRowM1[k] * pRowM2[k * N] + pRowM1[k + 1] * pRowM2[(k + 1) * N];
                        }

                        for (; k < N; k++)
                        {
                            sum += pRowM1[k] * pRowM2[k * N];
                        }

                        pRowC[j] = sum;
                    }

                    pRowM2 -= N;
                }
            }

            return reslt;
        }

        /*                   Helper Functions                                                  */

        // calculating a smaller matrix for divide step
        unsafe static private int[,] GetSubmatrix(int[,] M, int row, int col, int size)
        {
            int[,] submatrix = new int[size, size];

            fixed (int* pM = M, pSub = submatrix)
            {
                int* pRowM = pM + row * M.GetLength(1) + col;
                int* pRowSub = pSub;

                for (int i = 0; i < size; i++, pRowM += M.GetLength(1), pRowSub += size)
                {
                    for (int j = 0; j < size; j++)
                    {
                        pRowSub[j] = pRowM[j];
                    }
                }
            }

            return submatrix;
        }

        // Add 2 matrices, still faster than multiplication based on Strassens
        static unsafe private int[,] Add(int[,] M1, int[,] M2, int N)
        {
            int[,] result = new int[N, N];

            fixed (int* pM1 = M1, pM2 = M2, pRes = result)
            {
                int* pRowM1 = pM1, pRowM2 = pM2, pRowRes = pRes;

                for (int i = 0; i < N; i++, pRowM1 += N, pRowM2 += N, pRowRes += N)
                {
                    for (int j = 0; j < N; j++)
                    {
                        pRowRes[j] = pRowM1[j] + pRowM2[j];
                    }
                }
            }

            return result;
        }

        // Subtract 2 matrices, still faster than multiplication based on Strassens
        static unsafe private int[,] Subtract(int[,] M1, int[,] M2, int N)
        {
            int[,] result = new int[N, N];

            fixed (int* pM1 = M1, pM2 = M2, pRes = result)
            {
                int* pRowM1 = pM1, pRowM2 = pM2, pRowRes = pRes;

                for (int i = 0; i < N; i++, pRowM1 += N, pRowM2 += N, pRowRes += N)
                {
                    for (int j = 0; j < N; j++)
                    {
                        pRowRes[j] = pRowM1[j] - pRowM2[j];
                    }
                }
            }

            return result;
        }

    }
}