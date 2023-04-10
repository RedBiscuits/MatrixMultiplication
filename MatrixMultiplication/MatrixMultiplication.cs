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
            // Apply thresholding - faster than Strassen in small matrices based on Winograd hypothesis due to High Constant factor
            if (size <= WINOGRAD_THRESHOLD) return SMM(mat1, mat2, size);
            
            /*                                    DIVIDE                                  */
            int h = size / 2;

            //Strassen's parameters.
            //dividing into submatrices 
            // processing each mat on its own for better locality
            int[,] topLeft = DM(mat1, 0, 0, h);
            int[,] botLeft = DM(mat1, h, 0, h);

            int[,] topRight = DM(mat1, 0, h, h);
            int[,] botRight = DM(mat1, h, h, h);

            int[,] topLeft2 = DM(mat2, 0, 0, h);
            int[,] botLeft2 = DM(mat2, h, 0, h);
            
            int[,] topRight2 = DM(mat2, 0, h, h);
            int[,] botRight2 = DM(mat2, h, h, h);


            /*                                       CONQUER                     */
            var tasks = new[] {
                Task.Run(() => MatrixMultiply(topLeft, SM(topRight2, botRight2, h), h)),
                Task.Run(() => MatrixMultiply(AM(topLeft, topRight, h), botRight2, h)),
                Task.Run(() => MatrixMultiply(AM(botLeft, botRight, h), topLeft2, h)),
                Task.Run(() => MatrixMultiply(botRight, SM(botLeft2, topLeft2, h), h)),

                Task.Run(() => MatrixMultiply(AM(topLeft, botRight, h), AM(topLeft2, botRight2, h), h)),
                Task.Run(() => MatrixMultiply(SM(topRight, botRight, h), AM(botLeft2, botRight2, h), h)),
                Task.Run(() => MatrixMultiply(SM(topLeft, botLeft, h),  AM(topLeft2, topRight2, h), h))
            };

            Task.WhenAll(tasks);


            /*                                     COMBINE                     */
            int[,] resTopLeft = AM(SM(AM(tasks[4].Result, tasks[3].Result, h), tasks[1].Result, h), tasks[5].Result, h);
            int[,] resBotLeft = AM(tasks[2].Result, tasks[3].Result, h);

            int[,] resBotRight = SM(SM(AM(tasks[4].Result, tasks[0].Result, h), tasks[2].Result, h), tasks[6].Result, h);           
            int[,] resTopRight = AM(tasks[0].Result, tasks[1].Result, h);


            int[,] finalMat = new int[size, size];

            //loop to combine the outout
            Parallel.For(0, h, i =>
            {
                for (int j = 0; j < h; j++)
                {
                    finalMat[i, j] = resTopLeft[i, j];
                    finalMat[i, j + h] = resTopRight[i, j];
                    finalMat[i + h, j] = resBotLeft[i, j];
                    finalMat[i + h, j + h] = resBotRight[i, j];
                }
            });
            return finalMat;
        }


        /*                   Helpers                                                  */

        // Normal O(N^3), Faster in small matrices based on Winograd hypothesis.
        static unsafe private int[,] SMM(int[,] M1, int[,] M2, int N)
        {
            int[,] reslt = new int[N, N];

            //unsafe pointers
            fixed (int* pM1 = M1, pM2 = M2, pC = reslt)
            {

                int* pRowM1 = pM1, pRowM2 = pM2, pRowC = pC;

                for (int i = 0; i < N; i++, pRowM1 += N, pRowC += N)
                {
                    for (int j = 0; j < N; j++, pRowM2++)
                    {

                        //Unrolled loops, for better cache localization
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


        // calculating a smaller matrix for divide step
        unsafe static private int[,] DM(int[,] M, int row, int col, int size)
        {
            //res
            int[,] submatrix = new int[size, size];

            //unsafe divide with pointers
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

        // Add 2 matrices, O(N^2)
        static unsafe private int[,] AM(int[,] M1, int[,] M2, int N)
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
        static unsafe private int[,] SM(int[,] M1, int[,] M2, int N)
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