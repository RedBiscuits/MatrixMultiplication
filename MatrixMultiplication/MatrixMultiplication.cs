using System;
using System.Threading.Tasks;

namespace Problem
{
    public static class MatrixMultiplication
    {
        static public int[,] MatrixMultiply(int[,] M1, int[,] M2, int N)
        {
            // Apply thresholding - faster than Strassen in small matrices
            if (N <= 128)
            {
                return StandardMatrixMultiply(M1, M2, N);
            }
            /*                                    DIVIDE           */
            // Partitioning to submatices 
            int[,] A11 = GetSubmatrix(M1, 0, 0, N / 2);
            int[,] A12 = GetSubmatrix(M1, 0, N / 2, N / 2);
            int[,] A21 = GetSubmatrix(M1, N / 2, 0, N / 2);
            int[,] A22 = GetSubmatrix(M1, N / 2, N / 2, N / 2);

            int[,] B11 = GetSubmatrix(M2, 0, 0, N / 2);
            int[,] B12 = GetSubmatrix(M2, 0, N / 2, N / 2);
            int[,] B21 = GetSubmatrix(M2, N / 2, 0, N / 2);
            int[,] B22 = GetSubmatrix(M2, N / 2, N / 2, N / 2);

            // Calculate intermediate matrices
            int[,] S1 = Subtract(B12, B22, N / 2);
            int[,] S2 = Add(A11, A12, N / 2);
            int[,] S3 = Add(A21, A22, N / 2);
            int[,] S4 = Subtract(B21, B11, N / 2);
            int[,] S5 = Add(A11, A22, N / 2);
            int[,] S6 = Add(B11, B22, N / 2);
            int[,] S7 = Subtract(A12, A22, N / 2);
            int[,] S8 = Add(B21, B22, N / 2);
            int[,] S9 = Subtract(A11, A21, N / 2);
            int[,] S10 = Add(B11, B12, N / 2);

            /*                                       CONQUER                     */

            // Getting submatrix using parallel processing
            var tasks = new[] {
                Task.Run(() => MatrixMultiply(A11, S1, N / 2)),
                Task.Run(() => MatrixMultiply(S2, B22, N / 2)),
                Task.Run(() => MatrixMultiply(S3, B11, N / 2)),
                Task.Run(() => MatrixMultiply(A22, S4, N / 2)),
                Task.Run(() => MatrixMultiply(S5, S6, N / 2)),
                Task.Run(() => MatrixMultiply(S7, S8, N / 2)),
                Task.Run(() => MatrixMultiply(S9, S10, N / 2))
            };

            Task.WhenAll(tasks);

            // Get the results of the tasks
            int[,] P1 = tasks[0].Result;
            int[,] P2 = tasks[1].Result;
            int[,] P3 = tasks[2].Result;
            int[,] P4 = tasks[3].Result;
            int[,] P5 = tasks[4].Result;
            int[,] P6 = tasks[5].Result;
            int[,] P7 = tasks[6].Result;

            /*                                     COMBINE                     */

            // Generate output matrix
            int[,] C11 = Add(Subtract(Add(P5, P4, N / 2), P2, N / 2), P6, N / 2);
            int[,] C12 = Add(P1, P2, N / 2);
            int[,] C21 = Add(P3, P4, N / 2);
            int[,] C22 = Subtract(Subtract(Add(P5, P1, N / 2), P3, N / 2), P7, N / 2);

            // Combine output matrices
            int[,] C = new int[N, N];
            SetSubmatrix(C, 0, 0, C11);
            SetSubmatrix(C, 0, N / 2, C12);
            SetSubmatrix(C, N / 2, 0, C21);
            SetSubmatrix(C, N / 2, N / 2, C22);

            return C;
        }


        // Normal O(N^3), Faster in small matrices
        static unsafe private int[,] StandardMatrixMultiply(int[,] M1, int[,] M2, int N)
        {
            int[,] C = new int[N, N];

            fixed (int* pM1 = M1, pM2 = M2, pC = C)
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

            return C;
        }

        /*                   Helper Functions                    */

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

        // Calculate matrix for combine step
        static unsafe private void SetSubmatrix(int[,] M, int row, int col, int[,] submatrix)
        {
            int size = submatrix.GetLength(0);

            fixed (int* pM = M, pSub = submatrix)
            {
                int* pRowM = pM + row * M.GetLength(1) + col;
                int* pRowSub = pSub;

                for (int i = 0; i < size; i++, pRowM += M.GetLength(1), pRowSub += size)
                {
                    for (int j = 0; j < size; j++)
                    {
                        pRowM[j] = pRowSub[j];
                    }
                }
            }
        }

        // Add 2 matrices, still faster than multiplication
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

        // Subtract 2 matrices, still faster than multiplication
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