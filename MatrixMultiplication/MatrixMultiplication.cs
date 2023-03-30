using System;
using System.Threading.Tasks;

namespace Problem
{
    public static class MatrixMultiplication
    {
        static public int[,] MatrixMultiply(int[,] M1, int[,] M2, int N)
        {
            // Apply thresholding
            if (N <= 32)
            {
                return StandardMatrixMultiply(M1, M2, N);
            }

            // Block partitioning
            int[,] A11 = GetSubmatrix(M1, 0, 0, N / 2);
            int[,] A12 = GetSubmatrix(M1, 0, N / 2, N / 2);
            int[,] A21 = GetSubmatrix(M1, N / 2, 0, N / 2);
            int[,] A22 = GetSubmatrix(M1, N / 2, N / 2, N / 2);

            int[,] B11 = GetSubmatrix(M2, 0, 0, N / 2);
            int[,] B12 = GetSubmatrix(M2, 0, N / 2, N / 2);
            int[,] B21 = GetSubmatrix(M2, N / 2, 0, N / 2);
            int[,] B22 = GetSubmatrix(M2, N / 2, N / 2, N / 2);

            // Compute intermediate matrices
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

            // Compute submatrix products using tasks
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

            // Compute output matrix
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

        static private int[,] StandardMatrixMultiply(int[,] M1, int[,] M2, int N)
        {
            int[,] C = new int[N, N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    int sum = 0;
                    int k = 0;

                    for (; k < N - 1; k += 2)
                    {
                        sum += M1[i, k] * M2[k, j] + M1[i, k + 1] * M2[k + 1, j];
                    }

                    for (; k < N; k++)
                    {
                        sum += M1[i, k] * M2[k, j];
                    }

                    C[i, j] = sum;
                }
            }

            return C;
        }


        static private int[,] GetSubmatrix(int[,] M, int row, int col, int size)
        {
            int[,] submatrix = new int[size, size];

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    submatrix[i, j] = M[row + i, col + j];
                }
            }

        return submatrix;
        }

        static private void SetSubmatrix(int[,] M, int row, int col, int[,] submatrix)
        {
            int size = submatrix.GetLength(0);

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    M[row + i, col + j] = submatrix[i, j];
                }
            }
        }

        static private int[,] Add(int[,] M1, int[,] M2, int N)
        {
            int[,] result = new int[N, N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    result[i, j] = M1[i, j] + M2[i, j];
                }
            }

            return result;
        }

        static private int[,] Subtract(int[,] M1, int[,] M2, int N)
        {
            int[,] result = new int[N, N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    result[i, j] = M1[i, j] - M2[i, j];
                }
            }

            return result;
        }

    }
}
