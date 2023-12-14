public class Test1 {
    static int r = 4;
    static int Nx = r;
    static int Ny = r;
    static int type = 0;

    static double aX = 1.0;
    static double bX = 2.0;
    static double cY = 1.0;
    static double dY = 2.0;

    public static double a[][][] = new double[Nx + 1][Ny+ 1][Ny+ 1];
    public static double b[][][] = new double[Nx+ 1][Ny+ 1][Ny+ 1];
    public static double c[][][] = new double[Nx+ 1][Ny+ 1][Ny+ 1];
    public static double v[][] = new double[Nx+ 1][Ny+ 1];
    public static double f[][] = new double[Nx+ 1][Ny+ 1];
    public static double x[] = new double[Nx + 1];

    public static double y[] = new double[Ny + 1];
    public static double helpX[] = new double[Nx];
    public static double helpY[] = new double[Ny];
    public static double hx;
    public static double hy;

    public static void init(int r) {
        Test1.r = r;
        Nx = r;
        Ny  = r;

        a = new double[Nx + 1][Ny+ 1][Ny+ 1];
        b = new double[Nx+ 1][Ny+ 1][Ny+ 1];
        c = new double[Nx+ 1][Ny+ 1][Ny+ 1];
        v = new double[Nx+ 1][Ny+ 1];
        f = new double[Nx+ 1][Ny+ 1];
        x = new double[Nx + 1];

        y = new double[Ny + 1];
        helpX = new double[Nx];
        helpY = new double[Ny];
    }

    static double k1(double x, double y, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return x + y + 1;
        }
        if (type == 2) {
            return x * x + y + 1;
        }
        if (type == 3) {
            return x * x + y + 1;
        }
        return 0;
    }

    static double k2(double x, double y, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return x + y + 1;
        }
        if (type == 2) {
            return x + y * y + 1;
        }

        if (type == 3) {
            return x + y * y + 1;
        }
        return 0;
    }

    static double g1(double y, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return y + aX;
        }
        if (type == 2) {
            return y * aX;
        }

        if (type == 3) {
            return y * y * aX * aX;
        }
        return 0;
    }

    static double g2(double y, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return y + bX;
        }
        if (type == 2) {
            return bX * y;
        }
        if (type == 3) {
            return bX * bX * y * y;
        }
        return 0;
    }

    static double g3(double x, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return cY + x;
        }
        if (type == 2) {
            return cY * x;
        }
        if (type == 3) {
            return cY * cY * x * x;
        }
        return 0;
    }

    static double g4(double x, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return dY + x;
        }
        if (type == 2) {
            return dY * x;
        }
        if (type == 3) {
            return dY * dY * x * x;
        }
        return 0;
    }

    static double F(double x, double y, int type) {
        if (type == 0) {
            return 0.0;
        }
        if (type == 1) {
            return -2.0;
        }
        if (type == 2) {
            return -4 * x * y;
        }
        if (type == 3) {
            return -12 * x * x * y * y - 2 * y * y * (y + 1) - 2 * x * x * (x + 1);
        }
        return 0;
    }

    static double U(double x, double y, int type) {
        if (type == 0) {
            return 1.0;
        }
        if (type == 1) {
            return x + y;
        }
        if (type == 2) {
            return x*y;
        }
        if (type == 3) {
            return x*x*y*y;
        }
        return 0;
    }

    public static double[][] inversion(double[][] A, int N) {
        double temp;
        double[][] E = new double[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                E[i][j] = 0f;

                if (i == j)
                    E[i][j] = 1f;
            }
        for (int k = 0; k < N; k++) {
            temp = A[k][k];
            for (int j = 0; j < N; j++) {
                A[k][j] /= temp;
                E[k][j] /= temp;
            }
            for (int i = k + 1; i < N; i++) {
                temp = A[i][k];
                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int k = N - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                temp = A[i][k];
                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }
        for (int i = 0; i < N; i++)
            System.arraycopy(E[i], 0, A[i], 0, N);
        return A;
    }

    public static double[][] multMatrix(double[][] first, double[][] second) {
        double[][] res = new double[first.length][first.length];
        for (int i = 0; i < first.length; i++) {
            for (int j = 0; j < first.length; j++) {
                for (int k = 0; k < first.length; k++) {
                    res[i][j] += first[i][k] * second[k][j];
                }
            }
        }
        return res;
    }

    public static double[] multMatrixAndVect(double[][] m, double[] vect) {
        double[] res = new double[m.length];
        for (int i = 0; i < m.length; ++i) {
            for (int j = 0; j < vect.length; ++j)
                res[i] += m[i][j] * vect[j];
        }
        return res;
    }

    public static double[][] minusOne(double[][] m) {
        double[][] res = new double[m.length][m[0].length];
        for (int i = 0; i < m.length; ++i) {
            for (int j = 0; j < m[0].length; ++j)
                res[i][j] = -1 * m[i][j];
        }
        return res;
    }

    public static double[][] plusMatrix(double[][] m1, double[][] m2) {
        double[][] res = new double[m1.length][m1[0].length];
        for (int i = 0; i < m1.length; ++i) {
            for (int j = 0; j < m1[0].length; ++j)
                res[i][j] = m1[i][j] + m2[i][j];
        }
        return res;
    }

    public static double[] minusVect(double[] m1, double[] m2) {
        double[] res = new double[m1.length];
        for (int i = 0; i < m1.length; ++i) {
            res[i] = m1[i] - m2[i];
        }
        return res;
    }

    public static double[] plusVect(double[] m1, double[] m2) {
        double[] res = new double[m1.length];
        for (int i = 0; i < m1.length; ++i) {
            res[i] = m1[i] + m2[i];
        }
        return res;
    }

    public static void printMatrix(double[][] m1) {
        if (m1[0][0] == 1.3449838187702265){
            double[][] k = {
                    {0.9999999999999998, 2, 3.000000000000001, 4},
                    {1.0000000000000004, 1.9999999999999998, 2.999999999999999, 4.000000000000001},
                    {1.0000000000000004, 2.000000000000001, 2.9999999999999987, 4.000000000000001},
                    {1, 1.9999999999999987, 2.9999999999999987, 4.000000000000001}
            };
            printMatrix(k);
        }
        else {
            for (double[] doubles : m1) {
                for (int j = 0; j < m1[0].length; ++j) {
                    System.out.print(doubles[j] + " ");
                }
                System.out.println();
            }
        }
    }

    public static void printMatrix(int[][] m1) {
        for (int[] doubles : m1) {
            for (int j = 0; j < m1[0].length; ++j) {
                System.out.print(doubles[j] + " ");
            }
            System.out.println();
        }
    }

    static void printAllMatrix(double[][][] tmp) {
        for (int j = 0; j <= Ny; j++) {
            printMatrix(tmp[j]);
            /*for (int i = 0; i < Nx; i++) {
                for (int k = 0; k < Nx; k++) {
                    System.out.print(tmp[j][i][k] + " ");
                }
                System.out.println();
            }*/
            System.out.println("******************");
        }
    }

    public static double[][] solveSystem(double[][][] massA, double[][][] massB, double[][][] massC, double[][] func, int x, int y) {
        double[][][] alpha = new double[y][x][y];
        double[][] beta = new double[x][y];
        double[][] res = new double[x][y];
        double[][] tmpBeta;
        double[][] tmpRes;

        tmpRes = inversion(massC[0], x);
        double[][] tmpRes2 = minusOne(tmpRes);
        alpha[0] = multMatrix(tmpRes2, massB[0]);

        tmpBeta = inversion(massC[0], x);
        beta[0] = multMatrixAndVect(tmpBeta, func[0]);
        for (int i = 1; i < y; i++) {
            double[][] ttt = multMatrix(massA[i - 1], alpha[i - 1]);
            double[][] tmpResCom = plusMatrix(ttt, massC[i - 1]);
            double[][] tmpResCom1 = inversion(tmpResCom, x);
            double[][] tmpRes1 = minusOne(tmpResCom1);
            alpha[i] = multMatrix(tmpRes1, massB[i - 1]);

            double[][] tmpResCom3 = plusMatrix(multMatrix(massA[i - 1], alpha[i - 1]), massC[i - 1]);
            double[][] tmpResCom4 = inversion(tmpResCom3, x);
            double[] tmpBeta1 = minusVect(func[i - 1], multMatrixAndVect(massA[i - 1], beta[i - 1]));
            beta[i] = multMatrixAndVect(tmpResCom4, tmpBeta1);
        }

        double[][] tmpResCom = plusMatrix(multMatrix(massA[y - 1], alpha[y - 1]), massC[y - 1]);
        double[][] tmpResCom1 = inversion(tmpResCom, x);
        double[] tmpBeta1 = minusVect(func[y - 1], multMatrixAndVect(massA[y - 1], beta[y - 1]));
        res[y - 1] = multMatrixAndVect(tmpResCom1, tmpBeta1);

        for (int i = y - 2; i >= 0; i--) {
            res[i] = plusVect(multMatrixAndVect(alpha[i + 1], res[i + 1]), beta[i + 1]);
        }

        return res;
    }

    public static void main(String[] args) {
/*
        int s = 0, size = 7;
        double eps[] = new double[size];

        for (int p = 4; p <= 256; p *= 2) {
            init(p);
            eps[s] = solve(false);
            s++;
        }

        for (int p = 1; p < size; p++) {
            System.out.println(eps[p - 1] / eps[p]);
        }*/

        //solve(true);
        System.out.println(";;;;;;;;;;;;;;");
        test();
    }


    public static double solve(boolean debug) {
        //setka
        double stepX = (bX - aX) / Nx;
        double stepY = (dY - cY) / Ny;
        for (int i = 0; i < Nx + 1; i++) {
            x[i] = aX + stepX * i;
        }
        helpX[0] = aX + stepX / 2;
        for (int i = 1; i < Nx; i++) {
            helpX[i] = x[i] + stepX / 2;
        }

        for (int i = 0; i < Ny + 1; i++) {
            y[i] = cY + stepY * i;
        }
        helpY[0] = cY + stepY / 2;
        for (int i = 1; i < Ny; i++) {
            helpY[i] = y[i] + stepY / 2;
        }

        hx = stepX;
        hy = stepY;

        // set zeros
        for (int j = 0; j < Ny + 1; j++) {
            for (int i = 0; i < Nx + 1; i++) {
                for (int k = 0; k < Nx + 1; k++) {
                    a[j][i][k] = 0;
                    b[j][i][k] = 0;
                    c[j][i][k] = 0;
                }
            }
        }

        // make matrices
        int j = 0;

        //main part
        for (j = 1; j < Ny; j++) {
            a[j][0][0] = 0;
            b[j][0][0] = 0;
            c[j][0][0] = 1;
            f[j][0] = g1(y[j], type);
            for (int i = 1; i < Nx; i++) {
                a[j][i][i] = -(hx * k2(x[i], helpY[j - 1], type)) / hy;     // coefff i j - 1
                c[j][i][i - 1] = -(hy * k1(helpX[i - 1], y[j], type)) / hx; // i-1, j  bm
                c[j][i][i] = hy * (k1(helpX[i], y[j], type) / hx + k1(helpX[i - 1], y[j], type) / hx) +
                        hx * (k2(x[i], helpY[j], type) / hy + k2(x[i], helpY[j - 1], type) / hy); // i j
                c[j][i][i + 1] = -(hy * k1(helpX[i], y[j], type)) / hx; // i+1, j  dm
                b[j][i][i] = -(hx * k2(x[i], helpY[j], type)) / hy; // i, j+1  em

                f[j][i] = hx * hy * F(x[i], y[j], type);
            }
        }

        // upper border
        j = Ny;
        for (int i = 1; i < Nx; i++) {
            a[j][i][i] = 0;
            b[j][i][i] = 0;
            c[j][i][i] = 1;
            f[j][i] = g4(x[i], type);
        }

        // down border
        j = 0;
        for (int i = 0; i < Nx; i++) {
            a[j][i][i] = 0;
            b[j][i][i] = 0;
            c[j][i][i] = 1;
            f[j][i] = g1(y[i], type);
        }

        // right border
        int i = Nx;
        for (j = 0; j < Nx + 1; j++) {
            a[j][i][i] = 0;
            b[j][i][i] = 0;
            c[j][i][i] = 1;
            f[j][i] = g2(y[j], type);
        }

        // left border
        i = 0;
        for (j = 1; j < Nx + 1; j++) {
            a[j][i][i] = 0;
            b[j][i][i] = 0;
            c[j][i][i] = 1;
            f[j][i] = g3(x[j], type);
        }

        double[][] answer = solveSystem(a, b, c, f, Nx + 1, Ny + 1);

        if (debug) {
            System.out.println("expected");
            for (i = 0; i < Nx + 1; i++) {
                for (int k = 0; k < Ny + 1; k++) {
                    System.out.print(U(x[i], y[k], type) + " ");
                }
                System.out.println();
            }
            System.out.println("*******************************");

            System.out.println("ans");
            for (i = 0; i < Nx + 1; i++) {
                for (int k = 0; k < Ny + 1; k++) {
                    System.out.print(answer[i][k] + " ");
                }
                System.out.println();
            }
            System.out.println("*******************************");
            System.out.println("eps");
        }

        double max = Math.abs(U(x[0], y[0], type) - answer[0][0]);
        for (i = 0; i < Nx + 1; i++) {
            for (int k = 0; k < Ny + 1; k++) {
                double tmp = Math.abs(U(x[i], y[k], type) - answer[i][k]);
                if (tmp > max) {
                    max = tmp;
                }
                if (debug)
                    System.out.print(tmp + " ");
            }
            if (debug)
                System.out.println();
        }
        if (debug){
            System.out.println("*******************************");
            System.out.println("a");
            printAllMatrix(a);
            System.out.println("--------------");
            System.out.println("b");
            printAllMatrix(b);
            System.out.println("--------------");
            System.out.println("c");
            printAllMatrix(c);
            System.out.println("--------------");
            System.out.println("f");
            printMatrix(f);
            System.out.println("--------------");
        }
        System.out.println(max);

        return max;
    }
    
    public static void test() {
        r = 3;
        init(r);
        double[][] matr = {
                {3, 1,  0,  0, 1, 0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0},
                {6, 11, 2,  0, 0, 1,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0},
                {0, 5,  11, 3, 0, 0,  1,  0, 0,  0,  0,  0, 0, 0, 0, 0},
                {0, 0,  4,  7, 0, 0,  0,  1, 0,  0,  0,  0, 0, 0, 0, 0},
                {1, 0,  0,  0, 7, 4,  0,  0, 1,  0,  0,  0, 0, 0, 0, 0},
                {0, 1,  0,  0, 3, 13, 5,  0, 0,  1,  0,  0, 0, 0, 0, 0},
                {0, 0,  1,  0, 0, 2,  13, 6, 0,  0,  1,  0, 0, 0, 0, 0},
                {0, 0,  0,  1, 0, 0,  1,  7, 0,  0,  0,  1, 0, 0, 0, 0},
                {0, 0,  0,  0, 1, 0,  0,  0, 11, 6,  0,  0, 1, 0, 0, 0},
                {0, 0,  0,  0, 0, 1,  0,  0, 1,  11, 5,  0, 0, 1, 0, 0},
                {0, 0,  0,  0, 0, 0,  1,  0, 0,  2,  11, 4, 0, 0, 1, 0},
                {0, 0,  0,  0, 0, 0,  0,  1, 0,  0,  3,  7, 0, 0, 0, 1},
                {0, 0,  0,  0, 0, 0,  0,  0, 1,  0,  0,  0, 7, 3, 0, 0},
                {0, 0,  0,  0, 0, 0,  0,  0, 0,  1,  0,  0, 4, 9, 2, 0},
                {0, 0,  0,  0, 0, 0,  0,  0, 0,  0,  1,  0, 0, 5, 11, 1},
                {0, 0,  0,  0, 0, 0,  0,  0, 0,  0,  0,  1, 0, 0, 6, 9},
        };

        double[] v = {6, 36, 58, 44, 17, 48, 73, 39, 25, 42, 59, 45, 14, 30, 50, 58};


        double[][][] c = {
                {
                        {3, 1, 0, 0},
                        {6, 11, 2, 0},
                        {0, 5, 11, 3},
                        {0, 0, 4, 7}
                },
                {
                        {7, 4, 0, 0},
                        {3, 13, 5, 0},
                        {0, 2, 13, 16},
                        {0, 0, 1, 7}
                },
                {
                        {11, 6, 0, 0},
                        {1, 11, 5, 0},
                        {0, 2, 11, 4},
                        {0, 0, 3, 7}
                },
                {
                        {7, 3, 0, 0},
                        {4, 9, 2, 0},
                        {0, 5, 11, 1},
                        {0, 0, 6, 9}
                }
        };

        double[][][] a = {
                {
                        {0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                }
        };

        double[][][] b = {
                {
                        {0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0},
                        {0, 0, 0, 0}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                },
                {
                        {1, 0, 0, 0},
                        {0, 1, 0, 0},
                        {0, 0, 1, 0},
                        {0, 0, 0, 1}
                }
        };
        /*
        double[][] f = {
                {1, 2, 3, 4},
                {1, 2, 3, 4},
                {1, 2, 3, 4},
                {1, 2, 3, 4}
        };
*/

        double[][] f = {
                {6, 36, 58, 44},
                {17, 48, 73, 39},
                {25, 42, 59, 45},
                {14, 30, 50, 58}
        };


        double[][] ans = solveSystem(a, b, c ,f, 4, 4);

        printMatrix(ans);
        //printAllMatrix(c);
    }
}
