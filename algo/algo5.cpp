
#include "math.h"
#include "stdlib.h"
#include <vector>
#include <iostream>

void free_matrix(double **res, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(res[i]);
    }
    free(res);
}

double **allocate(int n, int m)
{
    double **res = (double**)calloc(n, sizeof(double*));
    if (!res)
        return NULL;
    for (int i = 0; i < n; i++)
    {
        res[i] = (double*)calloc(m, sizeof(double));
        if (!res)
        {
            free_matrix(res, n);
            return NULL;
        }
    }
    return res;
}


int not_zero_to_top(double **matr, double *vect, int col, int n)
{
    int not_zero_number = -1;
    for (int i = col; i < n; i++)
    {
        if (fabs(matr[i][col]) > 0)
        {
            not_zero_number = i;
            break;
        }
    }

    if (not_zero_number == -1)
    {

        return 0;
    }


    double *temp = matr[col];;
    matr[col] = matr[not_zero_number];
    matr[not_zero_number] = temp;

    double t = vect[col];
    vect[col] = vect[not_zero_number];
    vect[not_zero_number] = t;
    return 1;
}

void subtraction_mult_rows(double *row1, double *row2, double &num1, double num2, int n, double mult)
{
    for (int i = 0; i < n; i++)
    {
        row1[i] -= mult * row2[i];
    }
    num1 -= num2 * mult;
}
void devide_row(double *row, double &num, int col, int n)
{
    double div = row[col];
    if (div == 0) return;
    for (int i = 0; i < n; i++)
    {
        row[i] /= div;
    }
    num /= div;
}

int check_result(double **matr, double *vect, int n)
{
    double diap = 0.0001;
    for (int i = 0; i < n; i++)
    {
        if (fabs(matr[i][i]) < diap && fabs(vect[i]) > diap)
            return -1;
    }
    for (int i = 0; i < n; i++)
    {
        if (fabs(matr[i][i]) < diap && fabs(vect[i]) < diap)
            return -2;
    }
    return 0;
}

int gauss_method(double **matr, double *vect, int n)
{
    int err = 0;
    for (int i = 0; i < n; i++)
    {
        if (not_zero_to_top(matr, vect, i, n))
        {
            devide_row(matr[i], vect[i], i, n);
            for (int j = 0; j < n; j++)
            {

                if (j != i)
                {
                    subtraction_mult_rows(matr[j], matr[i], vect[j], vect[i], n, matr[j][i]);
                }
            }
        }
    }

    err = check_result(matr, vect, n);

    return err;
}


double f(double g, double T, double V, std::vector<double> x_vect, std::vector<int> z_vect)
{
    double sum = 0;
    for (int i = 2; i <= 5; i++)
    {
        sum += exp(x_vect[i - 1]) * pow(z_vect[i - 1], 2) / (1 + pow(z_vect[i - 1], 2)*(g / 2));
    }
    //std::cout << "sum " << sum << "\n";
    double y = pow(g, 2) - 5.87 * pow(10, 10) * (1 / pow(T, 3)) * (exp(V) / (1 + g / 2) + sum);
    //std::cout << "y " << y << "\n";
    return y;
}

void  get_d_E(std::vector<double>& d_E, double T, std::vector<int> z_vect, double g)
{
    d_E.clear();
    double zn;
    for (int i = 0; i < 4; i++)
    {
        zn = 8.61 * pow(10, -5) * T * log((1 + pow(z_vect[i + 1], 2) * (g / 2)) * (1 + g / 2) / (1 + pow(z_vect[i], 2)*(g / 2)));
        //std::cout << "zn "<< log((1 + pow(z_vect[i + 1], 2) * (g / 2)) * (1 + g / 2) / (1 + pow(z_vect[i], 2)*(g / 2))) << "\n";
        d_E.push_back(zn);
    }
}



std::vector<double> Q(5, 0);

void get_Q(double t)
{
    if (t <= 4000)
    {
        Q = {1, 4.05, 5.15, 6, 7};
    }
    else if (t <= 8000)
    {
        Q = {1, 4.3, 5.98, 6, 7};
    }
    else if (t <= 10000)
    {
        Q = {1.0025, 4.44, 6.47, 7, 8};
    }
    else if (t <= 12000)
    {
        Q = {1.020, 4.57, 6.96, 7.5, 8 };
    }
    else if (t <=14000)
    {
        Q = {1.0895, 4.65, 7.41, 8, 9};
    }
}

void  get_K(std::vector<double>& K, std::vector<double> E, std::vector<double> d_E, double T, std::vector<int> z_vect, double g)
{
    K.clear();
    double zn;

    get_Q(T);

    for (int i = 0; i < 4; i++)
    {
        //std::cout << "fg  " << get_Q(T, i + 1) << " "<<  get_Q(T, i) << "\n";
        //std::cout << -(E[i] - d_E[i]) * 11606 / T << "\n";
        //std::cout << "QQ "<< get_Q(T, i + 1) << " " << get_Q(T, i) << "\n";
        zn = 2 * 2.415 * pow(10,-3)* Q[i + 1] / Q[i] * pow(T, 1.5) * exp(-(E[i] - d_E[i]) * (11606 / T));
        //std::cout << "znnn  "<< zn << "\n";
        K.push_back(zn);
    }
}

double get_x_pol(double a, double b)
{
    return (a + b) / 2;
}
double polovDel(double a, double b, double eps, double T, double V, std::vector<double> x_vect, std::vector<int> z_vect)
{
    int max_iter = 1000;
    int k = 0;
    double c = 0;
    while (fabs(b - a) > eps)
    {
        if (k > max_iter)
            break;
        c = get_x_pol(a, b);
        //std::cout << "res =  " << c << "\n";
        //std::cout << "a b" << a << " " << b << "\n";
        if (f(b, T, V, x_vect, z_vect) * f(c, T, V, x_vect, z_vect) <= 0)
            a = c;
        else if (f(a, T, V, x_vect, z_vect) * f(c, T, V, x_vect, z_vect) < 0)
            b = c;
        k++;
    }
    return get_x_pol(a, b);
}



int check_epsilon(std::vector<double> err, double eps)
{
    double max = fabs(err[0]);
    double abs;

    for (size_t i = 1; i < err.size(); ++i)
    {
        abs = fabs(err[i]);
        if (max < abs)
            max = abs;
    }

    if (max > eps)
        return 1;


    return 0;
}

void create_matrix(double ***matr, double **vect, double V, std::vector<double> x_vect, std::vector<double> K, std::vector<int> z_vect, double P, double T, double g)
{
    int n = 6;
    *matr = allocate(n, n);
    *vect = (double*)calloc(n, sizeof(double));
    (*matr)[0][0] = 1; (*matr)[0][1] = -1; (*matr)[0][2] = 1; (*matr)[0][3] = 0; (*matr)[0][4] = 0; (*matr)[0][5] = 0;
    (*matr)[1][0] = 1; (*matr)[1][1] = 0; (*matr)[1][2] = -1; (*matr)[1][3] = 1; (*matr)[1][4] = 0; (*matr)[1][5] = 0;
    (*matr)[2][0] = 1; (*matr)[2][1] = 0; (*matr)[2][2] = 0; (*matr)[2][3] = -1; (*matr)[2][4] = 1; (*matr)[2][5] = 0;
    (*matr)[3][0] = 1; (*matr)[3][1] = 0; (*matr)[3][2] = 0; (*matr)[3][3] = 0; (*matr)[3][4] = -1; (*matr)[3][5] = 1;
    (*matr)[4][0] = exp(V); (*matr)[4][1] = -z_vect[0]*exp(x_vect[0]);  (*matr)[4][2] = -z_vect[1] * exp(x_vect[1]);  (*matr)[4][3] = -z_vect[2] * exp(x_vect[2]);  (*matr)[4][4] = -z_vect[3] * exp(x_vect[3]);  (*matr)[4][5] = -z_vect[4] * exp(x_vect[4]);
    (*matr)[5][0] = exp(V); (*matr)[5][1] = exp(x_vect[0]);  (*matr)[5][2] = exp(x_vect[1]);  (*matr)[5][3] = exp(x_vect[2]);  (*matr)[5][4] = exp(x_vect[3]);  (*matr)[5][5] = exp(x_vect[4]);
    (*vect)[0] = -(V + x_vect[1] - x_vect[0] - log(K[0]));
    //std::cout << "log(K_i(i)) " << log(K[0]) << "\n";


    (*vect)[1] = -(V + x_vect[2] - x_vect[1] - log(K[1]));
    (*vect)[2] = -(V + x_vect[3] - x_vect[2] - log(K[2]));
    (*vect)[3] = -(V + x_vect[4] - x_vect[3] - log(K[3]));
    (*vect)[4] = -(exp(V) - (z_vect[1] * exp(x_vect[1]) + z_vect[2] * exp(x_vect[2]) + z_vect[3] * exp(x_vect[3]) + z_vect[4] * exp(x_vect[4])));
    (*vect)[5] = -(exp(V) + exp(x_vect[0]) + exp(x_vect[1]) + exp(x_vect[2])  + exp(x_vect[3])  + exp(x_vect[4]) - 0.285 * 1e-11 * pow(g*T,3) - P * 7.242 * 1e3 / T);
}

int main()
{
    double V = 0.02;
    double d_V = 0.02;
    std::vector<double> x_vect{-0.5, 0.02, -7.5,  -12, -15};
    std::vector<double> d_x_vect{ -0.5, 0.02, -7.5,  -12, -15};
    std::vector<int> z_vect{ 0,1,2,3,4 };
    std::vector <double> err(6, 0);
    double g;
    double T;
    double P;
    std::vector<double> E{ 12.08, 20.98, 31.00, 45.00 };
    std::vector<double> d_E;
    std::vector<double> K;
    double eps = 1e-12; // Dichotomy
    double eps2 = 1e-5; // check eps
    //T = 8000;
    //P = 5;
    double **matr;
    double *vect;
    int k = 0;
    std::cout << "Enter T: ";
    std::cin >> T;
    std::cout << "Enter P: ";
    std::cin >> P;
    //std::cout << f(10, T, V, x_vect, z_vect) << "\n";
    do
    {
        g = polovDel(0, 3, eps, T, V, x_vect, z_vect);
        //std::cout << "g " << g << "\n";
        get_d_E(d_E, T, z_vect, g);
        //std::cout << "d_E " << d_E[0] << " " << d_E[1] << " " << d_E[2] << " " << d_E[3] << "\n";
        get_K(K, E, d_E, T, z_vect, g);
        //std::cout << "K "<< K[0] << " " << K[1] << " " << K[2] << " " << K[3] << "\n";

        create_matrix(&matr, &vect, V, x_vect, K, z_vect, P, T, g);

        int n = 6;
        /*for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                std::cout << matr[i][j] << " ";
            }
            std::cout << "\n";
        }
        printf("\n");*/
        /*for (int i = 0; i < n; i++)
        {
            std::cout << vect[i] << " ";
        }
        printf("\n");*/
        gauss_method(matr, vect, n);
        //std::cout << "right part ";
        /*for (int i = 0; i < n; i++)
        {
            printf("%lf ", vect[i]);
        }
        printf("\n");*/
        d_V = vect[0];
        err[0] = vect[0] / V;
        V += d_V;
        for (int i = 0; i < 5; i++)
        {
            d_x_vect[i] = vect[i + 1];
            err[i + 1] = vect[i + 1] / x_vect[i];
            x_vect[i] += d_x_vect[i];
        }
        //std::cout <<"v " << V << "\n";
        //std::cout <<"x_vect " << x_vect[0] << " " << x_vect[1] << " " << x_vect[2] << " " << x_vect[3] << " " << x_vect[4] << "\n";


    } while (check_epsilon(err, eps2));


    std::cout << "Ne = " << exp(V) << std::endl;

    for (size_t i = 0; i < 5; ++i)
    {
        std::cout << "N" << i + 1 << " = " <<  exp(x_vect[i]) << std::endl;
    }

    return 0;

}

