
// algo6.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <vector>

double get_polinom_leg(int n, double t)
{
	//std::vector<double> koef;

	double p0 = 1;
	double p1 = t;
	double p_tmp;
	if (n == 0)
		return p0;
	if (n == 1)
		return p1;
	int n_now = 1;
	do
	{
		//koef.clear();
		n_now++;
		p_tmp = p1;
		p1 = (1.0 / n_now) * ((2 * n_now - 1) * t * p_tmp - (n_now - 1) * p0);
		//koef.push_back(n_now);
		//koef.push_back(p_tmp);
		//koef.push_back(p0);
		//std::cout << p1 << std::endl;
		p0 = p_tmp;	
	} while (n_now != n);
	return p1;
}
/*
double put_in_polynom_leg(int n, double t, std::vector<double> koef)
{

	if (n == 0)
		return 0;
	if (n == 1)
		return t;

	return (1.0 / koef[0]) * ((2 * koef[0] - 1) * t * koef[1] - (koef[0] - 1) * koef[2]);
}*/

double get_x_pol(double a, double b)
{
	return (a + b) / 2;
}

void finds_roots_lagrange(std::vector<double>& roots_leg, int n)
{
	//std::vector<double> koef = get_polinom_leg(n, a);
	double a = -1;
	double b_end = 1;
	int max_iter = 10000;
	int k = 0;
	double c = 0;
	double eps = 1e-6;
	double delta = 2. / n;

	double b;
	while (roots_leg.size() < n)
	{
		roots_leg.clear();
		a = -1;
		b = a;
		while (b < b_end)
		{
			a = b;
			b += delta;
			if (get_polinom_leg(n, a) * get_polinom_leg(n, b) > 0)
			{
				//std::cout << "b a" << b << " " << a << "\n";
				continue;
			}
			
			k = 0;
			double t;
			while (fabs(b - a) > eps)
			{
				if (k > max_iter)
					break;
				c = get_x_pol(a, b);
				//std::cout << "res =  " << c << "\n";
				//std::cout << "b a" << b << " " << a << "\n";
				t = get_polinom_leg(n, c);
				if (get_polinom_leg(n, b) * t <= 0)
					a = c;
				else if (get_polinom_leg(n, a) * t < 0)
					b = c;
				k++;
			}
			roots_leg.push_back(get_x_pol(a, b));
			//std::cout << get_x_pol(a, b) << std::endl;
		}
		delta /= 2;
		/*for (int i = 0; i < roots_leg.size(); i++)
		{
			std::cout << roots_leg[i] << " ";
		}
		std::cout << std::endl;*/
	}
	
}

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

void gauss(double **a, double **y, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	//x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			std::cout << "Решение получить невозможно из-за нулевого столбца ";
			std::cout << index << " матрицы A" << std::endl;
			return;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = (*y)[k];
		(*y)[k] = (*y)[index];
		(*y)[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			(*y)[i] = (*y)[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			(*y)[i] = (*y)[i] - (*y)[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		//x[k] = y[k];
		for (int i = 0; i < k; i++)
			(*y)[i] = (*y)[i] - a[i][k] * (*y)[k];
	}
	//return x;
}

void create_table(int n, double ***matr, double **vect, std::vector<double> roots_leg)
{
	*matr = allocate(n, n);
	*vect = (double*)calloc(n, sizeof(double));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			(*matr)[i][j] = pow(roots_leg[j], i);
		}
	}
	for (int i = 0; i < n; i++)
	{
		if (i % 2 == 0)
			(*vect)[i] = 2. / (i + 1);
		else
			(*vect)[i] = 0;
	}

	/*for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << (*matr)[i][j] << " ";
		}
		std::cout << "\n";
	}
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		std::cout << (*vect)[i] << " ";
	}*/

}

#define _USE_MATH_DEFINES // for C++  
#include <cmath>  
const double MY_PI = 3.14159265358979323846;

double func(double x)
{
	return x * x;
}

double part_func_laplas(double x)
{
	return exp(-x * x / 2);
}

double from_t_to_x(double a, double b, double t)
{
	return ((b - a) / 2) * t + ((b + a) / 2);
}


double find_integral(double a, double b, std::vector<double> roots_leg, double *vect, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++)
	{
		res += vect[i] * part_func_laplas(from_t_to_x(a, b, roots_leg[i]));
	}
	res *= ((b - a) / 2);
	return res;
}



double func_laplas(double x, double alpha, int n)
{
	double a = 0;
	double b = x;

	std::vector<double> roots_leg;
	double **matr;
	double *vect;
	finds_roots_lagrange(roots_leg, n);
	create_table(n, &matr, &vect, roots_leg);
	double *ans;
	//ans = gauss(matr, vect, n);
	gauss(matr, &vect, n);
	double res = find_integral(a, b, roots_leg, vect, n);

	free_matrix(matr, n);
	free(vect);
	//free(ans);
	roots_leg.clear();

	res *= (1 / sqrt(2 * MY_PI));
	res -= alpha;
	return res;
}


double find_x_in_laplas(double alpha, int n)
{
	double c = 0;
	double eps = 1e-7;
	double a = 0;
	double b = 10;

	while (fabs(b - a) > eps)
	{
		c = get_x_pol(a, b);
		//std::cout << "res =  " << c << "\n";
		//std::cout << "a b" << a << " " << b << "\n";
		//std::cout << "a b" << " " << a << " " << b << "\n";
		//std::cout << "func_laplas(b, alpha) - func_laplas(c, alpha)" << " " << func_laplas(b, alpha) << " " << func_laplas(c, alpha) << "\n";
		if (func_laplas(b, alpha, n) * func_laplas(c, alpha, n) <= 0)
			a = c;
		else if (func_laplas(a, alpha, n) * func_laplas(c, alpha, n) < 0)
			b = c;
	}
	return get_x_pol(a, b);
}



int main(void)
{
	int n = 90;
	std::vector<double> roots_leg;
	double **matr;
	double *vect;
	//std::cout << "leg " << get_polinom_leg(n, 2) << std::endl;
	finds_roots_lagrange(roots_leg, n);
	for (int i = 0; i < roots_leg.size(); i++)
	{
		std::cout << roots_leg[i] << " ";
	}
	std::cout << std::endl;
	create_table(n, &matr, &vect, roots_leg);

	//gauss_method(matr, vect, n);
	//double *vect1;
	gauss(matr, &vect, n);
	/*for (int i = 0; i < n; i++)
	{
		std::cout << vect[i] << " ";
	}*/
	double a = 0;
	double b = 5;
	double res = find_integral(a, b, roots_leg, vect, n);
	//std::cout << "roots_leg size = " << roots_leg.size() << std::endl;
	std::cout << "res = " << res << std::endl;
	//std::cout << "lapl = " << func_laplas(1, 0) << "\n";
	double alpha;
	std::cout << "Enter a: ";
	std::cin >> alpha;
	std::cout << "Enter n: ";
	std::cin >> n;
	//double alpha = 0.341;
	std::cout << "lapl = " << find_x_in_laplas(alpha, n) << "\n";
	std::cout << "check = " << func_laplas(find_x_in_laplas(alpha,n), alpha, n) << "\n";
	//std::cout << "check = " << func_laplas(2.48, 0, n) << "\n";

	system("pause");
    return 0;
	
}

