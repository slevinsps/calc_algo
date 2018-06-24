#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define NUMBER_OF_POINTS 20
#define NUMBER_OF_DEGREE_POLINOM 6
#define NUMBER_OF_POINTS_EQUATION 4

double func(double x)
{
	double pi = 3.1415926;
	//double y = x*x - 2;
	//double y = sin(x) + cos(2 * x) * cos(2 * x) + x * x * x;
	double y = 2*x + 3;
	//double y = sin(pi*x/6);
	return y;	
}

void free_matrix_rows(double **data, int n)
{
    for (int i = 0;i < n; i++)
    {
        free((data)[i]);        
    }
    free(data);    
}
  
// Выделение памяти  
double **allocate_matrix_row(int n, int m)
{
    double **data = (double**)calloc(n,sizeof(double*));
    if (!data)
        return NULL;
    for (int i = 0;i < n;i++)
    {
        data[i] = (double*)calloc(m,sizeof(double));
        if (!data[i])
        {
            free_matrix_rows(data, n);
            return NULL;
        } 
    }
    return data;
}

void create_polinom_table(double ***matrix, double *x_num , double *y_num , int num_points, int max_degree_polinom, double (*func)(double))
{
	int col = max_degree_polinom + 2;
	int row = num_points;
	*matrix = allocate_matrix_row(row,col);

	for (int i = 0; i < row; i++)
	{
		(*matrix)[i][0] = x_num[i];
		(*matrix)[i][1] = y_num[i];
	}
	
	for (int j = 2; j < col; j++)
	{
		for (int i = 0; i < row - j + 1; i++)
		{
			(*matrix)[i][j] = ((*matrix)[i + 1][j - 1] - (*matrix)[i][j - 1])/((*matrix)[i+(j-1)][0] - (*matrix)[i][0]);
		}	
	}
}

void print_matrix(double **matrix, int row, int col)
{
	printf("      x         y       y(x^2)      y(x^3)      y(x^4)    y(x^5)      y(x^6)      y(x^7)\n");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			if (j <= 1 || i <= row - j )
				printf("%10.6f ", matrix[i][j]);
		}
		printf("\n");
	}
}

double find_y(double **matrix, int num_points, int max_degree_polinom, double x, int degree_polinom)
{
	int col = max_degree_polinom + 2;
	int row = num_points;
	
	int i_begin;
	if (x < matrix[0][0])
		i_begin = 0;
	else if (x > matrix[row-1][0])
		i_begin = row - degree_polinom - 1;
	else
	{
		for (int i = 0; i < row-1; i++)
		{
			if (x > matrix[i][0] && x <= matrix[i+1][0])
			{
				
				i_begin = i;
				break;
			}
		}
		if (i_begin + (degree_polinom + 1) > row)
			i_begin = row - degree_polinom - 1;
	}	
			
	double result = 0;
	double mnoj = 1;
	
	for (int j = 1; j <= (degree_polinom + 1); j++)
	{
		mnoj = matrix[i_begin][j];
		for (int i = i_begin; i < i_begin + (j - 1); i++)
		{	
			mnoj *= (x - matrix[i][0]);
		}
		result += mnoj;
		
	}
	return result;
}

double find_root(double (*func)(double), double x1, double x2, int *err)
{
	double eps = 0.1;
	int degree_polinom = 3;
	double **matrix;
	double begin_number;
	int is_find = 0;
	double i = x1;
	while(i <= x2)
	{
		if (func(i) * func(i + eps) <= 0)
		{
			begin_number = i;	
			is_find = 1;
			break;
		}
		i+=eps;
	}
	if (!is_find)
	{
		*err = -1;
		return 0;
	}
	
	begin_number -= eps;
	
	double *x_num = calloc(NUMBER_OF_POINTS_EQUATION, sizeof(double));
	double *y_num = calloc(NUMBER_OF_POINTS_EQUATION, sizeof(double));
	for (int i = 0; i < NUMBER_OF_POINTS_EQUATION; i++)
	{
		x_num[i] = begin_number;
		y_num[i] = func(x_num[i]);
		begin_number += eps;
	}
	create_polinom_table(&matrix, y_num, x_num, NUMBER_OF_POINTS_EQUATION, degree_polinom, func);
	double res = find_y(matrix, NUMBER_OF_POINTS_EQUATION, degree_polinom, 0, degree_polinom);
	free(x_num);
	free(y_num);
	free_matrix_rows(matrix,NUMBER_OF_POINTS_EQUATION);
	return res;
}





void fill_num(double *x_num, double *y_num, int kol, double (*func)(double))
{
	int begin = -10;
	for (int i = 0; i < kol; i++)
	{
		x_num[i] = begin;
		y_num[i] = func(x_num[i]);
		begin += 1;
	}
}


int main(void)
{
	setbuf(stdout,NULL); 
	int row = NUMBER_OF_POINTS;
	int col = NUMBER_OF_DEGREE_POLINOM + 2;
	double **matrix;
	double begin_number = 0, step = 1;
	double val;
	int n;
	int choose;
	while(1)
	{
		printf("1 - вычислить значение, используя интерполяционный полином\n");
		printf("2 - найти корень функции\n\n");
		printf("0 - выход\n");
		printf("Выбор: ");
		if (scanf("%d", &choose) != 1 || choose > 2 || choose < 0)
		{
			
			printf("Некорректно введен пункт меню\n");
		
			fflush(stdin);
			continue;
		}
		if (choose == 1)
		{
			printf("Программа вычисляет значение, используя интерполяционный полином\n");
			printf("Таблица точек:\n");
			double *x_num = calloc(NUMBER_OF_POINTS, sizeof(double));
			double *y_num = calloc(NUMBER_OF_POINTS, sizeof(double));
			fill_num(x_num, y_num, NUMBER_OF_POINTS, func);
			create_polinom_table(&matrix, x_num, y_num, NUMBER_OF_POINTS, NUMBER_OF_DEGREE_POLINOM, func);
			print_matrix(matrix, row, col);
			while (1)
			{
				printf("Введите x : ");	
				if (scanf("%lf", &val) != 1)
				{
					printf("Некорректно введен x\n");
					fflush(stdin);
					continue;
				}
				break;
			}
			while (1)
			{
				printf("Введите степень полинома : ");	
				if (scanf("%d", &n) != 1)
				{
					printf("Некорректно введена степень полинома\n");
					fflush(stdin);
					continue;
				}
				break;
			}
			
			double res = find_y(matrix, NUMBER_OF_POINTS, NUMBER_OF_DEGREE_POLINOM, val, n);
			printf("Полученное значение: y = %f\nРеальное значение: y = %f\n", res, func(val));
			free(x_num);
			free(y_num);
			free_matrix_rows(matrix, NUMBER_OF_POINTS);
		}
		if (choose == 2)
		{
			double x1, x2;
			while (1)
			{
				printf("Введите интервал поиска корня через пробел: ");	
				if (scanf("%lf %lf", &x1, &x2) != 2 || x1 >= x2)
				{
					printf("Некорректно введен интервал\n");
					fflush(stdin);
					continue;
				}
				break;
			}
			int err = 1;
			double res_eq = find_root(func, x1, x2, &err);
			if (err == -1)
				printf("Нет корней\n");
			else
				printf("Значение корня = %f\n", res_eq);
		}
		if (choose == 0)
		{
			break;
		}
	}
}
