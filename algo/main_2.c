#include "stdio.h"
#include "stdlib.h"
#include "math.h"

#define NUMBER_OF_POINTS 20
#define NUMBER_OF_POINTS_X 10
#define NUMBER_OF_POINTS_Y 10
#define NUMBER_OF_DEGREE_POLINOM 6
#define NUMBER_OF_POINTS_EQUATION 4

double func(double x, double y)
{
	double pi = 3.1415926;
	//double y = x*x - 2;
	//double y = sin(x) + cos(2 * x) * cos(2 * x) + x * x * x;
	//double z = x + y*y;
	double z = sin(pi*x/6) - y*y;
	return z;	
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
void create_polinom_table(double ***matrix, double *x_num , double *y_num , int num_points, int max_degree_polinom)
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

double find_inter_res(double **matrix, int num_points, int max_degree_polinom, double x, int degree_polinom)
{
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

/* double find_root(double (*func)(double), double x1, double x2, int *err)
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
} */
double **fill_num(int kol_x, int kol_y, double (*func)(double, double))
{
	double **matrix_res = allocate_matrix_row(kol_x + 1, kol_y + 1);
	int begin_x = -5;
	int begin_y = -5;
	for (int i = 1; i < kol_x + 1; i++)
	{
		matrix_res[i][0] = begin_x;
		for (int j = 1; j < kol_y + 1; j++)
		{
			if (i == 1)
				matrix_res[0][j] = begin_y;
			matrix_res[i][j] = func(matrix_res[i][0], matrix_res[0][j]);
			begin_y++;
		}
		begin_x++;
	}
	return matrix_res;
}

int find_begin(double val, double n, double **source_matrix, int kol_x, int kol_y)
{
	int begin;

	if (val < source_matrix[0][0])
		begin = 1;
	else if (val > source_matrix[0][kol_x])
		begin = kol_x + 1 - n - 1;
	else
	{
		for (int i = 1; i < kol_x; i++)
		{
			if (source_matrix[i][0] <= val && source_matrix[0][i + 1] > val)
			{
				begin = i;
				break;
			}
		}
		if (begin + (n + 1) > kol_x)
			begin = kol_x + 1 - n - 1;
	}
	// printf("!! %d\n", begin);
	return begin;
}

double find_res(double val_x, double val_y, double n_x, double n_y, double **source_matrix, int kol_x, int kol_y)
{
	int begin_x = find_begin(val_x, n_x, source_matrix, kol_x, kol_y);
	int begin_y = find_begin(val_y, n_y, source_matrix, kol_x, kol_y);
	double **matrix_inter;
	
	double *x_num = calloc(n_x + 1, sizeof(double));
	for (int i = 0; i < n_x + 1; i++)
	{
		x_num[i] = source_matrix[begin_x + i][0];
	}
	double *y_num = calloc(n_y + 1, sizeof(double));
	for (int i = 0; i < n_y + 1; i++)
	{
		y_num[i] = source_matrix[0][begin_y + i];
	}
	/* printf("\n////\n");
	for (int i = 0; i < n_x + 1; i++)
	{
		printf("%lf ", x_num[i]);
	}
	printf("\n////\n");
	for (int i = 0; i < n_y + 1; i++)
	{
		printf("%lf ", y_num[i]);
	}
	printf("\n////\n"); */
	
	double *z_y_num = calloc(n_y + 1, sizeof(double));
	double *z_x_num = calloc(n_x + 1, sizeof(double));
	
	for (int i = begin_x; i < begin_x + n_x + 1; i++)
	{
		for (int j = begin_y; j < begin_y + n_y + 1; j++)
		{
			z_y_num[j - begin_y] = source_matrix[i][j];
		}
		create_polinom_table(&matrix_inter, y_num, z_y_num, n_y + 1, n_y);
		z_x_num[i - begin_x] = find_inter_res(matrix_inter, n_y + 1, n_y, val_y, n_y);
		free_matrix_rows(matrix_inter, n_y + 1);
	}
	create_polinom_table(&matrix_inter, x_num, z_x_num, n_x + 1, n_x);
	double res = find_inter_res(matrix_inter, n_x + 1, n_x, val_x, n_x);

	return res;
}

void print_begin_table(double **source_matrix, int k_x, int k_y)
{
	for (int i = 0; i < k_x + 1; i++)
	{
		for (int j = 0; j < k_y + 1; j++)
		{
			if (i == 0 && j == 0)
				printf("   x\\y  ");
			else
				printf("%7.2f ", source_matrix[i][j]);
		}
		printf("\n");
	} 
}

int main(void)
{
	setbuf(stdout,NULL); 
	double val_x;
	double val_y;
	int n_x;
	int n_y;
	
	double **source_matrix = fill_num(NUMBER_OF_POINTS_X, NUMBER_OF_POINTS_Y, func);
	
	
	
	printf("Программа вычисляет значение, используя многомерную интерполяцию\n");
	printf("Таблица точек:\n");
	print_begin_table(source_matrix, NUMBER_OF_POINTS_X, NUMBER_OF_POINTS_Y);
	
	while (1)
	{
		printf("Введите x и степень полинома для x через пробел: ");	
		if (scanf("%lf %d", &val_x, &n_x) != 2)
		{
			printf("Некорректно введены значения\n");
			fflush(stdin);
			continue;
		}
		break;
	}
	while (1)
	{
		printf("Введите y и степень полинома для y через пробел: ");	
		if (scanf("%lf %d", &val_y, &n_y) != 2)
		{
			printf("Некорректно введены значения\n");
			fflush(stdin);
			continue;
		}
		break;
	}
	
	double res = find_res(val_x, val_y, n_x, n_y, source_matrix, NUMBER_OF_POINTS_X, NUMBER_OF_POINTS_Y);
	printf("результат = %lf\nверный результат = %lf", res, func(val_x, val_y));	
}
