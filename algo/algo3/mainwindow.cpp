#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "windows.h"
#include "math.h"
#include "time.h"

#include "QStandardItemModel"
#include "QStandardItem"
#include <QMessageBox>
#include <QColorDialog>

#include <stdio.h>    // Заголовочный файл для стандартной библиотеки ввода/вывода
#include <iostream>
#include <gl\gl.h>    // Заголовочный файл для библиотеки OpenGL32
#include <gl\glu.h>   // Заголовочный файл для библиотеки GLu32
//#include <gl\glaux.h> // Заголовочный файл для библиотеки GLaux
#include <math.h>     // Заголовочный файл для математической библиотеки ( НОВОЕ )
#include <stdarg.h>   // Заголовочный файл для функций для работы с переменным



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

}

MainWindow::~MainWindow()
{
    delete ui;
}


#include "qcustomplot.h"

#include <iostream>
#include <vector>
#include <cstdio>

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_DEPRECATE
#pragma warning(disable : 4996)

typedef struct data_s
{
    double x;
    double y;
    double p;
}data_points;


double func(double x, int k)
{

    double y = pow(x, k);
    return y;
}

double func2(double x, int k)
{

    double y = sin(x*k);
    return y;
}


using namespace std;

/*void getMatrixWithoutRowAndCol(double **matrix, int size, int row, int col, double **newMatrix) {
    int offsetRow = 0; //Смещение индекса строки в матрице
    int offsetCol = 0; //Смещение индекса столбца в матрице
    for(int i = 0; i < size-1; i++) {
        //Пропустить row-ую строку
        if(i == row) {
            offsetRow = 1; //Как только встретили строку, которую надо пропустить, делаем смещение для исходной матрицы
        }

        offsetCol = 0; //Обнулить смещение столбца
        for(int j = 0; j < size-1; j++) {
            //Пропустить col-ый столбец
            if(j == col) {
                offsetCol = 1; //Встретили нужный столбец, проускаем его смещением
            }

            newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
        }
    }
}

double Determinant(double  **matrix, int size) {
    double det = 0;
    int degree = 1;

    if(size == 1) {
        return matrix[0][0];
    }
    //Условие выхода из рекурсии
    else if(size == 2) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else {
        //Матрица без строки и столбца
        double **newMatrix = new double*[size-1];
        for(int i = 0; i < size-1; i++) {
            newMatrix[i] = new double[size-1];
        }

        for(int j = 0; j < size; j++) {
            //Удалить из матрицы i-ю строку и j-ый столбец
            //Результат в newMatrix
            getMatrixWithoutRowAndCol(matrix, size, 0, j, newMatrix);
            det = det + (degree * matrix[0][j] * Determinant(newMatrix, size-1));
            degree = -degree;
        }

        for(int i = 0; i < size-1; i++) {
            delete [] newMatrix[i];
        }
        delete [] newMatrix;
    }

    return det;
}
*/

#include <cstdlib>

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


void create_base_matrix(vector<data_points> data_array, double ***matr, double **vect, int n, double (*func)(double x, int k))
{
    *matr = allocate(n+1, n+1);
    *vect = (double*)calloc(n+1, sizeof(double));
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            for (int k = 0; k < data_array.size(); k++)
            {
                (*matr)[i][j] += func(data_array[k].x, i) * func(data_array[k].x, j) * data_array[k].p;
                //printf("!!!! %lf \n", (*matr)[i][j]);
            }

        }
        for (int k = 0; k < data_array.size(); k++)
        {
            (*vect)[i] += func(data_array[k].x, i) * data_array[k].y * data_array[k].p;
        }
    }

    /*for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            printf("%lf ",(*matr)[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    for (int i = 0; i <= n; i++)
    {
        printf("%lf ", (*vect)[i]);
    }
    printf("\n");*/
}




/*double determ_col(double **matr, double *vect, int n, int k)
{
    vector<double> copy;
    for (int i = 0; i <= n; i++)
    {
        copy.push_back(matr[i][k]);
        matr[i][k] = vect[i];
    }

    double res = Determinant(matr, n + 1);

    for (int i = 0; i <= n; i++)
    {
        matr[i][k] = copy[i];
    }
    return res;
}

vector<double> find_koef_kramer(double **matr, double *vect, int n, int &err)
{
    vector<double> koef;
    double det_begin = Determinant(matr, n + 1);
    //qDebug() << "! " << det_begin;
    if (fabs(det_begin) < 0.001)
    {
        err = -1;
        return koef;
    }
    for (int i = 0; i <= n; i++)
    {
        double delta_det = determ_col(matr, vect, n, i);
        koef.push_back(delta_det / det_begin);
        //qDebug() << "! " << delta_det << " "<< det_begin;
    }
    return koef;
}
*/
// Основная функция


int not_zero_to_top(double **matr, double *vect, int col, int n)
{
    int not_zero_number = -1;
    for(int i = col; i <= n; i++)
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
    for (int i = 0; i <= n; i++)
    {
        row1[i] -= mult * row2[i];
    }
    num1 -= num2 * mult;
}
void devide_row(double *row, double &num, int col, int n)
{
    double div = row[col];
    if (div == 0) return;
    for (int i = 0; i <= n; i++)
    {
        row[i] /= div;
    }
    num /= div;
}

int check_result(double **matr, double *vect, int n)
{
    double diap = 0.0001;
    for (int i = 0; i <= n; i++ )
    {
        if (fabs(matr[i][i]) < diap && fabs(vect[i]) > diap)
            return -1;
    }
    for (int i = 0; i <= n; i++ )
    {
        if (fabs(matr[i][i]) < diap && fabs(vect[i]) < diap)
            return -2;
    }
    return 0;
}

int gauss_method(double **matr, double *vect, int n)
{
    int err = 0;
    for (int i = 0; i <= n; i++)
    {
        if (not_zero_to_top(matr, vect, i, n))
        {
            devide_row(matr[i], vect[i], i, n);
            for (int j = 0; j <= n; j++)
            {

                if (j != i)
                {
                    subtraction_mult_rows(matr[j], matr[i], vect[j], vect[i], n, matr[j][i]);
                }
            }
        }
    }

    err =  check_result(matr, vect, n);

    /*for (int i = 0; i <= n; i++)
      {
          for (int j = 0; j <= n; j++)
          {
              printf("%lf ",matr[i][j]);
          }
          printf("\n");
      }
      printf("\n");
      for (int i = 0; i <= n; i++)
      {
          printf("%lf ", vect[i]);
      }
      printf("\ngggggg\n");*/
      return err;
}


typedef struct steps_angle_s
{
    int angle;
    int steps;
}steps_angle;



QVector<double> x_arr, y_arr;

QVector<double> x_arr2, y_arr2;


typedef struct max_min_s
{
    double x_max;
    double y_max;
    double x_min;
    double y_min;
}max_min;

max_min znach;
max_min y_border;
vector<data_points> data_array;

void points_for_graf(vector<double> koef, int n, double (*func)(double x, int k))
{
    y_border.y_max = znach.y_max;
    y_border.y_min = znach.y_min;
    double i = znach.x_min;

    while (i <= znach.x_max)
    {
        //qDebug() << "!!!!!  =  " << i;
        x_arr.push_back(i);
        double res = 0;
        for(int j = 0; j <= n; j++)
        {
            res += koef[j] * func(i,j);
        }
        if (res > y_border.y_max)
            y_border.y_max = res;
        if (res < y_border.y_min)
            y_border.y_min = res;

        y_arr.push_back(res);
        //qDebug() << i << " " << res;
        i += 0.001;
    }
}

void points_for_graf2(std::vector<data_points> points)
{
    for(int j = 0; j < points.size(); j++)
    {
        x_arr2.push_back(points[j].x);
        y_arr2.push_back(points[j].y);
    }
}


max_min find_max_min(std::vector<data_points> points)
{
    max_min coords;

    double min_x = points[0].x;
    double max_x = points[0].x;
    double min_y = points[0].y;
    double max_y = points[0].y;

    for (int i = 1; i < points.size(); i++)
    {
        if (min_x > points[i].x)
            min_x = points[i].x;
        if (max_x < points[i].x)
            max_x = points[i].x;

        if (min_y > points[i].y)
            min_y = points[i].y;
        if (max_y < points[i].y)
            max_y = points[i].y;
    }

    coords.x_max = max_x;
    coords.y_max = max_y;
    coords.x_min = min_x;
    coords.y_min = min_y;

    return coords;
}

Widget_graf::Widget_graf(QWidget *parent) :
    QWidget(parent)
{


    resize(500,400);
    setWindowTitle(QString::fromUtf8("График"));
    customPlot = new QCustomPlot(this);
    vbox = new QVBoxLayout(this);


    vbox->addWidget(customPlot);

    setLayout(vbox);


    customPlot->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    //customPlot->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    customPlot->legend->setVisible(true);   //Включаем Легенду графика
        // Устанавливаем Легенду в левый верхний угол графика
        customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);

        // Инициализируем график и привязываем его к Осям
        QCPGraph *graphic = new QCPGraph(customPlot->xAxis, customPlot->yAxis);
        customPlot->addPlottable(graphic);  // Устанавливаем график на полотно
        //graphic->setName("");       // Устанавливаем
        graphic->setPen(QPen(QColor(Qt::red))); // Устанавливаем цвет графика
        graphic->setAntialiased(false);         // Отключаем сглаживание, по умолчанию включено
        //graphic->setLineStyle(QCPGraph::lsImpulse); // График в виде импульсных тиков

    graphic->setData(x_arr, y_arr);


    QCPGraph *graphic2 = new QCPGraph(customPlot->xAxis, customPlot->yAxis);
    customPlot->addPlottable(graphic2);  // Устанавливаем график на полотно



    graphic2->setPen(QColor(50, 50, 50, 255));//задаем цвет точки
    graphic2->setLineStyle(QCPGraph::lsNone);//убираем линии
    graphic2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4)); //формируем вид точе

    graphic2->setData(x_arr2, y_arr2);

    //Подписываем оси Ox и Oy
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");

    //Установим область, которая будет показываться на графике



    customPlot->xAxis->setRange(znach.x_min - 2, znach.x_max + 2);//Для оси Ox

    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах

    customPlot->yAxis->setRange(y_border.y_min - 2, y_border.y_max + 2);//Для оси Oy
    customPlot->legend->setVisible(true);
    //И перерисуем график на нашем widget
    customPlot->replot();
}

#include <fstream>
int load_file(std::vector<data_points> &points)
{

    data_points coords;
    QString fileName;
    fileName = QFileDialog::getOpenFileName(0,
                                QString::fromUtf8("Загрузить файл"),
                                QDir::currentPath(),
                                "(*.txt *.doc);;All files (*.*)");
    if (!fileName.isEmpty())
    {
        //points.clear();

        QFile file(fileName);

        if (!file.open(QIODevice::Text | QIODevice::ReadOnly)) return -1; // пытаемся открыть файл
         // читаем первую строку

         // делим строку на кусочки

        while(!file.atEnd())
        {
            QByteArray ba = file.readLine();
            ba = ba.simplified();
            QList<QByteArray> baList = ba.split(' ');
            if (!baList[0].isEmpty())
            {
                coords.x = baList[0].toDouble();
                coords.y = baList[1].toDouble();
                coords.p = baList[2].toDouble();
                points.push_back(coords);
                //qDebug() << baList;
            }
        }
        file.close();
        return 0;
    }
    return -1;

}


int load_file_test(double ***matr, double **vector, int &n)
{
    QString fileName;
    fileName = QFileDialog::getOpenFileName(0,
                                QString::fromUtf8("Загрузить файл"),
                                QDir::currentPath(),
                                "(*.txt *.doc);;All files (*.*)");
    if (!fileName.isEmpty())
    {

        //points.clear();

        QFile file(fileName);

        if (!file.open(QIODevice::Text | QIODevice::ReadOnly)) return -1; // пытаемся открыть файл
         // читаем первую строку

         // делим строку на кусочки

        QString s = file.readLine();
        s = s.simplified();
        n = s.toInt();
        *matr = allocate(n, n);
        *vector = (double*)calloc(n, sizeof(double));

        for(int i = 0; i < n; i++)
        {
            QByteArray ba = file.readLine();
            ba = ba.simplified();
            QList<QByteArray> baList = ba.split(' ');
            for (int j = 0; j < n; j++)
            {
                (*matr)[i][j] = baList[j].toDouble();
            }

        }

            QByteArray ba = file.readLine();
            ba = ba.simplified();
            QList<QByteArray> baList = ba.split(' ');
            for(int i = 0; i < n; i++)
            {
                (*vector)[i] = baList[i].toDouble();
            }

        file.close();


        for (int i = 0; i < n; i++)
          {
              for (int j = 0; j < n; j++)
              {
                  printf("%lf ",(*matr)[i][j]);
              }
              printf("\n");
          }
          printf("\n");
          for (int i = 0; i < n; i++)
          {
              printf("%lf ", (*vector)[i]);
          }
          printf("\ngggggg\n");
        return 0;
    }



    return -1;
}

void MainWindow::on_pushButton_5_clicked()
{

    vector<double> koef;
    data_array.clear();
    x_arr.clear();
    y_arr.clear();
    x_arr2.clear();
    y_arr2.clear();


    int err = 0;

    err = load_file(data_array);

    if (err == 0)
    {
        QString n_s = ui->lineEdit_n->text();

        bool ok;
        int n = n_s.toInt(&ok);
        if (!ok || n < 0)
        {
            QMessageBox::information(this, "Ошибка", "Неверно введена степеь!", "Ok");
            return;
        }


        double **matr;
        double *vect;
        //if (ui->radioButton->isChecked())
        //{
                                create_base_matrix(data_array, &matr, &vect, n, func);
        //}
        /*if (ui->radioButton_2->isChecked())
        {
            create_base_matrix(data_array, &matr, &vect, n, func2);
        }*/



        /*for (int i = 0; i <= n; i++)
            {
                for (int j = 0; j <= n; j++)
                {
                    printf("%lf ",matr[i][j]);
                }
                printf("\n");
            }
            printf("\n");
            for (int i = 0; i <= n; i++)
            {
                printf("%lf ", vect[i]);
            }*/

        int err = 0;
        /*load_file_test(&matr, &vect, n);*/
        err = gauss_method(matr, vect, n);
        if (err == -1)
        {
            QMessageBox::information(this, "Ошибка", "Определитель равен нулю!", "Ok");
            return;
        }
        if (err == -2)
        {
            QMessageBox::information(this, "Ошибка", "Бесконечное множество решений!", "Ok");
            return;
        }
        koef.clear();
        for (int i = 0; i <= n; i++)
        {
            koef.push_back(vect[i]);
        }
        /*koef = find_koef_kramer(matr, vect, n, err);
        if (err == -1)
        {
            QMessageBox::information(this, "Ошибка", "Определитель равен нулю!", "Ok");
            return;
        }*/

        /*for (int i = 0; i < koef.size(); i++)
        {
            qDebug() << koef[i] << " ";
        }*/
                              znach = find_max_min(data_array);

        //if (ui->radioButton->isChecked())
        //{
                              points_for_graf(koef, n, func);
        //}
        //if (ui->radioButton_2->isChecked())
        //{
        //    points_for_graf(koef, n, func2);
        //}
                               points_for_graf2(data_array);

        Widget_graf *scene = new Widget_graf;
        scene->resize(600, 600);

        //scene->rndres();
        scene->show();
        free_matrix(matr, n+1);
        free(vect);
        koef.clear();
    }
}
