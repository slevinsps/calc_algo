#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QWidget>
#include <QMainWindow>
#include <QtGui>
#include <QtCore>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_pushButton_5_clicked();

private:
    Ui::MainWindow *ui;
};

#include <QVBoxLayout>
#include <QPushButton>
#include "qcustomplot.h"




class Widget_graf : public QWidget
{
    Q_OBJECT

public:
    explicit Widget_graf(QWidget *parent = 0);

private:
    QCustomPlot *customPlot;
    QVBoxLayout *vbox;
    QPushButton *res;
    QCPBars *fossil;

};

#endif // MAINWINDOW_H
