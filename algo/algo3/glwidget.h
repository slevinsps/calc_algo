#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QTimer>

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(QWidget *parent = 0);
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    double x_n;
    double y_n;
    double x_k;
    double y_k;
    int method;
    QColor color;

private:
    QTimer timer;
};

#endif // GLWIDGET_H
