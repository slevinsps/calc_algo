#include "glwidget.h"
#include <gl\gl.h>    // Заголовочный файл для библиотеки OpenGL32
#include <gl\glu.h>
#include <iostream>
#include "mainwindow.h"



GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
    //connect(&timer, SIGNAL(timeout()), this, SLOT(update()));
    // timer.start(16);
    //this->angle = 0;
}

typedef struct options_ss
{
    double x_n;
    double y_n;
    double x_k;
    double y_k;
    int method;
    int lines;
    double angle;
    double rad;
    QColor color;
}options_s;

extern options_s  options;
extern QColor color_fone;
extern std::vector<options_s> opt;


void GLWidget::initializeGL()
{

    /*glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);*/

    // glClearColor(1.0, 1.0, 1.0, 1.0);
    glShadeModel( GL_FLAT );
}

void rotate_func(double *x, double *y, double angle, double xc, double yc)
{

    double x_rot = xc + (*x - xc)*cos(angle) + (*y - yc)*sin(angle);
    double y_rot = yc - (*x - xc)*sin(angle) + (*y - yc)*cos(angle);
    *x = x_rot;
    *y = y_rot;
}

void scale_func(double *x, double *y, double k_x, double k_y, double x_scale, double y_scale)
{
    k_x *= 1;
    k_y *= 1;

    double x_sc = k_x * (*x) + (1 - k_x) * x_scale;
    double y_sc = k_y * (*y) + (1 - k_y) * y_scale;
    *x = x_sc;
    *y = y_sc;
}

int sign(double x)
{
    if (x < 0)
        return -1;
    else if (x == 0)
        return 0;
    else
        return 1;

}


void draw_line_cda(options_s  options, int paint)
{
    glColor3f(options.color.redF(),options.color.greenF(),options.color.blueF());
    double dx = fabs(options.x_k - options.x_n);
    double dy = fabs(options.y_k - options.y_n);
    double l;
    if (dx > dy)
        l = dx;
    else
        l = dy;
    double sx = (options.x_k - options.x_n)/l;
    double sy = (options.y_k - options.y_n)/l;
    double x = options.x_n;
    double y = options.y_n;
    glBegin( GL_POINTS ); // тело
        for (int i = 1; i <= l+1; i++)
        {
             if (paint)
             {
                 glVertex2i(round(x), round(y));
             }
             x += sx;
             y += sy;
        }
    glEnd();
}


void draw_brez_double(options_s  options, int paint)
{
    glColor3f(options.color.redF(),options.color.greenF(),options.color.blueF());
    int dx = options.x_k - options.x_n;
    int dy = options.y_k - options.y_n;
    int sx = sign(dx);
    int sy = sign(dy);
    dx = abs(dx);
    dy = abs(dy);
    int obmen;
    if (dy > dx)
    {
        int t = dy;
        dy = dx;
        dx = t;
        obmen = 1;
    }
    else
        obmen = 0;
    double m = dy/dx;
    double e = m - 0.5;
    int x = options.x_n;
    int y = options.y_n;
    glBegin( GL_POINTS ); // тело
        for (int i = 1; i <= dx+1; i++)
        {
             //glVertex2f(floor(x+0.5), floor(y+0.5));
             if (paint)
             {
                 glVertex2d(x, y);
             }
             if (e >= 0)
             {
                 if (obmen)
                     x += sx;
                 else
                     y += sy;
                 e -= 1;
             }

             if (!obmen)
                 x += sx;
             else
                 y += sy;
             e += m;
        }
   glEnd();
}

void draw_brez_steps(options_s  options, int paint)
{
    double dx = options.x_k - options.x_n;
    double dy = options.y_k - options.y_n;
    double sx = sign(dx);
    double sy = sign(dy);
    dx = fabs(dx);
    dy = fabs(dy);
    double obmen;
    if (dy > dx)
    {
        double t = dy;
        dy = dx;
        dx = t;
        obmen = 1;
    }
    else
        obmen = 0;
    double m = dy/dx;
    double I = 255;
    double e = 0.5 * I;
    m *= I;
    double w = I - m;

    double x = options.x_n;
    double y = options.y_n;

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin( GL_POINTS ); // тело
        for (int i = 0; i <= dx; i++)
        {
             //glVertex2f(floor(x+0.5), floor(y+0.5));
             if (paint)
             {
                 qDebug() << "e = " << e << "w = " << w << "I/e = " << e/I << endl;
                 //options.color.setAlphaF(e);
                 glColor4f(options.color.redF(),options.color.greenF(),options.color.blueF(), e/I );
                 glVertex2d(x, y);
             }
             if (e >= w)
             {
                 if (obmen)
                     x += sx;
                 else
                     y += sy;
                 e -= I;
             }
             if (!obmen)
                 x += sx;
             else
                 y += sy;
             e += m;
        }
   glEnd();
   glDisable(GL_BLEND);

}


void draw_brez_int(options_s  options, int paint, int &steps)
{
    steps = 0;
    glColor3f(options.color.redF(),options.color.greenF(),options.color.blueF());
    int dx = options.x_k - options.x_n;
    int dy = options.y_k - options.y_n;
    int sx = sign(dx);
    int sy = sign(dy);
    dx = abs(dx);
    dy = abs(dy);
    int obmen;
    if (dy > dx)
    {
        int t = dy;
        dy = dx;
        dx = t;
        obmen = 1;
    }
    else
        obmen = 0;
    //double m = dy/dx;
    int e = 2*dy - dx;
    int x = options.x_n;
    int y = options.y_n;
    glBegin( GL_POINTS ); // тело
        for (int i = 1; i <= dx+1; i++)
        {
             //glVertex2f(floor(x+0.5), floor(y+0.5));
            if (paint)
            {
                glVertex2i(x, y);
            }
            if (e >= 0)
             {
                 if (obmen)
                     x += sx;
                 else
                     y += sy;
                 e -= 2*dx;
                 steps++;
             }

             if (!obmen)
                 x += sx;
             else
                 y += sy;
             e += 2*dy;
        }
   glEnd();

}

void draw_standart(options_s  options)
{
    glColor3f(options.color.redF(),options.color.greenF(),options.color.blueF());
    glBegin( GL_LINES ); // тело
        glVertex2f(options.x_n, options.y_n);
        glVertex2f(options.x_k, options.y_k);
    glEnd();
}


void draw_lines(double angle, double rad, int i, QColor color, int paint)
{
    int steps;
    double angle_g = -angle * 3.1415926/180;
    double xc = 300;
    double yc = 300;
    double xk = xc + rad;
    double yk = 300;
    double ang = 0;
    while (ang < 360)
    {
        options.x_n = xc;
        options.y_n = yc;
        options.x_k = xk;
        options.y_k = yk;
        options.color = color;
        if (i == 1)
        {
            draw_line_cda(options, paint);
        }
        else if (i == 2)
        {
            draw_brez_double(options, paint);
        }
        else if (i == 3)
        {
            draw_brez_int(options, paint, steps);
        }
        else if (i == 4)
        {
            draw_brez_steps(options, paint);
        }
        else if (i == 5)
        {
            draw_standart(options);
        }
        rotate_func(&xk, &yk, angle_g, xc, yc);
        ang += fabs(angle);
    }
}


void GLWidget::paintGL()
{
    glClearColor(color_fone.redF(), color_fone.greenF(), color_fone.blueF(),1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    int steps;

    int paint = 1;
    for (int i = 0; i < opt.size(); i++)
    {
        if (opt[i].lines)
        {
            draw_lines(opt[i].angle, opt[i].rad, opt[i].method, opt[i].color, paint);
        }
        else
        {

            if (opt[i].method == 1)
            {
                draw_line_cda(opt[i], paint);
            }
            else if (opt[i].method == 2)
            {
                draw_brez_double(opt[i], paint);
            }
            else if (opt[i].method == 3)
            {
                draw_brez_int(opt[i], paint, steps);
            }
            else if (opt[i].method == 4)
            {
                draw_brez_steps(opt[i], paint);
            }
            else if (opt[i].method == 5)
            {
                draw_standart(opt[i]);
            }
        }
    }


}

void GLWidget::resizeGL(int w, int h)
{
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    glOrtho(0, 600, 0, 600, -1.0, 1.0);

    glViewport(0, 0, (GLint)w, (GLint)h);
}
