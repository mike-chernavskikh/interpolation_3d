#include "glwidget.h"

#include <stdio.h>
#include <math.h>

#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define MAX(a, b) ((a) < (b) ? (b) : (a))
#define sign(a) ((a) > (0) ? (1) : (-1))
double f(double x, double y);
double dfdx(double x, double y);
double dfdy(double x, double y);
double ddfdxdy(double x, double y);
double dfdx_2(double x, double y, double dx, double dy);
double dfdy_2(double x, double y, double dx, double dy);
double ddfdxdy_2(double x, double y, double dx, double dy);

void myGLWidget::initializeGL()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);

	setDefaultCamera();
}

void A_inversed(double h, double a[4][4])
{
	double h_1 = 1 / h;
	double h_2 = 1 / (h * h); 
	double h_3 = h_2 / h;
	a[0][0] = 1;
	a[0][1] = 0;
	a[0][2] = 0;
	a[0][3] = 0;

	a[1][0] = 0;
	a[1][1] = 1;
	a[1][2] = 0;
	a[1][3] = 0;

	a[2][0] = -3 * h_2;
	a[2][1] = -2 * h_1;
	a[2][2] = 3 * h_2;
	a[2][3] = -h_1;

	a[3][0] = 2 * h_3;
	a[3][1] = h_2;
	a[3][2] = -2 * h_3;
	a[3][3] = h_2;
}

void A_inversed_T(double h, double a[4][4])
{
	double h_1 = 1 / h;
	double h_2 = 1 / (h * h); 
	double h_3 = h_2 / h;
	a[0][0] = 1;
	a[1][0] = 0;
	a[2][0] = 0;
	a[3][0] = 0;

	a[0][1] = 0;
	a[1][1] = 1;
	a[2][1] = 0;
	a[3][1] = 0;

	a[0][2] = -3 * h_2;
	a[1][2] = -2 * h_1;
	a[2][2] = 3 * h_2;
	a[3][2] = -h_1;

	a[0][3] = 2 * h_3;
	a[1][3] = h_2;
	a[2][3] = -2 * h_3;
	a[3][3] = h_2;
}

void F_matrix(double x_i, double x_i_1, double y_i, double y_i_1, double F[4][4])
{
	F[0][0] = f(x_i, y_i);
	F[0][1] = dfdy(x_i, y_i);
	F[0][2] = f(x_i, y_i_1);
	F[0][3] = dfdy(x_i, y_i_1);

	F[1][0] = dfdx(x_i, y_i);
	F[1][1] = ddfdxdy(x_i, y_i);
	F[1][2] = dfdx(x_i, y_i_1);
	F[1][3] = ddfdxdy(x_i, y_i_1);

	F[2][0] = f(x_i_1, y_i);
	F[2][1] = dfdy(x_i_1, y_i);
	F[2][2] = f(x_i_1, y_i_1);
	F[2][3] = dfdy(x_i_1, y_i_1);

	F[3][0] = dfdx(x_i_1, y_i);
	F[3][1] = ddfdxdy(x_i_1, y_i);
	F[3][2] = dfdx(x_i_1, y_i_1);
	F[3][3] = ddfdxdy(x_i_1, y_i_1);
}

void F_matrix_2(double x_i, double x_i_1, double y_i, double y_i_1, double F[4][4])
{
	F[0][0] = f(x_i, y_i);
	F[0][1] = dfdy_2(x_i, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[0][2] = f(x_i, y_i_1);
	F[0][3] = dfdy_2(x_i, y_i_1, x_i_1 - x_i, y_i_1 - y_i);

	F[1][0] = dfdx_2(x_i, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[1][1] = ddfdxdy_2(x_i, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[1][2] = dfdx_2(x_i, y_i_1, x_i_1 - x_i, y_i_1 - y_i);
	F[1][3] = ddfdxdy_2(x_i, y_i_1, x_i_1 - x_i, y_i_1 - y_i);

	F[2][0] = f(x_i_1, y_i);
	F[2][1] = dfdy_2(x_i_1, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[2][2] = f(x_i_1, y_i_1);
	F[2][3] = dfdy_2(x_i_1, y_i_1, x_i_1 - x_i, y_i_1 - y_i);

	F[3][0] = dfdx_2(x_i_1, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[3][1] = ddfdxdy_2(x_i_1, y_i, x_i_1 - x_i, y_i_1 - y_i);
	F[3][2] = dfdx_2(x_i_1, y_i_1, x_i_1 - x_i, y_i_1 - y_i);
	F[3][3] = ddfdxdy_2(x_i_1, y_i_1, x_i_1 - x_i, y_i_1 - y_i);
}

double func(double x, double y)
{
	return sin(x) * cos(y);
}

double f(double x, double y)
{
	return sin(x) * cos(y);
}

double dfdx(double x, double y)
{
	(void)y;
	return cos(x) * cos(y);
}

double dfdy(double x, double y)
{
	(void)x;
	return -sin(x) * sin(y);
}

double ddfdxdy(double x, double y)
{
	(void)x;
	(void)y;
	return -cos(x) * sin(y);
}

double dfdx_2(double x, double y, double dx, double dy)
{
	(void)dy;
	double f2 = (f(x, y) - f(x - dx, y)) / dx;
	double f3 = (f(x + dx, y) - f(x, y)) / dx;
	return ((sign(f3) != sign(f2)) ? 0 : sign(f2)) * fmin(fabs(f3), fabs(f2));
}

double dfdy_2(double x, double y, double dx, double dy)
{
	(void)dx;
	double f2 = (f(x, y) - f(x, y - dy)) / dy;
	double f3 = (f(x, y + dy) - f(x, y)) / dy;
	return ((sign(f3) != sign(f2)) ? 0 : sign(f2)) * fmin(fabs(f3), fabs(f2));
}

double ddfdxdy_2(double x, double y, double dx, double dy)
{
	return (f(x + dx, y + dy) -  f(x, y + dy) - f(x + dx, y) + f(x, y)) / dx / dy;
}

void print_A(double A[4][4])
{
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		printf("%lf ", A[i][j]);
		printf("\n");
	}
}

double calculate(double x, double y, double F[4][4])
{
	double X = 1;
	double Y;
	double SUM = 0;
	for(int i = 0; i < 4; i++)
	{
		Y = 1;
		for(int j = 0; j < 4; j++)
		{
			SUM += F[i][j] * X * Y;
			Y *= y;
		}
		X *= x;
	}
	return SUM;
}

void mult_matr(double A[4][4], double B[4][4], double C[4][4])
{
	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 4; k++)
		C[i][j] += A[i][k] * B[k][j];
	}
}


void myGLWidget::paintGL()
{
	setProjectionMatrix();

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);

	glBegin(GL_QUADS);

double	x_min = -1, x_max = 1, y_min = -1, y_max = 1;
int n = 16;

if(variant == 0){
	int	x_n = 300, y_n = 300;

	glColor3d(1.0,0.0,0.0);

	for (int i = 0; i < x_n - 1; i++)
		for (int j = 0; j < y_n - 1; j++) {
			double	x, y, z;

			glColor3d(1.0 * (x_n - i) / x_n, 1.0 * j / y_n, 0.0);

			x = (x_max - x_min) * i / (x_n - 1) + x_min;
			y = (y_max - y_min) * j / (y_n - 1) + y_min;
			z = func(x, y);
			glVertex3d(x, y, z);
			x = (x_max - x_min) * (i + 1) / (x_n - 1) + x_min;
			y = (y_max - y_min) * j / (y_n - 1) + y_min;
			z = func(x, y);
			glVertex3d(x, y, z);
			x = (x_max - x_min) * (i + 1) / (x_n - 1) + x_min;
			y = (y_max - y_min) * (j + 1) / (y_n - 1) + y_min;
			z = func(x, y);
			glVertex3d(x, y, z);
			x = (x_max - x_min) * i / (x_n - 1) + x_min;
			y = (y_max - y_min) * (j + 1) / (y_n - 1) + y_min;
			z = func(x, y);
			glVertex3d(x, y, z);
		}
}

if(variant == 1)
{
	glColor3d(1.0,0.0,0.0);

	int i, j, k, l, M;
	double dx, dy, ddx, ddy, x, y, x_i, y_j;
	double A[4][4], B[4][4], F[4][4], C[4][4];
	dx = (x_max - x_min) / point_number; // шаг по х
	dy = (y_max - y_min) / point_number; // шаг по у
	A_inversed(dx, A); // вычислили А
	A_inversed_T(dy, B); // вычислили В
	for (i = 0; i < point_number; i++)
	{
		x_i = x_min + i * dx;
		for (j = 0; j < point_number; j++)
		{
			
			y_j = y_min + j * dy;
			F_matrix(x_i, x_i + dx, y_j, y_j + dy, F); // вычислили F
			mult_matr(A, F, C); // перемножили
			mult_matr(C, B, F); // перемножили
			M = MAX(n / point_number, 5); 
			ddx = dx / M;
			ddy = dy / M;
			for(k = 0; k < M; k++)
			{
				for(l = 0; l < M; l++)
				{
					x = x_i + k * ddx;
					y = y_j + l * ddy;
					glColor3d(1.0 * (x - x_min) / (x_max - x_min), 0.5 * (y - y_min) / (y_max - y_min),  0.1);
					// поставили цвет
					glVertex3d(x, y, calculate(k * ddx, l * ddy, F));
					//проставляем точки квадрата
					x = x_i + (k + 1) * ddx;
					y = y_j + l * ddy;
					glVertex3d(x, y, calculate((k + 1) * ddx, l * ddy, F));
					x = x_i + (k + 1) * ddx;
					y = y_j + (l + 1) * ddy;
					glVertex3d(x, y, calculate((k + 1) * ddx, (l + 1) * ddy, F));
					x = x_i + k * ddx;
					y = y_j + (l + 1) * ddy;
					glVertex3d(x, y, calculate(k * ddx, (l + 1) * ddy, F));
				}
			}
		}
	}


}

if(variant == 2) // смотри комментарии; отличается только заполнение матрицы F
{
	glColor3d(1.0,0.0,0.0);

	int i, j, k, l, M;
	double dx, dy, ddx, ddy, x, y, x_i, y_j;
	double A[4][4], B[4][4], F[4][4], C[4][4];
	dx = (x_max - x_min) / point_number;
	dy = (y_max - y_min) / point_number;
	A_inversed(dx, A);// problem
	A_inversed_T(dy, B); //problem
	for (i = 0; i < point_number; i++)
	{
		x_i = x_min + i * dx;
		for (j = 0; j < point_number; j++)
		{
			
			y_j = y_min + j * dy;
			F_matrix_2(x_i, x_i + dx, y_j, y_j + dy, F);
			mult_matr(A, F, C);
			mult_matr(C, B, F);
			M = MAX(n / point_number, 5);
			ddx = dx / M;
			ddy = dy / M;
			for(k = 0; k < M; k++)
			{
				for(l = 0; l < M; l++)
				{
					x = x_i + k * ddx;
					y = y_j + l * ddy;
					glColor3d(0.3 * (x - x_min) / (x_max - x_min), 0.5 * (y - y_min) / (y_max - y_min),  0.5);
					glVertex3d(x, y, calculate(k * ddx, l * ddy, F));
					x = x_i + (k + 1) * ddx;
					y = y_j + l * ddy;
					glVertex3d(x, y, calculate((k + 1) * ddx, l * ddy, F));
					x = x_i + (k + 1) * ddx;
					y = y_j + (l + 1) * ddy;
					glVertex3d(x, y, calculate((k + 1) * ddx, (l + 1) * ddy, F));
					x = x_i + k * ddx;
					y = y_j + (l + 1) * ddy;
					glVertex3d(x, y, calculate(k * ddx, (l + 1) * ddy, F));
				}
			}
		}
	}


}

	glEnd();

	glDisable(GL_DEPTH_TEST);
}

void myGLWidget::resizeGL(int nWidth, int nHeight)
{
	glViewport(0, 0, nWidth, nHeight);
	aspect = 1.0 * nWidth / nHeight;
	updateGL();
}

#define ANGLE_DIFF	(5)
#define POSITION_DIFF	(0.1)

void myGLWidget::keyPressEvent(QKeyEvent* e)
{
	switch (e->key()) {
	case Qt::Key_C:
		setDefaultCamera();
		break;
	case Qt::Key_Up:
		angle_v = MIN(angle_v + 5.0, 80);
		break;
	case Qt::Key_Down:
		angle_v = MAX(angle_v - 5.0, -80);
		break;
	case Qt::Key_Left:
		angle_h -= 5.0;
		break;
	case Qt::Key_Right:
		angle_h += 5.0;
		break;
	case Qt::Key_Plus:
		camera_p = MAX(camera_p - 0.1, 7);
		break;
	case Qt::Key_P: // увеличнение количества точек интерполяции
		point_number += 2;
		break;
	case Qt::Key_M: // уменьшение точек интерполяции
		point_number -= 2;
		break;
	case Qt::Key_F1: // переключение режима
		variant = 0;
		break;
	case Qt::Key_F2:
		variant = 1;
		break;
	case Qt::Key_F3:
		variant = 2;
		break;
	case Qt::Key_Minus:
		camera_p += 0.1;
		break;
	}

	updateGL();
}

void myGLWidget::setProjectionMatrix()
{
	GLfloat view[16] = {0}, projection[16] = {0}, tmp[16] = {0};

	static GLfloat	near = 5, top = 2, bottom, left, right;

	bottom = -top;
	right = top * aspect;
	left = -right;

	projection[0] = 2 * near / (right - left);
	projection[2] = (right + left) / (right - left);
	projection[5] = 2 * near / (top - bottom);
	projection[6] = (top + bottom) / (top - bottom);
	projection[10] = - 1;
	projection[11] = - 2 * near;
	projection[14] = -1;

	GLfloat	cam_x, cam_y, cam_z;
	cam_x = 0;
	cam_y = 0;
	cam_z = camera_p;

	view[0] = 1;
	view[6] = -1;
	view[9] = 1;
	view[15] = 1;

	view[12] = -cam_x;
	view[13] = -cam_y;
	view[14] = -cam_z;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glRotated(angle_h, 0, 0, 1);
	glGetFloatv(GL_PROJECTION_MATRIX, tmp);

	glLoadTransposeMatrixf(projection);
	glMultMatrixf(view);

	glRotated(angle_h, 0, 0, 1);
	glRotated(angle_v, tmp[0], tmp[4], tmp[8]);
}

void myGLWidget::setDefaultCamera()
{
	camera_p = 7;
	angle_h = 45;
	angle_v = 20;
	variant = 0;// вид интерполяции 0 - оригинал, 1 - кубический Эрмита 2 - кубический  с граничными узлами
	point_number = 4; 
	// количество узлов интерполяции на каждой стороне, 
	//если н = 4, то представляем квадрат 4х4
	aspect = 1.0 * width() / height();
}
