#ifndef _my_widget
#define _my_widget

#include <qgl.h>

#include <QKeyEvent>

#include <qnamespace.h>

class myGLWidget:public QGLWidget
{
  Q_OBJECT

  protected:
	virtual void paintGL();
	virtual void initializeGL();
	virtual void resizeGL(int nWidth, int nHeight);
	virtual void keyPressEvent(QKeyEvent* e);

	void setProjectionMatrix();
	void setDefaultCamera();

	float	angle_h, angle_v;
	float	camera_p;
	float	aspect;
	int variant = 1;
	int point_number;
};

#endif
