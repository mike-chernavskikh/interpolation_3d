#include <qapplication.h>
#include <qgl.h>

#include "glwidget.h"

int main(int argc, char** argv)
{
	QApplication app(argc, argv);
	myGLWidget myglw;

	app.setActiveWindow(&myglw);
	myglw.show();
	return app.exec();
}
