#include <QtGui/QMouseEvent>
#include <QtOpenGL>
#include "GLWidget.hpp"
#include <list>
#include "CParticle.hpp"
#include <iostream>
#include "GL/glu.h"

GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent) {
    setMouseTracking(true);
}

void GLWidget::initializeGL() {
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0, 0, 0, 0);
}

void GLWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    std::cout << "(w, h)= (" << w << ", " << h << ")\n";
}

void GLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(0,0,1);
    
    
    for (std::list<CParticle>::iterator p = vorticesList.begin();
				p!=vorticesList.end(); ++p) {
			std::cout << "Draw\n"; 
			int x=p->get_x();
			int y=p->get_y();
			glBegin(GL_POLYGON);
			glVertex2f(x+0,y+0);
			glVertex2f(x+10,y+50);
			glVertex2f(x+50,y+10);
			glEnd();
					
		}
    
    

}

void GLWidget::mousePressEvent(QMouseEvent *event) {

}
void GLWidget::mouseMoveEvent(QMouseEvent *event) {
    //printf("%d, %d\n", event->x(), event->y());
}

void GLWidget::keyPressEvent(QKeyEvent* event) {
    switch(event->key()) {
    case Qt::Key_Escape:
        close();
        break;
    default:
        event->ignore();
        break;
    }
}

void GLWidget::bindData(std::list<CParticle>* vorticesList_)
{
	vorticesList=*vorticesList_;
}



void GLWidget::draw()
{
		updateGL();
	
}
