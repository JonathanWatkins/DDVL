#ifndef GLWIDGET_HPP
#define GLWIDGET_HPP

#include <QtOpenGL/QGLWidget>
#include <list>
#include "CParticle.hpp"


class GLWidget : public QGLWidget {

    Q_OBJECT // must include this if you use Qt signals/slots

private:
		std::list<CParticle> vorticesList;

		
public:
    GLWidget(QWidget *parent = NULL);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    
public:
		void bindData(std::list<CParticle>* vorticesList_);
		void draw();

};

#endif  /* _GLWIDGET_H */
