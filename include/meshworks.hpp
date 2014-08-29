#ifndef _MESHWORKS_H_
	#define _MESHWORKS_H_
//#include <fstream>

using namespace std;

#include "myClass.h"
#include <list>
#include <vector>

#if defined (__WINDOWS__)
	#include <windows.h>
	#include "SDL.h"
	#include "SDL_gfxPrimitives.h"
	#include "SDL_ttf.h"
		#include <SDL_opengl.h>

GLuint GLpinsList;
GLuint GLvortex;
GLfloat drawUnit;
					
double GLxoffset;   
double GLyoffset;   
double GLzoffset;	
      
      
GLfloat drawnVortexSize;

void makePinsList( list<CVortex> pinsList) {
	GLpinsList = glGenLists(1);
    
      
		glNewList( GLpinsList, GL_COMPILE );
			glColor3f(.9f,0.9f,0.9f);
		for (list<CVortex>::iterator p = pinsList.begin();
				p != pinsList.end(); ++p) {
					 //if (p->get_x() < sourceWidth || p->get_x() > sourceWidth+channelWidth )
						//continue;
				
						double x = GLxoffset+drawUnit*p->get_x()/a0;
						double y2 = GLyoffset+drawUnit*p->get_y()/a0;
						double z2 =GLzoffset;
					
		
						
						GLUquadricObj *quadObj = gluNewQuadric();
		
						glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
						
						glTranslatef(x,y2,z2); 
										
						gluSphere(quadObj,drawnVortexSize/2.0 ,20 , 40); 
						gluDeleteQuadric(quadObj); 
						
						
		}
	
		glEndList();
	
	}

bool initialiseSDLandGL (SDL_Surface **screen,	TTF_Font **font,	TTF_Font **graphfont,	SDL_Color text_color, list<CVortex> pinsList) {
	
	

  
	//initialise Screen
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) return 1;
   
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
	
    //if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, DEPTH, SDL_FULLSCREEN|SDL_HWSURFACE)))
    if (!(*screen = SDL_SetVideoMode(WIDTH, HEIGHT, 0, SDL_HWSURFACE | SDL_GL_DOUBLEBUFFER | SDL_OPENGL)))
    {
        SDL_Quit();
        return 1;
    }
    
    glEnable(GL_MULTISAMPLE);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST );
 
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
    
    glClearColor(255, 255, 255, 0);
 
    glViewport(0, 0, WIDTH, HEIGHT);
 
    glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glTranslatef(-0.9f,-0.9f,0); 
   
 
    glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
 
    //glEnable(GL_TEXTURE_2D);
			
			glClearDepth(1.0f);
			
		 glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
 
   	
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);          // Really Nice Perspective Calculations


    
	glEnable(GL_LIGHTING);
    	
		glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
   glEnable ( GL_COLOR_MATERIAL ) ;
    
     glEnable(GL_LIGHT0);
     //GLfloat lightpos0[] = {-.5, .5, -1., 0.};
		GLfloat lightpos0[] = {1, 1, -1., 0.};
  glLightfv(GL_LIGHT0, GL_POSITION, lightpos0);
    
    //glEnable(GL_LIGHT1);
     //GLfloat lightpos1[] = {-.5, .5, -1., 0.};
  //glLightfv(GL_LIGHT1, GL_POSITION, lightpos1);
    
    
    
	SDL_WM_SetCaption( "meshworks (Tube Edition) 3.0", 0 );
	  if (runtype==0 || runtype==2) {
	// Initialize SDL_ttf library
   if (TTF_Init() != 0)
   {
      cerr << "TTF_Init() Failed: " << TTF_GetError() << endl;
      SDL_Quit();
      return 1;
   }

   // Load a font
  
   *font = TTF_OpenFont("FreeSans.ttf", 20);
   if (*font == NULL)
   {
      cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << endl;
      TTF_Quit();
      SDL_Quit();
      return 1;
   }
	
	*graphfont = TTF_OpenFont("FreeSans.ttf", 14);
   if (*font == NULL)
   {
      cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << endl;
      TTF_Quit();
      SDL_Quit();
      return 1;
   }

	}
		
		drawUnit= /*2.0**/(GLfloat)(1.5*a0/(sourceWidth+channelWidth+sinkWidth)); // lengths in units of a0 needs to be scaled by this for drawing
					
		GLxoffset =0.15;   
		GLyoffset =0;   
		GLzoffset = 0;	
      
      
		drawnVortexSize = (GLfloat)drawUnit*vortexSize/a0;
	
		makePinsList(pinsList);
			
	
	
	
	
	return 1;
}



bool initialiseSDLandGLTube (SDL_Surface **screen,	TTF_Font **font,	TTF_Font **graphfont,	SDL_Color text_color) {
	
	GLxoffset =.1;   
	drawUnit= (GLfloat)(1.5/((sourceWidth+channelWidth+sinkWidth)/a0)); // lengths in units of a0 needs to be scaled by this for drawing
		

  drawnVortexSize = (GLfloat)drawUnit*vortexSize/a0;
	
	//initialise Screen
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) return 1;
   
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
	
    //if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, DEPTH, SDL_FULLSCREEN|SDL_HWSURFACE)))
    if (!(*screen = SDL_SetVideoMode(WIDTH, HEIGHT, 0, SDL_HWSURFACE | SDL_GL_DOUBLEBUFFER | SDL_OPENGL)))
    {
        SDL_Quit();
        return 1;
    }
    
    glEnable(GL_MULTISAMPLE);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST );
 
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
    
    glClearColor(255, 255, 255, 0);
 
    glViewport(0, 0, WIDTH, HEIGHT);
 
    glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glTranslatef(-0.9f,-0.9f,0); 
   
 
    glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
 
    //glEnable(GL_TEXTURE_2D);
			
			glClearDepth(1.0f);
			
		 glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
 
   	
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);          // Really Nice Perspective Calculations


    
	glEnable(GL_LIGHTING);
    	
		glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
   glEnable ( GL_COLOR_MATERIAL ) ;
    
     glEnable(GL_LIGHT0);
     //GLfloat lightpos0[] = {-.5, .5, -1., 0.};
		GLfloat lightpos0[] = {1, 1, -1., 0.};
  glLightfv(GL_LIGHT0, GL_POSITION, lightpos0);
    
    //glEnable(GL_LIGHT1);
     //GLfloat lightpos1[] = {-.5, .5, -1., 0.};
  //glLightfv(GL_LIGHT1, GL_POSITION, lightpos1);
    
    
    
	SDL_WM_SetCaption( "meshworks (Tube Edition) 3.0", 0 );
	  if (runtype==0 || runtype==2) {
	// Initialize SDL_ttf library
   if (TTF_Init() != 0)
   {
      cerr << "TTF_Init() Failed: " << TTF_GetError() << endl;
      SDL_Quit();
      return 1;
   }

   // Load a font
  
   *font = TTF_OpenFont("FreeSans.ttf", 20);
   if (*font == NULL)
   {
      cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << endl;
      TTF_Quit();
      SDL_Quit();
      return 1;
   }
	
	*graphfont = TTF_OpenFont("FreeSans.ttf", 14);
   if (*font == NULL)
   {
      cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << endl;
      TTF_Quit();
      SDL_Quit();
      return 1;
   }

	}
	
		
		GLvortex = glGenLists(1);
    
      
		glNewList( GLvortex, GL_COMPILE );
		
						
		GLUquadricObj *quadObj = gluNewQuadric();
		gluSphere(quadObj,drawnVortexSize ,20 , 40); 
		
		gluDeleteQuadric(quadObj); 
						
						
		
		glEndList();
		
	
	
	return 1;
}

void currentRotation(GLfloat xRot, GLfloat yRot) {
			
			glMatrixMode(GL_MODELVIEW);
		  
			glLoadIdentity(); // Reset The ModelView Matrix
			glTranslatef(0.2f,-0.7f,0.0f); // use this when displaying plane aswell
			//glTranslatef(0.2f,0.0f,0.0f);
			glRotatef(yRot,1.0f,0.0f,0.0f); //rotate our camera on the x-axis (left and right)
			glRotatef(xRot,0.0f,1.0f,0.0f); //rotate our camera on the x-axis (left and right)
		
			//glRotatef(30,0.0f,1.0f,0.0f); //rotate our camera on the y-axis (up and down)
			//glRotatef(30,0.0f,0.0f,1.0f); //rotate our camera on the y-axis (up and down)
			
		
			 //translate the screen to the position of our camera

}

void glEnable2D()
{
	int vPort[4];
  
	glGetIntegerv(GL_VIEWPORT, vPort);
  
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
  
	glOrtho(0, vPort[2], 0, vPort[3], -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}

void glDisable2D()
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();   
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();	
}



float GLscaleY( double y) {

	return (float)(.4+y/(sourceHeight));
}

float GLscaleX( double x) {

	return (float)(.4+x/(channelWidth+sourceWidth+sinkWidth));
}


void drawSystem(list<CVortex> &delVortexList, list<CVortex> &pinsList, GLfloat &xRot, GLfloat &yRot, list<CDelLine> &delLinesList) {
    
    drawUnit=drawUnit*zoom;
    
    // with drawn source && sink
    //double GLxoffset =.15;   
		
		
		//glRotatef(90,0.0f,1.0f,0.0f); 
		//GLfloat num=t*0.01;
		
		// with drawn source && sink
		//GLfloat drawUnit= (GLfloat)(1.5/((sourceWidth+channelWidth+sinkWidth)/a0)); // lengths in units of a0 needs to be scaled by this for drawing
		
     //draw Delaunay points
    	
		 // draw remove vortices boxes
		if(DEBUG) {
		
		glColor3f(0.1f,0.1f,0.1f);
		
		double rx0 = GLxoffset+drawUnit*removechannelx0/a0;
		double rx1 = GLxoffset+drawUnit*removechannelx1/a0;
		double ry0 = GLyoffset+drawUnit*removetopchannely/a0;
		double ry1 = GLyoffset+drawUnit*removebottomchannely/a0;
		
		glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
		glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(rx0, ry0, 0.0f);              // Top Left
        glVertex3f(rx1, ry0, 0.0f);              // Top Right
        glVertex3f(rx1, ry1, 0.0f);              // Bottom Right
        glVertex3f(rx0, ry1, 0.0f);              // Bottom Left
    glEnd();  
	
	
		rx0 = GLxoffset+drawUnit*removesourcex/a0;
		rx1 = GLxoffset+drawUnit*(removesourcex+0.1*a0)/a0;
		ry0 = GLyoffset+drawUnit*removesourcey0/a0;
		ry1 = GLyoffset+drawUnit*removesourcey1/a0;
		
		glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
		glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(rx0, ry0, 0.0f);              // Top Left
        glVertex3f(rx1, ry0, 0.0f);              // Top Right
        glVertex3f(rx1, ry1, 0.0f);              // Bottom Right
        glVertex3f(rx0, ry1, 0.0f);              // Bottom Left
    glEnd();  
		rx0 = GLxoffset+drawUnit*(removesinkx-0.1*a0)/a0;
		rx1 = GLxoffset+drawUnit*removesinkx/a0;
		ry0 = GLyoffset+drawUnit*removesourcey0/a0;
		ry1 = GLyoffset+drawUnit*removesourcey1/a0;
		glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
		glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(rx0, ry0, 0.0f);              // Top Left
        glVertex3f(rx1, ry0, 0.0f);              // Top Right
        glVertex3f(rx1, ry1, 0.0f);              // Bottom Right
        glVertex3f(rx0, ry1, 0.0f);              // Bottom Left
    glEnd();  
		}
		
		
		//// Hack for pinned edges //////////////////
	
		double ax0 = GLxoffset+drawUnit*sourceWidth/a0;
		double ax1 = GLxoffset+drawUnit*(sourceWidth+channelWidth)/a0;
		double ay0 = GLyoffset+drawUnit*(extraPins-0.1*a0)/a0;
		double ay1 = GLyoffset+drawUnit*(extraPins)/a0;
		glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
		glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(ax0, ay0, 0.0f);              // Top Left
        glVertex3f(ax1, ay0, 0.0f);              // Top Right
        glVertex3f(ax1, ay1, 0.0f);              // Bottom Right
        glVertex3f(ax0, ay1, 0.0f);              // Bottom Left
    glEnd();  
		
		ax0 = GLxoffset+drawUnit*sourceWidth/a0;
		ax1 = GLxoffset+drawUnit*(sourceWidth+channelWidth)/a0;
		ay0 = GLyoffset+drawUnit*(channelHeight-extraPins)/a0;
		ay1 = GLyoffset+drawUnit*(channelHeight-extraPins+0.1*a0)/a0;
		
		glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
		glBegin(GL_QUADS);                      // Draw A Quad
        glVertex3f(ax0, ay0, 0.0f);              // Top Left
        glVertex3f(ax1, ay0, 0.0f);              // Top Right
        glVertex3f(ax1, ay1, 0.0f);              // Bottom Right
        glVertex3f(ax0, ay1, 0.0f);              // Bottom Left
    glEnd();  
		
		
		// end hack ////////////////////////
		
		
		
			/*GLuint vortex = glGenLists(1);
      GLUquadricObj *quadObj = gluNewQuadric();
      
			glNewList( vortex, GL_COMPILE );
				
				gluSphere(quadObj,drawnVortexSize ,20 , 40);
			glEndList();
			*/

      // planet1
      /*planet1 = gl.glGenLists(1);
      if( planet1 == 0 ) {
	System.out.println("GenList Error");
      } else {
	System.out.println("GenList planet1:"+planet1);//ddd
	gl.glNewList( planet1, GL.GL_COMPILE ); {
	  gl.glTranslatef ( 3.0f, 0.0f, 0.0f); 
	  gl.glColor4f( 0f, 1f, 0f, 1f );
	  //setSomeBlueMaterial( gl );
	  glu.gluSphere( qobj0, 0.3f, 20, 20);is 
	} gl.glEndList();
      }
		 */
		   
     for (list<CVortex>::iterator p = delVortexList.begin();
				p != delVortexList.end(); ++p) {
			  if (true==p->get_ghost())
					continue;
			  
			  //if (p->get_x() < sourceWidth || p->get_x() > sourceWidth+channelWidth )
				//	continue;
				if (p->get_coordNum()!=5 && p->get_coordNum()!=7 && p->get_coordNum()!=6)	
					continue;
				
				if (false==triangulateInRead && p->get_coordNum()==6 && runtype!=0)
					continue;
							
				if (true==triangulateInRead &&	p->get_coordNum()==6)
					continue;
				
			  		double x2 = GLxoffset+drawUnit*p->get_x()/a0;
				    
						//double y2= p->get_y()/(channelHeight*1.5);
						double y2= GLyoffset+drawUnit*p->get_y()/a0;
						
            double z2=GLzoffset;
					
					
						//currentRotation(xRot,yRot);
			
			      if (p->get_coordNum()==6) {
								//white sphere
						
								glColor3f(1.0f,1.0f,1.0f);
						}
						else if (p->get_coordNum()==5) {
								
								glColor3f(1.0f,0.0f,0.0f);
						}
						else if (p->get_coordNum()==7) {
								// blue sphere
			  	
								glColor3f(.0f,0.0f,1.0f);
						}
						else {
							//grey sphere
				
							glColor3f(0.0f,0.0f,0.0f);
						}
							
						
						GLUquadricObj *quadObj = gluNewQuadric();
						
						glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
				
						glTranslatef(x2,y2,z2); 
						
						//glCallList(vortex);
			
						gluSphere(quadObj,drawnVortexSize ,20 , 40);  
		
						gluDeleteQuadric(quadObj);
						
					
            	
      
			}
	
	glCallList(GLpinsList);
		
	/*for (list<CVortex>::iterator p = pinsList.begin();
				p != pinsList.end(); ++p) {
					 //if (p->get_x() < sourceWidth || p->get_x() > sourceWidth+channelWidth )
						//continue;
				
						double x = GLxoffset+drawUnit*p->get_x()/a0;
						double y2 = GLyoffset+drawUnit*p->get_y()/a0;
						double z2 =GLzoffset;
					
		
						/*GLUquadricObj *quadObj = gluNewQuadric();
						glColor3f(.9f,0.9f,0.9f);
		
						glMatrixMode(GL_MODELVIEW);
						glLoadIdentity();
						
						glTranslatef(x,y2,z2); 
										
						gluSphere(quadObj,drawnVortexSize/2.0 ,20 , 40); 
						gluDeleteQuadric(quadObj); 
						*/
	//}
		
		  
	// draw nn lines
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
				
				
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
			
			 if ((p->get_x1() < sourceWidth && p->get_y1() < 0  &&  p->get_x2() > sourceWidth)
				|| (p->get_x2() < sourceWidth && p->get_y2() < 0  &&  p->get_x1() > sourceWidth)) 
						continue;
		   if ((p->get_x1() < sourceWidth && p->get_y1() > channelHeight+b0  &&  p->get_x2() > sourceWidth)
				|| (p->get_x2() < sourceWidth && p->get_y2() > channelHeight+b0  &&  p->get_x1() > sourceWidth)) 
						continue;
			 if ((p->get_x1() > sourceWidth+channelWidth && p->get_y1() < 0  &&  p->get_x2() < sourceWidth+channelWidth)
				|| (p->get_x2() > sourceWidth+channelWidth && p->get_y2() < 0  &&  p->get_x1() < sourceWidth+channelWidth)) 
						continue;
		   if ((p->get_x1() > sourceWidth+channelWidth && p->get_y1() > channelHeight+b0  &&  p->get_x2() < sourceWidth+channelWidth)
				|| (p->get_x2() > sourceWidth+channelWidth && p->get_y2() > channelHeight+b0  &&  p->get_x1() < sourceWidth+channelWidth)) 
						continue;
			 
			  //cout  << "line: " << scaleX(p->get_x1()) << ", " << scaleY(p->get_y1()) << " -->  " << scaleX(p->get_x2()) << ", " << scaleY(p->get_y2()) << endl;
				GLfloat x1 = GLxoffset+drawUnit*p->get_x1()/a0;
				 GLfloat x2 = GLxoffset+drawUnit*p->get_x2()/a0;
				 GLfloat y1 = GLyoffset+drawUnit*p->get_y1()/a0;
				 GLfloat y2 = GLyoffset+drawUnit*p->get_y2()/a0;
				 GLfloat z = GLzoffset;
				 
					glColor3f(0.9f,0.9f,0.9f);
					glLineWidth(2); 
					glBegin(GL_LINES); 
					glVertex3f( x1, y1,z); 
					glVertex3f(x2, y2,z); 
					glEnd(); 
				
		
		}
    
    /*if(t==numTimeSteps) {
			oss.str("");
			oss << "[" << sourceDensity << "](" << sourceWidth<< "x" << sourceHeight << ")(" <<  channelWidth << "x" << channelHeight << ")["<< sinkDensity << "]-" << At << "(" << t<< ").bmp";
			
			renderStr = oss.str();
		
			// convert str stream to char
			strcpy(renderChar,renderStr.c_str());
			SDL_SaveBMP(screen, renderChar);	
     }*/  	  
		  
		  
		  
}	

void drawSystemTube(list<CVortex> &delVortexList, list<CVortex> &pinsList, GLfloat &xRot, GLfloat &yRot, list<CDelLine> &delLinesList) {
    
    // cylinder
		
		currentRotation(xRot,yRot);
		
		GLUquadricObj *quadObj = gluNewQuadric();
		
		//funnel code
		
		glRotatef(90,0.0f,1.0f,0.0f); 
		//GLfloat num=t*0.01;
		
		
		
				
		int numsteps = (int)(funnelWidth/(2*a0));
		/*glTranslatef(0.0f,0.0f,0.4-drawUnit*2.0);
		
		for (int i=0; i<= numsteps-1;i++) {
			if (i!=0) glTranslatef(0.0f,0.0f,-drawUnit*2.0); 
		
			if (i%2==1) glColor3f(.3f,0.3f,0.3f);
			if (i%2==0) glColor3f(.7f,0.7f,0.7f);
			double yrange1= channelHeight+b0+fabs(2.0*(2*a0)*i*tan(funnelAngle)); 
			double yrange2= channelHeight+b0+fabs(2.0*(2*a0)*(i+1)*tan(funnelAngle)); 
			
			double radius1 = ((yrange1)/(2.0*pi))*tubeRadiusScale;	
			double radius2 = ((yrange2)/(2.0*pi))*tubeRadiusScale;	
			
			gluCylinder(quadObj, radius2, radius1, drawUnit*2.0 ,30, 40);
		}
		*/
		
		//glLoadIdentity(); // Reset The ModelView Matrix
		
		//cylinder code
		currentRotation(xRot,yRot);			
		glRotatef(90,0.0f,1.0f,0.0f); 
		
		glTranslatef(0.0f,0.0f,GLxoffset);
		numsteps = (int)((sourceWidth+sinkWidth+channelWidth)/a0);
		for (int i=1; i<= numsteps/2.0;i++) {
			if (i!=1) glTranslatef(0.0f,0.0f,drawUnit*2.0); 
		
			if (i%2==1) glColor3f(.3f,0.3f,0.3f);
			if (i%2==0) glColor3f(.7f,0.7f,0.7f);
			
			gluCylinder(quadObj, 0.1, 0.1, drawUnit*2.0 ,30, 40);
		}
		
						
		
		
      
     //draw Delaunay points
    	
		 
		   
     for (list<CVortex>::iterator p = delVortexList.begin();
				p != delVortexList.end(); ++p) {
			  
			  
			   if (p->get_ghost()!=true) {
					
				
				    double x = GLxoffset+drawUnit*p->get_x()/a0;
				    double x2 = GLxoffset+drawUnit*p->get_x()/a0;
				    
           double y,y2,z,z2;
           if (p->get_x()>0) { 
						double radius = ((channelHeight+b0)/(2.0*pi))*tubeRadiusScale;
						//cout << "Radius : " << radius << endl;
            double theta = 2.0*pi*p->get_y()/(b0+channelHeight);
            if (theta<0) theta=theta+2*pi;
            if (theta>2*pi) theta=theta-2*pi;
            
            y = radius*sin(theta);
            y2= drawUnit*p->get_y()/a0;
            z = radius*cos(theta);
            z2=0;
					}
					else {
						
						double yrange= channelHeight+b0+fabs(2.0*p->get_x()*tan(funnelAngle)); 
						double radius = ((yrange)/(2.0*pi))*tubeRadiusScale;			
						
						double theta = 2.0*pi*p->get_y()/yrange;
            if (theta<0) theta=theta+2.0*pi;
            if (theta>2*pi) theta=theta-2.0*pi;
            //cout << "Y Mapping : (y->theta) " << p->get_y()/yrange << " -> " << theta << endl;
            
            y = radius*sin(theta);
            y2= drawUnit*p->get_y()/a0;
            z = radius*cos(theta);
						z2=0;
						
					}
            
         
					//	 cout << "del Vortex: " << x << ", " << y << endl; 
			
			
			
            if (p->get_coordNum()==6) {
				
				//white sphere
				currentRotation(xRot,yRot);
						  
		
		
		
				glTranslatef(x,y,z); 
				//glTranslatef(0.0f,0.5f,0); 
				glColor3f(1.0f,1.0f,1.0f);
				glCallList(GLvortex);
				
	
				/*glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x2,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				glCallList(GLvortex);
				*/
			
				
					}
					else if (p->get_coordNum()==5) {
					// triangle
					
        currentRotation(xRot,yRot);   
           	  
		
						
	
				glTranslatef(x,y,z); 
				//glTranslatef(0.5f,0.5f,0); 
				glColor3f(1.0f,0.0f,0.0f);
				
				glCallList(GLvortex);
	        	
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x2,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				glCallList(GLvortex);
	  	
             
			}
			else if (p->get_coordNum()==7) {
			  // blue sphere
			  	
				currentRotation(xRot,yRot);
		
		
				glTranslatef(x,y,z); 
				//glTranslatef(0.5f,0.5f,0); 
				glColor3f(.0f,0.0f,1.0f);
				glCallList(GLvortex);
		  
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x2,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				glCallList(GLvortex);
	
			  
			}
			else {
				//grey sphere
				
				currentRotation(xRot,yRot);
		
		
		
				glTranslatef(x,y,z); 
				//glTranslatef(0.5f,0.5f,0); 
				glColor3f(.3f,0.3f,0.3f);
				glCallList(GLvortex);
	
				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x2,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				glCallList(GLvortex);
	
			}
                 
         
		 }
		 /*else {  //draw ghosts
					  double x = 0.4+0.8*3.0*2.0*pi*0.1*(p->get_x()/(channelWidth+sourceWidth+sinkWidth));
            double theta = 2.0*pi*p->get_y()/(b0+channelHeight);
            if (theta<0) theta=theta+2*pi;
            if (theta>2*pi) theta=theta-2*pi;
            
            double y = 0.1*sin(theta);
            double y2= p->get_y()/channelHeight;
            double z = 0.1*cos(theta);
            double z2 = 0;
            
            //currentRotation(xRot,yRot);
		
		GLUquadricObj *quadObj = gluNewQuadric();
		
		
		//glTranslatef(x,y,z); 
		//glTranslatef(0.5f,0.5f,0); 
		glColor3f(.9f,0.9f,0.9f);
		//gluSphere(quadObj,0.01 ,20 , 40); 
		
				 glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				gluSphere(quadObj,0.01 ,20 , 40);  
            
            
		 
	 }*/
		 
	 }     
	
	
	/*for (list<CVortex>::iterator p = pinsList.begin();
				p != pinsList.end(); ++p) {
		 double x = GLxoffset+drawUnit*p->get_x()/a0;  
       
           double y2;
           
					  y2= p->get_y()/(channelHeight*1.5);
          
						double z2 =0;
					
		
		GLUquadricObj *quadObj = gluNewQuadric();
		
		
		
		
		glColor3f(.3f,0.3f,0.3f);
	
		
				 glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();
				
				glTranslatef(x,y2,z2); 
				glTranslatef(0.0f,-0.25f,0); 
			
				gluSphere(quadObj,0.01 ,20 , 40);  
		
		
		
	}*/
		
		  
	// draw nn lines
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
				
				
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
			  //cout  << "line: " << scaleX(p->get_x1()) << ", " << scaleY(p->get_y1()) << " -->  " << scaleX(p->get_x2()) << ", " << scaleY(p->get_y2()) << endl;
				 if (p->get_y1() <= removetopchannely || p->get_y2() <= removetopchannely || p->get_y1() >=removebottomchannely || p->get_y2() >= removebottomchannely)
					continue;
				 
				 GLfloat x1 = GLxoffset+drawUnit*p->get_x1()/a0;
				 GLfloat x2 = GLxoffset+drawUnit*p->get_x2()/a0;
				 GLfloat y1 = drawUnit*p->get_y1()/a0-.25f;
				 GLfloat y2 = drawUnit*p->get_y2()/a0-.25f;
				 
				 
				 glColor3f(0.9f,0.9f,0.9f);
				glLineWidth(2); 
    glBegin(GL_LINES); 
    glVertex2f( x1, y1); 
    glVertex2f(x2, y2); 
    glEnd(); 
		
		
		}
    gluDeleteQuadric(quadObj);
		    	  
		  
}	



void doEvents(SDL_Event event, int &t) {
	
	while (SDL_PollEvent(&event))
			{
				// Check for the quit message
				if (event.type == SDL_QUIT)
				{
					// Quit the program
					t=numTimesteps;
					cerr << "SDL_QUIT event"<< endl;
					break;
				}
				if (event.type == SDL_KEYDOWN)
				{
					if (event.key.keysym.sym == SDLK_ESCAPE) {
						// Handle mouse clicks here.
						t=numTimesteps;
						//break;
					}
					if (event.key.keysym.sym == SDLK_v) {
						// Handle mouse clicks here.
						outputvelfield=true;
						
						//break;
					}
					
					if (event.key.keysym.sym == SDLK_w) {
						rollup=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_s) {
						rolldown=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_a) {
						rollleft=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_d) {
						rollright=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_p) {
						zoom=zoom+0.1;
						
						
					}
					if (event.key.keysym.sym == SDLK_o) {
						zoom=zoom-0.1;
						
						
					}
					if (event.key.keysym.sym == SDLK_SPACE) {
						if (pauseActive==true) pauseActive=false;
						else if (pauseActive==false) pauseActive=true;
						
					} 
					if (event.key.keysym.sym == SDLK_l) {
						
						if (pinLayer==0) {
							pinLayer=1;
						  cout << "PIN" << endl;
						}
						
					} 
					
					
					
					
					
				}
				if (event.type == SDL_KEYUP)
				{
					if (event.key.keysym.sym == SDLK_w) {
						rollup=false;
						
						
					}
					if (event.key.keysym.sym == SDLK_s) {
						rolldown=false;
					}
					if (event.key.keysym.sym == SDLK_a) {
						rollleft=false;
						
						
					}
					if (event.key.keysym.sym == SDLK_d) {
						rollright=false;
					}
					
					if (event.key.keysym.sym == SDLK_p) {
						zoom=zoom+0.1;
						
						
					}
					if (event.key.keysym.sym == SDLK_o) {
						zoom=zoom-0.1;
						
						
					}
					
					
				}
				if(event.type == SDL_MOUSEBUTTONDOWN) {
				
					if( event.button.button == SDL_BUTTON_LEFT ){
						
					}
				}
				if(event.type == SDL_MOUSEBUTTONUP) {
				
					if( event.button.button == SDL_BUTTON_LEFT ){
						
						
						
					}
				}
				
					
				
				
				
				
			
			
				
			}
			
}



void doEventsTube(SDL_Event event, int &t) {
	
	while (SDL_PollEvent(&event))
			{
				// Check for the quit message
				if (event.type == SDL_QUIT)
				{
					// Quit the program
					t=numTimesteps;
					cerr << "SDL_QUIT event"<< endl;
					break;
				}
				if (event.type == SDL_KEYDOWN)
				{
					if (event.key.keysym.sym == SDLK_ESCAPE) {
						// Handle mouse clicks here.
						t=numTimesteps;
						//break;
					}
					if (event.key.keysym.sym == SDLK_v) {
						// Handle mouse clicks here.
						outputvelfield=true;
						
						//break;
					}
					
					if (event.key.keysym.sym == SDLK_w) {
						rollup=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_s) {
						rolldown=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_a) {
						rollleft=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_d) {
						rollright=true;
						
						
					}
					if (event.key.keysym.sym == SDLK_SPACE) {
						if (pauseActive==true) pauseActive=false;
						else if (pauseActive==false) pauseActive=true;
						
					}
					
					
				}
				if (event.type == SDL_KEYUP)
				{
					if (event.key.keysym.sym == SDLK_w) {
						rollup=false;
						
						
					}
					if (event.key.keysym.sym == SDLK_s) {
						rolldown=false;
					}
					if (event.key.keysym.sym == SDLK_a) {
						rollleft=false;
						
						
					}
					if (event.key.keysym.sym == SDLK_d) {
						rollright=false;
					}
					
					
					
				}
				if(event.type == SDL_MOUSEBUTTONDOWN) {
				
					if( event.button.button == SDL_BUTTON_LEFT ){
						
					}
				}
				if(event.type == SDL_MOUSEBUTTONUP) {
				
					if( event.button.button == SDL_BUTTON_LEFT ){
						
						
						
					}
				}
				
					
				
				
				
				
			
			
				
			}
			
}




#endif

bool CVortexSort (CVortex i,CVortex j) { return (i.get_x()<j.get_x()); }





/*void wrapVortices(list<CVortex>& vorticesList, list<CVortex>& wrappedVorticesList, int t, double forceRange) {
	wrappedVorticesList=vorticesList;
  for (list<CVortex>::iterator p = vorticesList.begin();
			p!=vorticesList.end(); ++p ) {
		// wrap vortices on tube
		if (p->get_x() >=0) {
			if (p->get_y() <= forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelHeight+b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= channelHeight-forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelHeight-b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
		}
		
		
		
		
	}
}
*/

/*
void wrapVorticesWithTSClass(MyThreadSafeList<CVortex>& newVorticesList, list<CVortex>& wrappedVorticesList, int t, double forceRange) {
	//wrappedVorticesList=vorticesList;
  wrappedVorticesList.insert(wrappedVorticesList.end(),newVorticesList.begin(),newVorticesList.end());
			
  for (list<CVortex>::iterator p = newVorticesList.begin();
			p!=newVorticesList.end(); ++p ) {
		// wrap vortices on tube
		if (p->get_x() >=0) {
			if (p->get_y() <= forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()+channelHeight+b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
			else if (p->get_y() >= channelHeight-forceRange) {
				CVortex newVortex;
				newVortex = (*p);
				newVortex.set_pos(newVortex.get_x(),newVortex.get_y()-channelHeight-b0);
				newVortex.set_ghost();
				wrappedVorticesList.push_back(newVortex);
			}
		}
		
		
		
		
	}
}

*/

void pause(int seconds) {
	int startpause=time(0);
	int endpause=time(0);
	while (endpause-startpause<seconds) endpause= time(0);
	
}



	
int scaleX( double x) {

	return (int)((x*SDLscaling)+SDLxOffset);
}

int scaleY( double y) {

	return (int)((y*SDLscaling) + SDLyOffset);
}

inline double gaussian(double v1[], double v2[],double R) {
	double modsquared = (v1[0]- v2[0])*(v1[0]- v2[0]) + (v1[1]- v2[1])*(v1[1]- v2[1]);
	return exp( - modsquared/(R*R));	
	//return 1+ x +(x*x)/2.0 +(x*x*x)/6.0+(x*x*x*x)/24.0+(x*x*x*x*x)/120.0+(x*x*x*x*x*x)/720.0;
	//return (R*R)/modsquared;
	 
	//return exp(x);
	
}
/*
inline double energyForm(int form, double dist, double R ) {
	
	
	if (form==0) {  // vortex-vortex Gaussian	
			//force = -(2/(R*R))*gaussian(v1,v2,R);
			
	
	}
	else if (form==1) { // vortex-pin Gaussian
			//force=-(2/(R*R))*Ap*gaussian(v1,v2,R);
		
	}
	else if (form==2 || form ==3) {
		return (Phi*Phi/(2*pi*mu0*lambda*lambda))*boost::math::cyl_bessel_k(0,  dist/lambda);
	}
	
	return false;
}
*/

void findCell(double xval, double yval, int& cellx, int& celly, double forceRange) {
	cellx = (int)ceil(xval / forceRange);
	celly= (int)ceil(yval /forceRange);
	
}

void normaliseSource(list<CDelLine>& delLinesList, list<CVortex>& vorticesList, int t, ofstream& newVortexfile, double Phi, string &normaliseSourceStr) {
	
	static int lastadded=0;
	double aaverage=0;
	int numa=0;
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if ( midx> sourceWidth-binsize/2.0 && midx < sourceWidth+binsize/2.0 &&
				(midy>0 && midy<channelHeight))
		{
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
											
			aaverage=aaverage+linelength;
			numa++;	
		}			
	}
	aaverage=aaverage/(double)numa;
	
	double Beff=2*Phi/(sqrt((double)3)*aaverage*aaverage);		
	ostringstream oss;
	oss.str("");
	
	oss << "Bl(x="<<sourceWidth/a0<<"a0)=" << Beff << "T";
	
	
	
	
	// calculate zone densities
	int sourceCount=0;
	int sinkCount=0;
	
	// count densities
	for (list<CVortex>::iterator p = vorticesList.begin();
		p != vorticesList.end(); ++p) {
		if (p->get_x() < sourceWidth && (p->get_y()<0 || p->get_y()>channelHeight))
			      sourceCount++;
		if (p->get_x() > sourceWidth+channelWidth && 
			     (p->get_y()<0 || p->get_y()>channelHeight)) sinkCount++;
	}
	
	if (Beff<sourceBfield && t-lastadded>=10) { // add a source vortex - at most every 10 steps
		lastadded=t;
	newVortexfile << setw(10) << t;
		
		
			CVortex newVortex;
			
			double xval = (sourceWidth-a0)*(rand() % 1000)/1000.0;
			double yval = channelOffset*(rand() % 1000)/1000.0;
			xval=xval+a0/2.0;
			// choose top || bottom
			if (rand()%2==0)
			{
				// top
				yval=yval-channelOffset;
				
			}
			else {
				//bottom
				yval=yval+channelHeight;
			}
			
			
			
			newVortex.set_pos
					( xval,yval);
			
			vorticesList.push_back(newVortex);
			
			//output newVortex added data
		
			
		
		newVortexfile << setw(10) << 1 << endl;
		
	}
	
	
	
	
	
	//cout << "Source: " << sourceCount << " Sink: " << sinkCount << endl;
	//add sourceVortex
	
	/*if ( sourceCount < sourceDensity ) { // add new source vortices
		newVortexfile << setw(10) << t;
		int numadded=0;
		for (int i = 0; i<(sourceDensity-sourceCount);i++) {
			CVortex newVortex;
			
			double xval = (sourceWidth-a0)*(rand() % 1000)/1000.0;
			double yval = channelOffset*(rand() % 1000)/1000.0;
			xval=xval+a0/2.0;
			// choose top || bottom
			if (rand()%2==0)
			{
				// top
				yval=yval-channelOffset;
				
			}
			else {
				//bottom
				yval=yval+channelHeight;
			}
			
			
			
			newVortex.set_pos
					( xval,yval);
			
			vorticesList.push_back(newVortex);
			
			//output newVortex added data
			numadded++;
			
		}
		newVortexfile << setw(10) << numadded << endl;
		
	}	*/
	
	
	
	
	//int removecount=0;
	
	/*
	if ( sinkCount > sinkDensity ) { // remove sink vortices from very start && very end of channel
		
		list<CVortex>::iterator itSink = vorticesList.begin();
		while (sinkCount > sinkDensity) {
		//cerr << "Remove count: " << removecount << endl;
			if ( (itSink->get_x()>sourceWidth+channelWidth) && 
			     (itSink->get_y()<0 || itSink->get_y()>channelHeight) )		
			 {
			  cout << "Removed by Normalise at " << t << " pos (" <<  itSink->get_x() << ", " << itSink->get_y() << ") vel (" << itSink->get_velx() << ", " << itSink->get_vely() << ")" << endl;  
				itSink=vorticesList.erase(itSink);
				sinkCount--;
			}
			else {
				++itSink;
			}
		}
		
	
	} 
		
	
	
	// calculate sink a0
	
	
	// 
	
	 */	
	normaliseSourceStr = oss.str();	
	
}

void normaliseSink(list<CDelLine>& delLinesList, list<CVortex> &vorticesList, double Phi, string &returnStr, int t) {
	// This routine calculates the effect B field of the source (not in the wings)
	// It then removes sink vortices from the wings until the density is low enough.
	static int lastremove=0;
	
	
	double aaverage=0;
	int numa=0;
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if ( (midx > sourceWidth+channelWidth-binsize/2.0 && midx<sourceWidth+channelWidth+binsize/2.0) &&
				(midy>0 && midy<channelHeight)) {
			//continue;
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
											
			aaverage=aaverage+linelength;
			numa++;	
		}				
	}
	aaverage=aaverage/(double)numa;
	
	double Beff=2*Phi/(sqrt((double)3)*aaverage*aaverage);		
	ostringstream oss;
	oss.str("");
	
	oss << "Br(x="<<(sourceWidth+channelWidth)/a0<<"a0)=" << Beff << "T";
	
	
	
	// calculate number of vortices in the wings
	int sinkCount=0;
	for (list<CVortex>::iterator p = vorticesList.begin();
		p != vorticesList.end(); ++p) {
		if (p->get_x() > sourceWidth+channelWidth && 
			     (p->get_y()<0 || p->get_y()>channelHeight)) sinkCount++;
	}
	
	
	if (Beff>sinkBfield && sinkCount!=0 && t-lastremove>=50  ) { // remove a sinkVortex
		bool removed=false;
		lastremove=t;
		list<CVortex>::iterator itSink = vorticesList.begin();
		while (false==removed) {
		//cerr << "Remove count: " << removecount << endl;
			if ( (itSink->get_x()>sourceWidth+channelWidth) && 
			     (itSink->get_y()<0 || itSink->get_y()>channelHeight) )		
			 {
			  cout << "Removed by NormaliseSink at " << t << " pos (" <<  itSink->get_x() << ", " << itSink->get_y() << ") vel (" << itSink->get_velx() << ", " << itSink->get_vely() << ")" << endl;  
				itSink=vorticesList.erase(itSink);
				removed=true;
				oss << " Sink Count: " << sinkCount-1;
			}
			itSink++;
		}
		
	
	}
	
		returnStr = oss.str();	
}



void normaliseDensitiesTube(list<CVortex>& vorticesList, int t, double Phi, ofstream& newVortexfile) {
	// calculate zone densities
	int sourceCount=0;
	int sinkCount=0;
	
	// count densities
	for (list<CVortex>::iterator p = vorticesList.begin();
		p != vorticesList.end(); ++p) {
		if (p->get_x() < -funnelWidth/2.0)
			      sourceCount++;
		if (p->get_x() > sourceWidth+channelWidth)
			      sinkCount++;
	}
	//cout << "Source: " << sourceCount << "of (" << sourceDensity << ") Sink: " << sinkCount << endl;
	//add sourceVortex
	
	if ( sourceCount < sourceDensity ) { // add new source vortices
		//cout << "Source: " << sourceCount << " density required: " << sourceDensity << endl;
	
		newVortexfile << setw(10) << t;
		int numadded=0;
		for (int i = 0; i<(sourceDensity-sourceCount);i++) {
			CVortex newVortex;
			
			double xval = (funnelWidth/2.0)*(rand() % 1000)/1000.0-funnelWidth; // generates xpos in first half of funnel
			double yrange= channelHeight+b0+2.0*fabs(xval)*tan(funnelAngle); 
			double yval = yrange*(rand() % 1000)/1000.0-fabs(xval)*tan(funnelAngle);
			
			//cout << "New funnel vortex added at pos: ( " << xval << " , " << yval << " )" << endl; 
			//cout << "yrange: " << yrange << endl;
			//cout << xval << ", " << yval << endl;
			// choose top or bottom
			/*if (rand()%2==0)
			{
				// top
				yval=yval-channelOffset;
				
			}
			else {
				//bottom
				yval=yval+channelHeight;
			}*/
			
			
			
			
			newVortex.set_pos
					( xval,yval);
			
			vorticesList.push_back(newVortex);
			
			//output newVortex added data
			numadded++;
			
		}
		newVortexfile << setw(10) << numadded << endl;
	}	
	//int removecount=0;
	if ( sinkCount > sinkDensity ) { // remove sink vortices from very start and very end of channel
		
		list<CVortex>::iterator itSink = vorticesList.begin();
		while (sinkCount > sinkDensity) {
			
		//cerr << "Remove count: " << removecount << endl;
		//cout << "Sink count: " << sinkCount << "of " << sinkDensity << endl;
			if ( (itSink->get_x()>sourceWidth+channelWidth)  )		
			 {
			  cout << "Removed by Normalise at " << t << " pos (" <<  itSink->get_x() << ", " << itSink->get_y() << ") vel (" << itSink->get_velx() << ", " << itSink->get_vely() << ")" << endl;  
				itSink=vorticesList.erase(itSink);
				sinkCount--;
			}
			else {
				++itSink;
			}
		}
		
	
	} 
		
	//cout << "end Normalise" << endl;
}


void vortexCounter(list<CVortex>& vorticesList, string &vortexCounterStr){
	
	ostringstream oss;
	oss.str("");
	
	int srcCount=0;
	int chnCount=0;
	int sinkCount=0;
	

	for (list<CVortex>::iterator p=vorticesList.begin();
			p!=vorticesList.end();++p) {
		if (p->get_x() <=sourceWidth) srcCount++;
		else if (p->get_x() >sourceWidth && p->get_x() < sourceWidth+channelWidth) chnCount++;
		else if (p->get_x() >=sourceWidth+channelWidth) sinkCount++;
		
		
	}
	oss << "Total Vortices: " << vorticesList.size() << " Src: "<< srcCount << " Chn: " << chnCount << " Sink: " << sinkCount<< " Sum: " << srcCount+chnCount+sinkCount;
	
	vortexCounterStr=oss.str();
	
}

void initialisePins(list<CVortex>& pinsList) {
		//double channelOffset=channelHeight/2.0;
		double locala0=a0;
		double localb0=b0;
	    double xPos;
        double yPos=-6*localb0-channelOffset;
	    
	    while (yPos<sourceHeight+6*localb0-channelOffset) {
	        xPos=-3*a0;
	        
			while (xPos<sourceWidth+channelWidth+sinkWidth+3*locala0) {
		
			
				CVortex newPin;
				newPin.set_pos(xPos,yPos+localb0/2.0);
				pinsList.push_back(newPin);
				
				newPin.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
				pinsList.push_back(newPin);
				
				
				
				xPos=xPos+locala0;
			
			}
			yPos=yPos+2*localb0;
		}
		
		//etch source, sink && channel
		bool removed;
		list<CVortex>::iterator p= pinsList.begin();
		
		while (p!=pinsList.end()) {
			
			removed=false;
			//etch source
			if (  p->get_x() > etchsourcex0 && p->get_x() < etchsourcex1
			  && p->get_y() > etchsourcey0 && p->get_y() < etchsourcey1 )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			else if (  p->get_x() > etchchannelx0 && p->get_x() < etchchannelx1
			  && p->get_y() > etchchannely0 && p->get_y() < etchchannely1 )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			else if (  p->get_x() > etchsinkx0 && p->get_x() < etchsinkx1
			  && p->get_y() > etchsinky0 && p->get_y() < etchsinky1 )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			
			
			if (removed==false) { ++p; }
		
				
		}
		
		cout << "Num of Pins: " << pinsList.size() << endl;

}


void initialisePinsTube(list<CVortex>& pinsList) {
		//double channelOffset=channelHeight/2.0;
		double locala0=a0;
		double localb0=b0;
	    double xPos;
        double yPos=-20*localb0-channelOffset;
	    
	    while (yPos<sourceHeight+20*localb0-channelOffset) {
	        xPos=-3*a0-funnelWidth;
	        
			while (xPos<sourceWidth+channelWidth+sinkWidth+3*locala0) {
		
			
				CVortex newPin;
				newPin.set_pos(xPos,yPos+localb0/2.0);
				pinsList.push_back(newPin);
				
				newPin.set_pos(xPos+locala0/2.0,yPos+3*localb0/2.0);
				pinsList.push_back(newPin);
				
				
				
				xPos=xPos+locala0;
			
			}
			yPos=yPos+2*localb0;
		}
		
		//etch source, sink and channel
		bool removed;
		list<CVortex>::iterator p= pinsList.begin();
		
		while (p!=pinsList.end()) {
			
			removed=false;
			//etch source
			/*
			if (  p->get_x() > etchsourcex0 and p->get_x() < etchsourcex1
			  and p->get_y() > etchsourcey0 and p->get_y() < etchsourcey1 )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			else if (  p->get_x() > etchchannelx0 and p->get_x() < etchchannelx1
			  //and p->get_y() > etchchannely0 and p->get_y() < etchchannely1
			   )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			else if (  p->get_x() > etchsinkx0 and p->get_x() < etchsinkx1
			 and p->get_y() > etchsinky0 and p->get_y() < etchsinky1 )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			*/
			if (  p->get_x() > etchsourcex0 && p->get_x() < etchsinkx1
			  )
			{
				p=pinsList.erase(p);
			
				removed=true;
			}
			
			if (removed==false) { ++p; }
		
				
		}
		
		cout << "Num of Pins: " << pinsList.size() << endl;

}



double calculatea0( double x ) {
		double funcp = 1.0 - (x-sourceWidth)/(channelWidth)*0.5;
		double funca0 =(0.00000000984)*pow((funcp),(-0.415)); 
	
		return funca0;
} 

bool leftslipplane(CCoord AP, CCoord AB) {
		
	if ((AB.get_x() * AP.get_y() - AB.get_y() * AP.get_x()) > 0 )	return true;
	else return false;
}

void loadVorticesFromFile(list<CVortex>& vorticesList, ifstream myfile){
double xval;
double yval;

if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      myfile >> xval;
      myfile >> yval;
      CVortex newVortex;
      newVortex.set_pos
					(xval,yval);
			vorticesList.push_back(newVortex);
    }
    myfile.close();
  }
	
			
		
	

}

void initialiseVortices(list<CVortex>& vorticesList, bool& file, double forceRange){
	
	ostringstream oss;
	char renderChar[100];
	string renderStr;
	oss.str("");
	oss << "startdataMesh//" << sourceBfield << "-" << sourceWidth/a0<< "x" << sourceHeight/b0 << "-" <<  channelWidth/a0 << "x" << channelHeight/b0 << "-"<< sinkBfield << ".txt";
	cout <<"oss: " << oss << endl;
	renderStr = oss.str();
	cout << renderStr << endl;
	strcpy(renderChar,renderStr.c_str());
			
  ifstream myfile (renderChar);
  cout << "initialise Vortices" << endl;
  cout << "sourceDensity: " << sourceDensity << endl;
  cout << "sinkDensity: " << sinkDensity << endl;
  
	if (myfile.is_open()) {
		cout << "Initial Vortex Positions From File" << endl;
		file = true;
		double xval;
		double yval;
    
    while ( myfile.good() )
    {
      myfile >> xval;
      myfile >> yval;
      // calculate cell containing vortex.
      int cellx,celly;
      findCell(xval,yval,cellx,celly, forceRange);
      
      CVortex newVortex;
      newVortex.set_pos
					(xval,yval);
			newVortex.set_cell(cellx,celly);
			vorticesList.push_back(newVortex);
			
			
    }
    myfile.close();

	}
	else {
		cout << "no start data" << endl;
		file = false;
	for (int i = 0; i<(sourceDensity);i++) {
			int cellx,celly;
			double xval,yval;
			
			xval = sourceWidth*(rand() % 1000)/1000.0;
			yval = sourceHeight*(rand() % 1000)/1000.0;
			yval = yval-channelOffset;
			//cout << "source vortex: " << xval << ", " << yval << endl;
			
			findCell(xval,yval,cellx,celly, forceRange);		
			
			CVortex newVortex;
				
			newVortex.set_pos(xval,yval);
			newVortex.set_cell(cellx,celly);
      		
			//newVortex.set_pos
			//		( ((sourceWidth/2.0)*10)/10.0
			//		,(sourceHeight*10)/10.0);
		
			vorticesList.push_back(newVortex);
		
	}
	for (int i = 0; i<(sinkDensity);i++) {
			CVortex newVortex;
			double xval= sinkWidth*(rand() % 1000)/1000.0+sourceWidth+channelWidth;
			double yval = sinkHeight*(rand() % 1000)/1000.0;
			yval = yval-channelOffset;
			//cout << "sink vortex: " << xval << ", " << yval << endl;
			newVortex.set_pos(
					xval,yval);
			vorticesList.push_back(newVortex);
			
		
	}	
	//cout << "sink done" << endl;
	
	
	
	for (int i = 0; i<(channelDensity);i++) {
			double xval = channelWidth*(rand() % 1000)/1000.0+sourceWidth;
			double yval = channelHeight*(rand() % 1000)/1000.0;
			//cout << "channel vortex: " << xval << ", " << yval << endl;
			CVortex newVortex;
			newVortex.set_pos(
					xval,yval);
			vorticesList.push_back(newVortex);
		
	}	
	
	
		
		
	}
   
  cout << "initialised" << endl;
}



void initialiseVorticesTube(list<CVortex>& vorticesList, bool& file, double forceRange){
	
	ostringstream oss;
	char renderChar[100];
	string renderStr;
	oss.str("");
	oss << "startdataMesh//TE-" << sourceBfield << "-" << sourceWidth/a0<< "x" << sourceHeight/b0 << "-" <<  channelWidth/a0 << "x" << channelHeight/b0 << "-"<< sinkBfield << ".txt";
	cout <<"oss: " << oss << endl;
	renderStr = oss.str();
	cout << renderStr << endl;
	strcpy(renderChar,renderStr.c_str());
			
  ifstream myfile (renderChar);
  cout << "initialise Vortices" << endl;
  cout << "sourceDensity: " << sourceDensity << endl;
  cout << "sinkDensity: " << sinkDensity << endl;
  
	if (myfile.is_open()) {
		cout << "Initial Vortex Positions From File" << endl;
		file = true;
		double xval;
		double yval;
    
    while ( myfile.good() )
    {
      myfile >> xval;
      myfile >> yval;
      // calculate cell containing vortex.
      
      CVortex newVortex;
      newVortex.set_pos
					(xval,yval);
			vorticesList.push_back(newVortex);
			
			
    }
    myfile.close();

	}
	else {
		cout << "no start data" << endl;
		file = false;
	for (int i = 0; i<(sourceDensity);i++) {
			int cellx,celly;
			double xval,yval;
			
			xval = sourceWidth*(rand() % 1000)/1000.0;
			yval = sourceHeight*(rand() % 1000)/1000.0;
			yval = yval-channelOffset;
			//cout << "source vortex: " << xval << ", " << yval << endl;
			
			
			CVortex newVortex;
				
			newVortex.set_pos(xval,yval);
					
			//newVortex.set_pos
			//		( ((sourceWidth/2.0)*10)/10.0
			//		,(sourceHeight*10)/10.0);
		
			vorticesList.push_back(newVortex);
		
	}
	for (int i = 0; i<(sinkDensity);i++) {
			CVortex newVortex;
			double xval= sinkWidth*(rand() % 1000)/1000.0+sourceWidth+channelWidth;
			double yval = sinkHeight*(rand() % 1000)/1000.0;
			yval = yval-channelOffset;
			//cout << "sink vortex: " << xval << ", " << yval << endl;
			newVortex.set_pos(
					xval,yval);
			vorticesList.push_back(newVortex);
			
		
	}	
	//cout << "sink done" << endl;
	
	
	
	for (int i = 0; i<(channelDensity);i++) {
			double xval = channelWidth*(rand() % 1000)/1000.0+sourceWidth;
			double yval = channelHeight*(rand() % 1000)/1000.0;
			//cout << "channel vortex: " << xval << ", " << yval << endl;
			CVortex newVortex;
			newVortex.set_pos(
					xval,yval);
			vorticesList.push_back(newVortex);
		
	}	
	
	
		
		
	}
   
  cout << "initialised" << endl;
}




void initialiseStaticVortices(list<CVortex>& vorticesList, bool& file, double slipPoint){
	
	
	
	CCoord A;
	A.set_coords(slipPoint,0);
	
	CCoord AB;
	AB.set_coords(calculatea0(slipPoint)/2,b0);
	
	double newa0=calculatea0(sourceWidth);
	double lasthighestx=0;
	for (double x= sourceWidth;x<=sourceWidth+channelWidth; x=x+newa0 ){
		
		for (double y=b0/2;y<=channelHeight-b0/2;y=y+2*b0) {
			CCoord AP;
			AP.set_coords(x-A.get_x(),y-A.get_y());
			
			if (leftslipplane(AP,AB)==true) {
			
				CVortex newVortex;
				
				newVortex.set_pos(x,y);
				vorticesList.push_back(newVortex);
				
				newVortex.set_pos(x+newa0/2.0,y+b0);
				vorticesList.push_back(newVortex);

			}
			
			newa0=calculatea0(x); //*(x-sourceWidth)/(sourceWidth+channelWidth);
	
		}
		
	}
	
	for (list<CVortex>::iterator p=vorticesList.begin();
			p!= vorticesList.end(); ++p ) {

		if (p->get_y()==b0/2 && p->get_x() >lasthighestx ) lasthighestx=p->get_x();
		
	}
	//cout << "lasthighestx" << lasthighestx << endl;
	
	
	// This code uses the cross product of two vectors to find out which of the
	// slip plane a given point in this channel is.
	// The vector uses CCoord type && defines the AB vector at 60 degrees across the channel
	// The vector AP, where P is the point, is calculated as follows (using the iterator p)
	// P(x,y) = (p->get_x()-Ax,p->get_y()-Ay)
	
	
			
	newa0=calculatea0(sourceWidth);
	double newb0=b0*1.1;
	//=newa0*(sqrt((double)3)/2.0);
	lasthighestx=lasthighestx+newa0*1.5;
	for (double x= lasthighestx;x<=sourceWidth+channelWidth; x=x+newa0 ){
		
		for (double y=b0/2;y<=channelHeight-newb0/2.0;y=y+2*newb0) {
			CCoord AP;
			AP.set_coords(x-A.get_x(),y-A.get_y());
			
			if (leftslipplane(AP,AB)==false) {
			
				CVortex newVortex;
				
				newVortex.set_pos(x,y);
				vorticesList.push_back(newVortex);
				if ((y+newb0)<channelHeight-newb0/2) {
					newVortex.set_pos(x+newa0/2.0,y+newb0);
					vorticesList.push_back(newVortex);
				}
			}
			
			newa0=calculatea0(x);
	
		}
	}
	
	
	
}
/*
void calculateForces(list <CVortex>:: iterator &q, vector <list <CVortex> > &cellLinkedList, list<CVortex>& pinsList, int t) {
	
	
	list<CVortex> wrappedVorticesList;
	
	
	//determine cell containing q particle
	int currentCell = (int)(floor)(q->get_x()/cellSize)+1;
	
	// calculate nearest neighbour cell list for particle q
	list <CVortex> nnList;
	nnList= cellLinkedList[currentCell];
	if (currentCell == 0) nnList.insert(nnList.end(),cellLinkedList[1].begin(),cellLinkedList[1].end());
	else if (currentCell == numcells-1) nnList.insert(nnList.end(),cellLinkedList[numcells-2].begin(),cellLinkedList[numcells-2].end());
	else {
		nnList.insert(nnList.end(),cellLinkedList[currentCell-1].begin(),cellLinkedList[currentCell-1].end());
		nnList.insert(nnList.end(),cellLinkedList[currentCell+1].begin(),cellLinkedList[currentCell+1].end());
	} 
	
	
	double forcex=0;
	double forcey=0;
	
	// temp solution
	double v1[2];
	double v2[2]; 
	double vortexSum[2]={0,0};
	double pinningSum[2]={0,0};
	double directionVector[2];
	double tempSum[2]={0,0};
	
	
	
	// vortex sum
	for (list<CVortex>::iterator p = nnList.begin();
			p != nnList.end(); ++p) {

		if (!eqtest(q->get_x(),p->get_x()) && !eqtest(q->get_y(),p->get_y())) { // do not include the same vortex
				
				v1[0]=q->get_x();
				v1[1]=q->get_y();
				v2[0]=p->get_x();
				v2[1]=p->get_y();
				
				double distcal= sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]));
			  
			  if (distcal!=distcal) { // error checking
					cout << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cout << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					cerr << "distcal error in routine: calculateForces.  Error: nan" << endl; 
					cerr << "t: " << t << "(" << v1[0] << ", " << v1[1] << ") (" << v2[0] << ", " << v2[1] << ")" << endl;
					
					
					}
			if (distcal < forceRange) { //only include close vortices
				directionVector[0]=v2[0]-v1[0];
				directionVector[1]=v2[1]-v1[1];
			
				if (directionVector[0] ==0 && directionVector[1]==0) {
					vortexSum[0]=vortexSum[0]-0.0000000001*fabs(gaussianRand());
					vortexSum[1]=vortexSum[1]-0.0000000001*fabs(gaussianRand());
				} else {
					vortexSum[0]=vortexSum[0]+forceForm(vvpot,v1,v2,Rv)*directionVector[0];
					vortexSum[1]=vortexSum[1]+forceForm(vvpot,v1,v2,Rv)*directionVector[1];
				}
				
			}
			
		}
			
		
		
	}
	
	
	if (vortexSum[0]!=vortexSum[0] || vortexSum[1] != vortexSum[1])
	cout << "t: " << t << "vort nan" << "(" << vortexSum[0] << ", " << vortexSum[1] << ")" << endl;
	
	
	
	// pin sum
	
	for (list<CVortex>::iterator p = pinsList.begin();
			p != pinsList.end(); ++p) {
	      v1[0]=q->get_x();
				v1[1]=q->get_y();
				v2[0]=p->get_x();
				v2[1]=p->get_y();
			
			if ( sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]) ) < forceRange) { //only include close vortices
			
				directionVector[0]=v2[0]-v1[0];
				directionVector[1]=v2[1]-v1[1];
			    //double forceFactor= forceForm(2,v1,v2,Rp);
						
				pinningSum[0]=pinningSum[0]+forceForm(vppot,v1,v2,Rp)*directionVector[0];
				pinningSum[1]=pinningSum[1]+forceForm(vppot,v1,v2,Rp)*directionVector[1];
				//pinningSum[0]=pinningSum[0]+forceFactor*directionVector[0];
				//pinningSum[1]=pinningSum[1]+forceFactor*directionVector[1];
				
			}
		}
	if (pinningSum[0]!=pinningSum[0] || pinningSum[1] != pinningSum[1])
	cout << "t: " << t << "pin nan" << "(" << pinningSum[0] << ", " << pinningSum[1] << ")" << endl;
	
	if ( rand()%100 >80) { // 1/5 of the time get temp force
			//tempSum[0]=directionVector[0]*At*(rand()%100)/100.0;
			//tempSum[1]=directionVector[1]*At*(rand()%100)/100.0;
			
			// At already scaled by a0 early in code
			// This form gives At as fraction of a0 fluctuated per timestep
			
			tempSum[0]=At*gaussianRand()/dt;
			tempSum[1]=At*gaussianRand()/dt;
			
			// This code ensures that vortices do not jump out of the channel
			// This gives a maximum effective At= 0.433
			//if (tempSum[0]*dt>0.5*b0) tempSum[0]=0.5*b0/dt;
			//if (tempSum[1]*dt>0.5*b0) tempSum[1]=0.5*b0/dt;
			
			//boost::math::cyl_bessel_k(1,  dist/lambda)
			
			//cout << "temp: " << "(" << tempSum[0] << ", " << tempSum[1] << ")"<< endl;
	}
	if (tempSum[0]!=tempSum[0] || tempSum[1] != tempSum[1]) {
	  tempSum[0]=0;
	  tempSum[1]=0;
		cout << "t: " << t << "temp nan" << "(" << tempSum[0] << ", " << tempSum[1] << ")" << endl;
	}
	//cout  << "t: " << t << "final temp sum: " << "(" << tempSum[0] << ", " << tempSum[1] << ")" << endl;
	//tempSum[0]=0;
	//tempSum[1]=0;
	
	
	//forcex= 1.1*(vortexSum[0]+pinningSum[0]+tempSum[0]);
	//forcey= 1.1*(vortexSum[1]+pinningSum[1]+tempSum[1]);
	
	forcex= vortexSum[0]+pinningSum[0];
	forcey= vortexSum[1]+pinningSum[1];
	//forcex=1e-8;
	//forcey=1e-8;
		
	//forcex= vortexSum[0]+pinningSum[0];
	//forcey= vortexSum[1]+pinningSum[1];
	//int lastcellx,lastcelly,newcellx,newcelly;
	//lastcellx = q->get_cellx();
	//lastcelly = q->get_celly();
	 
	if (forcex!=forcex || forcey != forcey || tempSum[0] != tempSum[0]
	 || tempSum[1]!=tempSum[1])
	cout << "t: " << t << "force nan" << "(" << forcex << ", " << forcey << ")" << endl;
	
	
	
	q->set_vel(forcex/eta+tempSum[0],forcey/eta+tempSum[1]);
	
	// Velocities are allowed to be a maximum 0.5*b0 /dt to avoid escaping vortices
	if(q->get_velx()*dt>0.5*b0) {
		cout << "vx velocity rectified " << q->get_velx()*dt/b0; 
		q->set_velx(0.5*b0/dt);
		cout << setw(15) << q->get_velx() << endl;
	}
	
	if(q->get_vely()*dt>0.5*b0) {
		cout << "vy velocity rectified " << q->get_velx()*dt/b0; 
	  q->set_vely(0.5*b0/dt);
	  cout << setw(15) << q->get_vely() << endl;
	}
	q->set_pos(q->get_x()+q->get_velx()*dt,q->get_y()+q->get_vely()*dt);

   
  
	
}
*/





/*void calculateCell(int i, vector <list <CVortex> > &cellLinkedList, list<CVortex>& pinsList, int t, MyThreadSafeList<CVortex>& newVorticesList) {
	
	
		
	for ( list <CVortex>:: iterator z= cellLinkedList[i].begin();
					z!=cellLinkedList[i].end(); ++z) {
						
			calculateForces(z,cellLinkedList,pinsList,t);
			CVortex newVortex;
			newVortex.set_pos(z->get_x(),z->get_y());
			newVortex.set_vel(z->get_velx(),z->get_vely());
			newVorticesList.push_back(newVortex);
	}
	
	
}
*/
/*
double calculateEnergy(list<CVortex> vorticesList, list<CVortex> pinsList, list<CVortex>::iterator q) {
	
	double Energy;
	double v1[2];
	double v2[2]; 
	double vortexSum=0;
	double pinningSum=0;
	
	
	// vortex sum
	for (list<CVortex>::iterator p = vorticesList.begin();
			p != vorticesList.end(); ++p) {
		
		v1[0]=q->get_x();
		v1[1]=q->get_y();
		v2[0]=p->get_x();
		v2[1]=p->get_y();
				
		double dist = sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]) );
		
		if ( dist < forceRange && dist!=0) { //only include close vortices
			
			vortexSum=vortexSum+energyForm(vvpot,dist,Rv);
				
				
		}
			
		
			
		
		
	}
	
	// pin sum
	
	for (list<CVortex>::iterator p = pinsList.begin();
			p != pinsList.end(); ++p) {
	     
		v1[0]=q->get_x();
		v1[1]=q->get_y();
		v2[0]=p->get_x();
		v2[1]=p->get_y();
		
		double dist = sqrt((double)(v1[0]- v2[0])*(v1[0]- v2[0])+ (v1[1]- v2[1])*(v1[1]- v2[1]) );
		
		    
		if ( dist < forceRange && dist!=0) {	
				
			pinningSum=pinningSum+energyForm(vppot,dist,Rp);
			
		}
	}
	
	Energy = vortexSum+pinningSum;
	
	
	
	return Energy;
	
}
*/
inline double dotProduct(CCoord a, CCoord b) {
	return a.get_x() * b.get_x() + a.get_y() * b.get_y();
	
}

inline double det(double a, double b, double c, double d, double e, double f, double g, double h, double i){
	/*
	  a  b  c 
	  d  e  f
	  g  h  i
	*/

	return a*( e*i - h*f) -b*(d*i-g*f) +c*(d*h -g*e);
}

inline double circum(list<CDelTriangle>::iterator p,list<CVortex>::iterator q) {
	double Ax, Bx, Cx, Dx, Ay, By, Cy, Dy;
	Ax=p->get_Ax();
	Ay=p->get_Ay();
	Bx=p->get_Bx();
	By=p->get_By();
	Cx=p->get_Cx();
	Cy=p->get_Cy();
	Dx=q->get_x();
	Dy=q->get_y();
	return det(	Ax-Dx,Ay-Dy,Ax*Ax-Dx*Dx+Ay*Ay-Dy*Dy,
				Bx-Dx,By-Dy,Bx*Bx-Dx*Dx+By*By-Dy*Dy,
				Cx-Dx,Cy-Dy,Cx*Cx-Dx*Dx+Cy*Cy-Dy*Dy);
	
}

bool pointInTriangle( CDelTriangle& p,  list<CVortex>::iterator& q) {
	
	CCoord Q;
	
	CCoord A;
	CCoord B;
	CCoord C;
	CCoord v0;
	CCoord v1;
	CCoord v2;
	
	double dot00, dot01, dot02, dot11, dot12, invDenom,u,v;
	bool inTriangle=false;
	
	A.set_coords(p.get_Ax(),p.get_Ay());
	B.set_coords(p.get_Bx(),p.get_By());
	C.set_coords(p.get_Cx(),p.get_Cy());
	Q.set_coords(q->get_x(),q->get_y());
	
	
	
	v0 = C - A;
	v1 = B - A;
	v2 = Q - A;
	//cout << "v0: ("  << C.get_x() << ", " << C.get_y() << ") - ("  << A.get_x() << ", " << A.get_y() << ") = (" << v0.get_x() << ", " << v0.get_y() << ")" << endl;
	//cout << "v1: ("  << B.get_x() << ", " << B.get_y() << ") - ("  << A.get_x() << ", " << A.get_y() << ") = (" << v1.get_x() << ", " << v1.get_y() << ")" << endl;
	//cout << "v2: ("  << Q.get_x() << ", " << Q.get_y() << ") - ("  << A.get_x() << ", " << A.get_y() << ") = (" << v2.get_x() << ", " << v2.get_y() << ")" << endl;
	
	// calculate dot products
	dot00 = dotProduct(v0, v0);
	dot01 = dotProduct(v0, v1);
	dot02 = dotProduct(v0, v2);
	dot11 = dotProduct(v1, v1);
	dot12 = dotProduct(v1, v2);

	// Compute barycentric coordinates
	invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	
	if ((u >= 0) && (v >= 0) && (u + v < 1)) {inTriangle =true;}
	return inTriangle;
	 /*
	
	// Compute vectors        
	v0 = C - A
	v1 = B - A
	v2 = Q - A

	
	dot00 = dot(v0, v0)
	dot01 = dot(v0, v1)
	dot02 = dot(v0, v2)
	dot11 = dot(v1, v1)
	dot12 = dot(v1, v2)

	// Compute barycentric coordinates
invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
u = (dot11 * dot02 - dot01 * dot12) * invDenom
v = (dot00 * dot12 - dot01 * dot02) * invDenom

// Check if point is in triangle
return (u >= 0) && (v >= 0) && (u + v < 1)
	*/
	
	
}

void delaunayTriangulation(list<CVortex>::iterator& corner, list<CVortex> vorticesList, list<CVortex>& delVortexList, list<CVortex>& pinsList, list<CDelLine>& delLinesList, list<CDelTriangle>& delTrianglesList) {
	
	//add edge of pins list to vortices list
	for(list<CVortex>::iterator p=pinsList.begin();
			p!=pinsList.end();p++){
		if (  p->get_x() > etchchannelx0 && p->get_x() < etchchannelx1
			  && p->get_y() > etchchannely0-a0 && p->get_y() < etchchannely1+a0 ){
				(*p).set_ghost();	
			vorticesList.push_back(*p);
		}
	}
			
	
	
	list<CDelLine> allEdges;
	list<CDelLine> uniqueEdges;	
	set<CCoord> uniquePoints;
		
	// add p_{-1} && p_{-2}
	delVortexList.clear();
	delLinesList.clear();
	delTrianglesList.clear();
	CVortex pminus2;
	
	CVortex pminus1;
	CCoord t0,t1, t2;  //corners of large triangle
	
	t0.set_coords(-5*a0,-sourceHeight);
	t1.set_coords(-5*a0,100*channelHeight);
	t2.set_coords((sourceWidth+channelWidth+sinkWidth)*100,-sourceHeight);
	
	
	pminus2.set_pos(t2.get_x(),t2.get_y());
	pminus1.set_pos(t1.get_x(),t1.get_y());
	
	delVortexList.push_back(pminus2);
	delVortexList.push_back(pminus1);
	
	CVortex p0;
	
	p0.set_pos(t0.get_x(),t0.get_y());
		
	delVortexList.push_back(p0);
	corner= delVortexList.end();
	--corner;
	
	//drawlines list
	CDelLine newDelLine;
	newDelLine.set_points(t0.get_x(),t0.get_y(),t2.get_x(),t2.get_y());
	delLinesList.push_back(newDelLine);	
	
		
	newDelLine.set_points(t2.get_x(),t2.get_y(),t1.get_x(),t1.get_y());
	delLinesList.push_back(newDelLine);	
	
	
	newDelLine.set_points(t1.get_x(),t1.get_y(),t0.get_x(),t0.get_y());
	delLinesList.push_back(newDelLine);	
	
	
	
	// add new delTriangle
	CDelTriangle newDelTriangle;
	newDelTriangle.set_vertices(t0.get_x(),t0.get_y(),
															t2.get_x(),t2.get_y(),
															t1.get_x(),t1.get_y());
	newDelTriangle.set_finalDaughter(false);
	delTrianglesList.push_back(newDelTriangle);
	
	
	// choose next vortex from vorticesList to add to the Delaunay Triangulation
	
	
	
	for (list<CVortex>::iterator q = vorticesList.begin();
			q!=vorticesList.end();++q) {	
		//if (q->get_x()>sourceWidth+latticeSpacing/2.0 && q->get_x() <sourceWidth+channelWidth-latticeSpacing/2.0 ) {
	
		list<CDelTriangle>::iterator p=delTrianglesList.begin();
		
		
		allEdges.clear();
		uniqueEdges.clear();
		
		while (p!=delTrianglesList.end()) {
			if (circum(p,q)>0) 
			{
			
				// add three edges to the list of edges affected
				CDelLine addedge;
				
				addedge.set_points(p->get_Ax(),p->get_Ay(),p->get_Bx(),p->get_By());
				allEdges.push_back(addedge);
				addedge.set_points(p->get_Bx(),p->get_By(),p->get_Cx(),p->get_Cy());
				allEdges.push_back(addedge);
				addedge.set_points(p->get_Cx(),p->get_Cy(),p->get_Ax(),p->get_Ay());
				allEdges.push_back(addedge);
				
				// remove illegal lines
				
				list<CDelLine>::iterator r=delLinesList.begin();
				
				while(r!=delLinesList.end()){
					if( 	(	r->get_x1()==p->get_Ax() && r->get_y1()==p->get_Ay() 
							&& r->get_x2()==p->get_Bx() && r->get_y2()==p->get_By() )
							
							||
							(	r->get_x1()==p->get_Bx() && r->get_y1()==p->get_By() 
							&& r->get_x2()==p->get_Ax() && r->get_y2()==p->get_Ay() )
							)
							
							 {
						r=delLinesList.erase(r);
						
							
					}
					else if(( 	r->get_x1()==p->get_Bx() && r->get_y1()==p->get_By() 
							&& r->get_x2()==p->get_Cx() && r->get_y2()==p->get_Cy() )
							||
							( 	r->get_x1()==p->get_Cx() && r->get_y1()==p->get_Cy() 
							&& r->get_x2()==p->get_Bx() && r->get_y2()==p->get_By() )
							)
							 {
								
						r=delLinesList.erase(r);
						
							
					}
					else if(( 	r->get_x1()==p->get_Cx() && r->get_y1()==p->get_Cy() 
							&& r->get_x2()==p->get_Ax() && r->get_y2()==p->get_Ay() )
							||
							(  	r->get_x1()==p->get_Ax() && r->get_y1()==p->get_Ay() 
							&& r->get_x2()==p->get_Cx() && r->get_y2()==p->get_Cy()  )
												
							) {
						r=delLinesList.erase(r);
						
					}
					else {
						++r;
					}
							
					
				}
				
				// remove triangle
				
				p=delTrianglesList.erase(p);
				
				
				
			}
			else {
				++p;
			}
			
			
		
		}
		//remove all double edges from edgebuffer,keeping only unique ones
		bool doubleEdge;
		for (list<CDelLine>::iterator r= allEdges.begin();
				r!=allEdges.end();++r) {
			doubleEdge=false;
			for (list<CDelLine>::iterator s= allEdges.begin();
					s!=allEdges.end();++s) {
				if (r!=s) {
					if ( 		(eqtest(r->get_x1(),s->get_x1()) && eqtest(r->get_y1(),s->get_y1()) && eqtest(r->get_x2(),s->get_x2()) && eqtest(r->get_y2(),s->get_y2()))
							||	(eqtest(r->get_x1(),s->get_x2()) && eqtest(r->get_y1(),s->get_y2()) && eqtest(r->get_x2(),s->get_x1()) && eqtest(r->get_y2(),s->get_y1()))	) {
						// duplicated line
						doubleEdge=true;
					}
				}	
			
			} 
			if (doubleEdge==false) {
				uniqueEdges.push_back(*r);
			}
			
			
		}
		//cout << "uniqueEdges Size: " << uniqueEdges.size() << endl;
		
		/*uniquePoints.clear();
		CCoord newCoord;*/
		
		// count unique points
		//double uniquePointsArray[50][2];
		//double firstInArray=false;
		//double secondInArray=false;
		
		
		
		// form new triangles && lines between edges && vertices and calculate coord num
		for (list<CDelLine>::iterator p = uniqueEdges.begin();
				p!=uniqueEdges.end();++p)
			{
					
				
				//form a new triangle between edge and vertex
				newDelTriangle.set_vertices(p->get_x1(),p->get_y1(),p->get_x2(),p->get_y2(),q->get_x(),q->get_y());
				delTrianglesList.push_back(newDelTriangle);
				
				//add new lines for sdl
				
				newDelLine.set_points(p->get_x1(),p->get_y1(),p->get_x2(),p->get_y2());
				delLinesList.push_back(newDelLine);	
	
				newDelLine.set_points(p->get_x2(),p->get_y2(),q->get_x(),q->get_y());
				delLinesList.push_back(newDelLine);	
		
				newDelLine.set_points(q->get_x(),q->get_y(),p->get_x1(),p->get_y1());
				delLinesList.push_back(newDelLine);	
				
			
				
				
			}
		
	
		delVortexList.push_back(*q);
		//}
		
	  
	}
	
	
	// check delLinesList for dupilcates
	list<CDelLine> uniqueLinesList;
	//CDelLine testLine;
	bool addLine;
	
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end();++p) {
		//testLine=*p;
		addLine=true;
		for (list<CDelLine>::iterator q = uniqueLinesList.begin();
				q!=uniqueLinesList.end();++q) {
			if (  
				(eqtest(p->get_x1(),q->get_x1()) && eqtest(p->get_y1(),q->get_y1()) && eqtest(p->get_x2(),q->get_x2()) && eqtest(p->get_y2(),q->get_y2())) ||
				(eqtest(p->get_x2(),q->get_x1()) && eqtest(p->get_y2(),q->get_y1()) && eqtest(p->get_x1(),q->get_x2()) && eqtest(p->get_y1(),q->get_y2()))
		
			
			  ) {
				addLine=false;
			}
		}
		if (addLine==true) {
			uniqueLinesList.push_back(*p);
		}
			
	}
	
	int numNeighbours=0;
	for (list<CVortex>::iterator p = delVortexList.begin();
			p!=delVortexList.end();++p) {
		
		numNeighbours=0;
		double a0Sum=0;
		for (list<CDelLine>::iterator q = uniqueLinesList.begin();
				q!=uniqueLinesList.end();++q) {
			
			if (   (eqtest( p->get_x(),q->get_x1()) && eqtest(p->get_y(),q->get_y1()))
							|| (eqtest(p->get_x(),q->get_x2()) && eqtest(p->get_y(),q->get_y2() )) ){
				  a0Sum += q->get_length();
					numNeighbours++;
			}
			
			
			
			
			
					
		}
		
		
		//cout << "Num Neighbours: " << numNeighbours << endl;
		p->set_coordNum(numNeighbours);
		p->set_a0(a0Sum/((double)numNeighbours));
	}

	
	
	// remove initial traingle
	list<CDelLine>::iterator itLines = delLinesList.begin();
	while (itLines!=delLinesList.end()) {
		if (	   (  eqtest(itLines->get_x1(),t0.get_x()) && eqtest(itLines->get_y1(),t0.get_y()))
				|| ( eqtest(itLines->get_x2(),t0.get_x()) && eqtest(itLines->get_y2(),t0.get_y()) )
				
				|| ( eqtest(itLines->get_x1(),t2.get_x()) && eqtest(itLines->get_y1(),t2.get_y()) )
				|| ( eqtest(itLines->get_x2(),t2.get_x()) && eqtest(itLines->get_y2(),t2.get_y()) ) 
				
				|| ( eqtest(itLines->get_x1(),t1.get_x())  && eqtest(itLines->get_y1(),t1.get_y()))
				|| ( eqtest(itLines->get_x2(),t1.get_x())  && eqtest(itLines->get_y2(),t1.get_y()))
				
			) {
	
			itLines=delLinesList.erase(itLines);
		} else	{
			//cout  <<"line: " << itLines->get_x1() << ", " << itLines->get_y1() << " -->  " << itLines->get_x2() << ", " << itLines->get_y2() << endl;
		
			++itLines;
		}
			
	}
	
	
	
	list<CVortex>::iterator d=delVortexList.begin();
	delVortexList.erase(d);
	
	d=delVortexList.begin();
	delVortexList.erase(d);
	
	d=delVortexList.begin();
	delVortexList.erase(d);
	
	
}

double systemCoordNum(list<CVortex>& vorticesList){
	double sum=0;
	double count=0;
	for (list<CVortex>::iterator p = vorticesList.begin();
			p!=vorticesList.end();++p) {
		if (p->get_x()>sourceWidth+2*latticeSpacing && p->get_x() < sourceWidth+channelWidth-2*latticeSpacing) {
			if ( p->get_y()>sourceHeight/2.0-channelHeight/2.0+2*latticeSpacing && p->get_y() < sourceHeight/2.0+channelHeight/2.0-2*latticeSpacing) {
				sum=sum+p->get_coordNum();
				count++;
			}
		}
	}
	return sum/count;
	
}

void calculateDislocationPaths(vector <vector<CTemporalCoord> >& pathsVector, list<CVortex>& delVortexList, int t) {
	
	//vector<CCoord> latestDislocationPosVector;
	
	
	
	for (list<CVortex>::iterator p = delVortexList.begin();
			p!=delVortexList.end();++p ) {
				
		if (p->get_coordNum()==5) {
		
		if (p->get_y()> bulky0 && p->get_y() < bulky1 && 
				p->get_x()> bulkx0 && p->get_x() < bulkx1) { // inside the white calculation zone 
				
			bool newpath=true;	
			for (vector< vector<CTemporalCoord> >::iterator iter_ii = pathsVector.begin();
					iter_ii != pathsVector.end(); ++iter_ii) {
				//iterates through the previous time step dislocation positions to
				// check if any of the new positions are related to the last
				
				int last = (*iter_ii).size();
				
				if ( (*iter_ii)[last-1].get_time()==(t-1)) {
					if  ( sqrt((double) ( p->get_x()-(*iter_ii)[last-1].get_x() )*(  p->get_x()-(*iter_ii)[last-1].get_x() )
								+  (p->get_y()-(*iter_ii)[last-1].get_y() ) * (  p->get_y()-(*iter_ii)[last-1].get_y() ) )  < 0.02*a0 ) {  //i,e dist per step = 0.14*.05 = 0.007
						// add dislocation to this current path
						
						newpath=false;
						CTemporalCoord newTemporalCoord;
						newTemporalCoord.add_TemporalCoord(p->get_x(), p->get_y(), t); 
						(*iter_ii).push_back(newTemporalCoord);
			
			
					}
					
					
					
				}
			}
			if (newpath==true) {
						// add new path
						
						CTemporalCoord newTemporalCoord;
						newTemporalCoord.add_TemporalCoord(p->get_x(), p->get_y(), t); 
						
						vector<CTemporalCoord> path;
						path.push_back(newTemporalCoord);
						pathsVector.push_back(path);
			}
		
		
		}
		}
	}
	
	
	
	
	
}

void readData(vector < list<CVortex> >& delVortexListVector){
	
	ifstream myfile ("guidata.txt");
	list<CVortex> delVortexList;
	
	if (myfile.is_open()) {
		
		
		string input;
		double xval;
		double yval;
		int coord_num;
		string a0str;
		double a0double;
		string timeline ("timestep");
		string nan ("nan");
		size_t found;
		bool first=true;
		
		while ( myfile.good() &&  (int)delVortexListVector.size()!=numTimesteps)
		{
			
			myfile >> input;
			found = input.find(timeline);
		    
			if (found==string::npos) {
				
				// add new points to delVortexList
				// converts first to double
				char charinput[256];
				strcpy(charinput, input.c_str());
				xval = atof (charinput); 
				
				// reads next three directly to correct variable type
				myfile >> yval;
				myfile >> coord_num;
				myfile >> a0str;
				
				if (a0str.find(nan)==string::npos) {
					strcpy(charinput, a0str.c_str());
					
					a0double = atof (charinput); 
					// add new vortex
					CVortex newVortex;
					newVortex.set_pos(xval,yval);
					newVortex.set_coordNum(coord_num);
					newVortex.set_a0(a0double);
					delVortexList.push_back(newVortex);
					//cout  << setw(20) << newVortex.get_x() << setw(20) << newVortex.get_y() << setw(20) << newVortex.get_coordNum() << setw(20) << newVortex.get_a0() << endl;
		
					
				} else {
					cout << "nan" << endl;
				}	
				
				
			} else if (found!=string::npos) {  // found line containing timestep
				
				if (first==true) {
					first=false;
				}
				else{
					//cout << "t:" << delVortexListVector.size() << endl;
					//cout << delVortexList.size() << endl;
					delVortexListVector.push_back(delVortexList);
					delVortexList.clear();
				}
				
			}
			
			
     /* myfile >> xval;
      myfile >> yval;
     
      
      
      CVortex newVortex;
      newVortex.set_pos
					(xval,yval);
			vorticesList.push_back(newVortex);
			
			*/
    }
    myfile.close();
		//cout << delVortexList.size() << endl;
		delVortexListVector.push_back(delVortexList);
		//cout << delVortexListVector.size() << endl;
					
	}
	
}

void readSingleDataStep(list<CVortex>& delVortexList, int t){
	
	if (!readStepDataFile.is_open()) {
	 
	 ostringstream jobstr;
		jobstr.str("");
   jobstr << ini.readDataLocation << "job" << jobnum << "//guidata.txt"; 
  
  char renderChar[100];
		string renderStr;
		renderStr = jobstr.str();
		strcpy(renderChar,renderStr.c_str());
   cout << jobstr.str();
	 
	 readStepDataFile.open(renderChar);
	 cout << "loadfile";	
	}
	
	 
	
	if (readStepDataFile.is_open()) {
		//cout << "Guidata.txt opened";
		
		string input;
		string input2;
		double xval;
		double yval;
		int coord_num;
		
		double a0double;
		string timeline ("timestep");
		string nan ("nan");
		string inf ("INF");
		size_t found;
		static bool first = true;
		bool endofstep=false;
		
		while ( readStepDataFile.good() && endofstep==false)
		{
			
			readStepDataFile >> input;
			found = input.find(timeline);
		    
			if (found==string::npos) {
				bool nanFound= false;
				// add new points to delVortexList
				// converts first to double
				//if (t>4000) cout << input << "  ";
				if (input.find("INF")==string::npos || input2.find("nan")==string::npos) { 
					char charinput[256];
					strcpy(charinput, input.c_str());
					xval = atof (charinput); 
				} else nanFound=true;
				
				// reads next three into strings, check if numnber
				// then converts to double
				readStepDataFile >> input2;
				//if (t>4000) cout << input2 << "  ";
				if (input2.find("INF")==string::npos || input2.find("nan")==string::npos ) { 
					char charinput[256];
					strcpy(charinput, input2.c_str());
					yval = atof (charinput);
				} else nanFound=true;
				
				readStepDataFile >> input2;
				//if (t>4000) cout << input2 << "  ";
				if (input2.find("INF")==string::npos || input2.find("nan")==string::npos) { 
					char charinput[256];
					strcpy(charinput, input2.c_str());
					coord_num = (int)atof (charinput);
				} else nanFound=true;
				
				readStepDataFile >> input2;
				//if (t>4000) cout << input2 << "  " << endl;
				if (input2.find("INF")==string::npos || input2.find("nan")==string::npos) { 
					char charinput[256];
					strcpy(charinput, input2.c_str());
					a0double = atof (charinput);
				} else nanFound=true;
				
				
				//cout << yval << endl;
				
				if (nanFound==false) {							
					//cout << "Found coordinates" << endl;
					// add new vortex
					CVortex newVortex;
					newVortex.set_pos(xval,yval);
					newVortex.set_coordNum(coord_num);
					newVortex.set_a0(a0double);
					delVortexList.push_back(newVortex);

					//cout  << setw(20) << newVortex.get_x() << setw(20) << newVortex.get_y() << setw(20) << newVortex.get_coordNum() << setw(20) << newVortex.get_a0() << endl;
		
					
				} else {
					cout << "nan" << endl;
				}	
				
				
			} else if (found!=string::npos) {  // found line containing timestep
				//cout << "input: " << input << endl;
				if (first==true) {
					first=false;
				}
				else {
					//cout << "End of read step" << endl;
					//cout << delVortexList.size() << endl;
					endofstep=true;
					
				}
				
			}
			
			
     /* myfile >> xval;
      myfile >> yval;
     
      
      
      CVortex newVortex;
      newVortex.set_pos
					(xval,yval);
			vorticesList.push_back(newVortex);
			
			*/
    
    }
    //cout << "finished read: return to main loop" << endl;
   // readStepDataFile.close();
		//cout << "Num vortices read in for timestep " << t << ": " << delVortexList.size() << endl;
	  
					
	}
	
}



void readIniFile(){
	
	ifstream read;
	read.open("config.ini");
	string input;
	
	 
	
	if (read.is_open()) {
		read >> input;
		ini.readDataLocation= input;
		
	}
}


void removeEscapedVortices(list<CVortex>& vorticesList, int t) {
  list<CVortex>::iterator p = vorticesList.begin();
  

  
  
	while (p != vorticesList.end()) {
		bool removed=false;
		if (p->get_x() <= removesourcex) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} else if (p->get_y() <= removesourcey0) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
		}
		else if (p->get_x() >= removesinkx) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} 
		else if (p->get_y() >= removesourcey1) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
		} 
		else if ( (p->get_x() >= removechannelx0 && p->get_x()<= removechannelx1) 
				&& ( p->get_y() <= removetopchannely  || p->get_y() >= removebottomchannely ) ) {
			cout << "Channel Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} 
		
		
		
	
	
	
		if (removed==false) { ++p; }
	
	}
}


void removeEscapedVorticesTube(list<CVortex>& vorticesList, int t) {
  list<CVortex>::iterator p = vorticesList.begin();
  

  
  
	while (p != vorticesList.end()) {
		bool removed=false;
		if (p->get_x() <= removefunnelx) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		}
		else if (p->get_y() <= -funnelWidth*tan(funnelAngle)-2*b0 || p->get_y() >= funnelWidth*tan(funnelAngle)+channelHeight+2*b0) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
		
		}
		else if (p->get_x() >= removesinkx) {
			cout << "Removed at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p=vorticesList.erase(p);
			//p--;
			removed=true;
			
		} 
		
		/* else if ((p->get_y() <= removesourcey0) and  (p->get_x() < removechannelx0 ) )  {
			cout << "wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			//p=vorticesList.erase(p);
			p->set_pos(p->get_x(),p->get_y()+channelHeight+b0);
			//p--;
			//removed=true;
		
		} else if ((p->get_y() > removesourcey1) and  (p->get_x() < removechannelx0) )  {
			cout << "wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			//p=vorticesList.erase(p);
			p->set_pos(p->get_x(),p->get_y()-channelHeight-b0);
			//p--;
			//removed=true;
		
		}*/
		 else if ((p->get_y() <= removetopchannely) &&  (p->get_x() > removesourcex) )  {
			cout << "wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p->set_pos(p->get_x(),p->get_y()+channelHeight+b0);
			//p=vorticesList.erase(p);
			//p--;
			//removed=true;
		}  else if ((p->get_y() >= removebottomchannely) &&  (p->get_x() > removesourcex) )  {
			cout << "wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p->set_pos(p->get_x(),p->get_y()-channelHeight-b0);
			//p=vorticesList.erase(p);
			//p--;
			//removed=true;
		}
		
		
		
		else if ((p->get_y() <= -fabs(p->get_x())*tan(funnelAngle) ) &&  (p->get_x() < removesourcex) )  {
			cout << "funnel wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p->set_pos(p->get_x(),p->get_y()+channelHeight+b0+2.0*p->get_x()*tan(funnelAngle));
			//p=vorticesList.erase(p);
			//p--;
			//removed=true;
		}  else if ((p->get_y() >= channelHeight+b0+fabs(p->get_x())*tan(funnelAngle)) &&  (p->get_x() < removesourcex) )  {
			cout << "funnel wrapped at " << t << " pos (" <<  p->get_x() << ", " << p->get_y() << ") vel (" << p->get_velx() << ", " << p->get_vely() << ")" << endl;  
			p->set_pos(p->get_x(),p->get_y()-channelHeight-b0-2.0*p->get_x()*tan(funnelAngle));
			//p=vorticesList.erase(p);
			//p--;
			//removed=true;
		}
		
		
		
		
		
		
	
	
	
		if (removed==false) { ++p; }
	
	}
}




void calculateCellNeighboursList(vector < list<int> > &cellNeighboursList, vector <list <CVortex> > &cellLinkedList, double forceRange) {
  //define cells using single index i
	
	cout << "Nearest Neighour cells list - number of cells: " << numcells << endl;
	// add first list
	list<int> newList;
	int newInt=1;
	newList.push_back(newInt);	
	cellNeighboursList.push_back(newList);
	
	for (int i =1;i<=numcells-2;i++) {
			newList.clear();
			list<int> newList;
			newInt=i-1;
			newList.push_back(newInt);	
			newInt=i+1;
			newList.push_back(newInt);	
			cellNeighboursList.push_back(newList);
			
		
		
	}
	
	//add last list
	
	newList.clear();
	newInt=numcells-2;
	newList.push_back(newInt);	
	cellNeighboursList.push_back(newList);
	
	
	
	// output vector of lists
	
	for (vector < list<int> > :: iterator p=cellNeighboursList.begin();
			p!=cellNeighboursList.end();++p) {
				cout << p- cellNeighboursList.begin() << ":  ";
				for (list<int> ::iterator q=(*p).begin();
						q!=(*p).end();++q) {
					cout << *q << "  ";
					
					
							
				}
				cout << "(" << forceRange*(p- cellNeighboursList.begin()-3)/a0  << " a0 -> "  << forceRange*(p- cellNeighboursList.begin()-2)/a0  << " a0)" << endl;
	}
	
	//initialise cellLinkedLIst
	
	for (int i =0;i<=numcells-1;i++) {
		list <CVortex> newList;
		cellLinkedList.push_back(newList);	
		cout << i << ": " << cellLinkedList[i].size() << endl;	
		
		
	}
	
	
}

void initialiseFiles(double tau) {
	// position data
	fileOss.str("");
	fileOss << dirOss.str() << "//posdata.txt";
	fileStr= fileOss.str();
	
	cout << "Initialied new Files" << endl;
  posfile.open(fileStr.c_str());
  posfile.precision(5);
  
  // density results
 
	fileOss.str("");
	fileOss << dirOss.str() << "//densitydata.txt";
	fileStr= fileOss.str();
  
  densityfile.open(fileStr.c_str());
  densityfile.precision(5);
  
  // density results
 
	fileOss.str("");
	fileOss << dirOss.str() << "//newVortexData.txt";
	fileStr= fileOss.str();
  
	newVortexfile.open(fileStr.c_str());
	newVortexfile.precision(5);
  
  
   
  // a0 results
 
	fileOss.str("");
	fileOss << dirOss.str() << "//a0data.txt";
	fileStr= fileOss.str();
  
  
  
  a0file.open(fileStr.c_str());
  a0file.precision(5);
  //output a0 file header
  /*a0file << setw(15) << "Time";
	for (int i=0; i<=numa0Bins-1;i++ ) {  // USE THIS IF A0 DATA OUTPUT IN AVERGAED BINS
		a0file <<setw(15) << (urectx1-urectx0)/numa0Bins*(i+.5);
	}	
	a0file << endl;	*/
  
  // averaged line lengths results
 
	fileOss.str("");
	fileOss << dirOss.str() << "//avlinelengthdata.txt";
	fileStr= fileOss.str();
  
  
  
  avlinelengthfile.open(fileStr.c_str());
  avlinelengthfile.precision(5);
  
  
  
   // b0 results
 
	fileOss.str("");
	fileOss << dirOss.str() << "//b0data.txt";
	fileStr= fileOss.str();
  
  
  
  b0file.open(fileStr.c_str());
  b0file.precision(5);
  
  
  // dislocation results
  fileOss.str("");
	fileOss << dirOss.str() << "//disdata.txt";
	fileStr= fileOss.str();
  

  disfile.open(fileStr.c_str());
  disfile.precision(5);
  disfile << "Modified output range (not whole channel)" <<endl;
  
  // raw velocity results
  fileOss.str("");
	fileOss << dirOss.str() << "//rawvdata.txt";
	fileStr= fileOss.str();
  

  rawvfile.open(fileStr.c_str());
  rawvfile.precision(5);
  
  
  
  
  // v(x) results
  fileOss.str("");
	fileOss << dirOss.str() << "//vxdata.txt";
	fileStr= fileOss.str();
	
 
  vxfile.open(fileStr.c_str());
  vxfile.precision(5);
  
  //vfile << setw(15) << "Time";
	/*for (int i=0; i<=numa0Bins-1;i++ ) {
		vfile <<setw(15) << (urectx1-urectx0)/numa0Bins*(i+.5);
	}	*/
	//vfile << endl;	
  
  // v(y) results
  fileOss.str("");
	fileOss << dirOss.str() << "//vydata.txt";
	fileStr= fileOss.str();
	
 
  vyfile.open(fileStr.c_str());
  vyfile.precision(5);
  
  /*for (int i=0; i<=numvyBins-1;i++ ) {
		vfile <<setw(15) << (urecty1-urecty0)/numvyBins*(i+.5);
	}	
	vfile << endl;	
  */
  
  // file needed for displaying gui data with read program
  fileOss.str("");
	fileOss << dirOss.str() << "//jobheader.txt";
	fileStr= fileOss.str();
	
 
  guiheader.open(fileStr.c_str());
  guiheader.precision(5);
  
  guiheader << setw(20) << seedtime << endl
  << setw(20) << sourceBfield << endl
  << setw(20) << sinkBfield << endl
  << setw(20) << sourceWidth/a0 << endl
  << setw(20) << sourceHeight/b0 << endl
  << setw(20) << channelWidth/a0 << endl
  << setw(20) << channelHeight/b0 << endl
  << setw(20) << numTimesteps << endl
  << setw(20) << tau << endl;
  
  
  fileOss.str("");
	fileOss << dirOss.str() << "//guidata.txt";
	fileStr= fileOss.str();
	 
  
  guifile.open(fileStr.c_str());
  guifile.precision(5);
  
  
  // velocity field data
	fileOss.str("");
	fileOss << dirOss.str() << "//velfielddata.txt";
	fileStr= fileOss.str();
	
 
  velfieldfile.open(fileStr.c_str());
  velfieldfile.precision(5);
  
  // averageda0data.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//averageda0data.csv";
	fileStr= fileOss.str();
	
 
  averageda0data.open(fileStr.c_str());
  averageda0data.precision(5);
  
  // averagedvxdata.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//averagedvxdata.csv";
	fileStr= fileOss.str();
	
 
  averagedvxdata.open(fileStr.c_str());
  averagedvxdata.precision(5);
  
  // averagedBfielddata.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//averagedBfielddata.csv";
	fileStr= fileOss.str();
	
 
  averagedBfielddata.open(fileStr.c_str());
  averagedBfielddata.precision(5);
  
  // Jyy.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//Jyy.dat";
	fileStr= fileOss.str();
	
 
  Jyydata.open(fileStr.c_str());
  Jyydata.precision(18);
  //Jyydata << "{";
  
  // Jxx.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//Jxx.dat";
	fileStr= fileOss.str();
	
 
  Jxxdata.open(fileStr.c_str());
  Jxxdata.precision(18);
  //Jyydata << "{";
  
  // Jxy.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//Jxy.dat";
	fileStr= fileOss.str();
	
 
  Jxydata.open(fileStr.c_str());
  Jxydata.precision(18);
  
  // Jyx.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//Jyx.dat";
	fileStr= fileOss.str();
	
 
  Jyxdata.open(fileStr.c_str());
  Jyxdata.precision(18);
  
  
  
  
  // mathematicaGuiFile.csv
	fileOss.str("");
	fileOss << dirOss.str() << "//mathematicaGuiFile.dat";
	fileStr= fileOss.str();
	
 
  mathematicaGuiFile.open(fileStr.c_str());
  mathematicaGuiFile.precision(18);
  //Jyydata << "{";
  
  
  
}

void calculateCellLinkedLists(list<CVortex> vorticesList, vector< list <CVortex> > &cellLinkedList, double cellSize) {
	// updates the vector of lists of pointers to each vortex in each cell

	//clear all lists in the vector
  for (int i =0;i<=numcells-1;i++) {
		cellLinkedList[i].clear();	
	}

	// sort the list into ascending order on x
	vorticesList.sort (CVortexSort);
	
	int currentcellnum=0;
	for (list<CVortex>::iterator p=vorticesList.begin();
			p!=vorticesList.end();++p) {
		if (p->get_x() > cellSize*(double)(currentcellnum-1) && p->get_x() <= cellSize*(double)(currentcellnum) ) {
			// add vortex to cell list
			cellLinkedList[currentcellnum].push_back(*p);
			
		} 
		else {currentcellnum++; --p;}
			 
				
	}
	
	
	
	
	
}
#if defined (__WINDOWS__)

void SDL_GL_RenderText(char *text, 
                      TTF_Font *font,
                      SDL_Color color,
                      SDL_Rect *location)
{
	SDL_Surface *initial;
	SDL_Surface *intermediary;
	
	int w,h;
	GLuint texture;
	//SDL_FillRect(initial, NULL, SDL_MapRGB(initial->format, 255, 255, 255));	
	/* Use SDL_TTF to render our text */
	initial = TTF_RenderText_Blended(font, text, color);
	
	/* Convert the rendered text to a known format */
	//w = nextpoweroftwo(initial->w);
	//h = nextpoweroftwo(initial->h);
	w=initial->w;
	h=initial->h;
	intermediary = SDL_CreateRGBSurface(0, w, h, 32, 
			0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
	SDL_FillRect(intermediary, NULL, SDL_MapRGB(intermediary->format, 255, 255, 255));	
	
	SDL_BlitSurface(initial, 0, intermediary, 0);
	
	/* Tell GL about our new texture */
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, 4, w, h, 0, GL_BGRA, 
			GL_UNSIGNED_BYTE, intermediary->pixels );
	
	/* GL_NEAREST looks horrible, if scaled... */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);	

	/* prepare to render our texture */
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture);
	glColor3f(1.0f, 1.0f, 1.0f);
	
	/* Draw a quad at location */
	glBegin(GL_QUADS);
		/* Recall that the origin is in the lower-left corner
		   That is why the TexCoords specify different corners
		   than the Vertex coors seem to. */
		glTexCoord2f(0.0f, 1.0f); 
			glVertex2f(location->x    , location->y);
		glTexCoord2f(1.0f, 1.0f); 
			glVertex2f(location->x + w, location->y);
		glTexCoord2f(1.0f, 0.0f); 
			glVertex2f(location->x + w, location->y + h);
		glTexCoord2f(0.0f, 0.0f); 
			glVertex2f(location->x    , location->y + h);
	glEnd();
	
	/* Bad things happen if we delete the texture before it finishes */
	glFinish();
	
	/* return the deltas in the unused w,h part of the rect */
	location->w = initial->w;
	location->h = initial->h;
	
	/* Clean up */
	SDL_FreeSurface(initial);
	SDL_FreeSurface(intermediary);
	glDeleteTextures(1, &texture);
}


bool writeTextToSurface(string renderStr, SDL_Surface *screen, int x ,int y , TTF_Font *font, SDL_Color text_color) {
		bool error =false;
		
		char renderChar[400];
		
		// convert str stream to char
		strcpy(renderChar,renderStr.c_str());
		
		SDL_Rect DestR;
		
		
		SDL_Surface *initial;
		initial = TTF_RenderText_Blended(font, renderChar, text_color);
		
		
		DestR.x = x;
		DestR.y = HEIGHT-initial->h-y;
		
		
	
		
		
		/*	
		// Write text to surface
		SDL_Surface *text;
		  
		  
		text = TTF_RenderText_Solid(font,	renderChar,	text_color);
		
		if (text == NULL)
		{
			cerr << "TTF_RenderText_Solid() Failed: " << TTF_GetError() << endl;
			TTF_Quit();
			SDL_Quit();
			error =true;
		}
		
		    // Apply the text to the display
		if (SDL_BlitSurface(text, NULL, screen, &DestR) != 0)
			{
			cerr << "SDL_BlitSurface() Failed: " << SDL_GetError() << endl;
			error=true;
		}
	
		// SDL_FreeSurface(text);
		*/
		
		glEnable2D();
		SDL_GL_RenderText(renderChar, font, text_color, &DestR);
		glDisable2D();
		
		SDL_FreeSurface(initial);
		
		return error;
	
}

#endif

void analysisData(int t, list <CVortex> vorticesList, list<CVortex> pinsList, list<CVortex> delVortexList, list<CDelLine> delLinesList, string &returnStr, double Phi) {
	
	
	
	static vector<Ca0Bin> vxBinVectorSum (numdensityBins);
	static vector<Ca0Bin> densityBinVectorSum (numdensityBins);
	
	
	ostringstream oss;
	oss.precision(5);
	oss.str("");
	
	vector<Ca0Bin> vxBinVector (numdensityBins);
	
		
		vector<Ca0Bin> densityBinVector (numdensityBins);
		int bin;
		
			
				// calculate the v(x) velocity profile  ( using Ca0Bin class)
			for (list<CVortex>::iterator p = vorticesList.begin();
					p!=vorticesList.end();++p) {
				if (p->get_x() >= 0 && p->get_x() <= sourceWidth+channelWidth+sinkWidth 
						&& p->get_y() >= 0 && p->get_y() <= channelHeight) {
					//cout << "bin" << endl;
					
					
					bin = (int)floor((p->get_x()+binsize/2.0)* (numdensityBins/(urectx1-urectx0)));
					//vxbin = bin;   // for vxfile
					
					 
					 
					vxBinVector[bin].set_bin(p->get_velx());  //for vxfile
					
			}
		}
				vxfile << "timestep: " << t << endl;
				int bincount=0;
			
				for (vector<Ca0Bin>::iterator p = vxBinVector.begin();
						p!= vxBinVector.end(); ++p ) {
				  
					vxfile<< setw(15) << binsize/2.0+urectx0+bincount*binsize  << setw(15) << p->get_a0() << endl;
					bincount++;
					
					
					
  				vxBinVectorSum[p-vxBinVector.begin()].set_bin(p->get_a0());
					
  
	
						
				}
				vxfile  << endl;
		
			
			// CALCULATION && FILE OUTPUT OF A0 USING AVERAGED BIN METHOD
			if (0==t%triangulationInterval) {
			for (list<CVortex>::iterator p = delVortexList.begin();
					p!=delVortexList.end();++p) {
				if (p->get_x() >= 0 && p->get_x() <= sourceWidth+channelWidth+sinkWidth 
						&& p->get_y() >= 0 && p->get_y() <= channelHeight) {
					//cout << "bin" << endl;
					
					bin = (int)floor((p->get_x()+binsize/2.0)* (numdensityBins/(urectx1-urectx0)));
					
					//vxbin = bin;   // for vxfile
					densityBinVector[bin].set_bin(p->get_a0());
					//vxBinVector[bin].set_bin(p->get_velx());  //for vxfile
					
				}
			}
				
			
			bincount=0;
			
				densityfile << "timestep: " << t << endl;
				
			
				for (vector<Ca0Bin>::iterator p = densityBinVector.begin();
						p!= densityBinVector.end(); ++p ) {
				  
					densityfile<< setw(15) << binsize/2.0+urectx0+bincount*binsize  << setw(15) << p->get_a0() << setw(15) << p->get_count() << endl;
					bincount++;
				  densityBinVectorSum[p-densityBinVector.begin()].set_bin(p->get_a0());
				
				}
				densityfile  << endl;
			
			  
				// OUTPUT A0 DATA FOR EACH VORTEX SEPARATELY
				avlinelengthfile << "Time: " << t << endl;
			
				for (list<CVortex>::iterator p = delVortexList.begin();
						p!= delVortexList.end(); ++p ) {
							if ( p->get_y() >= bulky0 &&
							p->get_y() <= bulky1) {
								avlinelengthfile<< setw(15) << p->get_x() << setw(15) << p->get_a0() << endl;
							}
				
				
				}
			
				avlinelengthfile << endl;
				
				for (vector<Ca0Bin>::iterator p =densityBinVector.begin();
					p!= densityBinVector.end(); ++p ) {
					double Beff=2*Phi/(sqrt((double)3)*p->get_a0()*p->get_a0());		
					if (boost::math::isinf((double)(Beff))) oss << setw(10) << "-";
					else oss<<  setw(10) << Beff;
			
				
					//cout << Beff << endl;
			
				
				}
	

			returnStr=oss.str();
				
				
				
				}
			
			
		
				
			
			if (t==numTimesteps) {
				vxfile << "Time average of binned vx values" << endl;
				bincount=0;
				for (vector<Ca0Bin>::iterator p = vxBinVectorSum.begin();
						p!= vxBinVectorSum.end(); ++p ) {
				  
					vxfile<< setw(15) << binsize/2.0+urectx0+bincount*binsize  << setw(15) << p->get_a0() << setw(15) << p->get_count() <<  endl;
					bincount++;
					
				}
				vxfile  << endl;
				
				densityfile << "Time average of binned density values" << endl;
				bincount=0;
				for (vector<Ca0Bin>::iterator p = densityBinVectorSum.begin();
						p!= densityBinVectorSum.end(); ++p ) {
				  
					densityfile<< setw(15) << binsize/2.0+urectx0+bincount*binsize  << setw(15) << p->get_a0() << setw(15) << p->get_count() <<  endl;
					bincount++;
					
				}
				densityfile  << endl;
				
				// output in csv format
				
				bincount=0;
				for (vector<Ca0Bin>::iterator p = vxBinVectorSum.begin();
						p!= vxBinVectorSum.end(); ++p ) {
				  
					averagedvxdata<<  binsize/2.0+urectx0+bincount*binsize  << "," << p->get_a0() << endl;
					bincount++;
					
				}
				
				bincount=0;
				for (vector<Ca0Bin>::iterator p = densityBinVectorSum.begin();
						p!= densityBinVectorSum.end(); ++p ) {
				  
					averageda0data<<  binsize/2.0+urectx0+bincount*binsize  << "," << p->get_a0() << endl;
					bincount++;
					
				}
				
				
				
				
				
				
				
			}
			
			
			/*
			vxfile << "  Time:"<< t << endl;
			for (list<CVortex>::iterator p = vorticesList.begin();
					p!=vorticesList.end();++p) {
				if (p->get_x() >= urectx0 && p->get_x() <= urectx1) {
						vxfile << setw(15) << p->get_x()  << setw(15) << p->get_y() << setw(15) << p->get_velx() << setw(15) << p->get_vely() <<endl;;
					
					
				}
			}
			vxfile << endl;
			*/
			
			/*
			//ux
			for (list<CVortex>::iterator p = vorticesList.begin();
					p!=vorticesList.end();++p) {
				if (p->get_x() >= urectx0 && p->get_x() <= urectx1
						&& p->get_y() >= urecty0 && p->get_y() <= urecty1) {
					//cout << "bin" << endl;
					double adjustedx = p->get_x() - urectx0;
					vbin = (int)floor(adjustedx * (numa0Bins/(urectx1-urectx0)));
					vBinVector[vbin].set_bin(p->get_velx());
		
				}
			}
			vfile << setw(15) << t;
			for (vector<Ca0Bin>::iterator p = vBinVector.begin();
					p!= vBinVector.end(); ++p ) {
				vfile << setw(15) << p->get_a0();
				
			}
			vfile << endl;
			
			//uy
			for (list<CVortex>::iterator p = vorticesList.begin();
					p!=vorticesList.end();++p) {
				if (p->get_x() >= urectx0+(urectx1-urectx0)/2.0-latticeSpacing && p->get_x() <= urectx0+(urectx1-urectx0)/2.0+latticeSpacing
						&& p->get_y() >= urecty0 && p->get_y() <= urecty1) {
					//cout << "bin" << endl;
					double adjustedy = p->get_y() - urecty0;
					vybin = (int)floor(adjustedy * (numvyBins/(urecty1-urecty0)));
					vyBinVector[vybin].set_bin(p->get_velx());
		
				}
			}
			vyfile << setw(15) << t;
			for (vector<Ca0Bin>::iterator p = vyBinVector.begin();
					p!= vyBinVector.end(); ++p ) {
				vyfile << setw(15) << p->get_a0();
				
			}
			vyfile << endl;
			
			
			*/
			
			// output a0 && bo data to file
			
			// This result comes from actual line lengths plotted at the mid point of the
		    // line
		    
		  
		  if (0==t%triangulationInterval) {  
			a0file << "  time:" << t << " " << endl;
				b0file << "  time:" << t << " " << endl;
				 
				for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
						double midy = (p->get_y1() + p->get_y2())/2.0;
						
						if ( midy >= bulky0 && midy <= bulky1) {
						//cout << "a0scaled: " << graphy1-a0scaling*p->get_a0() << endl;
						
						// find line length
						double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
						double rowspacing= fabs(p->get_y2()-p->get_y1());												
						
						//find mid point
											
					
						
							// draw horizontal lines as different colours to diagonal
							if (fabs(p->get_y2()-p->get_y1())< linelength/4.0) // horizontal
							{ 
								a0file<< setw(15) << (p->get_x1() +p->get_x2())/2.0 << setw(15) << linelength << "   horizontal a0"<< endl;
								
							}
							else {
								a0file<< setw(15) << (p->get_x1() +p->get_x2())/2.0 << setw(15) << linelength << "  diagonal a0" << endl;
							
								
							}
							
							if (rowspacing > linelength/4.0) {  // calculate b0 from a0 lines only 
								b0file << setw(15) << (p->get_x1() +p->get_x2())/2.0 << setw(15) << rowspacing << endl;
								
							}
					    
					    					
					
					
				
						}
				}
				
				a0file << endl;
			
			
	  }
	
	
	

			
			
			
			
			// output dislocation data
			
			
		if (t%25==0) { // set to 50 after test
			vector<CVortex> unsortedVector;
			disfile << setw(10) << t;
		for (list<CVortex>::iterator p = delVortexList.begin();
				p != delVortexList.end(); ++p) {
			   
			if (p->get_x() >= dislocationx0 && p->get_x() <= dislocationx1
						&& p->get_y() >= dislocationy0 && p->get_y() <= dislocationy1) {
			
				if (p->get_coordNum()== 5) {
							// add to unsorted vector
						CVortex newVector;
						newVector.set_pos(p->get_x(),p->get_y());
						unsortedVector.push_back(newVector);
						
							
					
						}
		 	
				}
			}
		//sort
		sort (unsortedVector.begin(),unsortedVector.end(),CVortexSort);
	
		//output sortedvector
		for (vector<CVortex>::iterator p=unsortedVector.begin();
			p!=unsortedVector.end(); ++p ) {
			disfile << setw(15) << p->get_x() << setw(15) << p->get_y();
			
			
		}
		disfile << endl;
		}
		
			
		
		
			
					
			// output dislocation path data
	  /*
		int pathnum=1;
		for (vector< vector<CTemporalCoord> >::iterator iter_ii = pathsVector.begin();
				iter_ii != pathsVector.end(); ++iter_ii) {
			if (iter_ii->size()>500) {		
				disfile << setw(3) << "path: " << setw(5) << pathnum << endl; 
				for (vector<CTemporalCoord>::iterator iter_jj= (*iter_ii).begin();
						iter_jj!=(*iter_ii).end(); ++iter_jj) {
					//disfile << "  (" << setw(7)<< iter_jj->get_x() << "," << setw(7) << iter_jj->get_y() << "," << iter_jj->get_time() <<")";	
					disfile << setw(7)<< iter_jj->get_x() <<  setw(15) << iter_jj->get_y() << setw(15) << iter_jj->get_time() << endl;	
				}	
				disfile << endl;
				pathnum++;
			}
		}
		*/
		
		
		
		// output velocity field data
		
	
		if (outputvelfield==true) {
			velfieldfile << "Time: " << t << endl;
			velfieldfile.setf(ios::fixed);
			velfieldfile.precision(20);
		
			velfieldfile << "{";
			for(list<CVortex>::iterator p = vorticesList.begin();
					p != vorticesList.end(); ++p) {
				velfieldfile << "{{" << p->get_x() << "," << p->get_y() << "},{"
							<< p->get_velx() << "," << p->get_vely() << "}},";
							
			}
			velfieldfile << "}"<< endl;
			
		

			
			}
		
		
			if (0==t%25) { 
			// Output Jyy as list for mathematica
			//Jyydata << "{";
			for (list<CVortex>::iterator p = delVortexList.begin();
					p!=delVortexList.end(); ++p) {
				
				Jyydata << fixed << "{" << p->get_x()/a0 <<"," << p->get_y()/a0 << "," <<p->get_Jyy()*10e11 <<"," << p->get_coordNum() << "}";
				Jxxdata << fixed << "{" << p->get_x()/a0 <<"," << p->get_y()/a0 << "," <<p->get_Jxx()*10e11 <<"," << p->get_coordNum() << "}";
				Jxydata << fixed << "{" << p->get_x()/a0 <<"," << p->get_y()/a0 << "," <<p->get_Jxy()*10e11 <<"," << p->get_coordNum() << "}";
				Jyxdata << fixed << "{" << p->get_x()/a0 <<"," << p->get_y()/a0 << "," <<p->get_Jyx()*10e11 <<"," << p->get_coordNum() << "}";
				
				//if (p==vorticesList.end()) Jyydata << "end";
				//if (p!=vorticesList.end()) Jyydata << ",";
				if ( std::distance(p,delVortexList.end()) != 1 ) {
					Jyydata << "  ";
					Jxxdata << "  ";
					Jxydata << "  ";
					Jyxdata << "  ";
				}
				
			}
			//Jyydata << "}";
			
			
			if (t!=numTimesteps) {
				Jyydata << endl; 
				Jxxdata << endl; 
				Jxydata << endl; 
				Jyxdata << endl; 
			}
			
			
			}
		  
		

		
}

void calculateBinnedBfield(list<CDelLine>& delLinesList, double Phi, string &returnStr1, string &returnStr2, int t) {
	// This routine calculates the effective B field of channel && source && sink using
	// the binned method (does not include the wings)
	
	ostringstream oss;
	oss.str("");
	oss.precision(5);
	vector<Ca0Bin> BfieldBinVector (numdensityBins);
	int bin;
	
	static vector<Ca0Bin> BfieldBinVectorSum (numvxBins);
	int bincount=0;
		
	if (0==t%triangulationInterval) {
	
	for (list<CDelLine>::iterator p = delLinesList.begin();
				p!=delLinesList.end(); ++p) {
		double midy = (p->get_y1() + p->get_y2())/2.0;
		double midx = (p->get_x1() + p->get_x2())/2.0;
						
		if (midy>0 && midy<channelHeight) {
			//continue;
			double linelength=sqrt((double) (p->get_x1()-p->get_x2())*(p->get_x1()-p->get_x2())
																		+ (p->get_y1()-p->get_y2())*(p->get_y1()-p->get_y2()));
			
			bin = (int)floor((midx+binsize/2.0)* (numdensityBins/(urectx1-urectx0)));
					//vxbin = bin;   // for vxfile
			BfieldBinVector[bin].set_bin(linelength);
					
	
		}

	}
	
			
	for (vector<Ca0Bin>::iterator p = BfieldBinVector.begin();
			p!= BfieldBinVector.end(); ++p ) {

			oss<< setw(10) << binsize/2.0+urectx0+bincount*binsize;
			bincount++;
			
				
	}
	returnStr1=oss.str();
	
	oss.str("");
	
	bincount=0;
		for (vector<Ca0Bin>::iterator p = BfieldBinVector.begin();
			p!= BfieldBinVector.end(); ++p ) {
			double Beff=2*Phi/(sqrt((double)3)*p->get_a0()*p->get_a0());		
			if (boost::math::isinf((double)(Beff))) oss << setw(10) << "-";
			else oss<<  setw(10) << Beff;
			bincount++;
			BfieldBinVectorSum[p-BfieldBinVector.begin()].set_bin(Beff);
		  
				
	}
	

	returnStr2=oss.str();
	}
	if (t==numTimesteps) {
	// output in csv format
				
		bincount=0;
		for (vector<Ca0Bin>::iterator p = BfieldBinVectorSum.begin();
				p!= BfieldBinVectorSum.end(); ++p ) {
				  
			averagedBfielddata<<  binsize/2.0+urectx0+bincount*binsize  << "," << p->get_a0() << endl;
			bincount++;
					
		}
	}
	
	
	
	
	
}

#endif
