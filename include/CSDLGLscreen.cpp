#include "CSDLGLscreen.hpp"

#include <list>
#include <vector>
#include <iostream>
#include <SDL.h>
#include <SDL_ttf.h>
#include <sstream>

#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CSimulation.hpp"
#include "CRowCount.hpp"

CSDLGLscreen::CSDLGLscreen()
{
	
	drawUnit=0;
	
	GLxoffset =0;   
	GLyoffset =0;   
	GLzoffset = 0;	
	
	WIDTH=900;
	HEIGHT=900;
	BPP=4;
	DEPTH=32;
	
	zoom=1;
	viewpointx=0;
	viewpointy=0;
	
}

int CSDLGLscreen::initialiseSDLandGL(CSimulation* sim_)
{
	std::cout << "Initialise OpenGL/SDL..." << std::endl;
	
	sim=sim_;
	
	//initialise Screen
	if (SDL_Init(SDL_INIT_VIDEO) < 0 ) return 1;
	
	//SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);

	/*SDL_GL_SetAttribute(SDL_GL_RED_SIZE,            8);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,          8);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,           8);
	SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,          8);

	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,          16);
	SDL_GL_SetAttribute(SDL_GL_BUFFER_SIZE,            32);

	SDL_GL_SetAttribute(SDL_GL_ACCUM_RED_SIZE,        8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_GREEN_SIZE,    8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_BLUE_SIZE,        8);
	SDL_GL_SetAttribute(SDL_GL_ACCUM_ALPHA_SIZE,    8);

	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS,  1);

	SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,  2);

	*/

	//if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, DEPTH, SDL_FULLSCREEN|SDL_HWSURFACE)))
	if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, 32, SDL_HWSURFACE | SDL_GL_DOUBLEBUFFER | SDL_OPENGL)))
	{
		SDL_Quit();
		return 1;
	}
	
	glClearColor(1.0f, 1.0f, 1.0f, 0);
	glClearDepth(1.0f);
	
	glViewport(0, 0, WIDTH, HEIGHT);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	glOrtho(0, WIDTH, HEIGHT, 0, 1, -1);
	
	glMatrixMode(GL_MODELVIEW);
	
	glEnable(GL_TEXTURE_2D);
	
	glLoadIdentity();
	
	/* ////////////////////////////////////////
	 
	glEnable(GL_MULTISAMPLE);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST );
	
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	
	glClearColor(255, 255, 255, 0);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (640.0/480.0), 0.1, 100.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	*////////////////////////////
	
	
	//glTranslatef(-0.9f,-0.9f,0); 
	
	
	/*glMatrixMode(GL_MODELVIEW);
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
	*/
	
	
	SDL_WM_SetCaption( "meshworks", 0 );
	
	
	// Initialize SDL_ttf library
	text_color.r=255;
	text_color.g=0;
	text_color.b=0;
		
	if (TTF_Init() != 0)
	{
		std::cerr << "TTF_Init() Failed: " << TTF_GetError() << std::endl;
		SDL_Quit();
		sim->end_simulation();
		return 1;
	}
	
	// Load a font
  
	largefont = TTF_OpenFont("FreeSans.ttf", 32);
	if (largefont == NULL)
	{
		std::cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << std::endl;
		TTF_Quit();
		SDL_Quit();
		sim->end_simulation();
		return 1;
	}
	
	smallfont = TTF_OpenFont("FreeSans.ttf", 20);
	if (smallfont == NULL)
	{
		std::cerr << "TTF_OpenFont() Failed: " << TTF_GetError() << std::endl;
		TTF_Quit();
		SDL_Quit();
		sim->end_simulation();
		return 1;
	}
	
	
	
	// set drawing variables
	/*
	 *  We want a border of 0.1 on the left and right of the system.
	 * 	We want the vertical to be centered along the y=0 line. 
	 * 	
	 * 	This is the coordinate system
	 * 
	 *    -1,1================1,1
	 *		|                  |
	 *		|                  |
	 *		|                  |
	 *		|                  |
	 *		|                  |
	 *	 -1,-1===============1,-1
	 * 
	 *  x offset
	 * 
	 *  
	 */
	
	
	drawUnit=1.8/sim->get_systemLength();
	
	GLxoffset =-0.9-sim->get_viewpoint().get_x()*drawUnit;    
	GLyoffset =-sim->get_systemWidth()*drawUnit/2.0-sim->get_viewpoint().get_y()*drawUnit;   
	GLzoffset =0.0;	
          
	drawnVortexSize = drawUnit*sim->get_vortexSize();
	drawnChannelPinSize = drawUnit*sim->get_vortexSize()/2.0;
	
	zoom = sim->get_zoom();
	
	x1Crop=(sim->get_viewpoint().get_x()+sim->get_systemLength())/2.0-1.2/zoom*sim->get_systemLength()/2.0;
	x2Crop=(sim->get_viewpoint().get_x()+sim->get_systemLength())/2.0+1.2/zoom*sim->get_systemLength()/2.0;
	
	//std::cout << "x1Crop: " << x1Crop << " x2Crop: " << x2Crop << std::endl;
	//std::cout << "FirstPin (x,y): (" << sim->get_firstPin().get_x() << ", " << sim->get_firstPin().get_y() << ")\n"; 
	std::cout << "   GL scalings\n"
			  << "   ----------------\n"
			  << "   Draw Unit: " << drawUnit << std::endl
			  << "   offsets (x,y,z): (" << GLxoffset << ", " << GLyoffset << ", " << GLzoffset << ")\n"
			  << "   Drawn vortex size: " << drawnVortexSize << std::endl;
	std::cout << "   OpenGL/SDL initialised.\n\n";
			 
		
	// make lists
	
	makeGLLists();
	
	
	return 1;
	
}


void CSDLGLscreen::drawSystem()
{
	
	if (zoom != sim->get_zoom() || viewpointx!= sim->get_viewpoint().get_x() || viewpointy!= sim->get_viewpoint().get_y())
	{
		zoom = sim->get_zoom();
		viewpointx = sim->get_viewpoint().get_x();
		viewpointy = sim->get_viewpoint().get_y();
		
		drawUnit=zoom*1.8/sim->get_systemLength();
		GLxoffset =-0.9*zoom-sim->get_viewpoint().get_x()*drawUnit;    
		GLyoffset =-sim->get_systemWidth()*drawUnit/2.0-sim->get_viewpoint().get_y()*drawUnit;   
		GLzoffset =0.0;	
    
		x1Crop=(sim->get_viewpoint().get_x()+sim->get_systemLength())/2.0-1.2/zoom*sim->get_systemLength()/2.0;
		x2Crop=(sim->get_viewpoint().get_x()+sim->get_systemLength())/2.0+1.2/zoom*sim->get_systemLength()/2.0;
    
    //std::cout << "x1Crop: " << x1Crop << " x2Crop: " << x2Crop << std::endl;
	 
          
		drawnVortexSize = drawUnit*sim->get_vortexSize();
		drawnChannelPinSize = drawUnit*sim->get_vortexSize()/2.0;
		
		// make lists
	
		makeGLLists();
		
		//std::cout << "Zoom: " << zoom << std::endl;
	
	}
	
	
	
	//std::cout << "Start draw" << std::endl;
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	///// SDL OUTPUT TO SCREEN ////////////////
	//SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));	
	
	
	glLoadIdentity();
	
	
	glCallList(GLpinsList);
	
	glCallList(GLdisorderList);
	
	// draw nn lines
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	drawCoordinateGrid();
	
	
		
	//if periodic draw boundary of periodic cell
	if (sim->get_geometry()==periodic)
	{
		glColor3f(0.0f,0.0f,1.0f);
		glLineWidth(2); 
		
		GLfloat x1 = GLxoffset+drawUnit*0;
		GLfloat y1 = GLyoffset+drawUnit*0;
		
		GLfloat x2 = GLxoffset+drawUnit*0;
		GLfloat y2 = GLyoffset+drawUnit*sim->get_channelWidth();
		
		GLfloat x3 = GLxoffset+drawUnit*sim->get_channelLength();
		GLfloat y3 = GLyoffset+drawUnit*sim->get_channelWidth();
		
		GLfloat x4 = GLxoffset+drawUnit*sim->get_channelLength();
		GLfloat y4 = GLyoffset+drawUnit*0;
		
		
		GLfloat z = GLzoffset;
		
		
		glBegin(GL_LINES); 
		glVertex3f(x1, y1,z); 
		glVertex3f(x2, y2,z); 
		glEnd(); 
		
		glBegin(GL_LINES); 
		glVertex3f( x2, y2,z); 
		glVertex3f(x3, y3,z); 
		glEnd(); 
		
		glBegin(GL_LINES); 
		glVertex3f( x3, y3,z); 
		glVertex3f(x4, y4,z); 
		glEnd(); 
		
		glBegin(GL_LINES); 
		glVertex3f( x4, y4,z); 
		glVertex3f(x1, y1,z); 
		glEnd(); 
		
		
	}			
				
	
	for (std::list<CDelLine>::iterator p = sim->get_delLinesList()->begin();
			p!=sim->get_delLinesList()->end(); ++p) {
		
		GLfloat x1 = GLxoffset+drawUnit*p->get_x1();
		GLfloat x2 = GLxoffset+drawUnit*p->get_x2();
		GLfloat y1 = GLyoffset+drawUnit*p->get_y1();
		GLfloat y2 = GLyoffset+drawUnit*p->get_y2();
		GLfloat z = GLzoffset;
		
		if (sim->get_geometry()==periodic &&  
			(		 
			    p->get_x1()>sim->get_channelLength() && p->get_x2()>sim->get_channelLength()
				|| p->get_x1()<0 && p->get_x2()<0
				|| p->get_y1()>sim->get_channelWidth() && p->get_y2()>sim->get_channelWidth()
				|| p->get_y1()<0 && p->get_y2()<0
			)
		)
		continue;
		
		
		glColor3f(0.9f,0.9f,0.9f);
		glLineWidth(2); 
		glBegin(GL_LINES); 
		glVertex3f( x1, y1,z); 
		glVertex3f(x2, y2,z); 
		glEnd(); 
		
	
	}
	
	
	
	
	//draw vortices
	
	for (std::list<CParticle>::iterator p = sim->get_delVortexList()->begin();
			p != sim->get_delVortexList()->end(); ++p) {
		//std::cout << "#: " <<  p->get_coord_num() << std::endl;
		if (true==p->get_ghost())
			continue;
		
		double x = GLxoffset+drawUnit*p->get_x();
		double y2 = GLyoffset+drawUnit*p->get_y();
		double z2 =GLzoffset;
		
		// if in bubble draw large circle around vortex
		if (p->get_in_bubble()==true)
		{	
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
		
			glTranslatef(x,y2,z2); 
			
			glColor3f(1.0f,1.0f,0.0f);
			
			GLUquadricObj *quadObj = gluNewQuadric();
			gluSphere(quadObj,2*drawnVortexSize ,20 , 40); 
	
			
			gluDeleteQuadric(quadObj); 
		}
		
		
		if (p->get_coord_num()==6 && !sim->get_drawSixFold() && !sim->get_showParticleTracker())
			continue;
	
		if (p->get_x()<x1Crop || p->get_x() > x2Crop)
			continue;
					  
		
		
		
	
		if(sim->get_showParticleTracker()==true)
		{
			if(p->get_Tracked()==true) glColor3f(0.0f,0.0f,1.0f);
			else glColor3f(0.9f,0.9f,0.9f);
		}
		else
		{
			if (p->get_coord_num()==6)
			{
				//yellow sphere
				glColor3f(1.0f,1.0f,0.0f);
			}
			else if (p->get_coord_num()==5)
			{
				glColor3f(1.0f,0.0f,0.0f);
			}
			else if (p->get_coord_num()==7) 
			{
				// blue sphere
				glColor3f(.0f,0.0f,1.0f);
			}
			else
			{
				//black sphere
				glColor3f(0.0f,0.0f,0.0f);
			}
		}
		//GLUquadricObj *quadObj = gluNewQuadric();
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glTranslatef(x,y2,z2); 
		
		//gluSphere(quadObj,drawnVortexSize ,20 , 40); 
		//gluDeleteQuadric(quadObj); */
		
		glCallList(GLvortexList);

		
		
	
	}
	
	if (sim->get_geometry()==channel || sim->get_geometry()==tube || sim->get_geometry()==wedge || sim->get_geometry()==BSCCO) drawBathEdges();
	
	if (sim->get_geometry()==channel || sim->get_geometry()==wedge || sim->get_geometry()==BSCCO) drawReboundWalls();
	
	if (sim->get_geometry()==channel || sim->get_geometry()==wedge || sim->get_geometry()==BSCCO) drawRemoveWalls();
	
	if (sim->get_geometry()==wedge) drawRemoveWallsWedge();
		
	if (sim->get_applyBounceBack()==true) drawBounceBackWalls();
	
	//drawBurgersCircuitCenter();
	
	writeText();
	
		
	glMatrixMode(GL_PROJECTION);
	
	glLoadIdentity(); // Reset Matrix
	//glTranslated(-1,0,0);
	
	SDL_GL_SwapBuffers();	
	
	SDL_FreeSurface(screen);	  
	
	//std::cout << "Drawn" << std::endl;
}	


void CSDLGLscreen::makeGLLists()
{
	
	// All pins list
	
	GLpinsList = glGenLists(1);
      
	glNewList( GLpinsList, GL_COMPILE );
	
	glColor3f(.9f,0.9f,0.9f);
	
	
	for (std::list<CParticle>::iterator p = sim->get_pinsList()->begin();
			p != sim->get_pinsList()->end(); ++p) {
		//if (p->get_x() < sourceWidth || p->get_x() > sourceWidth+channelWidth )
		//continue;
		if (p->get_x()<x1Crop || p->get_x() > x2Crop)
			continue;
		
		
		double x = GLxoffset+drawUnit*p->get_x();
		double y2 = GLyoffset+drawUnit*p->get_y();
		double z2 =GLzoffset;
		
		
		
		GLUquadricObj *quadObj = gluNewQuadric();
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glTranslatef(x,y2,z2); 
		
		gluSphere(quadObj,drawnVortexSize ,20 , 40); 
		gluDeleteQuadric(quadObj); 
		
	
	}
	
	glEndList();
	
		
	// Single vortex list
	
	GLvortexList = glGenLists(1);
    
  glNewList( GLvortexList, GL_COMPILE );
							
	GLUquadricObj *quadObj = gluNewQuadric();
	gluSphere(quadObj,drawnVortexSize ,20 , 40); 
	
			
	gluDeleteQuadric(quadObj); 
	
	glEndList();
	
	// Channel Disorder List
	
	GLdisorderList = glGenLists(1);
      
	glNewList( GLdisorderList, GL_COMPILE );
	
	glColor3f(.2f,0.2f,0.2f);
	
	for (std::list<CParticle>::iterator p = sim->get_disorderList()->begin();
			p != sim->get_disorderList()->end(); ++p) {
		//if (p->get_x() < sourceWidth || p->get_x() > sourceWidth+channelWidth )
		//continue;
		
		double x = GLxoffset+drawUnit*p->get_x();
		double y2 = GLyoffset+drawUnit*p->get_y();
		double z2 =GLzoffset;
		
		
		
		GLUquadricObj *quadObj = gluNewQuadric();
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glTranslatef(x,y2,z2); 
		
		gluSphere(quadObj,drawnChannelPinSize ,20 , 40); 
		gluDeleteQuadric(quadObj); 
		
	
	}
	
	glEndList();
	
}


void CSDLGLscreen::SDL_GL_RenderText(char *text_, TTF_Font *font_, SDL_Color color_, SDL_Rect *location_)
{
	SDL_Surface *initial;
	SDL_Surface *intermediary;
	
	int w,h;
	GLuint texture;
	//SDL_FillRect(initial, NULL, SDL_MapRGB(initial->format, 255, 255, 255));	
	/* Use SDL_TTF to render our text */
	initial = TTF_RenderText_Blended(font_, text_, color_);
	
	/* Convert the rendered text to a known format */
	//w = nextpoweroftwo(initial->w);
	//h = nextpoweroftwo(initial->h);
	w=initial->w;
	h=initial->h;
	//std::cout << w << std::endl;
	//std::cout << h << std::endl;
	
	intermediary = SDL_CreateRGBSurface(0, w, h, 32, 
			0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
	SDL_FillRect(intermediary, NULL, SDL_MapRGB(intermediary->format, 0, 255, 255));	
	
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
	glColor3f(0.00f, 1.0f, 1.0f);
	
	/* Draw a quad at location */
	glBegin(GL_QUADS);
		/* Recall that the origin is in the lower-left corner
		   That is why the TexCoords specify different corners
		   than the Vertex coors seem to. */
		glTexCoord2f(0.0f, 1.0f); 
			glVertex2f(location_->x    , location_->y);
		glTexCoord2f(1.0f, 1.0f); 
			glVertex2f(location_->x + w, location_->y);
		glTexCoord2f(1.0f, 0.0f); 
			glVertex2f(location_->x + w, location_->y + h);
		glTexCoord2f(0.0f, 0.0f); 
			glVertex2f(location_->x    , location_->y + h);
	glEnd();
	
	/* Bad things happen if we delete the texture before it finishes */
	glFinish();
	
	/* return the deltas in the unused w,h part of the rect */
	location_->w = initial->w;
	location_->h = initial->h;
	
	/* Clean up */
	SDL_FreeSurface(initial);
	SDL_FreeSurface(intermediary);
	glDeleteTextures(1, &texture);
}	
	
bool CSDLGLscreen::writeTextToSurface(std::string renderStr_, GLfloat x_ , GLfloat y_)
{
		
		bool error = false;
		
		char renderChar[1000];
		
		// convert str stream to char
		strcpy(renderChar,renderStr_.c_str());
		
		SDL_Rect DestR;
				
		SDL_Surface *initial;
		initial = TTF_RenderText_Blended(smallfont, renderChar, text_color);
		
		
		DestR.x = x_;
		DestR.y = y_;//HEIGHT-initial->h-y_;
		
		
	
		
		
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
		
		//glEnable2D();
		std::cout << renderStr_.c_str() << std::endl;
		
		SDL_GL_RenderText(renderChar, smallfont, text_color, &DestR);
		//glDisable2D();
		
		SDL_FreeSurface(initial);
		
		return error;
	
}

void CSDLGLscreen::writeText()
{
	std::ostringstream oss;

	oss.str("");
	oss << sim->get_t() << "/" << sim->get_simulation_time();
	if (oss.str()!="") RenderText(largefont,255,0,0,-.9,.95, 0,oss.str());


	
	if (sim->get_geometry()==channel || sim->get_geometry()==tube || sim->get_geometry()==wedge || sim->get_geometry()==BSCCO)
	{
		oss.str("");
		oss << sim->get_normaliseSourceStr();
		
		if (oss.str()!="") RenderText(largefont,255,0,0,-.9,.9, 0,oss.str());

		oss.str("");
		oss << sim->get_normaliseSinkStr();
		if (oss.str()!="") RenderText(largefont,255,0,0,-.9,.85, 0,oss.str());
		
	}
	else if (sim->get_geometry()==periodic)
	{
		oss.str("");
		oss << sim->get_BfieldStr();
		if (oss.str()!="") RenderText(smallfont,255,0,0,-.9,.9, 0,oss.str());	
	}
	
	oss.str("");
	oss << sim->get_thermostat() << " TS - Temperature = " << sim->get_temp();
	if (oss.str()!="") RenderText(smallfont,255,0,0,-.9,.70, 0,oss.str());
	
	oss.str("");
	oss << "A = " << sim->get_A() << "  M2av = " << sim->get_M2Average()
	    << "  M2Fullav = " << sim->get_M2FullAverage();
	if (oss.str()!="") RenderText(largefont,255,0,0,-.9,.65, 0,oss.str());

	
	oss.str("");
	oss << sim->get_finishTimeStr();
	if (oss.str()!="") RenderText(smallfont,255,0,0,-.9,-.9, 0,oss.str());

	oss.str("");
	oss << "Time (dt*t) : " << sim->get_time();
	if (oss.str()!="") RenderText(largefont,255,0,0,.5,.95, 0,oss.str());

	oss.str("");
	oss << "frame_force_d = " << sim->get_frame_force_d() << "  frame_force_t = " << sim->get_frame_force_t();
	if (oss.str()!="") RenderText(smallfont,255,0,0,-.3,.6, 0,oss.str());
	
	oss.str("");
	oss << "av force_d = " << sim->get_av_force_d() << "  av force_t = " << sim->get_av_force_t();
	if (oss.str()!="") RenderText(smallfont,255,0,0,-.3,.55, 0,oss.str());
	



}




void CSDLGLscreen::RenderText(const TTF_Font *Font_, const GLubyte& R_, const GLubyte& G_, const GLubyte& B_,
                const double& X_, const double& Y_, const double& Z_,  const std::string& Text_)
{
	
		
	/*Create some variables.*/
	SDL_Color Color = {255, 0, 0};

	SDL_Surface *Message = TTF_RenderText_Blended(const_cast<TTF_Font*>(Font_), Text_.c_str(), Color);
	
	SDL_Surface *blank = SDL_CreateRGBSurface(0, Message->w, Message->h, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);
	
	//SDL_Surface *blank = SDL_CreateRGBSurface(0, Message->w, Message->h, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0);
		
	//SDL_Surface *blank = SDL_CreateRGBSurface(0, Message->w, Message->h, 32, 255, 255, 255, 0);
		
	SDL_FillRect(blank, NULL, SDL_MapRGB(blank->format, 255, 255, 255));
	
	SDL_BlitSurface(Message, NULL, blank, NULL);
	
		
	unsigned Texture = 0;
 
	/*Generate an OpenGL 2D texture from the SDL_Surface*.*/
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &Texture);
	glBindTexture(GL_TEXTURE_2D, Texture);
	
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	
	glColor3f(1.0f,1.0f,1.0f);
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, blank->w, blank->h, 0, GL_RGBA,
	             GL_UNSIGNED_BYTE, blank->pixels);
	
	
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, Message->w, blank->h, 0, GL_RGBA,
	 //            GL_UNSIGNED_BYTE, Message->pixels);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	
	/*Draw this texture on a quad with the given xyz coordinates.*/
	glBegin(GL_QUADS);
		glTexCoord2d(0, 1); glVertex3f(X_, Y_, Z_);
		glTexCoord2d(1, 1); glVertex3f(X_+Message->w/(double)WIDTH, Y_, Z_);
		glTexCoord2d(1, 0); glVertex3f(X_+Message->w/(double)WIDTH, Y_+Message->h/(double)HEIGHT, Z_);
		glTexCoord2d(0, 0); glVertex3f(X_, Y_+Message->h/(double)HEIGHT, Z_);
	glEnd();
 
	/*Clean up.*/
	glDeleteTextures(1, &Texture);
	SDL_FreeSurface(Message);
	SDL_FreeSurface(blank);

	
}

void CSDLGLscreen::drawCoordinateGrid()
{
	/* Draw a grid of lines b0 appart starting from the bottom of the channel
	*  up to the channelwidth+b0.
	* 
	*  Note: the coorinates shown are relative to the bottom of the channel not 0.
	*  The channel is offset -b0/2 due to the etch process
	* 
	*  weff= w+b0
	* 
	*  w=channelWidth and is specified in the job batch file
	* 
	*/ 
	
	
	//if periodic draw boundary of periodic cell
	if (!sim->get_drawCoordinateGrid()) return;

	
	GLfloat x1= GLxoffset+drawUnit*sim->get_viewpoint().get_x();
	GLfloat x2=GLxoffset+drawUnit*(sim->get_systemLength()+sim->get_viewpoint().get_x());
	GLfloat z = GLzoffset;
	
	std::ostringstream oss;
	oss.str("");
	oss << "b0";
	if (oss.str()!="") RenderText(smallfont,255,0,0,x1-drawUnit*2.5*sim->get_a0(),
				GLyoffset+drawUnit*(sim->get_channelWidth()/2.0-sim->get_b0()/2.0)-drawUnit*sim->get_b0()/1.5,0,oss.str());
	
	for (double y = 0-sim->get_b0()/2.0; y<=sim->get_channelWidth()+sim->get_b0()*1.01-sim->get_b0()/2.0; y+=sim->get_b0())
	{
		
		GLfloat y1 = GLyoffset+drawUnit*y;
		
		
		GLfloat y2 = GLyoffset+drawUnit*y;
		
	
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		glColor3f(0.7f,0.7f,1.0f);
		glLineWidth(2); 
	
		glBegin(GL_LINES); 
		glVertex3f(x1, y1,z); 
		glVertex3f(x2, y2,z); 
		glEnd(); 
		
		std::ostringstream oss;
		oss.str("");
		oss << (y+sim->get_b0()/2.0)/sim->get_b0();
		if (oss.str()!="") RenderText(smallfont,255,0,0,x1-drawUnit*sim->get_a0(),y1-drawUnit*sim->get_b0()/1.5,0,oss.str());
	
	}
	
}

void CSDLGLscreen::drawBathEdges()
{
	// source
	GLfloat x1= GLxoffset+drawUnit*sim->get_bathLength();
	
	// sink
	GLfloat x2=GLxoffset+drawUnit*(sim->get_bathLength()+sim->get_channelLength());
	GLfloat z = GLzoffset;
	
	
	GLfloat y1 = GLyoffset+drawUnit*0;
		
		
	GLfloat y2 = GLyoffset+drawUnit*sim->get_channelWidth();
		
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	glColor3f(0.7f,0.7f,1.0f);
	
	glLineWidth(2); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x1, y2,z); 
	glEnd(); 
	
	glBegin(GL_LINES); 
	glVertex3f(x2, y1,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 


}


void CSDLGLscreen::drawReboundWalls()
{

	GLfloat x1 = GLxoffset+drawUnit*sim->get_reboundx0();
	GLfloat x2 = GLxoffset+drawUnit*sim->get_reboundx1();
	GLfloat y1 = GLyoffset+drawUnit*sim->get_reboundy0();
	GLfloat y2 = GLyoffset+drawUnit*sim->get_reboundy1();

	GLfloat z = GLzoffset;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	glColor3f(1.0f,0.0f,0.0f);
	
	glLineWidth(1); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x2, y1,z); 
	glEnd(); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y2,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 

	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x1, y2,z); 
	glEnd(); 

	glBegin(GL_LINES); 
	glVertex3f(x2, y1,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 

	
	
}


void CSDLGLscreen::drawBounceBackWalls()
{

	GLfloat x1 = GLxoffset+drawUnit*sim->get_bouncebackx0();
	GLfloat x2 = GLxoffset+drawUnit*sim->get_bouncebackx1();
	GLfloat y1 = GLyoffset+drawUnit*sim->get_bouncebacky0();
	GLfloat y2 = GLyoffset+drawUnit*sim->get_bouncebacky1();

	GLfloat z = GLzoffset;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	glColor3f(0.0f,0.0f,1.0f);
	
	glLineWidth(1); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x2, y1,z); 
	glEnd(); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y2,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 

	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x1, y2,z); 
	glEnd(); 

	glBegin(GL_LINES); 
	glVertex3f(x2, y1,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 

	
	
}



void CSDLGLscreen::drawRemoveWalls()
{
	// source
	GLfloat x1= GLxoffset+drawUnit*sim->get_bathLength();
	
	// sink
	GLfloat x2=GLxoffset+drawUnit*(sim->get_bathLength()+sim->get_channelLength());
	GLfloat z = GLzoffset;
	
	
	GLfloat y1 = GLyoffset+drawUnit*sim->get_removesourcey0();
		
		
	GLfloat y2 = GLyoffset+drawUnit*sim->get_removesourcey1();
		
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	glColor3f(0.0f,1.0f,0.5f);
	
	glLineWidth(2); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y1,z); 
	glVertex3f(x2, y1,z); 
	glEnd(); 
	
	glBegin(GL_LINES); 
	glVertex3f(x1, y2,z); 
	glVertex3f(x2, y2,z); 
	glEnd(); 

	
}


void CSDLGLscreen::drawRemoveWallsWedge()
{

	GLfloat x0 = GLxoffset+drawUnit*sim->get_removewedgex0();
	GLfloat y0 = GLyoffset+drawUnit*sim->get_removewedgey0();
	GLfloat x1 = GLxoffset+drawUnit*sim->get_removewedgex1();
	GLfloat y1 = GLyoffset+drawUnit*sim->get_removewedgey1();
	GLfloat y2 = GLyoffset+drawUnit*sim->get_removewedgey2();

	GLfloat z = GLzoffset;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
		
	glColor3f(0.0f,1.0f,0.5f);
	
	glLineWidth(2); 
	
	glBegin(GL_LINES); 
	glVertex3f(x0, y0,z); 
	glVertex3f(x1, y1,z); 
	glEnd(); 
	
	glBegin(GL_LINES); 
	glVertex3f(x0, y0,z); 
	glVertex3f(x1, y2,z); 
	glEnd(); 

	
	
}
