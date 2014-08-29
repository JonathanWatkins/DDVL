#ifndef CSDLGLSCREEN_HPP
#define CSDLGLSCREEN_HPP

#include <SDL.h>
//#include <SDL_gfxPrimitives.h>
#include <SDL_opengl.h>
#include <SDL_ttf.h>
#include <list>

#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CSimulation.hpp"

class CSDLGLscreen {

	SDL_Surface *screen;

	CSimulation* sim;

	double zoom;
	
	double x1Crop;
	double x2Crop;
	
	double viewpointx;
	double viewpointy;
	

	GLfloat drawUnit;
	
	GLfloat drawnVortexSize;
	GLfloat drawnChannelPinSize;
					
	GLfloat GLxoffset;   
	GLfloat GLyoffset;   
	GLfloat GLzoffset;	
	
	int WIDTH;
	int HEIGHT;
	int BPP;
	int DEPTH;

	double systemLength;
	double systemWidth;

	TTF_Font *largefont;
	TTF_Font *smallfont;
	SDL_Color text_color; 

	// gl Lists
	GLuint GLpinsList;
	GLuint GLvortexList;
	GLuint GLdisorderList;
	
	std::vector<CParticle> *burgers_circuit;

public:	
	CSDLGLscreen();
	~CSDLGLscreen(){};
	
	
	int initialiseSDLandGL(CSimulation* sim_);
	
	void drawSystem();

	void SDL_GL_RenderText(char *text_, TTF_Font *font_, SDL_Color color_, SDL_Rect *location_);


private:

	void makeGLLists();
	
	void writeText();
	
	bool writeTextToSurface(std::string renderStr_, GLfloat x_ , GLfloat y_);
	
	void RenderText(const TTF_Font *Font_, const GLubyte& R_, const GLubyte& G_, const GLubyte& B_,
                const double& X_, const double& Y_, const double& Z_,  const std::string& Text_);
	
	void drawCoordinateGrid();
	
	void drawBathEdges();
	
	void drawReboundWalls();
	
	void drawRemoveWalls();
	
	void drawBounceBackWalls();
	
	void drawRemoveWallsWedge();
	
};

#endif
