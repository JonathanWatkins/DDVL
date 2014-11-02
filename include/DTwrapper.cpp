#include "DTwrapper.hpp"

#include "delaunay.h"
//#include "delaunay.c"


void ComputationalGeometry::DelaunayTriangulation(const std::list<CParticle> &vorticesList_, std::list<CParticle> * delVortexList_, std::list<CDelLine> * delLinesList_)
{
	delVortexList_->clear();
	delLinesList_->clear();
	
	std::list<CLineIDs> lines;
	
	/* Define input points. */

	int numberofpoints = vorticesList_.size();
	if (numberofpoints==0) return;
	
	del_point2d_t* pointlist= new del_point2d_t[numberofpoints];
	
	std::vector<CParticle> vorticesVector;
	std::copy( vorticesList_.begin(), vorticesList_.end(), std::back_inserter( vorticesVector ) );
	
	int count=0;
	for (std::vector<CParticle>::iterator p = vorticesVector.begin();
			p!=vorticesVector.end(); ++p )
	{
		p->set_coord_num(0);
		pointlist[count].x = p->get_x();
		pointlist[count].y = p->get_y();
		count++;
	}
	
	
	int offset = 0;
	delaunay2d_t*	res = delaunay2d_from(pointlist,numberofpoints, NULL);
	
	int num_faces = res->num_faces;
	
	std::vector<int> poly;
	for(int i = 0; i < num_faces; i++ )
	{
	poly.clear();
	
		
		int num_verts = res->faces[offset];
		
		offset++;
		for( int j = 0; j < num_verts; j++ )
		{
			int p0 = res->faces[offset + j];
			int p1 = res->faces[offset + (j+1) % num_verts];
			
			poly.push_back(p0);
			
			
		}
		
		offset += num_verts;
	
		
		
		for (int p=0;p<poly.size()-1;p++)
		{
			CLineIDs newline(poly[p],poly[p+1]);
			lines.push_back(newline);
		}
	}
	
	delaunay2d_release(res);
	delete [] pointlist;
	// check delLinesList for dupilcates
	lines.sort();
	lines.unique();
	
	for (std::list<CLineIDs>::iterator p = lines.begin();
		p!=lines.end(); ++p)
	{
		//std::cout << vorticesVector[p->id1].get_x() << " " << vorticesVector[p->id1].get_y() << " " << vorticesVector[p->id2].get_x() << " " << vorticesVector[p->id2].get_y() << std::endl;
		
		CDelLine newDelLine;
		newDelLine.set_points(vorticesVector[p->id1].get_x(),vorticesVector[p->id1].get_y(),vorticesVector[p->id2].get_x(),vorticesVector[p->id2].get_y());
		delLinesList_->push_back(newDelLine);	
	
		vorticesVector[p->id1].coordPlusOne();
		vorticesVector[p->id2].coordPlusOne();
	}
	
	std::copy( vorticesVector.begin(), vorticesVector.end(), std::back_inserter( *delVortexList_ ) );
    
}
