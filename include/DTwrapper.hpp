#ifndef DTWRAPPER_HPP
#define DTWRAPPER_HPP

#include "CParticle.hpp"
#include "CDelLine.hpp"
#include "CLineIDs.hpp"

#include <list>
#include <vector>
#include <iterator>

namespace ComputationalGeometry
{
	void DelaunayTriangulation(const std::list<CParticle> &vorticesList_, std::list<CParticle> * delVortexList_, std::list<CDelLine> * delLinesList_);
}

#endif
