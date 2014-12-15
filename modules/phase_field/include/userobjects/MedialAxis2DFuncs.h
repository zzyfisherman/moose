#ifndef MEDIALAXIS2DFUNCS_H
#define MEDIALAXIS2DFUNCS_H

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/elem.h"
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

static Real minSegmLength(const std::vector<Point> &segm1, const std::vector<Point> &segm2)
{
  Real h1 = std::sqrt((segm1[0]-segm1[1]).size_sq());
  Real h2 = std::sqrt((segm2[0]-segm2[1]).size_sq());
  return std::min(h1, h2);
}

static bool twoPointsOverlap(Point p1, Point p2, Real ref_len)
{
  Real tol = 1.0e-8;
  Real dist = std::sqrt((p1-p2).size_sq());
  if ((dist/ref_len) < tol)
    return true;
  else
    return false;
}

static Real isLeft(Point P0, Point P1, Point P2)
{
// see http://geomalgorithms.com/a03-_inclusion.html
// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
  return ((P1(0) - P0(0))*(P2(1) - P0(1)) - (P2(0) -  P0(0))*(P1(1) - P0(1)));
}

static int wn_PnPoly(Point P, std::vector<std::vector<Point> > &poly_segms)
{
// see http://geomalgorithms.com/a03-_inclusion.html
// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
  int wn = 0;  // the  winding number counter

  // loop through all edges of the polygon
  for (int i = 0; i < poly_segms.size(); ++i)
  { 
    Point v0 = poly_segms[i][0];
    Point v1 = poly_segms[i][1];
    if (v0(1) <= P(1))
    {
      if (v1(1)  > P(1))      // an upward crossing
        if (isLeft(v0, v1, P) > 0.0)  // P left of  edge
          ++wn;            // have  a valid up intersect
    }
    else
    {                        // start y > P.y (no test needed)
      if (v1(1)  <= P(1))     // a downward crossing
        if (isLeft(v0, v1, P) < 0.0)  // P right of  edge
          --wn;            // have  a valid down intersect
    }
  }
  return wn;
}

#endif
