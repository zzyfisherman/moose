/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "Isocontour.h"
#include "MedialAxis2DFuncs.h"

Isocontour::Isocontour(Real iso_width, const std::vector<std::vector<Point> > & segms) :
    _iso_width(iso_width),
    _iso_segms(segms)
{
}

Point
Isocontour::getFirstPoint() const
{
  return _iso_segms[0][0];
}

Point
Isocontour::getLastPoint() const
{
  return _iso_segms[_iso_segms.size()-1][1];
}

bool
Isocontour::isClockWise() const
{
  // get strictly closed isocontour
  std::vector<std::vector<Point> > new_segms;
  for (unsigned int i = 0; i < _iso_segms.size(); ++i)
  {
    new_segms.push_back(_iso_segms[i]);
    int iplus1(i < _iso_segms.size()-1 ? i+1 : 0);
    Real ref_len = minSegmLength(_iso_segms[i], _iso_segms[iplus1]);
    Point this_p1 = _iso_segms[i][1];
    Point next_p0 = _iso_segms[iplus1][0];
    if (!twoPointsOverlap(this_p1, next_p0, ref_len))
    {
      std::vector<Point> temp_segm(2, Point(0.0,0.0,0.0));
      temp_segm[0] = this_p1;
      temp_segm[1] = next_p0;
      new_segms.push_back(temp_segm);
    }
  }

  // check the enclosed area
  Real area = 0.0;
  for (unsigned int j = 0; j < new_segms.size(); ++j)
    area += (new_segms[j][1](0) - new_segms[j][0](0))*(new_segms[j][1](1) + new_segms[j][0](1));
  if (area > 0.0)
    return true;
  else
    return false;
}

bool
Isocontour::isClosed() const
{
  // we say an isocontour is closed only if its two end points are close enough
  Point first_p = getFirstPoint();
  Point last_p = getLastPoint();
  Real dist = std::sqrt((first_p - last_p).size_sq());
  Real search_radius = 1.1*_iso_width; // changeable
  if (dist < search_radius)
    return true;
  else
    return false;
}

bool
Isocontour::isUnderdeveloped(Real ref_len) const
{
  // check if the isocontour is like a circle
  Point p0 = getFirstPoint();
  Real d_max = -1.0;
  for (unsigned int i = 0; i < _iso_segms.size(); ++i)
  {
    Real dist = std::sqrt((p0 - _iso_segms[i][1]).size_sq());
    if (dist > d_max) d_max = dist;
  }
  if (d_max <= ref_len)
    return true;
  else
    return false;
}

void
Isocontour::reverse()
{
  unsigned int n_segm = _iso_segms.size();
  if (n_segm > 1)
  {
    std::vector<Point> segm_temp(2, Point(0.0,0.0,0.0));
    for (unsigned int j = 0; j <= floor(0.5*(n_segm-2)); ++j) // reverse
    {
      segm_temp = _iso_segms[j];
      _iso_segms[j] = _iso_segms[n_segm-1-j];
      _iso_segms[n_segm-1-j] = segm_temp;
    } // j
  }
  for (unsigned int j = 0; j < n_segm; ++j) // swap the two end points
  {
    Point p_temp = _iso_segms[j][0];
    _iso_segms[j][0] = _iso_segms[j][1];
    _iso_segms[j][1] = p_temp;
  }
}

void
Isocontour::getIsoSegms(std::vector<std::vector<Point> > & segms) const
{
  segms = _iso_segms;
}

void
Isocontour::addIsoSegmsFrom(const Isocontour & other_contour)
{
  for (unsigned int i = 0; i < other_contour._iso_segms.size(); ++i)
    _iso_segms.push_back(other_contour._iso_segms[i]);
}

void
Isocontour::printSegms() const
{
  std::cout << "***** print isosegms in a isocontour *****" << std::endl;
  for (unsigned int i = 0; i < _iso_segms.size(); ++i)
  {
    std::cout << _iso_segms[i][0](0) << ", " << _iso_segms[i][0](1) << ", "
              << _iso_segms[i][1](0) << ", " << _iso_segms[i][1](1) << std::endl;
  }
}
