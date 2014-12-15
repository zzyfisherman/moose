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

#ifndef ISOCONTOUR_H
#define ISOCONTOUR_H

#include "ElementUserObject.h"

class Isocontour
{
public:
  Isocontour(Real iso_width, const std::vector<std::vector<Point> > & segms);
  Point getFirstPoint() const;
  Point getLastPoint() const;
  bool isClockWise() const;
  bool isClosed() const;
  bool isUnderdeveloped(Real ref_len) const;
  void reverse();
  void getIsoSegms(std::vector<std::vector<Point> > & segms) const;
  void addIsoSegmsFrom(const Isocontour & other_contour);
  void printSegms() const;

private:
  Real _iso_width;
  std::vector<std::vector<Point> > _iso_segms;
};

#endif
