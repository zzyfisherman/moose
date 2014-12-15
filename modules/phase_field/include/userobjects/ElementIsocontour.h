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

#ifndef ELEMENTISOCONTOUR_H
#define ELEMENTISOCONTOUR_H

#include "ElementUserObject.h"
#include "Isocontour.h"

//Forward Declarations
class ElementIsocontour;

template<>
InputParameters validParams<ElementIsocontour>();

class ElementIsocontour : public ElementUserObject
{
public:

  ElementIsocontour(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
//  virtual Real getValue();

  virtual void finalize();
  Real getIsoVal() const;
  void getIsoContours(std::vector<Isocontour> & isocontours) const;

protected:

  Real _iso_val;
  Real _iso_width;
  VariableValue *_p_pf;
  std::vector<Isocontour> _isocontours;
  void sort_segms(std::vector<std::vector<Point> > &iso_segms, std::vector<std::vector<int> > &iso_index);
  void close_isocontours();
  XFEM* _xfem;

private:

  std::vector<std::vector<Point> > _iso_segms;
  std::vector<std::vector<int> > _iso_index; 
  void combine_isocontours(unsigned int iso_id);
  void flattenPointVector(std::vector<std::vector<Point> > &segms, std::vector<Real> &single_v);
  void reshapePointVector(std::vector<Real> &single_v, std::vector<std::vector<Point> > &segms);
  bool truncateSegm(std::vector<Point> &iso_segm);
  bool intersectTwoSegms(std::vector<Point> &segm1, std::vector<Point> &segm2, Point &p_int);
  bool pointExistInSegm(Point p, std::vector<Point> &segm, Real h_ref);
};

#endif
