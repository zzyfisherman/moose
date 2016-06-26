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

#ifndef PLATEHOLESIG12_H
#define PLATEHOLESIG12_H

#include "Function.h"

class PlateHoleSig12;

template<>
InputParameters validParams<PlateHoleSig12>();

class PlateHoleSig12 : public Function
{
public:
  PlateHoleSig12(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _sig0;
  Real _a; // hole radius
};

#endif //PLATEHOLESIG12_H
