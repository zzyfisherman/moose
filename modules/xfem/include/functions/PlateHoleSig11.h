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

#ifndef PLATEHOLESIG11_H
#define PLATEHOLESIG11_H

#include "Function.h"

class PlateHoleSig11;

template<>
InputParameters validParams<PlateHoleSig11>();

class PlateHoleSig11 : public Function
{
public:
  PlateHoleSig11(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _sig0;
  Real _a; // hole radius
};

#endif //PLATEHOLESIG11_H
