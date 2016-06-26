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

#ifndef PLATEHOLEDISP1_H
#define PLATEHOLEDISP1_H

#include "Function.h"

class PlateHoleDisp1;

template<>
InputParameters validParams<PlateHoleDisp1>();

class PlateHoleDisp1 : public Function
{
public:
  PlateHoleDisp1(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _sig0;
  Real _a; // hole radius
  Real _E; // young's modulus
  Real _nu; // possion ratio
};

#endif //PLATEHOLEDISP1_H
