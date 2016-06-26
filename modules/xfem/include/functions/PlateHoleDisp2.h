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

#ifndef PLATEHOLEDISP2_H
#define PLATEHOLEDISP2_H

#include "Function.h"

class PlateHoleDisp2;

template<>
InputParameters validParams<PlateHoleDisp2>();

class PlateHoleDisp2 : public Function
{
public:
  PlateHoleDisp2(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p);

protected:
  Real _sig0;
  Real _a; // hole radius
  Real _E; // young's modulus
  Real _nu; // possion ratio
};

#endif //PLATEHOLEDISP2_H
