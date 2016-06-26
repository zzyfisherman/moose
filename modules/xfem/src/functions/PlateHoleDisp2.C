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

#include "PlateHoleDisp2.h"
#include <math.h>

template<>
InputParameters validParams<PlateHoleDisp2>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("sig0", 1.0, "The value of tensile traction in x");
  params.addParam<Real>("a", 1.0, "The radius of the hole");
  params.addParam<Real>("E", 1.0, "Young's modulus");
  params.addParam<Real>("nu", 1.0, "Poission ratio");
  return params;
}

PlateHoleDisp2::PlateHoleDisp2(const InputParameters & parameters) :
    Function(parameters),
    _sig0(getParam<Real>("sig0")),
    _a(getParam<Real>("a")),
    _E(getParam<Real>("E")),
    _nu(getParam<Real>("nu"))
{}

Real
PlateHoleDisp2::value(Real /*t*/, const Point & p)
{
  Real r = 0.0;
  Real theta = 0.0;
  Real k = 0.0; // kappa - Kolosov constant
  Real mu = 0.0; // Shear modulus
  Real disp2 = 0.0;

  r = sqrt(p(0)*p(0) + p(1)*p(1));
  theta = atan2(p(1),p(0));
  k = 3.0 - 4.0*_nu;
  mu = 0.5*_E/(1+_nu);
  if (r != 0.0)
    disp2 = _sig0*(_a/(8.0*mu))*((r/_a)*(k-3.0)*sin(theta) + 2.0*(_a/r)*((1.0-k)*sin(theta) + sin(3.0*theta)) - 2.0*((_a*_a*_a)/(r*r*r))*sin(3.0*theta));

  return disp2;
}
