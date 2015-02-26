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

#include "DiffusionNitschesBC.h"
#include "Function.h"
#include "libmesh/numeric_vector.h"
#include <cmath>

template<>
InputParameters validParams<DiffusionNitschesBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<FunctionName>("function", "The function.");
  params.addRequiredParam<Real>("alpha", "Penalty");
  return params;
}

DiffusionNitschesBC::DiffusionNitschesBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    _func(getFunction("function")),
    _alpha(getParam<Real>("alpha"))
{
}

Real
DiffusionNitschesBC::computeQpResidual()
{ 
  Real r = 0.0;
  r -= _grad_u[_qp] * _normals[_qp] * _test[_i][_qp];
  r -= _grad_test[_i][_qp] *_normals[_qp] * ( _u[_qp] - _func.value(_t, _q_point[_qp]));
  r += _alpha * _test[_i][_qp] * (_u[_qp] - _func.value(_t, _q_point[_qp]));
  return r;
}

Real
DiffusionNitschesBC::computeQpJacobian()
{
  Real r = 0.0;
  r -= (_grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp]);
  r -= (_grad_test[_j][_qp] * _normals[_qp] * _phi[_i][_qp]);
  r += _alpha * _phi[_j][_qp] * _test[_i][_qp];

  return r;
}
