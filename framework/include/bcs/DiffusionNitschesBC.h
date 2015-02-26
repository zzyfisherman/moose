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

#ifndef DIFFUSIONNITSCHESBC_H
#define DIFFUSIONNITSCHESBC_H

#include "IntegratedBC.h"


class DiffusionNitschesBC;

template<>
InputParameters validParams<DiffusionNitschesBC>();

/**
 * Implements a simple constant Neumann BC where grad(u)=value on the boundary.
 * Uses the term produced from integrating the diffusion operator by parts.
 */
class DiffusionNitschesBC : public IntegratedBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  DiffusionNitschesBC(const std::string & name, InputParameters parameters);


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  Function & _func; // g
  Real _alpha;
};


#endif //NEUMANNBC_H
