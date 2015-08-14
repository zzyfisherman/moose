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

#ifndef ELEMENTINTEGRALPOSTPROCESSOR_H
#define ELEMENTINTEGRALPOSTPROCESSOR_H

#include "ElementPostprocessor.h"

//Forward Declarations
class ElementIntegralPostprocessor;

template<>
InputParameters validParams<ElementIntegralPostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementIntegralPostprocessor : public ElementPostprocessor
{
public:
  ElementIntegralPostprocessor(const InputParameters & parameters);
  ElementIntegralPostprocessor(const std::string & deprecated_name, InputParameters parameters); // DEPRECATED CONSTRUCTOR

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual Real getValue();

protected:
  virtual Real computeQpIntegral() = 0;
  virtual Real computeIntegral();

  unsigned int _qp;

  Real _integral_value;

  // ZZY brings XFEM stuff in
  XFEM *_xfem;

  /// XFEM quadrature rule WJ
  std::string _xfem_qrule;
  std::vector<Real> _xfem_weights;
  void get_xfem_weights(std::vector<Real> & _xfem_weights);
};

#endif
