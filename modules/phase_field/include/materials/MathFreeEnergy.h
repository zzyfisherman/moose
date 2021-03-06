#ifndef MATHFREEENERGY_H
#define MATHFREEENERGY_H

#include "DerivativeBaseMaterial.h"

//Forward Declarations
class MathFreeEnergy;

template<>
InputParameters validParams<MathFreeEnergy>();

/**
 * Material class that creates the math free energy and its derivatives
 * for use with CHParsed and SplitCHParsed. F = 1/4*(1 + c)^2*(1 - c)^2.
 */
class MathFreeEnergy : public DerivativeBaseMaterial
{
public:
  MathFreeEnergy(const std::string & name,
             InputParameters parameters);

protected:
  virtual Real computeF();
  virtual Real computeDF(unsigned int j_var);
  virtual Real computeD2F(unsigned int j_var, unsigned int k_var);
  virtual Real computeD3F(unsigned int j_var, unsigned int k_var, unsigned int l_var);

private:
  /// Coupled variable value for the concentration \f$ \c \f$.
  VariableValue & _c;
  unsigned int _c_var;
};

#endif //MATHFREEENERGY_H
