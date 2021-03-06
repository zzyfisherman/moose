#ifndef TENSORMECHANICSPLASTICWEAKPLANESHEAR_H
#define TENSORMECHANICSPLASTICWEAKPLANESHEAR_H

#include "TensorMechanicsPlasticModel.h"
#include "TensorMechanicsHardeningModel.h"


class TensorMechanicsPlasticWeakPlaneShear;


template<>
InputParameters validParams<TensorMechanicsPlasticWeakPlaneShear>();

/**
 * Rate-independent associative weak-plane tensile failure
 * with hardening/softening.  The cone's tip is smoothed.
 */
class TensorMechanicsPlasticWeakPlaneShear : public TensorMechanicsPlasticModel
{
 public:
  TensorMechanicsPlasticWeakPlaneShear(const std::string & name, InputParameters parameters);


 protected:

  /// Hardening model for cohesion
  const TensorMechanicsHardeningModel & _cohesion;

  /// Hardening model for tan(phi)
  const TensorMechanicsHardeningModel & _tan_phi;

  /// Hardening model for tan(psi)
  const TensorMechanicsHardeningModel & _tan_psi;


  /**
   * The yield function
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the yield function
   */
  Real yieldFunction(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to stress
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return df_dstress(i, j) = dyieldFunction/dstress(i, j)
   */
  RankTwoTensor dyieldFunction_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to the internal parameter
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the derivative
   */
  Real dyieldFunction_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The flow potential
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return the flow potential
   */
  RankTwoTensor flowPotential(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to stress
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dstress(i, j, k, l) = dr(i, j)/dstress(k, l)
   */
  RankFourTensor dflowPotential_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to the internal parameter
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dintnl(i, j) = dr(i, j)/dintnl
   */
  RankTwoTensor dflowPotential_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The yield function is modified to
   * f = sqrt(s_xz^2 + s_yz^2 + a) + s_zz*_tan_phi - _cohesion
   * where "a" depends on the tip_scheme.  Currently _tip_scheme is
   * 'hyperbolic', where a = _small_smoother2
   * 'cap' where a = _small_smoother2 + (p(stress(2,2) - _cap_start))^2
   *       with the function p(x)=x(1-exp(-_cap_rate*x)) for x>0, and p=0 otherwise
   */
  MooseEnum _tip_scheme;

  /// smoothing parameter for the cone's tip - see doco for _tip_scheme
  Real _small_smoother2;

  /// smoothing parameter dictating when the 'cap' will start - see doco for _tip_scheme
  Real _cap_start;

  /// dictates how quickly the 'cap' degenerates to a hemisphere - see doco for _tip_scheme
  Real _cap_rate;

  /// Function that's used in dyieldFunction_dstress and flowPotential
  RankTwoTensor df_dsig(const RankTwoTensor & stress, const Real & _tan_phi_or_psi) const;

  /// returns the 'a' parameter - see doco for _tip_scheme
  virtual Real smooth(const RankTwoTensor & stress) const;

  /// returns the da/dstress(2,2) - see doco for _tip_scheme
  virtual Real dsmooth(const RankTwoTensor & stress) const;

  /// returns the d^2a/dstress(2,2)^2 - see doco for _tip_scheme
  virtual Real d2smooth(const RankTwoTensor & stress) const;

  /// cohesion as a function of internal parameter
  virtual Real cohesion(const Real internal_param) const;

  /// d(cohesion)/d(internal_param)
  virtual Real dcohesion(const Real internal_param) const;

  /// tan_phi as a function of internal parameter
  virtual Real tan_phi(const Real internal_param) const;

  /// d(tan_phi)/d(internal_param);
  virtual Real dtan_phi(const Real internal_param) const;

  /// tan_psi as a function of internal parameter
  virtual Real tan_psi(const Real internal_param) const;

  /// d(tan_psi)/d(internal_param);
  virtual Real dtan_psi(const Real internal_param) const;

};

#endif // TENSORMECHANICSPLASTICWEAKPLANESHEAR_H
