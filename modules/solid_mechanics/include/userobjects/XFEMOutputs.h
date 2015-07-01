/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef XFEMOUTPUTS_H
#define XFEMOUTPUTS_H

#include "ElementUserObject.h"
#include <cmath>
#include <set>
#include <algorithm>
#include <map>
#include <vector>

class XFEMOutputs;

template<>
InputParameters validParams<XFEMOutputs>();

class XFEMOutputs : public ElementUserObject
{
public:
  XFEMOutputs(const InputParameters & parameters);

  ~XFEMOutputs(); // the destructor closes the output file

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & u );
  virtual void finalize();

protected:

  VariableValue * _disp_x;
  VariableValue * _disp_y;
  VariableValue * _disp_z;
  std::vector<VariableValue *> _disp;
  std::map<unsigned int, Point> _id_to_nodes;
  std::map<unsigned int, Real> _id_to_values;
  std::vector<std::vector<unsigned int> > _connectivity;
  unsigned int _num_emb_nodes;
  unsigned int _num_nodes;
  unsigned int _mesh_dim;
  int _disp_dir;
  XFEM* _xfem;

private:

  template <typename T> std::string numberToString(T number);
};

#endif
