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

#ifndef XFEMCUTELEM_H
#define XFEMCUTELEM_H

// Local Includes -----------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh.h" // libMesh::invalid_uint
#include "libmesh/location_maps.h"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
#include "libmesh/node.h"
#include "libmesh/vector_value.h"
#include "libmesh/point_locator_base.h"

// C++ Includes   -----------------------------------
#include <vector>

#include "MaterialData.h"
#include "ElementFragmentAlgorithm.h"

using namespace libMesh;

class XFEMCutElem
{
public:

  XFEMCutElem(Elem* elem, unsigned int n_qpoints);
  virtual ~XFEMCutElem();

protected:

  unsigned int _n_nodes;
  unsigned int _n_qpoints;
  std::vector<Node*> _nodes;
  Real _elem_volume;
  Real _physical_volfrac;
  std::vector<Real> _new_weights; // quadrature weights from moment fitting
  virtual Point get_node_coords(EFAnode* node, MeshBase* displaced_mesh = NULL) const = 0;

public:

  virtual void calc_physical_volfrac() = 0;
  Real get_physical_volfrac() const;
  virtual void calc_mf_weights() = 0; // ZZY
  Real get_mf_weights(unsigned int i_qp) const; // ZZY
  virtual Point get_origin(unsigned int plane_id, MeshBase* displaced_mesh=NULL) const = 0;
  virtual Point get_normal(unsigned int plane_id, MeshBase* displaced_mesh=NULL) const = 0;
  virtual void get_frag_faces(std::vector<std::vector<Point> > &frag_faces, MeshBase* displaced_mesh=NULL) const = 0;
  virtual const EFAelement * get_efa_elem() const = 0;
  virtual unsigned int num_cut_planes() const = 0;
  virtual void get_frag_values(std::vector<Real> &elem_node_vals, std::vector<Real> &frag_node_vals) const = 0;
  virtual void get_frag_node_ids(std::vector<int> &frag_node_ids) const = 0;
};

#endif
