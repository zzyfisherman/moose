
#ifndef XFEM_CIRCLE_CUT_H
#define XFEM_CIRCLE_CUT_H

#include "XFEM_geometric_cut.h"

class XFEM_circle_cut : public XFEM_geometric_cut
{
public:

  XFEM_circle_cut(std::vector<Real> square_nodes);
  ~XFEM_circle_cut();

  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_elem_by_geometry(const Elem* elem, std::vector<cutFace> & cutFaces, Real time);

  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_edges, std::vector<cutEdge> & cutEdges, Real time);
  virtual bool cut_frag_by_geometry(std::vector<std::vector<Point> > & frag_faces, std::vector<cutFace> & cutFaces, Real time);

private:

  std::vector<Point> _vertices;
  Point _center;
  Point _normal;
  Real _radius;
  Real _angle;

private:

  bool intersect_with_edge(Point p1, Point p2, Point &pint);
  int plane_normal_line_exp_int_3d(double pp[3], double normal[3], double p1[3], double p2[3], double pint[3]);
  bool line_exp_is_degenerate_nd(int dim_num, double p1[], double p2[]);
  double r8vec_norm(int n, double a[]);
  void r8vec_copy(int n, double a1[], double a2[]);
  bool r8vec_eq(int n, double a1[], double a2[]);
  double r8vec_dot_product(int n, double a1[], double a2[]);
  bool isInsideCutPlane(Point p);
  bool isInsideEdge(Point p1, Point p2, Point p);
  Real getRelativePosition(Point p1, Point p2, Point p);
};

#endif

