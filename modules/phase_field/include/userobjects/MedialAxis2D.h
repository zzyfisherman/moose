#ifndef MEDIALAXIS2D_H
#define MEDIALAXIS2D_H

#include "GeneralUserObject.h"
#include "ElementIsocontour.h"

class CutNode
{
public:

  CutNode(Point pt0, Point tang, int prop);
  ~CutNode();

private:

  Point _point;
  Point _tang;
  int _prop;

public:

  bool isInSegm;
  bool isNewCut;
  Point getPoint() {return _point;}
  Point getTang() {return _tang;}
  int getProp() {return _prop;}
  void switchPoint(Point new_p) {_point = new_p;}
  bool hasSimilarPropWith(CutNode &cut_node2);
};

class ElemSegms
{
public:

  ElemSegms(int elem_id);
  ~ElemSegms();

private:

  int _elem_id;

public:

  std::vector<std::vector<Point> > medial_segms;
  std::vector<CutNode> cut_nodes;
  void add_segm(std::vector<Point> & line_segm) {medial_segms.push_back(line_segm);}
  void add_cut_node(CutNode p) {cut_nodes.push_back(p);}
  int num_segms() {return medial_segms.size();}
  int num_cut_nodes() {return cut_nodes.size();}
  int get_elem_id() {return _elem_id;}
  void addCutNodesFrom(ElemSegms & other_elem);
  void addSegmsFrom(ElemSegms & other_elem);
  bool cutNodeExist(CutNode & cut_node, Real tol);
};

class MedialAxis2D : public GeneralUserObject
{
public:

  MedialAxis2D(const InputParameters & parameters);

  virtual ~MedialAxis2D(){}

  virtual void initialSetup();

  virtual void residualSetup() {}

  virtual void timestepSetup() {}

  virtual void execute();

  virtual void initialize() {}
  virtual void finalize();
  void getAllMedialAxis(std::vector<std::vector<Point> > &segms);
  void getNewMedialAxis(std::vector<std::vector<Point> > &segms);

protected:

  Real  _l;
  std::string _aux_mesh;
  int _nx;
  int _ny;
  Real _xmin;
  Real _xmax;
  Real _ymin;
  Real _ymax;
  const ElementIsocontour & _elem_isocontour;
  Real _iso_width;
  bool _cut_mesh;
  bool _single_crack;

  std::vector<std::vector<Point> > _all_medial_axis;
  std::vector<std::vector<Point> > _new_medial_axis;
  std::vector<ElemSegms> _new_elem_segms;
  std::vector<ElemSegms> _all_elem_segms;
  std::vector<Point> _aux_nodes;
  std::vector<std::vector<int> > _aux_ix;
  std::vector<std::vector<int> > _nodal_elem_ix;
  std::vector<int> _nodal_elems;
  std::vector<std::vector<int> > _nodal_node_ix;
  std::vector<int> _nodal_nodes;
  std::vector<bool> _boundary;

  void uniformMesh();
  void readNodes();
  void readElements();
  void readIsoSegms(std::vector<std::vector<Point> > &iso_segms); // DEBUG ONLY
  void meshTopology();
  void computeMedialAxis2D(Isocontour & isocontour);

private:

  void closestNodalPatchSearch(Point &one_inter, std::vector<Point> &side_segm, std::vector<std::vector<Point> > &isosegm,
                               int elem_id, int side_id, std::vector<bool> &checked_node, std::vector<CutNode> &cut_nodes);
  bool getIsoTip(std::vector<std::vector<Point> > &iso_segms, int elem_id, CutNode cut_node, Point &tip_node);
  bool isInsideElem2D(Point &p, int elem_id);
  void correctNewCutNode(CutNode &cut_node, int elem_id);
  bool getOldSegmInterNewCut(CutNode new_cut, int elem_id, Point &p_int);
  bool getBranchPoint(std::vector<std::vector<Point> > &iso_segms, std::vector<CutNode> &cut_nodes,
                      Point &branch, int elem_id);
  void tangentMatchedSegms(std::vector<CutNode> &cut_nodes, std::vector<std::vector<Point> > &elem_segms);
  void removeInvalidSegms(std::vector<CutNode> &cut_nodes, std::vector<std::vector<int> > &segms_jcut, 
                          std::vector<std::vector<int> > &all_segm_id);
  bool checkSegmOverlap(std::vector<Point> &segm1, std::vector<Point> &segm2);

  bool intersectTwoLines(std::vector<Point> &segm1, std::vector<Point> &segm2, bool finite, Point &p_int);
  void intersectLines(std::vector<Point> &segm1, std::vector<std::vector<Point> > &segms2,
                      std::vector<Point> &p_inters, std::vector<int> &index2);
  bool intersectLineAndSegm(std::vector<Point> &segm1, std::vector<Point> &segm2, Point &p_int);

  void sortAscend(std::vector<Real> &v, std::vector<int> &index);
  template <class T> int getMinIndex(std::vector<T> &v);
  template <class T> int getMaxIndex(std::vector<T> &v);
  template <class T> void reArrange(std::vector<T> &v, std::vector<int> &index);
  template <class T> void shrinkArray(std::vector<T> &v, std::vector<int> &index);
  Real getVectorAngle(Point v1, Point v2, bool acute_angl);
  Real getMinSideLength(int elem_id);
  Real getMaxSideLength(int elem_id);
  Real safeAcos(Real x);
  void printMedialAxis();
  void writeMedialAxis(std::vector<Isocontour> & isocontours);
  template <typename T> std::string numberToString(T number);
  bool isValidTwoCutConnection(CutNode &cut_node1, CutNode &cut_node2);
  bool isInsideContour(Point p, std::vector<std::vector<Point> > &iso_segms);
};

template<>
InputParameters validParams<MedialAxis2D>();

#endif
