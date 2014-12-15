#include "MedialAxis2D.h"
#include "EFAfuncs.h"
#include "XFEM_geometric_cut_2d.h"
#include "MedialAxis2DFuncs.h"
#include <cmath>
#include <stdio.h>

CutNode::CutNode(Point pt0, Point tang, int prop):
  _point(pt0),
  _tang(tang),
  _prop(prop)
{
  // prop = 0 -> ingoing,  prop = 1 -> outgoing
  // prop = 2 -> sweeping, prop = 3 -> tip
  isInSegm = false;
  isNewCut = true;
}

CutNode::~CutNode()
{
}

bool
CutNode::hasSimilarPropWith(CutNode &cut_node2)
{
  bool is_similar = false;
  if (_prop == 3 && cut_node2.getProp() == 3)
    is_similar = true;
  else if (_prop != 3 && cut_node2.getProp() != 3)
    is_similar = true;
  return is_similar;
}

ElemSegms::ElemSegms(int elem_id):
  _elem_id(elem_id)
{
  medial_segms.clear();
  cut_nodes.clear();
}

void
ElemSegms::addCutNodesFrom(ElemSegms & other_elem)
{
  for (unsigned int i = 0; i < other_elem.num_cut_nodes(); ++i)
    cut_nodes.push_back(other_elem.cut_nodes[i]);
}

void
ElemSegms::addSegmsFrom(ElemSegms & other_elem)
{
  for (unsigned int i = 0; i < other_elem.num_segms(); ++i)
    medial_segms.push_back(other_elem.medial_segms[i]);
}

bool
ElemSegms::cutNodeExist(CutNode & cut_node, Real tol)
{
  bool node_exist = false;
  for (unsigned int j = 0; j < cut_nodes.size(); ++j)
  {
    Real diff = std::sqrt((cut_node.getPoint() - cut_nodes[j].getPoint()).size_sq());
    if (diff < tol && cut_node.hasSimilarPropWith(cut_nodes[j]))
      node_exist = true;
  } // j
  return node_exist;
}

ElemSegms::~ElemSegms()
{
}

template<>
InputParameters validParams<MedialAxis2D>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<Real>("l", 0.0, "phase field characteristic length");
  params.addParam<std::string>("aux_mesh", "default", "aux mesh name" );
  params.addParam<int>("nx", 0, "number of elements in x direction of the aux mesh");
  params.addParam<int>("ny", 0, "number of elements in y direction of the aux mesh");
  params.addParam<Real>("xmin", 0.0, "min x of the aux mesh");
  params.addParam<Real>("xmax", 0.0, "max x of the aux mesh");
  params.addParam<Real>("ymin", 0.0, "min y of the aux mesh");
  params.addParam<Real>("ymax", 0.0, "max y of the aux mesh");
  params.addRequiredParam<UserObjectName>("isocontour","The phase field isocontour user object name");
  params.addParam<bool>("cut_mesh", false, "if you want to actually cut the mesh");
  params.addParam<bool>("single_crack", false, "if the test has only one crack");
  return params;
}

MedialAxis2D::MedialAxis2D(const InputParameters & params)
  :GeneralUserObject(params),
   _l(getParam<Real>("l")),
   _aux_mesh(getParam<std::string>("aux_mesh")),
   _nx(getParam<int>("nx")),
   _ny(getParam<int>("ny")),
   _xmin(getParam<Real>("xmin")),
   _xmax(getParam<Real>("xmax")),
   _ymin(getParam<Real>("ymin")),
   _ymax(getParam<Real>("ymax")),
   _elem_isocontour(getUserObject<ElementIsocontour>("isocontour")),
   _iso_width(-2.0*_l*log(_elem_isocontour.getIsoVal())),
   _cut_mesh(getParam<bool>("cut_mesh")),
   _single_crack(getParam<bool>("single_crack"))
{
  _iso_width *= 3.0; // actual _iso_width -- 3.0 may not be very accurate
}

void 
MedialAxis2D::initialSetup()
{
  // Warning message
  std::cout << "WARNING: MedialAxis2D is not 100% robust. We need some computer graphics experts to develop a better version!" << std::endl;

  // Get aux mesh
  if (_aux_mesh == "default")
    uniformMesh();
  else
  {
    readNodes();
    readElements();
  }

  // Get mesh topology arrays
  meshTopology();

  // Intialize _all_new_segms
  _all_medial_axis.clear();
  _all_elem_segms.clear();
  for (unsigned int e = 0; e < _aux_ix.size(); ++e) 
    _all_elem_segms.push_back(ElemSegms(e));
}

void
MedialAxis2D::execute()
{
  // retrieve all isocontours
  std::vector<Isocontour> isocontours;
  _elem_isocontour.getIsoContours(isocontours);

  // DEBUG
  for (unsigned int i = 0; i < isocontours.size(); ++i)
    if (processor_id() == 0)
      isocontours[i].printSegms();

  // DEBUG, clear iso_segms and read iso_segms from outside
//  std::vector<std::vector<Point> > iso_segms;
//  readIsoSegms(iso_segms);
//  isocontours.clear();
//  isocontours.push_back(Isocontour(0.1, iso_segms));

  _new_elem_segms.clear(); // stores all new cut nodes and medial segments for all isocontours
  _new_medial_axis.clear();
  if (!isocontours.empty())
  {
    // initialize _new_elem_segms
    for (unsigned int e = 0; e < _aux_ix.size(); ++e)
      _new_elem_segms.push_back(ElemSegms(e));

    // compute medial axis segments
    for (unsigned int i = 0; i < isocontours.size(); ++i) // loop over separate isocontours
      computeMedialAxis2D(isocontours[i]);

    // save _new_elem_segms to _all_elem_segms
    for (unsigned int e = 0; e < _aux_ix.size(); ++e)
    {
      _all_elem_segms[e].addCutNodesFrom(_new_elem_segms[e]);
      _all_elem_segms[e].addSegmsFrom(_new_elem_segms[e]);
    }  
  } // if isocontour segments exist  

  // DEBUG
  if (processor_id() == 0)
  {
    printMedialAxis();
    writeMedialAxis(isocontours);
  }
}

void
MedialAxis2D::finalize()
{
  if (_cut_mesh)
  {
    XFEM* xfem = _fe_problem.get_xfem();
    for (unsigned int i = 0; i < _new_medial_axis.size(); ++i)
    {
      Real x0 = _new_medial_axis[i][0](0);
      Real y0 = _new_medial_axis[i][0](1);
      Real x1 = _new_medial_axis[i][1](0);
      Real y1 = _new_medial_axis[i][1](1);
      xfem->addGeometricCut(new XFEM_geometric_cut_2d(x0, y0, x1, y1, 0.0, 0.0));
    }
  }
}

void
MedialAxis2D::getAllMedialAxis(std::vector<std::vector<Point> > &segms)
{
  segms.clear();
  for (unsigned int i = 0; i < _all_medial_axis.size(); ++i)
    segms.push_back(_all_medial_axis[i]);
}

void
MedialAxis2D::getNewMedialAxis(std::vector<std::vector<Point> > &segms)
{
  segms.clear();
  for (unsigned int i = 0; i < _new_medial_axis.size(); ++i)
    segms.push_back(_new_medial_axis[i]);
}

void
MedialAxis2D::uniformMesh()
{
  if ((_nx == 0) || (_ny == 0))
    mooseError("no element in aux mesh defined!");  

  Real hx = (_xmax - _xmin)/_nx;
  Real hy = (_ymax - _ymin)/_ny;
  int numnp_x = _nx + 1;
  int numnp_y = _ny + 1;

  if ((hx <= _iso_width) || (hy <= _iso_width))
    mooseError("the default aux mesh is too fine!");

  // get nodal coords
  for (unsigned int j = 0; j < numnp_y; ++j)
    for (unsigned int i = 0; i < numnp_x; ++i)
      _aux_nodes.push_back(Point(i*hx, j*hy, 0.0));

  // get mesh connectivity
  for (unsigned int j = 0; j < _ny; ++j)
  {
    for (unsigned int i = 0; i < _nx; ++i)
    {
      std::vector<int> line_ix(4,0);
      int first_node = numnp_x*j + i;
      line_ix[0] = first_node;
      line_ix[1] = first_node + 1;
      line_ix[2] = first_node + 1 + numnp_x;
      line_ix[3] = first_node + numnp_x;
      _aux_ix.push_back(line_ix);
    }
  }
}

void
MedialAxis2D::readNodes()
{
  std::ifstream node_file;
  std::string filename = _aux_mesh + ".nod";

  node_file.open(filename.c_str());
  if (node_file.is_open())
  {
    while (1)
    {
      std::string line;
      std::vector<Real> point_coord;
      if (getline(node_file,line) == 0)
        break;
      std::istringstream is(line);
      Real coord;
      while(is >> coord)
        point_coord.push_back(coord);
      if (point_coord.size() == 3)
      {
        Point p = Point(point_coord[0], point_coord[1], point_coord[2]);
        _aux_nodes.push_back(p);
      }
      else
        mooseError("Need 3 coords for a mesh node");
    }
    node_file.close();
  }
  else
    mooseError("Mesh node file cannot be opened");
}

void
MedialAxis2D::readElements()
{
  std::ifstream elem_file;
  std::string filename = _aux_mesh + ".ele";

  elem_file.open(filename.c_str());
  if (elem_file.is_open())
  {
    while (1)
    {
      std::string line;
      std::vector<int> elem_ix;
      if (getline(elem_file,line) == 0)
        break;
      std::istringstream is(line);
      int n;
      while (is >> n)
        elem_ix.push_back(n);
      _aux_ix.push_back(elem_ix);
    }
    elem_file.close();
  }
  else
    mooseError("Mesh element file cannot be opened");
}

void
MedialAxis2D::readIsoSegms(std::vector<std::vector<Point> > &iso_segms) // only for debugging
{
  std::ifstream iso_file;
//  std::string filename = _aux_mesh + ".iso";
  std::string filename = "test.iso";

  iso_file.open(filename.c_str());
  if (iso_file.is_open())
  {
    while (1)
    {
      std::string line;
      std::vector<Real> coords;
      if (getline(iso_file,line) == 0)
        break;
      std::istringstream is(line);
      Real x;
      while (is >> x)
        coords.push_back(x);
      if (coords.size() == 4)
      {
        std::vector<Point> segm_line(2,Point(0.0,0.0,0.0));
        segm_line[0] = Point(coords[0],coords[1],0.0);
        segm_line[1] = Point(coords[2],coords[3],0.0);
        iso_segms.push_back(segm_line);
      }
      else
        mooseError("each line in the .iso file should have 4 values");
    }
    iso_file.close();
  }
  else
    mooseError("Isocontour file cannot be opened");
}

void
MedialAxis2D::meshTopology()
{
  // Initialize variables
  int nen = _aux_ix[0].size();
  int numnp = _aux_nodes.size();
  int numel = _aux_ix.size();
  std::vector<int> ib(numnp,0);
  _nodal_elem_ix.resize(numnp); // nifd
  for (unsigned int i = 0; i < numnp; ++i) _nodal_elem_ix[i].resize(2,0);
  _nodal_elems.resize(numel*nen,0); // nif
  _nodal_node_ix.resize(numnp); // nfsd
  for (unsigned int i = 0; i < numnp; ++i) _nodal_node_ix[i].resize(2,0);
  
  // Get patch element information 
  for (unsigned int j = 0; j < numel; ++j)
    for (unsigned int k = 0; k < nen; ++k)
        ib[_aux_ix[j][k]] += 1;

  int ip = -1;
  int iq = -1;
  for (unsigned int i = 0; i < numnp; ++i)
  {
    _nodal_elem_ix[i][0] = ip + 1;
    ip += ib[i];
    _nodal_elem_ix[i][1] = ip;
    ib[i] = 0;
  }

  for (unsigned j = 0; j < numel; ++j)
  {
    for (unsigned k = 0; k < nen; ++k)
    {
      ib[_aux_ix[j][k]] += 1;
      _nodal_elems[_nodal_elem_ix[_aux_ix[j][k]][0] + ib[_aux_ix[j][k]] - 1] = j;
    }
  }

  // Get patch node information
  ip = -1;
  iq = -1;
  for (unsigned int i = 0; i < numnp; ++i) // loop over all nodal points
  {
    _nodal_node_ix[i][0] = iq + 1; // iq counts the total patch nodes
    int n_patch_elem = _nodal_elem_ix[i][1] - _nodal_elem_ix[i][0] + 1; // number of patch elements of this node

    // loop over nodal elements
    for (unsigned int m = 1; m <= n_patch_elem; ++m)
    {
      int j = _nodal_elems[_nodal_elem_ix[i][0] + m - 1];
      ip = _nodal_elem_ix[i][0] + m - 1;
      for (unsigned int k = 0; k < nen; ++k)
      {
        if (_aux_ix[j][k] == i)
        {
          int node1(k < nen-1 ? k+1 : 0);
          int node2(k > 0 ? k-1 : nen-1);

          if (_nodal_node_ix[i][0] > iq)
          {
            iq += 1; // save node1
            _nodal_nodes.push_back(_aux_ix[j][node1]);
            iq += 1; // save node2
            _nodal_nodes.push_back(_aux_ix[j][node2]);
          }
          else
          {
            bool node1_exist = false; // save node1
            for (unsigned int l = _nodal_node_ix[i][0]; l <= iq; ++l) // check if already counted
            {
              if (_aux_ix[j][node1] == _nodal_nodes[l])
              {
                node1_exist = true; // node1 has been stored
                break;
              }
            }
            if (!node1_exist)
            {
              iq += 1;
              _nodal_nodes.push_back(_aux_ix[j][node1]);
            }
                   
            bool node2_exist = false; // save node2
            for (unsigned int l = _nodal_node_ix[i][0]; l <= iq; ++l)
            {
              if (_aux_ix[j][node2] == _nodal_nodes[l])
              {
                node2_exist = true;
                break;
              }
            }
            if (!node2_exist)
            {
              iq += 1;
              _nodal_nodes.push_back(_aux_ix[j][node2]);
            }
          }
        }
      } // k, loop over patch element nodes
    } // m, loop over patch elements
    _nodal_node_ix[i][1] = iq;
  } // i, loop over nodes

  // get boundary code
  _boundary.resize(numnp,false);
  for (unsigned int i = 0; i < numnp; ++i)
  {
    int n_patch_elem = _nodal_elem_ix[i][1] - _nodal_elem_ix[i][0] + 1;
    int n_patch_node = _nodal_node_ix[i][1] - _nodal_node_ix[i][0] + 1;
    if (n_patch_node != n_patch_elem) _boundary[i] = true;
  }
}

void
MedialAxis2D::computeMedialAxis2D(Isocontour & isocontour)
{
  // get iso-segms of this isocontour
  std::vector<std::vector<Point> > isosegm;
  isocontour.getIsoSegms(isosegm);

  // loop all aux elements to compute medial axis segments
  for (unsigned int e = 0; e < _aux_ix.size(); ++e)
  {
    // =============== PART I: find cut ndoes ===============
    // check if the isocontour is under-developed
    if (isocontour.isUnderdeveloped(1.1*getMaxSideLength(e))) // changeable
      continue;

    // loop over elem sides to compute cut nodes
    std::vector<CutNode> cut_nodes;
    int nen = _aux_ix[e].size();
    std::vector<bool> checked_node(nen,false); // a flag
    Real char_h = getMinSideLength(e);

    for (unsigned int s = 0; s < nen; ++s)
    {
      int s_node1 = s;
      int s_node2(s < nen-1 ? s+1 : 0);
      std::vector<Point> side_segm(2,Point(0.0,0.0,0.0));
      side_segm[0] = _aux_nodes[_aux_ix[e][s_node1]];
      side_segm[1] = _aux_nodes[_aux_ix[e][s_node2]];
      Point side_normal(side_segm[1](1)-side_segm[0](1), side_segm[0](0)-side_segm[1](0), 0.0);

      std::vector<Point> side_inters;
      std::vector<int> inter_index;
      intersectLines(side_segm,isosegm,side_inters,inter_index);

      if (side_inters.size() == 1) // this side only has one intersection with isosegms
        closestNodalPatchSearch(side_inters[0],side_segm,isosegm,e,s,checked_node,cut_nodes);
      else if (side_inters.size() > 1) // this side has multiple intersections with isosegms
      {
        // get intersections' tangents and normals
        std::vector<Point> inter_tangs;
        std::vector<Point> inter_norms;
        for (int i = 0; i < side_inters.size(); ++i)
        {
          Point tang = isosegm[inter_index[i]][1] - isosegm[inter_index[i]][0];
          normalize(tang);
          inter_tangs.push_back(tang);
          inter_norms.push_back(Point(tang(1),-tang(0),0.0));
        }

        // sort intersections according to their distance to the first side node
        std::vector<Real> dist(side_inters.size(),0.0);
        std::vector<int> sorted_ix;
        for (unsigned int i = 0; i < side_inters.size(); ++i) dist[i] = sqrt((side_inters[i]-side_segm[0]).size_sq());
        sortAscend(dist,sorted_ix);
        reArrange(side_inters,sorted_ix);
        reArrange(inter_tangs,sorted_ix);
        reArrange(inter_norms,sorted_ix);
        reArrange(inter_index,sorted_ix);

        // form intersection pairs and get a cut node for each pair
        std::vector<bool> paired_inter(side_inters.size(),false);
        for (unsigned int k = 0; k < side_inters.size()-1; ++k)
        {
          unsigned int l = k + 1;
          Point p_cut = 0.5*(side_inters[k] + side_inters[l]);
          Real dotp1 = (p_cut - side_inters[k])*inter_norms[k];
          Real dotp2 = (p_cut - side_inters[l])*inter_norms[l];
          if (dotp1 < 0.0 && dotp2 < 0.0) // the cut node is inside the isocontour
          {
            paired_inter[k] = true; // mark the paired intersection
            paired_inter[l] = true;
            Point tang(inter_index[k] < inter_index[l] ? 0.5*(inter_tangs[k]-inter_tangs[l]) : 0.5*(inter_tangs[l]-inter_tangs[k]));
            normalize(tang);
            int cut_prop(tang*side_normal > 0.0 ? 1 : 0); // 1: outgoing; 0: ingoing
            cut_nodes.push_back(CutNode(p_cut,tang,cut_prop)); // save to element data
          }
        } // k

        // for the lonely inersection, use patch search
        for (unsigned int k = 0; k < side_inters.size(); ++k)
          if (!paired_inter[k])
            closestNodalPatchSearch(side_inters[k],side_segm,isosegm,e,s,checked_node,cut_nodes);
      }
    } // s, loop over sides
    if (cut_nodes.size() == 0)
      continue; // jump to next element

    // Find the nucleation cut node
    if (cut_nodes.size() == 1 && cut_nodes[0].getProp() != 2)
    {
      Point tip_node(0.0,0.0,0.0);
      if (getIsoTip(isosegm, e, cut_nodes[0], tip_node))
      {
        Point tang = cut_nodes[0].getTang(); // the same tangent as the single cut node
        cut_nodes.push_back(CutNode(tip_node, tang, 3));
        std::swap(cut_nodes[0], cut_nodes[1]); // put the nucleation node in the 1st position
      }
    }

    // get NEW cut nodes for this element and correct new cut nodes if necessary
    for (unsigned int i = 0; i < cut_nodes.size(); ++i)
    {
      if (_all_elem_segms[e].cutNodeExist(cut_nodes[i], _iso_width))
        cut_nodes[i].isNewCut = false;
      if (cut_nodes[i].isNewCut)
        correctNewCutNode(cut_nodes[i], e); // cut_nodes[i] changed here
    } // i

    // =============== PART II: find medial segments ===============
    // find all old-segm-intersected new cut nodes (and reset is_new_cut) and get related medial segms
    // also remember to save the cut nodes and segms
    std::vector<std::vector<Point> > elem_segms;
    if (_all_elem_segms[e].num_segms() > 0)
    {
      for (unsigned int i = 0; i < cut_nodes.size(); ++i)
      {
        bool isInterFound = false;
        Point p_int(0.0,0.0,0.0);
        if ((cut_nodes[i].isNewCut) && (cut_nodes[i].getProp() != 2))
          isInterFound = getOldSegmInterNewCut(cut_nodes[i], e, p_int);
        if (isInterFound && isInsideContour(p_int, isosegm))
        {
          std::vector<Point> line_segm(2,Point(0.0,0.0,0.0));
          line_segm[0] = p_int;
          line_segm[1] = cut_nodes[i].getPoint();
          elem_segms.push_back(line_segm);
          cut_nodes[i].isInSegm = true;
          _new_elem_segms[e].add_cut_node(cut_nodes[i]);
          cut_nodes[i].isNewCut = false;
        }
      } // i
    }

    // delete the the cut nodes which already exist
    std::vector<int> rm_ix;
    for (unsigned int i = 0; i < cut_nodes.size(); ++i)
      if (!cut_nodes[i].isNewCut) rm_ix.push_back(i);
    shrinkArray(cut_nodes, rm_ix);

    // get element medial axis segments from the rest new cut nodes
    if (cut_nodes.size() == 2)
    {
      if (isValidTwoCutConnection(cut_nodes[0], cut_nodes[1]))
      {
        std::vector<Point> line_segm(2,Point(0.0,0.0,0.0));
        line_segm[0] = cut_nodes[0].getPoint();
        line_segm[1] = cut_nodes[1].getPoint();
        cut_nodes[0].isInSegm = true;
        cut_nodes[1].isInSegm = true;
        elem_segms.push_back(line_segm);
      }
    }
    else if (cut_nodes.size() == 3)
    {
      Point branch(0.0,0.0,0.0);
      if (getBranchPoint(isosegm, cut_nodes, branch, e))
      {
        for (unsigned int i = 0; i < cut_nodes.size(); ++i)
        {
          std::vector<Point> line_segm(2,Point(0.0,0.0,0.0));
          line_segm[0] = branch;
          line_segm[1] = cut_nodes[i].getPoint();
          cut_nodes[i].isInSegm = true;
          if (!twoPointsOverlap(line_segm[0], line_segm[1], char_h))
            elem_segms.push_back(line_segm); // only save line_segm with length > 0
        }
      }
      else
        tangentMatchedSegms(cut_nodes,elem_segms);
    }
    else if (cut_nodes.size() >= 4)
      tangentMatchedSegms(cut_nodes,elem_segms);

    // save new cut nodes to _new_elem_segms[e] (only save the cut nodes that compose segms)
    for (unsigned int i = 0; i < cut_nodes.size(); ++i)
      if (cut_nodes[i].isInSegm) _new_elem_segms[e].add_cut_node(cut_nodes[i]);

    // save elem_segms to _new_elem_segms[e] and medial_axis
    for (unsigned int i = 0; i < elem_segms.size(); ++i)
    {
      _new_elem_segms[e].add_segm(elem_segms[i]); // save elem new segms to _new_elem_segms[e]
      bool elem_segm_exist = false;
      for (unsigned int j = 0; j < _all_medial_axis.size(); ++j)
      {
        elem_segm_exist = checkSegmOverlap(elem_segms[i],_all_medial_axis[j]);
        if (elem_segm_exist) break;
      } // j
      if (!elem_segm_exist)
      {
        _all_medial_axis.push_back(elem_segms[i]);
        _new_medial_axis.push_back(elem_segms[i]);
      }
    } // i
  } //e, loop over elements
}

void
MedialAxis2D::closestNodalPatchSearch(Point &one_inter, std::vector<Point> &side_segm, std::vector<std::vector<Point> > &isosegm,
                                      int elem_id, int side_id, std::vector<bool> &checked_node, std::vector<CutNode> &cut_nodes)
{
  // output: checked_nodes, cut_nodes, cut_tangs, cut_props
  // Find the side node closer to one_inter
  Real dist1 = sqrt((one_inter - side_segm[0]).size_sq());
  Real dist2 = sqrt((one_inter - side_segm[1]).size_sq());
  int center_loc(dist1 < dist2 ? side_id : side_id + 1); // dummy value
  if (center_loc > _aux_ix[elem_id].size()-1) center_loc = 0;
  int center_glb = _aux_ix[elem_id][center_loc]; // global node ID of the patch center

  // Compute the cut node
  if (!checked_node[center_loc])
  {
    int n_patch_node = _nodal_node_ix[center_glb][1] - _nodal_node_ix[center_glb][0] + 1;
    std::vector<Point> patch_inters;
    std::vector<Point> patch_tangs;
    std::vector<int> patch_index;

    for (unsigned int i = 0; i < n_patch_node; ++i) // loop over all patch sides
    {
      int patch_node_id = _nodal_nodes[_nodal_node_ix[center_glb][0] + i]; // global node ID
      if ((!_boundary[center_glb]) || (_boundary[center_glb] && _boundary[patch_node_id]))
      {
        std::vector<Point> patch_segm(2,Point(0.0,0.0,0.0));
        patch_segm[0] = _aux_nodes[center_glb];
        patch_segm[1] = _aux_nodes[patch_node_id];
        std::vector<Point> inters;
        std::vector<int> index; // saves isosegm IDs
        intersectLines(patch_segm,isosegm,inters,index);
        
        if (inters.size() > 0)
        {
          int imin = 0; // default value for one intersection
          if (inters.size() > 1)
          {
            std::vector<Real> dist(inters.size(),0.0);
            for (unsigned int j = 0; j < inters.size(); ++j)
              dist[j] = sqrt((inters[j]-_aux_nodes[center_glb]).size_sq());
            imin = getMinIndex(dist);
          }
          Point tang = isosegm[index[imin]][1] - isosegm[index[imin]][0];
          normalize(tang);
          patch_inters.push_back(inters[imin]);
          patch_index.push_back(index[imin]);
          patch_tangs.push_back(tang);
        } // this patch side has intersection with isosegms
      }
    } // i, loop over all patch sides

    // check if the center node is inside isocontour
    if (patch_inters.size() == 0) mooseError("patch search must find at least one intersection!");
    bool validCenter = true; // if the center node is inside the isocontour, it is valid
    int temp_count = 0;
    for (unsigned int i = 0; i < patch_inters.size(); ++i)
    {
      Point inter_normal(patch_tangs[i](1), -patch_tangs[i](0), 0.0);
      if ((_aux_nodes[center_glb] - patch_inters[i])*inter_normal < 0.0)
        temp_count += 1;
    } // i
    if (temp_count != patch_inters.size()) validCenter = false;

    // get cut nodes, their tangents and their properties
    if (validCenter)
    {
      Point cut_p = Point(0.0,0.0,0.0); // cut node coord
      if (!_boundary[center_glb]) // non-boundary center node
      {
        for (unsigned int i = 0; i < patch_inters.size(); ++i) cut_p += patch_inters[i];
        cut_p *= (1.0/patch_inters.size());
      }
      else // boundary center node
      {
        if (patch_inters.size() != 2) mooseError("why not two intersections for a boundary patch?");
        dist1 = sqrt((_aux_nodes[center_glb] - patch_inters[0]).size_sq());
        dist2 = sqrt((_aux_nodes[center_glb] - patch_inters[1]).size_sq());
        Point xp2 = _aux_nodes[center_glb];
        Point xp1(dist1 > dist2 ? patch_inters[0] : patch_inters[1]);
        Real xi(dist1 > dist2 ? dist2/dist1 : dist1/dist2);
        cut_p = 0.5*(1.0-xi)*xp1 + 0.5*(1.0+xi)*xp2;
      }

      // get cut node's tangent
      Point tang = 0.5*(patch_tangs[getMinIndex(patch_index)] - patch_tangs[getMaxIndex(patch_index)]);
      normalize(tang);

      // get cut node property
      int center_plus1(center_loc < _aux_ix[elem_id].size()-1 ? _aux_ix[elem_id][center_loc+1] : _aux_ix[elem_id][0]);
      int center_minus1(center_loc > 0 ? _aux_ix[elem_id][center_loc-1] : _aux_ix[elem_id][_aux_ix[elem_id].size()-1]);
      Point ez(0.0,0.0,1.0);

      Point normal_plus1 = (_aux_nodes[center_plus1] - _aux_nodes[center_glb]).cross(ez);
      normalize(normal_plus1);
      Real dotp_plus1 = tang*normal_plus1;

      Point normal_minus1 = (_aux_nodes[center_glb] - _aux_nodes[center_minus1]).cross(ez);
      normalize(normal_minus1);
      Real dotp_minus1 = tang*normal_minus1;

      int cut_prop = -1; // dummy value
      if (dotp_plus1*dotp_minus1 < 0.0)
        cut_prop = 2; // sweeping cut node
      else
      {
        if (dotp_plus1 > 0.0 || dotp_minus1 > 0.0)
          cut_prop = 1; // outgoing cut node
        else if (dotp_plus1 < 0.0 || dotp_minus1 < 0.0)
          cut_prop = 0; // ingoing cut node
        else
          mooseError("unable to determine the status of the cut node");
      }

      cut_nodes.push_back(CutNode(cut_p,tang,cut_prop));// save to cut_nodes
      checked_node[center_loc] = true;// mark this elem node as processed
    } // if validCenter
  } // if this center node has not been processed
}

bool
MedialAxis2D::getIsoTip(std::vector<std::vector<Point> > &iso_segms, int elem_id, CutNode cut_node, Point &tip_node)
{
  // Purpose: find any nucleation point in this element
  // get all "door segments" - the virtual segments that close the isocontour
  std::vector<std::vector<Point> > door_segms;
  for (unsigned int i = 0; i < iso_segms.size(); ++i)
  {
    int iplus1(i < iso_segms.size()-1 ? i+1 : 0);
    Real ref_len = minSegmLength(iso_segms[i], iso_segms[iplus1]);
    Point this_p1 = iso_segms[i][1];
    Point next_p0 = iso_segms[iplus1][0];
    if (!twoPointsOverlap(this_p1, next_p0, ref_len))
    {
      std::vector<Point> temp_segm(2, Point(0.0,0.0,0.0));
      temp_segm[0] = this_p1;
      temp_segm[1] = next_p0;
      door_segms.push_back(temp_segm);
    }
  }

  bool tip_found = false;
  for (unsigned int i = 0; i < door_segms.size(); ++i)
  {
    if (isInsideElem2D(door_segms[i][0], elem_id) || isInsideElem2D(door_segms[i][1], elem_id))
    {
      tip_node = 0.5*(door_segms[i][0] + door_segms[i][1]);
      tip_found = true;
      break; // N.B. we are assuming there is only one door segm in this element
    }
  }
  return tip_found;
}

bool
MedialAxis2D::isInsideElem2D(Point &p, int elem_id)
{
  int nen = _aux_ix[elem_id].size();
  int count = 0;
  for (unsigned int i = 0; i < nen; ++i)
  {
    int node1 = _aux_ix[elem_id][i];
    int node2(i < nen-1 ? _aux_ix[elem_id][i+1] : _aux_ix[elem_id][0]);
    Point middle2p = p - 0.5*(_aux_nodes[node1] + _aux_nodes[node2]);
    Point side_tang = _aux_nodes[node2] - _aux_nodes[node1];
    Point side_norm = Point(side_tang(1),-side_tang(0),0.0);

    normalize(middle2p);
    normalize(side_norm);
    if (middle2p*side_norm < 0.0) count += 1;
  } // i

  bool isInside = false;
  if (count == nen) isInside = true;
  return isInside;
}

void
MedialAxis2D::correctNewCutNode(CutNode &cut_node, int elem_id)
{
  // Purpose: For a new cut_node in aux elem elem_id, check if it has a counterpart in any
  //          generalized neighbor element. If so, move this new cut_node to the counterpart

  // get all generalized neighbors' ID
  std::set<unsigned int> neighbor_ids;
  for (unsigned int i = 0; i < _aux_ix[elem_id].size(); ++i)
  {
    unsigned int start_index = _nodal_elem_ix[_aux_ix[elem_id][i]][0];
    unsigned int end_index = _nodal_elem_ix[_aux_ix[elem_id][i]][1];
    for (unsigned int j = start_index; j <= end_index; ++j)
    {
      if (_nodal_elems[j] != elem_id)
        neighbor_ids.insert(_nodal_elems[j]);
    } // j
  } // i

  // loop over all generalized neighbors
  std::set<unsigned int>::iterator sit;
  for (sit = neighbor_ids.begin(); sit != neighbor_ids.end(); ++sit)
  {
    int neigh_elem_id = (*sit);
    std::vector<CutNode> &neigh_cut_nodes = _all_elem_segms[neigh_elem_id].cut_nodes;
    bool point_corrected = false;
    for (unsigned int i = 0; i < neigh_cut_nodes.size(); ++i)
    {
      Real dist = std::sqrt((cut_node.getPoint()-neigh_cut_nodes[i].getPoint()).size_sq());
      if ((dist > 1.0e-12 && dist < _iso_width) && cut_node.hasSimilarPropWith(neigh_cut_nodes[i]))
      {
        cut_node.switchPoint(neigh_cut_nodes[i].getPoint());
        point_corrected = true;
        break;
      }
    } // i
    if (point_corrected) break;
  } // sit
}

bool
MedialAxis2D::getOldSegmInterNewCut(CutNode new_cut, int elem_id, Point &p_int)
{
  bool interFound = false;
  std::vector<Point> temp_line(2, Point(0.0,0.0,0.0));
  temp_line[0] = new_cut.getPoint();
  temp_line[1] = new_cut.getPoint() + new_cut.getTang();
  std::vector<std::vector<Point> > &old_segms = _all_elem_segms[elem_id].medial_segms;

  // check if this new_cut can intersection with any old_segms
  Real d_min = 1.0e15; // dummy value
  Point temp_p(0.0,0.0,0.0);
  for (unsigned int i = 0; i < old_segms.size(); ++i)
  {
    if (intersectLineAndSegm(temp_line, old_segms[i], temp_p))
    {
      interFound = true;
      Real d = sqrt((new_cut.getPoint()-temp_p).size_sq());
      if (d < d_min)
      {
        d_min = d;
        p_int = temp_p;
      }
    }
  } // i
  return interFound;
}

bool
MedialAxis2D::getBranchPoint(std::vector<std::vector<Point> > &iso_segms, std::vector<CutNode> &cut_nodes,
                             Point &branch, int elem_id)
{
  if (cut_nodes.size() != 3) mooseError("getBranchPoint can only be called with three cut nodes");
  bool branchFound = false;
  branch = Point(0.0,0.0,0.0);

  int n_ingoing = 0;
  int n_outgoing = 0;
  int id_ingoing = -1;
  for (unsigned int i = 0; i < cut_nodes.size(); ++i)
  {
    if (cut_nodes[i].getProp() == 0)
    {
      n_ingoing += 1;
      id_ingoing = i;
    }
    else if (cut_nodes[i].getProp() == 1)
      n_outgoing += 1;
  }

  if ((n_ingoing == 1) && (n_outgoing == 2)) // 1 ingoing cut node + 2 outgoing cut nodes
  {
    branchFound = true;
    std::vector<std::vector<Point> > three_vecs;
    for (unsigned int i = 0; i < cut_nodes.size(); ++i)
    {
      std::vector<Point> line_vector(2,Point(0.0,0.0,0.0));
      line_vector[0] = cut_nodes[i].getPoint();
      line_vector[1] = cut_nodes[i].getPoint() + cut_nodes[i].getTang();
      three_vecs.push_back(line_vector);
    }

    int n_inter = 0;
    for (unsigned int i = 0; i < cut_nodes.size()-1; ++i)
    {
      for (unsigned int j = i+1; j < cut_nodes.size(); ++j)
      {
        Point p_int(0.0,0.0,0.0);
        bool isIntersect = intersectTwoLines(three_vecs[i],three_vecs[j],false,p_int);
        if (isIntersect && isInsideElem2D(p_int, elem_id) && isInsideContour(p_int, iso_segms))
        {
          branch += p_int;
          n_inter += 1;
        }
      } // j
    } // i

    if (n_inter > 0)
      branch *= (1.0/n_inter);
    else
      branch = cut_nodes[id_ingoing].getPoint();
  } // 1 ingoing and 2 outgoings

  return branchFound;
}

void 
MedialAxis2D::tangentMatchedSegms(std::vector<CutNode> &cut_nodes, std::vector<std::vector<Point> > &elem_segms)
{
  if (cut_nodes.size() < 3) mooseError("tangentMatchedSegms cannot be implemented for fewer than 3 cut nodes");
  
  std::vector<std::vector<int> > all_segm_id;
  for (unsigned int j = 0; j < cut_nodes.size(); ++j)
  {
    std::vector<int> the_segm_id(2,0);
    bool the_segm_found = false;
    Real alpha_min = 1000.0; // dummy value
    for (unsigned int k = 0; k < cut_nodes.size(); ++k)
    {
      if ((k != j) && isValidTwoCutConnection(cut_nodes[j], cut_nodes[k]))
      {
        Point segm_dir = cut_nodes[k].getPoint() - cut_nodes[j].getPoint();
        Real alpha = std::max(getVectorAngle(segm_dir, cut_nodes[j].getTang(), true),
                              getVectorAngle(segm_dir, cut_nodes[k].getTang(), true));
        if (alpha < alpha_min)
        {
          alpha_min = alpha;
          the_segm_id[0] = j;
          the_segm_id[1] = k;
          the_segm_found = true;
        }
      }
    } // k

    if (the_segm_found)
    {
      bool the_segm_exist = false;
      for (unsigned int k = 0; k < all_segm_id.size(); ++k)
      {
        if (((the_segm_id[0] == all_segm_id[k][0]) && (the_segm_id[1] == all_segm_id[k][1])) || 
            ((the_segm_id[0] == all_segm_id[k][1]) && (the_segm_id[1] == all_segm_id[k][0])))
        {
          the_segm_exist = true;
          break;
        }
      }
      if (!the_segm_exist) all_segm_id.push_back(the_segm_id);
    }
  } // j
 
  // check if any cut node forms more than one segms. If so, remove invalid segms
  for (unsigned int j = 0; j < cut_nodes.size(); ++j)
  {
    std::vector<std::vector<int> > segms_jcut; // segms that contains the j-th cut node
    for (unsigned int k = 0; k < all_segm_id.size(); ++k)
    {
      if ((all_segm_id[k][0] == j) || (all_segm_id[k][1] == j))
        segms_jcut.push_back(all_segm_id[k]);
    }
    if (segms_jcut.size() > 1)
      removeInvalidSegms(cut_nodes, segms_jcut, all_segm_id);
  }

  // get medial segments 
  for (unsigned int j = 0; j < all_segm_id.size(); ++j)
  {
    std::vector<Point> line_segm(2,Point(0.0,0.0,0.0));
    line_segm[0] = cut_nodes[all_segm_id[j][0]].getPoint();
    line_segm[1] = cut_nodes[all_segm_id[j][1]].getPoint();
    cut_nodes[all_segm_id[j][0]].isInSegm = true; // important
    cut_nodes[all_segm_id[j][1]].isInSegm = true;
    elem_segms.push_back(line_segm);
  }
}

void
MedialAxis2D::removeInvalidSegms(std::vector<CutNode> &cut_nodes, std::vector<std::vector<int> > &segms_jcut,
                                 std::vector<std::vector<int> > &all_segm_id)
{
  // only pick the most aligned pair from segms_jcut
  std::vector<Real> alpha;
  for (unsigned int j = 0; j < segms_jcut.size(); ++j)
  {
    Point segm_dir = cut_nodes[segms_jcut[j][1]].getPoint() - cut_nodes[segms_jcut[j][0]].getPoint();
    Real a = std::max(getVectorAngle(segm_dir, cut_nodes[segms_jcut[j][0]].getTang(), true),
                      getVectorAngle(segm_dir, cut_nodes[segms_jcut[j][1]].getTang(), true));
    alpha.push_back(a);
  }
  int imin = getMinIndex(alpha);
  segms_jcut.erase(segms_jcut.begin() + imin); // now it contains all invalid segments
  std::vector<int> rm_ix;
  for (unsigned int j = 0; j < all_segm_id.size(); ++j)
  {
    for (unsigned int k = 0; k < segms_jcut.size(); ++k)
    {
      if ((all_segm_id[j][0] == segms_jcut[k][0]) && (all_segm_id[j][1] == segms_jcut[k][1]))
      {
        rm_ix.push_back(j);
        break;
      }
    } // k
  } // j
  shrinkArray(all_segm_id,rm_ix);
}

bool
MedialAxis2D::checkSegmOverlap(std::vector<Point> &segm1, std::vector<Point> &segm2)
{
  Real len = minSegmLength(segm1, segm2); // reference length
  bool isOverlap = false;
  if ((twoPointsOverlap(segm1[0],segm2[0],len) && twoPointsOverlap(segm1[1],segm2[1],len)) ||
      (twoPointsOverlap(segm1[0],segm2[1],len) && twoPointsOverlap(segm1[1],segm2[0],len)))
    isOverlap = true;
  return isOverlap;
}

bool
MedialAxis2D::intersectTwoLines(std::vector<Point> &segm1, std::vector<Point> &segm2, bool finite, Point &p_int)
{
  // see stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
  // Get p, q, r, s
  Point p = segm1[0];
  Point q = segm2[0];
  Point r = segm1[1] - segm1[0];
  Point s = segm2[1] - segm2[0];

  // Get the intersection
  bool haveIntersection = false;
  p_int = Point(0.0,0.0,0.0);
  Real rxs = (r.cross(s))(2);
  Real tol = 1.0e-14;
  if (std::abs(rxs) > tol) // r x s != 0
  {
    Real q_pxr = ((q - p).cross(r))(2);
    Real q_pxs = ((q - p).cross(s))(2);
    Real t = q_pxs/rxs;
    Real u = q_pxr/rxs;
    p_int = p + t*r;
    if (finite)
    {  
      if ((t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0))
        haveIntersection = true;
    }
    else
      haveIntersection = true;
  }
  return haveIntersection;
}

bool
MedialAxis2D::intersectLineAndSegm(std::vector<Point> &segm1, std::vector<Point> &segm2, Point &p_int)
{
  // see stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
  // N.B. segm1 is an infinite line and segm2 is a finite segment
  // Get p, q, r, s
  Point p = segm1[0];
  Point q = segm2[0];
  Point r = segm1[1] - segm1[0];
  Point s = segm2[1] - segm2[0];

  // Get the intersection
  bool haveIntersection = false;
  p_int = Point(0.0,0.0,0.0);
  Real rxs = (r.cross(s))(2);
  Real tol = 1.0e-14;
  if (std::abs(rxs) > tol) // r x s != 0
  {
    Real q_pxr = ((q - p).cross(r))(2);
    Real q_pxs = ((q - p).cross(s))(2);
    Real t = q_pxs/rxs; // ratio of segm1
    Real u = q_pxr/rxs; // ratio of segm2
    p_int = p + t*r;
    if (u >= 0.0 && u <= 1.0)
      haveIntersection = true;
  }
  return haveIntersection;
}

void 
MedialAxis2D::intersectLines(std::vector<Point> &segm1, std::vector<std::vector<Point> > &segms2, 
                             std::vector<Point> &p_inters, std::vector<int> &index2)
{
  // TODO: this function is very slow. We need some comp sci experts to optimize it
  // get intersections
  index2.clear();
  p_inters.clear();

  for (unsigned int j = 0; j < segms2.size(); ++j)
  {
    Point p_int = Point(0.0,0.0,0.0);
    bool haveIntersect = intersectTwoLines(segm1, segms2[j], true, p_int);
    if (haveIntersect)
    {
      bool point_exist = false;
      Real tol = 1.0e-14;
      for (unsigned int k = 0; k < p_inters.size(); ++k)
      {
        Real dist = sqrt((p_inters[k] - p_int).size_sq());
        if (dist < tol)
        {
          point_exist = true;
          break;
        }
      }
      if (!point_exist)
      {
        index2.push_back(j);
        p_inters.push_back(p_int);
      }
    }
  } // j
}

void 
MedialAxis2D::sortAscend(std::vector<Real> &v, std::vector<int> &index)
{
  if (v.size() == 0)
    mooseError("The input vector is empty!");

  index.clear();
  index.resize(v.size(),0);
  for (unsigned int i = 0; i < index.size(); ++i) 
    index[i] = i;
  
  if (v.size() > 1)
  {
    for (int j = v.size()-2; j >= 0; --j)
    {
      bool swapped = false;
      for (unsigned int i = 0; i <= j; ++i)
      {
        if (v[i] > v[i+1])
        {
          Real temp = v[i]; // swap v
          v[i] = v[i+1];
          v[i+1] = temp;
          int itemp = index[i]; // swap index
          index[i] = index[i+1];
          index[i+1] = itemp;
          swapped = true;
        }
      }
      if (!swapped) break;
    }
  }
}

template <class T> int
MedialAxis2D::getMinIndex(std::vector<T> &v)
{
  int imin = -1;
  if (v.size() > 0)
  {
    imin = 0;
    T min_val = v[0];
    for (unsigned int i = 1; i < v.size(); ++i)
    {
      if (v[i] < min_val)
      {
        min_val = v[i];
        imin = i;
      }
    }
  }
  return imin;
}

template <class T> int
MedialAxis2D::getMaxIndex(std::vector<T> &v)
{
  int imax = -1;
  if (v.size() > 0)
  {
    imax = 0;
    T max_val = v[0];
    for (unsigned int i = 1; i < v.size(); ++i)
    {
      if (v[i] > max_val)
      {
        max_val = v[i];
        imax = i;
      }
    }
  }
  return imax;
}

template <class T>
void
MedialAxis2D::reArrange(std::vector<T> &v, std::vector<int> &index)
{
  // rearrange the elements in v according to index
  if (v.size() > 0)
  {
    if (v.size() != index.size())
      mooseError("v's size should be equal to index's size");

    std::vector<T> v_copy;
    for (unsigned int i = 0; i < v.size(); ++i)
      v_copy.push_back(v[index[i]]);

    for (unsigned int i = 0; i < v_copy.size(); ++i)
      v[i] = v_copy[i];
  }
}

template <class T>
void
MedialAxis2D::shrinkArray(std::vector<T> &v, std::vector<int> &rm_ix)
{
  // delete the indexed array elements
  if ((v.size() > 0) && (rm_ix.size() > 0))
  {
    std::vector<T> v_copy;
    for (unsigned int j = 0; j < v.size(); ++j)
    {
      bool needToBeRemoved = false;
      for (unsigned int k = 0; k < rm_ix.size(); ++k)
        if (j == rm_ix[k])
        {
          needToBeRemoved = true;
          break;
        }
      if (!needToBeRemoved) v_copy.push_back(v[j]);
    } // j
    v.clear();
    for (unsigned int j = 0; j < v_copy.size(); ++j)
      v.push_back(v_copy[j]);
  }
}

Real
MedialAxis2D::getVectorAngle(Point v1, Point v2, bool acute_angl)
{
  Real dotp = v1*v2;
  Real norm1 = std::sqrt(v1.size_sq());
  Real norm2 = std::sqrt(v2.size_sq());
  if (norm1 < 1.0e-14 || norm2 < 1.0e-14)
    mooseError("vector norm is zero");
  else
    dotp /= (norm1*norm2);
  Real alpha = (180.0/3.14159265359)*safeAcos(dotp);
  if (acute_angl && alpha > 90.0) alpha = 180.0 - alpha;
  return alpha;
}

Real
MedialAxis2D::getMinSideLength(int elem_id)
{
  int nen = _aux_ix[elem_id].size();
  Real len_min = 1.0e20;
  for (unsigned int s = 0; s < nen; ++s)
  {
    int s_node1 = s;
    int s_node2(s < nen-1 ? s+1 : 0);
    Point side_p1 = _aux_nodes[_aux_ix[elem_id][s_node1]];
    Point side_p2 = _aux_nodes[_aux_ix[elem_id][s_node2]];
    Real len = sqrt((side_p1-side_p2).size_sq());
    if (len < len_min) len_min = len;
  }
  return len_min;
}

Real
MedialAxis2D::getMaxSideLength(int elem_id)
{
  int nen = _aux_ix[elem_id].size();
  Real len_max = -1.0;
  for (unsigned int s = 0; s < nen; ++s)
  {
    int s_node1 = s;
    int s_node2(s < nen-1 ? s+1 : 0);
    Point side_p1 = _aux_nodes[_aux_ix[elem_id][s_node1]];
    Point side_p2 = _aux_nodes[_aux_ix[elem_id][s_node2]];
    Real len = sqrt((side_p1-side_p2).size_sq());
    if (len > len_max) len_max = len;
  }
  return len_max;
}

Real
MedialAxis2D::safeAcos(Real x)
{
  if (x < -1.0) x = -1.0;
  else if (x > 1.0) x = 1.0;
  return std::acos(x);
}

void
MedialAxis2D::printMedialAxis()
{
  std::cout << "****** All medial axis ******" << std::endl;
  for (unsigned int i = 0; i < _all_medial_axis.size(); ++i)
    std::cout << _all_medial_axis[i][0](0) << ", " <<  _all_medial_axis[i][0](1) 
              << ", " << _all_medial_axis[i][1](0) << ", " << _all_medial_axis[i][1](1) << std::endl;
  std::cout << "****** New medial axis ******" << std::endl;
  for (unsigned int i = 0; i < _new_medial_axis.size(); ++i)
    std::cout << _new_medial_axis[i][0](0) << ", " <<  _new_medial_axis[i][0](1) 
              << ", " << _new_medial_axis[i][1](0) << ", " << _new_medial_axis[i][1](1) << std::endl;
}

void
MedialAxis2D::writeMedialAxis(std::vector<Isocontour> & isocontours)
{
  if (isocontours.size() == 0)
    return;

  // get the file name for iscontour and medial axis
  std::string file_id = numberToString(_t_step);
  std::string iso_filename = "iso_contour" + file_id + ".txt";
  std::string axis_filename = "medial_axis" + file_id + ".txt";

  // write isocontour
  std::ofstream contour_file;
  contour_file.open(iso_filename.c_str());
  for (unsigned int i = 0; i < isocontours.size(); ++i)
  {
    std::vector<std::vector<Point> > segms;
    isocontours[i].getIsoSegms(segms);
    for (unsigned int j = 0; j < segms.size(); ++j)
    {
      contour_file << segms[j][0](0) << ", " << segms[j][0](1) << ", "
                   << segms[j][1](0) << ", " << segms[j][1](1) << "\n";
    }
  } 
  contour_file.close();

  // write medial axis
  std::ofstream axis_file;
  axis_file.open(axis_filename.c_str());
  for (unsigned int i = 0; i < _all_medial_axis.size(); ++i)
  {
    axis_file << _all_medial_axis[i][0](0) << ", " <<  _all_medial_axis[i][0](1) 
              << ", " << _all_medial_axis[i][1](0) << ", " << _all_medial_axis[i][1](1) << "\n";
  }
  axis_file.close();
}

template <typename T> std::string
MedialAxis2D::numberToString(T number)
{
  std::ostringstream ss;
  ss << std::setfill('0') << std::setw(4) << number;
  return ss.str();
}

bool
MedialAxis2D::isValidTwoCutConnection(CutNode &cut_node1, CutNode &cut_node2)
{
  bool isValid = true;
  Real angle_c(_single_crack ? 90.0 : 30.0); // if only one single crack, no critical angle restriction
  Point segm_dir = cut_node2.getPoint() - cut_node1.getPoint();
  Real alpha = std::max(getVectorAngle(segm_dir, cut_node1.getTang(), true),
                        getVectorAngle(segm_dir, cut_node2.getTang(), true));
  if (cut_node1.getProp() == cut_node2.getProp()) // can't be both ingoing, outgoing or sweeping
    isValid = false;
  else if (cut_node1.getTang()*cut_node2.getTang() <= 0.0) // can't have opposite tangent vectors
    isValid = false;
  else if (alpha >= angle_c) // alignment measurement can't be too large -- 30.0 is changeable
    isValid = false;
  return isValid;
}

bool
MedialAxis2D::isInsideContour(Point p, std::vector<std::vector<Point> > &iso_segms)
{
  // get strictly closed isocontour
  std::vector<std::vector<Point> > new_segms;
  for (unsigned int i = 0; i < iso_segms.size(); ++i)
  {
    new_segms.push_back(iso_segms[i]);
    int iplus1(i < iso_segms.size()-1 ? i+1 : 0);
    Real ref_len = minSegmLength(iso_segms[i], iso_segms[iplus1]);
    Point this_p1 = iso_segms[i][1];
    Point next_p0 = iso_segms[iplus1][0];
    if (!twoPointsOverlap(this_p1, next_p0, ref_len))
    {
      std::vector<Point> temp_segm(2, Point(0.0,0.0,0.0));
      temp_segm[0] = this_p1;
      temp_segm[1] = next_p0;
      new_segms.push_back(temp_segm);
    }
  }

  // use winding number test to determine if p is inside iso_segms
  if (wn_PnPoly(p, new_segms) == 0)
    return false;
  else
    return true;
}
