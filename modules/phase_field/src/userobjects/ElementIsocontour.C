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

#include "ElementIsocontour.h"
#include "MedialAxis2DFuncs.h"
#include <math.h>

template<>
InputParameters validParams<ElementIsocontour>()
{
  InputParameters params = validParams<ElementUserObject>();

  params.addRequiredCoupledVar("phase_field", "Phase field variable");
  params.addParam<Real>("iso_val",0.0,"The iso-value required to generate the isocontour of phase field");
  params.addParam<Real>("l", 0.0, "phase field characteristic length");
  return params;
}

ElementIsocontour::ElementIsocontour(const InputParameters & parameters) :
    ElementUserObject(parameters),
    _iso_val(getParam<Real>("iso_val")),
    _iso_width(2.0*getParam<Real>("l")),
    _p_pf(&coupledNodalValue("phase_field", 0)), // pointer of field field variable
    _xfem(_fe_problem.get_xfem())
{
}

Real
ElementIsocontour::getIsoVal() const
{
  return _iso_val;
}

void
ElementIsocontour::getIsoContours(std::vector<Isocontour> & isocontours) const
{
  isocontours = _isocontours;
}

void
ElementIsocontour::initialize()
{
  _iso_segms.clear();
  _iso_index.clear();
  _isocontours.clear();
}

void
ElementIsocontour::execute()
{
  // Get element nodes
  unsigned int _n_nodes = _current_elem->n_nodes();
  std::vector<Node*> _nodes(_n_nodes,NULL);
  for (unsigned int i = 0; i < _n_nodes; ++i)
    _nodes[i] = _current_elem->get_node(i);

  // Get the isocontour cut points on all element sides
  std::vector<Point> iso_cuts;
  for (unsigned int i = 0; i < _n_nodes; ++i)
  {
    int node_id1 = i;
    int node_id2(i < _n_nodes-1 ? i+1 : 0);
    if (((*_p_pf)[node_id1] - _iso_val)*((*_p_pf)[node_id2] - _iso_val) < 0.0)
    {
      Real xi = (2.0*_iso_val - (*_p_pf)[node_id1] - (*_p_pf)[node_id2])/((*_p_pf)[node_id2] - (*_p_pf)[node_id1]);
      Real p_x = 0.5*(1.0-xi)*(*_nodes[node_id1])(0) + 0.5*(1.0+xi)*(*_nodes[node_id2])(0);
      Real p_y = 0.5*(1.0-xi)*(*_nodes[node_id1])(1) + 0.5*(1.0+xi)*(*_nodes[node_id2])(1);
      iso_cuts.push_back(Point(p_x, p_y, 0.0));
    }
  }

  // Get the isocontour segments in this element
  std::vector<std::vector<Point> > elem_segms;
  if (!iso_cuts.empty())
  {
    if (iso_cuts.size() == 2)
      elem_segms.push_back(iso_cuts);
    else if (iso_cuts.size() == 4)
    {
      std::vector<Point> iso_segm_line(2,Point(0.0,0.0,0.0));
      if (0.5*((*_p_pf)[0] + (*_p_pf)[2]) > 0.5*((*_p_pf)[1] + (*_p_pf)[3]))
      {
        iso_segm_line[0] = iso_cuts[0];
        iso_segm_line[1] = iso_cuts[1];
        elem_segms.push_back(iso_segm_line); // first iso-segm
        iso_segm_line[0] = iso_cuts[2];
        iso_segm_line[1] = iso_cuts[3];
        elem_segms.push_back(iso_segm_line); // second iso-segm
      }
      else
      {
        iso_segm_line[0] = iso_cuts[0];
        iso_segm_line[1] = iso_cuts[3];
        elem_segms.push_back(iso_segm_line); // first iso-segm
        iso_segm_line[0] = iso_cuts[1];
        iso_segm_line[1] = iso_cuts[2];
        elem_segms.push_back(iso_segm_line); // second iso-segm
      }
    }
    else
      mooseError("The number of isocontour cut points is neither 2 nor 4!");
  }

  // Truncate each elem segm if current elem is a partial elem, and save elem_segms to _iso_segms
  for (unsigned int i = 0; i < elem_segms.size(); ++i)
  {
    if (truncateSegm(elem_segms[i]))
      _iso_segms.push_back(elem_segms[i]);
  }
}


void
ElementIsocontour::threadJoin(const UserObject & y)
{
  const ElementIsocontour & uo = static_cast<const ElementIsocontour &>(y);
  const std::vector<std::vector<Point> > & this_segms = uo._iso_segms;
  for (unsigned int i = 0; i < this_segms.size(); ++i)
    _iso_segms.push_back(this_segms[i]);
}

void
ElementIsocontour::finalize()
{
  // gather _iso_segms from different processors
  std::vector<Real> isosegms_tmp;
  flattenPointVector(_iso_segms, isosegms_tmp);
  _communicator.allgather(isosegms_tmp, false);
  reshapePointVector(isosegms_tmp, _iso_segms);

  // connect adjacent iso-segments
  if (!_iso_segms.empty())
    sort_segms(_iso_segms, _iso_index);

  // generate vector of isocontour object
  for (unsigned int i = 0; i < _iso_index.size(); ++i) // loop over separate isocontours
  {
    std::vector<std::vector<Point> > one_isosegms;
    for (unsigned int j = _iso_index[i][0]; j <= _iso_index[i][1]; ++j)
      one_isosegms.push_back(_iso_segms[j]);
    _isocontours.push_back(Isocontour(_iso_width, one_isosegms));
  }

  // close open isocontours if necessary
  close_isocontours(); // _isocontours changed here

  // force each isocontour be counterclockwise
  for (unsigned int i = 0; i < _isocontours.size(); ++i)
  {
    if (_isocontours[i].isClockWise())
      _isocontours[i].reverse();
  }
}

void
ElementIsocontour::sort_segms(std::vector<std::vector<Point> > &segms, std::vector<std::vector<int> > &index)
{ 
  index.push_back(std::vector<int> (2,0));
  int n_contour = 0; // number of independent contour - 1
  bool reverse = false;
  int last_iter = 0;

  // sort the iso-segments
  for (unsigned int i = 0; i < segms.size()-1; ++i)
  {
    Point p = segms[i][1];
    int Q = i + 1;
    int R = -1; // dummy value
    bool isSamePoint1 = false;
    bool isSamePoint2 = false;
    for (unsigned int j = Q; j < segms.size(); ++j) // search from Q to the end
    {
      Real len = minSegmLength(segms[i], segms[j]);
      isSamePoint1 = twoPointsOverlap(p, segms[j][0], len);
      isSamePoint2 = twoPointsOverlap(p, segms[j][1], len);
      if (isSamePoint1 || isSamePoint2)
      {
        R = j;
        break;
      }
    } // j
    if (R >= 0) // a connecting segment found
    {
      if (R != Q) // swap R-segm and Q-segm
      {
        std::vector<Point> segm_temp(2, Point(0.0,0.0,0.0));
        segm_temp[0] = segms[Q][0];
        segm_temp[1] = segms[Q][1];
        if (isSamePoint1) // R's first end point is the connection point
        {
          segms[Q][0] = segms[R][0];
          segms[Q][1] = segms[R][1];
        }
        else if (isSamePoint2) // R's second end point is the connection point
        {
          segms[Q][0] = segms[R][1];
          segms[Q][1] = segms[R][0];
        }
        else
          mooseError("No connecting segment found");

        segms[R][0] = segm_temp[0];
        segms[R][1] = segm_temp[1];        
      }
      else // if R == Q
      {
        if (isSamePoint2)
        {
          Point p_temp = segms[R][0];
          segms[R][0] = segms[R][1];
          segms[R][1] = p_temp;
        }
      }
    }
    else // if R == -1
    {
      int i0 = index[n_contour][0]; // first index of this contour
      p = segms[i0][0];
      for (unsigned int j = Q; j < segms.size(); ++j) // search connecting segm in the other direction
      {
        Real len = minSegmLength(segms[i0], segms[j]);
        isSamePoint1 = twoPointsOverlap(p, segms[j][0], len);
        isSamePoint2 = twoPointsOverlap(p, segms[j][1], len);
        if (isSamePoint1 || isSamePoint2)
        {
          R = j;
          break;
        }
      } // j
      if (R >= 0) // connecting segment found in the other direction
      {
        reverse = true;
        std::vector<Point> segm_temp(2, Point(0.0,0.0,0.0));
        if (i > i0) // reverse the segms from i0 to i in current contour
        {
          for (unsigned int j = 0; j <= floor(0.5*(i-i0-1)); ++j) // reverse
          {
            segm_temp[0] = segms[j+i0][0];
            segm_temp[1] = segms[j+i0][1];
            segms[j+i0][0] = segms[i-j][0];
            segms[j+i0][1] = segms[i-j][1];
            segms[i-j][0] = segm_temp[0];
            segms[i-j][1] = segm_temp[1];          
          } // j
        }
        for (unsigned int j = i0; j <= i; ++j) // swap two end points
        {
          Point p_temp = segms[j][0];
          segms[j][0] = segms[j][1];
          segms[j][1] = p_temp;
        } // j
        if (R != Q) // swap R-segm and Q-segm
        {
          segm_temp[0] = segms[Q][0];
          segm_temp[1] = segms[Q][1];
          if (isSamePoint1) // R's first end point is the connection point
          {
            segms[Q][0] = segms[R][0];
            segms[Q][1] = segms[R][1];
          }
          else if (isSamePoint2) // R's second end point is the connection point
          {
            segms[Q][0] = segms[R][1];
            segms[Q][1] = segms[R][0];
          }
          else
            mooseError("No connecting segment found");
        
          segms[R][0] = segm_temp[0];
          segms[R][1] = segm_temp[1];        
        }
        else // if R == Q
        {
          if (isSamePoint2)
          {
            Point p_temp = segms[R][0];
            segms[R][0] = segms[R][1];
            segms[R][1] = p_temp;
          }
        }
      }
      else // R == -1, connecting segm can't be found anyway
      {
        if (reverse)
        {
          if (i > i0) // reverse again the segms from i0 to i
          {
            std::vector<Point> segm_temp(2, Point(0.0,0.0,0.0));
            for (unsigned int j = 0; j <= floor(0.5*(i-i0-1)); ++j) // reverse
            {
              segm_temp[0] = segms[j+i0][0];
              segm_temp[1] = segms[j+i0][1];
              segms[j+i0][0] = segms[i-j][0];
              segms[j+i0][1] = segms[i-j][1];
              segms[i-j][0] = segm_temp[0];
              segms[i-j][1] = segm_temp[1];
            } // j
          }
          for (unsigned int j = i0; j <= i; ++j) // swap two end points
          {
            Point p_temp = segms[j][0];
            segms[j][0] = segms[j][1];
            segms[j][1] = p_temp;
          } // j
          reverse = false;
        }
        index[n_contour][1] = i;
        n_contour = n_contour + 1;
        index.push_back(std::vector<int> (2,0));
        index[n_contour][0] = i + 1;
      }
    }
    last_iter = i;
  } // i
  index[n_contour][1] = last_iter + 1;
}

void
ElementIsocontour::close_isocontours()
{
  // Purpose: find another isocontour that can pair with the current one
  bool check_contours = true;
  while (check_contours)
  {
    bool has_non_closed = false;
    unsigned int non_closed_id = 0;
    for (unsigned int i = 0; i < _isocontours.size(); ++i)
    {
      if (!_isocontours[i].isClosed())
      {
        non_closed_id = i;
        has_non_closed = true;
        break;
      }
    } // i
    if (has_non_closed)
      combine_isocontours(non_closed_id); // _isocontours changed here
    else
      check_contours = false;
  } // while check_contours = true
}

void
ElementIsocontour::combine_isocontours(unsigned int iso_id)
{
  // Purpose: find a sibling isocontour that can pair with _isocontours[iso_id], and combine
  //          these two isocontours together
  bool iso_sibling_found = false;
  unsigned int sibling_id = 0;
  Point last_p = _isocontours[iso_id].getLastPoint();
  Real dist_min = 1.0e20; // dummy big value
  Real search_radius = 1.1*_iso_width; // changeable
  for (unsigned int i = 0; i < _isocontours.size(); ++i)
  {
    if (i != iso_id)
    {
      Real dist1 = std::sqrt((last_p - _isocontours[i].getFirstPoint()).size_sq());
      Real dist2 = std::sqrt((last_p - _isocontours[i].getLastPoint()).size_sq());
      Real dist = std::min(dist1, dist2);
      if (dist < dist_min && dist < search_radius)
      {
        dist_min = dist;
        sibling_id = i;
        iso_sibling_found = true;
      }
    }
  }
  if (iso_sibling_found)
  {
    Real dist1 = std::sqrt((last_p - _isocontours[sibling_id].getFirstPoint()).size_sq());
    Real dist2 = std::sqrt((last_p - _isocontours[sibling_id].getLastPoint()).size_sq());
    if (dist2 < dist1) _isocontours[sibling_id].reverse(); // re-order sibling isocontour
    _isocontours[iso_id].addIsoSegmsFrom(_isocontours[sibling_id]); // combine the two isocontours
    _isocontours.erase(_isocontours.begin() + sibling_id); // delete the sibling isocontour
  }
  else
    mooseError("Isocontour " << iso_id << " is open but we can't find a sibling isocontour to close it!");
}

void
ElementIsocontour::flattenPointVector(std::vector<std::vector<Point> > &segms, std::vector<Real> &single_v)
{
  single_v.clear();
  for (unsigned int i = 0; i < segms.size(); ++i)
    for (unsigned int j = 0; j < segms[i].size(); ++j)
      for (unsigned int k = 0; k < 3; ++k)
        single_v.push_back(segms[i][j](k));
}

void
ElementIsocontour::reshapePointVector(std::vector<Real> &single_v, std::vector<std::vector<Point> > &segms)
{
  if (single_v.size() % 6 != 0)
    mooseError("Incorrect size of input single vector");

  segms.clear();
  unsigned int n_row = single_v.size()/6;
  for (unsigned int i = 0; i < n_row; ++i)
  {
    std::vector<Point> line_v(2, Point(0.0,0.0,0.0));
    line_v[0] = Point(single_v[6*i],   single_v[6*i+1], single_v[6*i+2]);
    line_v[1] = Point(single_v[6*i+3], single_v[6*i+4], single_v[6*i+5]);
    segms.push_back(line_v);
  }
}

bool
ElementIsocontour::truncateSegm(std::vector<Point> &iso_segm)
{
  bool valid_segm_obtained = false;
  if (_xfem->is_elem_cut(_current_elem))
  {
    std::vector<std::vector<Point> > frag_edges;
    _xfem->get_frag_faces(_current_elem, frag_edges); // get fragment edges

    std::vector<Point> truncated_segm;
    for (unsigned int i = 0; i < frag_edges.size(); ++i)
    {
      Point p_int(0.0,0.0,0.0);
      Real h_char = minSegmLength(iso_segm, frag_edges[i]);
      if (intersectTwoSegms(iso_segm, frag_edges[i], p_int) && (!pointExistInSegm(p_int, truncated_segm, h_char)))
        truncated_segm.push_back(p_int); // get truncated isocontour segment
    }

    iso_segm.clear();
    if (truncated_segm.size() == 2)  // save truncated isocontour segment to iso_segm
      iso_segm = truncated_segm;
  }
  if (!iso_segm.empty())
    valid_segm_obtained = true;
  return valid_segm_obtained;
}

bool
ElementIsocontour::intersectTwoSegms(std::vector<Point> &segm1, std::vector<Point> &segm2, Point &p_int)
{
  // see stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
  // Get p, q, r, s
  Point p = segm1[0];
  Point q = segm2[0];
  Point r = segm1[1] - segm1[0];
  Point s = segm2[1] - segm2[0];

  // Get the intersection
  Real tol = 1.0e-8;
  bool haveIntersection = false;
  p_int = Point(0.0,0.0,0.0);
  Real rxs = (r.cross(s))(2);
  if (std::abs(rxs) > tol) // r x s != 0
  {
    Real q_pxr = ((q - p).cross(r))(2);
    Real q_pxs = ((q - p).cross(s))(2);
    Real t = q_pxs/rxs;
    Real u = q_pxr/rxs;
    p_int = p + t*r;
    if ((t > -tol && t < 1.0+tol) && (u > -tol && u < 1.0+tol))
      haveIntersection = true;
  }
  return haveIntersection;
}

bool
ElementIsocontour::pointExistInSegm(Point p, std::vector<Point> &segm, Real h_ref)
{
  for (unsigned int i = 0; i < segm.size(); ++i)
    if (twoPointsOverlap(p, segm[i], h_ref))
      return true;
  return false;
}
