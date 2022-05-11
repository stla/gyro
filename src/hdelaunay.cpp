#ifndef _GYROHEADER_
#include "gyro.h"
#endif

template <typename HDtT, typename HPointT>
Rcpp::List hdelaunay(const Rcpp::NumericMatrix points) {
  std::vector<HPointT> hpts;
  const unsigned npoints = points.ncol();
  hpts.reserve(npoints);
  for(unsigned i = 0; i != npoints; i++) {
    const Rcpp::NumericVector pt = points(Rcpp::_, i);
    hpts.emplace_back(HPointT(pt(0), pt(1)));
  }
  HDtT hdt;
  hdt.insert(hpts.begin(), hpts.end());
  Rcpp::NumericMatrix Vertices(2, hdt.number_of_vertices());
  {
    int index = 0;
    for(typename HDtT::All_vertices_iterator vd = hdt.all_vertices_begin();
        vd != hdt.all_vertices_end(); ++vd) {
      HPointT pt = vd->point();
      Vertices(0, index) = CGAL::to_double(pt.x());
      Vertices(1, index) = CGAL::to_double(pt.y());
      index++;
      vd->id() = index;
    }
  }
  const size_t nedges = hdt.number_of_hyperbolic_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  {
    size_t i = 0;
    for(typename HDtT::All_edges_iterator ed = hdt.all_edges_begin();
        ed != hdt.all_edges_end(); ++ed) {
      Rcpp::IntegerVector edge_i(2);
      typename HDtT::Vertex_handle sVertex =
          ed->first->vertex(HDtT::cw(ed->second));
      edge_i(0) = sVertex->id();
      typename HDtT::Vertex_handle tVertex =
          ed->first->vertex(HDtT::ccw(ed->second));
      edge_i(1) = tVertex->id();
      Edges(Rcpp::_, i) = edge_i;
      i++;
    }
  }
  const size_t nfaces = hdt.number_of_hyperbolic_faces();
  Rcpp::IntegerMatrix Faces(3, nfaces);
  {
    size_t i = 0;
    for(typename HDtT::All_faces_iterator fd = hdt.all_faces_begin();
        fd != hdt.all_faces_end(); ++fd) {
      Rcpp::IntegerVector face_i(3);
      face_i(0) = fd->vertex(0)->id();
      face_i(1) = fd->vertex(1)->id();
      face_i(2) = fd->vertex(2)->id();
      Faces(Rcpp::_, i) = face_i;
      i++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                            Rcpp::Named("edges") = Edges,
                            Rcpp::Named("faces") = Faces);
}

// [[Rcpp::export]]
Rcpp::List hdelaunay_K(const Rcpp::NumericMatrix points) {
  return hdelaunay<HDt, HPoint>(points);
}

// [[Rcpp::export]]
Rcpp::List hdelaunay_EK(const Rcpp::NumericMatrix points) {
  return hdelaunay<EHDt, EHPoint>(points);
}
