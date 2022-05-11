#ifndef _GYROHEADER_
#define _GYROHEADER_
#endif

// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]

#define CGAL_EIGEN3_ENABLED 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/number_utils.h>

#include <Rcpp.h>
#include <RcppEigen.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<K> HDtt;
typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<EK> EHDtt;
typedef HDtt::Point_2 HPoint;
typedef EHDtt::Point_2 EHPoint;
typedef CGAL::Triangulation_data_structure_2<
  CGAL::Triangulation_vertex_base_with_id_2<HDtt>,
  CGAL::Hyperbolic_triangulation_face_base_2<HDtt>>
    HTds;
typedef CGAL::Triangulation_data_structure_2<
  CGAL::Triangulation_vertex_base_with_id_2<EHDtt>,
  CGAL::Hyperbolic_triangulation_face_base_2<EHDtt>>
    EHTds;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<HDtt, HTds> HDt;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<EHDtt, EHTds> EHDt;
