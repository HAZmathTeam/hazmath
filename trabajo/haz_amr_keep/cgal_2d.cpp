/*! \file src/amr/haz_cgal.c
 *
 *  Created by James Adler, Xiaozhe Hu, and Ludmil Zikatanov on 20210313.
 *  Copyright 2019__HAZMATH__. All rights reserved.
 *
 *  \note routines doing meshes using cgal library
 *
 *  \note: by ltz on 20210313
 */
#if WITH_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <iostream>
#include <fstream>
//#include <CGAL/IO/write_vtu.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;
//
typedef CDT::Vertex_handle Vertex_handle;
typedef typename CDT::Face_handle Face_handle;
typedef typename CDT::Finite_vertices_iterator  Finite_vertices_iterator;
typedef typename CDT::Finite_faces_iterator  Finite_faces_iterator;
typedef typename CDT::Finite_edges_iterator  Finite_edges_iterator;
typedef CDT::Point Point;
typedef CDT::Edge  Edge;
typedef CDT::Face  Face;
#endif
//
extern "C" {
#include "hazmath.h"
}
static void do_tri_2d(REAL *x_in, INT *edges,			\
		      REAL *seedpt,				\
		      INT np, INT ne, INT dim,			\
		      REAL shape_bound,				\
		      REAL hmax,				\
		      INT **nodes, REAL **xv,			\
		      INT *nvo, INT *nso)
{
  CDT cdt;
  Vertex_handle va,vb,v0;
  INT ii,jj,i,j,dim1=dim+1;
  std::vector<Vertex_handle> vh;
  for(j=0;j<np;++j){
    jj=dim*j;
    va=cdt.insert(Point(x_in[jj],x_in[jj+1]));
    vh.push_back(va);
  }
  std::size_t nv0=i;
  //  std::cout <<std::endl <<std::endl <<std::endl;
  for (i=0;i<ne;++i){
    ii=edges[2*i];
    jj=edges[2*i+1];
    cdt.insert_constraint(vh[ii], vh[jj]);
    //    std::cout <<"!!!e(" << i << ")=(" << ii << "," << jj << ")" <<std::endl;
  }
  vh.clear();
  //
  std::list<Point> list_of_seeds;
  list_of_seeds.push_back(Point(seedpt[0],seedpt[1]));
  //  std::cout << "Meshing the domain..." << std::endl;
  // 0.125 is the default shape bound. It corresponds to abound 20.6 degree.
  // 0.5 is the upper bound on the length of the longuest edge.
  // See reference manual for Delaunay_mesh_size_traits_2<K>.
  CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),Criteria(shape_bound,hmax));
  //simplices:
  std::size_t ns;
  std::map<Face_handle, std::size_t > tall;
  Finite_faces_iterator fit;
  ns=0;
  for(fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end(); ++fit){
    if(!fit->is_in_domain()) continue;
    tall[fit] = ns++;
  }
  //vertices
  std::size_t nv;
  std::map<Vertex_handle, std::size_t> vall;
  Finite_vertices_iterator vit;
  nv=0;
  for (vit=cdt.finite_vertices_begin();vit != cdt.finite_vertices_end(); ++vit)
    vall[vit] = nv++;
  *nvo=nv;
  *nso=ns;
  // num simplices:
#ifdef __cpluspplus
  extern "C" {
#endif
    xv[0]=(REAL *)calloc(nv*dim,sizeof(REAL));
    nodes[0]=(INT *)calloc(ns*dim1,sizeof(INT));
#ifdef __cpluspplus
  }
#endif
  //vertices
  i=0;
  for (vit=cdt.finite_vertices_begin();vit != cdt.finite_vertices_end(); ++vit){
    for(j = 0; j < dim; j++)
      xv[0][i*dim+j]=vit->point()[j];
    i++;
  }
  //simplices
  i=0;
  for(fit = cdt.finite_faces_begin();fit != cdt.finite_faces_end(); ++fit) {
    if(!fit->is_in_domain()) continue; // grab only faces interior for the domain.
    for(j = 0; j < dim1; j++)
      nodes[0][i*dim1+j]=vall[fit->vertex(j)];
    i++;
  }
  std::cout << "%%Number of vertices: " << nv << std::endl;
  std::cout << "%%Number of triangles: " << ns << std::endl;
  int count = 0;
  // writes vtu in wrong way:infinite triangles and also places have
  //  two double quotes '""' instead of 1 double quote'"'
  //  std::ofstream filevtu("2dout.vtu");
  //  CGAL::write_vtu(filevtu,cdt,CGAL::IO::ASCII); filevtu.close();
  for(CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
      eit != cdt.finite_edges_end(); ++eit){
    //  for (const Edge& e : cdt.finite_edges())
        if (cdt.is_constrained(*eit))    ++count;
	//	std::cout << "Info:" << e.first->info() << "\n";
	//	std::cout << "Info:" << e.first->neighbor(e.second) << "\n";
  }
  std::cout << "%%The number of resulting constrained edges is: ";
  std::cout <<  count << std::endl;
}
extern "C" {
void c_do_tri_2d(REAL *x_in, INT *edges,			\
		 REAL *seedpt,					\
		 INT np, INT ne, INT dim,			\
		 REAL shape_bound,				\
		 REAL hmax,					\
		 INT **nodes, REAL **xv,			\
		 INT *nvo, INT *nso)
{
#if WITH_CGAL
       do_tri_2d(x_in, edges,					\
		 seedpt,					\
		 np, ne, dim,					\
		 shape_bound,					\
		 hmax,						\
		 nodes, xv,					\
		 nvo,nso);
#else
       error_extlib(250, __FUNCTION__, "CGAL");
       return;
#endif
     }
}
