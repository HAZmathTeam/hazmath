/*ltz
 Assembles -div (a(x) grad(u) + (b . u)); 
*/
{
  void SchurProduct(const Mesh &mesh, const int dim, SparseMatrix &Ain, Vector &dmass) {
    int i,j,jk,nrow=Ain.Size();
    int *I=Ain.GetI(),*J=Ain.GetJ();
    double *A=Ain.GetData();
    int nv=mesh.GetNV();
    Vector diag0(nrow);
    diag0=0.; //overloaded "=" for vector class. 
    Vector advcoeff(dim), xy(nv*dim),xmid(dim),te(dim);
    double xi,yi,zi,xj,yj,zj,bte,bern,alpe;
    int nv2=nv+nv; //;
    mesh.GetVertices(xy);
    for (i = 0; i < nrow; i++) {
      xi=xy[i];
      yi=xy[nv+i];
      zi=xy[nv2+i];
      for (jk=I[i]; jk<I[i+1]; jk++){
	j=J[jk]; 
	if(i != j){
	  xj=xy[j];
	  yj=xy[nv+j];
	  zj=xy[nv2+j];
	  // compute the advection field at the middle of the edge and
	  // then the bernoulli function
	  te[0]  = xi - xj;
	  te[1]  =  yi - yj;
	  te[2]  = zi - zj;
	  xmid[0] = (xi + xj)*0.5e+0;
	  xmid[1] = (yi + yj)*0.5e+0;
	  xmid[2] = (zi + zj)*0.5e+0;
	  advectEAFE(xmid,advcoeff);
	  //	  bte = (beta*t_e); //c++ overloaded "*" s..t we can multiply these...
	  //	  I do it the old way, unrolled it is faster
	  bte = advcoeff[0]*te[0]+ advcoeff[1]*te[1]+ advcoeff[2]*te[2];
	  ///
	  // diffusion coefficient for the flux J = a(x)\nabla u + \beta u;
	  alpe=diffuse(xmid);
	  //	  xmid.Print(std::cout,3);
	  //	  std::cout << "alpe = "<<alpe<<"; bte="<< bte<<std::endl<<std::flush;      
	  // alpe=a(xmid)\approx harmonic_average=|e|/(int_e 1/a);
	  // should be computed by quadrature in general for 1/a(x).
	  // a_{ij}=B(beta.t_e/harmonic_average)*harmonic_average*omega_e;
	  //	  for (i,j):    B(beta.t_e/harmonic_average);   
	  //	  for (j,i):    B(-beta.t_e/harmonic_average), the minus comes from
	  //     is because t_e=i-j;
	  A[jk] *= alpe*(bernoulli(bte/alpe)); 
    // the diagonal is equal to the negative column sum;
	  diag0[j]-=A[jk];
	}
      } 
    }
    // Another loop to set up the diagonal equal to the negative column sum;
    for (i = 0; i < nrow; i++) 
      for (jk=I[i]; jk<I[i+1]; jk++){
     	j = J[jk]; 
     	if(i == j) A[jk]=diag0[i]+dmass[i];
      }
  }
};
///
void LumpMassBndry(const int dim, const int bndry_robin,
		   DVECTOR dmass, dvector rhs) 
{
  //calculates boundary integral using lumped mass. 
  int i, j,k,nve,attri;
  Element *el=NULL;
  int *node=NULL;
  double v[3],w[3],sixth=1./6.,voltri=1e10,distz=0e0;
  Vector bnormal[dim],xyz[dim];
  dmass=0.;
  fprintf(stdout,"NumOfBdrElements =%i\n",NumOfBdrElements);
  for (i = 0; i < NumOfBdrElements; i++)   {
    el=boundary[i];
    nve = el->GetNVertices();
    node = el->GetVertices();
    attri = el->GetAttribute(); //top=6.
    // we only look at the top bounddary
    //    distz=0e0;
    //    for (j=0; j<nve; ++j){
    //      distz+=fabs(vertices[node[j]](2)-tN);
    //    }
    //    if(distz>1e-12)continue;
       if(attri != bndry_robin)continue;
    for(j=0;j<dim;j++) {
      v[j]=vertices[node[1]](j)-vertices[node[0]](j);
      w[j]=vertices[node[2]](j)-vertices[node[0]](j);
    }
    voltri = (v[1] * w[2]-v[2]*w[1])*(v[1] * w[2]-v[2]*w[1]);
    voltri += (v[2]*w[0]-v[0]*w[2])*(v[2]*w[0]-v[0]*w[2]);
    voltri+=  (v[0] * w[1]-v[1] * w[0])* (v[0] * w[1]-v[1] * w[0]);
    voltri=sixth*pow(voltri,0.5);
     // std::cout << " ==================area= "<< voltri*3.<<std::endl;
     // for (j=0; j<nve; ++j){
     //   std::cout <<" vert"<<j+1<<"="<<node[j]+1<<"(";
     //   for (k=0; k<dim; ++k){
     // 	std::cout <<  vertices[node[j]](k);
     // 	if(k<dim-1) std::cout<<",";  else std::cout<<")"<<std::endl;
     //   }
     // }
    // this should only be added for vertices interior to the boundary at t=T or z=T;
    for (j=0;j<nve;j++){
      xyz[0]=vertices[node[j]](0);
      xyz[1]=vertices[node[j]](1);
      xyz[2]=vertices[node[j]](2);
      advectEAFE(xyz,bnormal);
      // bnormal is b*normal=b\cdot (0,0,1); 
      dmass[node[j]]-=voltri * (bnormal[2]*eps0);      
      rhs[node[j]] += voltri * standard_n_rhs(xyz); 
    }
  }
}
///
int main(int argc, char *argv[])
{
  int i,j;//ltz
   // 1. Parse command-line options.
   EAFE doEAFE;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh;
   ifstream imesh(mesh_file);
   if (!imesh)
   {
      cerr << "\nCan not open mesh file: " << mesh_file << '\n' << endl;
      return 2;
   }
   mesh = new Mesh(imesh, 1, 1, 0);
   imesh.close();
   int dim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   {
     // ref_levels is in constnts.hpp
      for (int l = 0; l < ref_levels; l++)
         mesh->UniformRefinement();
   }
   // 4. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   if (order > 0)
      fec = new H1_FECollection(order, dim);
   else if (mesh->GetNodes())
      fec = mesh->GetNodes()->OwnFEC();
   else
      fec = new H1_FECollection(order = 1, dim);
   // scalar.
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
   // For the vector field.
   std::cout << "Number of unknowns: " << fespace->GetVSize() << std::endl;

   // 5. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm *b = new LinearForm(fespace);
   //   ConstantCoefficient one(1.0);
   FunctionCoefficient fabc(f_rhs);
   b->AddDomainIntegrator(new DomainLFIntegrator(fabc));
   b->Assemble();
   // 6. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;
   FunctionCoefficient uabc(u_exact);

   // 7. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator and imposing homogeneous Dirichlet boundary
   //    conditions. The boundary conditions are implemented by marking all the
   //    boundary attributes from the mesh as essential (Dirichlet). After
   //    assembly and finalizing we extract the corresponding sparse matrix A.
   BilinearForm *a = new BilinearForm(fespace);
   MatrixFunctionCoefficient KDiff(dim,DiffusionTensorEAFE);
   a->AddDomainIntegrator(new DiffusionIntegrator(KDiff));
   //a->AddDomainIntegrator(new DiffusionIntegrator(one));
   a->Assemble(0);
   Array<int> ess_bdr(mesh->bdr_attributes.Max());
   ess_bdr = 1;
   ess_bdr[6-1]=0; //make anything with attribute 6 to be like
		   //interior pt.
   a->Finalize(1);
   //   No constants; original: const SparseMatrix &A = a->SpMat();
   SparseMatrix &A = a->SpMat(); 
   Vector dmass(mesh->GetNV());
   mesh->LumpMassBndry(dim, 6, dmass, *b);
   doEAFE.SchurProduct(*mesh, dim,  A, dmass);
   x.ProjectBdrCoefficient(uabc,ess_bdr);
   a->EliminateEssentialBC(ess_bdr, x, *b); // *b is then a vector
  ////////////////////////////////////////////////////
  // if size is small, then open a stream and write the matrix, the
  // rhs and the umfpack (exact) solution to the linear system.
  /////////////////////
   //  int myn=A.Size();
  // if(myn < 5000) {
  //   std::fstream fp;
  //   fp.open ("a_mfem.dat",std::fstream::out);
  //   (a->SpMat()).zPrintCSR(fp,1);
  //   fp.close();
  //   fp.open ("bx_mfem.dat",std::fstream::out);
  //   b->zPrint(fp);
  //   x.zPrint(fp);
  //   fp.close();
  // }
  ///////////////////////
  int myn=A.Size();
#ifndef MFEM_USE_SUITESPARSE
  int *ia=A.GetI();
  int *ja=A.GetJ(); 
  double *a01=A.GetData();
  double *zrhs=b->GetData();
  // use now multigraph
  std::cout<<" *** Starting multigraph Solver" <<std::endl;
  /* 
  //uncomment if you want to print the matrix as (i,j,v) on a file file:  
     zwrijv(ia, ja, a01, &myn,  &myn,"symm_ref5_sldi.coo");
  */
  mgraph(ia, ja, a01, &myn,zrhs, x, NULL);
#else
   // 8. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   //   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
   umf_solver.SetOperator(A);
   umf_solver.Mult(*b, x);
#endif

   // 9. Save the refined mesh and the solution. This output can be viewed later
   //    using GLVis: "glvis -m refined.mesh -g sol.gf".
   std::cout << std::setw(10) << "EAFE: ref_levels= " << ref_levels	\
	     << "; Ndofs= " 						\
	     << fespace->GetVSize()				\
	     << "; E(L2)= "<< std::setprecision(3)		\
	     << std::scientific					\
	     << x.ComputeL2Error(uabc) << std::endl;
   if(myn<50000) {
    ofstream mesh_ofs("refined.mesh");
    mesh_ofs.precision(8);
    mesh->Print(mesh_ofs);
    ofstream sol_ofs("sol.gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);
    // 10. Send the solution by socket to a GLVis server.
    if (visualization)
      {
	char vishost[] = "localhost";
	int  visport   = 19916;
	socketstream sol_sock(vishost, visport);
	sol_sock.precision(8);
	sol_sock << "solution\n" << *mesh << x << flush;
      }
  }
   // 11. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   if (order > 0)
      delete fec;
   delete mesh;
   return 0;
}
///////////////////////////////
