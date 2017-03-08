
/* MULTIGRAPH PARAMETERS for the call to mginit_()*/
INT ispd=0; // is this an SPD problem ?
INT ncfact=4; // ?factor for coarsening
INT maxlvl=20; // max number of levels
INT maxfil=32; // max fill in
INT method=0; // ilu method (check the multigraph-2.0 user guide for
	      // details)
REAL dtol=1e-5; // drop tollerance
INT iflag=-16;  // error flag, set to an unusual value
INT lvl=-16; // levels counter, set to an unusual value

/* Multigrpah parameters for the call to mg_() */
REAL eps1=1e-8; // tolerance 
INT mxcg=10000; // max cg iterations
REAL relerr=1e0; // check the user manual

/* from multigraph-2.0 source
  subroutine mginit(n,ispd,nblock,ib,maxja,ja,maxa,a,ncfact,
	  +      maxlvl,maxfil,ka,lvl,dtol,method,iflag)
*/
void  mginit_(INT *n, INT *ispd, INT *nblock, INT *iblock,		\
	      INT *maxja, INT *jareb, INT *maxa, REAL *areb,		\
	      INT *ncfact, INT *maxlvl,	INT *maxfil, INT *ka,
	      INT *lvl, REAL *dtol, INT *method, 
	      INT *iflag );

  /*   from multigraph-2.0 source
       subroutine mg(ispd,lvl,mxcg,eps1,ja,a,dr,br,ka,relerr,
       +      iflag,hist)
  */
void mg_(INT *ispd, INT *lvl, INT *mxcg, REAL *eps1,			\
	 INT *jareb, REAL *areb,		\
	 REAL *sol, REAL *rhs, INT *ka, REAL *relerr,		\
	 INT *iflag, REAL *hist );

