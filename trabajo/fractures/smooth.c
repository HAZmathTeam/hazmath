#include "hazmath.h"
/****************************************/
/**
 * \fn void smoother_dcsr_Schwarz_forward0 (Schwarz_data  *Schwarz,
 *                                         Schwarz_param *param,
 *                                         dvector *x, dvector *b)
 *
 * \brief Schwarz smoother: forward sweep
 *
 * \param Schwarz Pointer to the Schwarz data
 * \param param   Pointer to the Schwarz parameter
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \note Needs improvment -- Xiaozhe
 * \note Improved (Ludmil)
 */
void smoother_dcsr_Schwarz_forward0 (Schwarz_data  *Schwarz,
				     Schwarz_param *param,
				     dvector       *x,
				     dvector       *b)
{
  INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

  // Schwarz partition
  INT  nblk = Schwarz->nblk;
  dCSRmat *blk = Schwarz->blk_data;
  INT  *iblock = Schwarz->iblock;
  INT  *jblock = Schwarz->jblock;
  INT  *mask   = Schwarz->mask;
  INT  block_solver = param->Schwarz_blksolver;


  // Schwarz data
  dCSRmat A = Schwarz->A;
  INT *ia = A.IA;
  INT *ja = A.JA;
  REAL *val = A.val;

  // Local solution and right hand vectors
  dvector rhs = Schwarz->rhsloc1;
  dvector u   = Schwarz->xloc1;

#if WITH_SUITESPARSE
  void **numeric = Schwarz->numeric;
#endif

  for (is=0; is<nblk; ++is) {
    // Form the right hand of eack block
    ibl0 = iblock[is];
    ibl1 = iblock[is+1];
    nloc = ibl1-ibl0;
    for (i=0; i<nloc; ++i ) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = i+1;
    }

    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki = jblock[iblk];
      rhs.val[i] = b->val[ki];
      iaa = ia[ki];
      iab = ia[ki+1];
      for (kij = iaa; kij<iab; ++kij) {
	kj = ja[kij];
	j  = mask[kj];
	if(j == 0) {
	  rhs.val[i] -= val[kij]*x->val[kj];
	}
      }
    }
    // Solve each block
    switch (block_solver) {
#if WITH_SUITESPARSE
    case SOLVER_UMFPACK: {
      /* use UMFPACK direct solver on each block */
      umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
      break;
    }
#endif
    default:
      /* use iterative solver on each block */
      u.row = blk[is].row;
      rhs.row = blk[is].row;
      dvec_set(u.row, &u, 0);
      dcsr_pvgmres(&blk[is], &rhs, &u, NULL, 1e-8, 20, 20, 1, 0);
    }

    //zero the mask so that everyting is as it was
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = 0;
      x->val[ki] = u.val[i];
    }
  }
}

/**
 * \fn void smoother_dcsr_Schwarz_backward0 (Schwarz_data  *Schwarz,
 *                                          Schwarz_param *param,
 *                                          dvector *x, dvector *b)
 *
 * \brief Schwarz smoother: backward sweep
 *
 * \param Schwarz Pointer to the Schwarz data
 * \param param   Pointer to the Schwarz parameter
 * \param x       Pointer to solution vector
 * \param b       Pointer to right hand
 *
 * \note Needs improvment -- Xiaozhe
 */
void smoother_dcsr_Schwarz_backward0 (Schwarz_data *Schwarz,
                                     Schwarz_param *param,
                                     dvector *x,
                                     dvector *b)
{
  INT i, j, iblk, ki, kj, kij, is, ibl0, ibl1, nloc, iaa, iab;

  // Schwarz partition
  INT  nblk = Schwarz->nblk;
  dCSRmat *blk = Schwarz->blk_data;
  INT  *iblock = Schwarz->iblock;
  INT  *jblock = Schwarz->jblock;
  INT  *mask   = Schwarz->mask;
  INT  block_solver = param->Schwarz_blksolver;


  // Schwarz data
  dCSRmat A = Schwarz->A;
  INT *ia = A.IA;
  INT *ja = A.JA;
  REAL *val = A.val;

  // Local solution and right hand vectors
  dvector rhs = Schwarz->rhsloc1;
  dvector u   = Schwarz->xloc1;

#if WITH_SUITESPARSE
  void **numeric = Schwarz->numeric;
#endif

  for (is=nblk-1; is>=0; --is) {
    // Form the right hand of eack block
    ibl0 = iblock[is];
    ibl1 = iblock[is+1];
    nloc = ibl1-ibl0;
    for (i=0; i<nloc; ++i ) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = i+1;
    }

    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki = jblock[iblk];
      rhs.val[i] = b->val[ki];
      iaa = ia[ki]-1;
      iab = ia[ki+1]-1;
      for (kij = iaa; kij<iab; ++kij) {
	kj = ja[kij]-1;
	j  = mask[kj];
	if(j == 0) {
	  rhs.val[i] -= val[kij]*x->val[kj];
	}
      }
    }

    // Solve each block
    switch (block_solver) {

#if WITH_SUITESPARSE
    case SOLVER_UMFPACK: {
      /* use UMFPACK direct solver on each block */
      umfpack_solve(&blk[is], &rhs, &u, numeric[is], 0);
      break;
    }
#endif
    default:
      /* use iterative solver on each block */
      rhs.row = blk[is].row;
      u.row   = blk[is].row;
      dvec_set(u.row, &u, 0);
      dcsr_pvgmres (&blk[is], &rhs, &u, NULL, 1e-8, 20, 20, 1, 0);
    }

    //zero the mask so that everyting is as it was
    for (i=0; i<nloc; ++i) {
      iblk = ibl0 + i;
      ki   = jblock[iblk];
      mask[ki] = 0;
      x->val[ki] = u.val[i];
    }
  }
}
/**END**/
/****************************************/
