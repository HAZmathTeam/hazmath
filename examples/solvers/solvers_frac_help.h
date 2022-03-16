
#ifndef READ_TO_EOF
#define READ_TO_EOF 0
#endif

#ifndef MAX_NUMBER_POLES
#define MAX_NUMBER_POLES 26
#endif
#ifndef TOL_POWER
#define TOL_POWER -46e0
#endif
#ifndef PRINT_LEVEL
#define PRINT_LEVEL 0
#endif

#ifndef FUNC_PARAM_S
#define FUNC_PARAM_S  0.5e0
#endif
#ifndef FUNC_PARAM_T
#define FUNC_PARAM_T  -0.65e0
#endif
#ifndef FUNC_PARAM_ALPHA
#define FUNC_PARAM_ALPHA  1e0
#endif
#ifndef FUNC_PARAM_BETA
#define FUNC_PARAM_BETA  10e0
#endif

/***************************************************************************/
/* ---------------------- */
/*  setup test problem: the exact solutions us.dat and u3.dat are for the following values: */
/* ---------------------- */
/* alpha = 1; s = 0.5; beta = 10; t = -0.65; */
/************************************************************************/
/*function to be approximated*/
/************************************************************************/
static REAL16 f_to_approx_local(REAL16 x,			\
				REAL16 s, REAL16 t,		\
				REAL16 alpha, REAL16 beta)
{  
  return 1./(alpha*powl(x,s)+beta*powl(x,t)); 
}
/**/
