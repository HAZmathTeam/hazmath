/************************************************************************/
/************************************************************************/
/*function to be approximated*/
REAL16 f_to_approx_l(REAL16 x,REAL16 s, REAL16 t, REAL16 alpha, REAL16 beta)
{  
  //  fprintf(stdout,"\ns1=%Lf; s2=%Lf; alpha=%Lf; beta=%Lf;",s,t,alpha,beta);
  //return 1./(alpha*powl(x,s)+beta*powl(x,t)); 
   return (alpha*powl(x,s)+beta*powl(x,t)); 
}
/**/
/*OLD:function to be approximated*/
REAL16 old_f_to_approx_l(REAL16 x,void *param)
{
  REAL16 *s,s1,s2,alpha,beta;
  if(param!=NULL){
    s=(REAL16 *)param;
    s1=s[0];
    s2=s[1];
    alpha=s[2];
    beta=s[3];
  } else {
    s1=0.5e0;
    s2=-0.5e0;
    alpha=1e0;
    beta=2e0;
  }
  //  fprintf(stdout,"\ns1=%Lf; s2=%Lf; alpha=%Lf; beta=%Lf;",s1,s2,alpha,beta);
  return f_to_approx_l(x,s1,s2,alpha,beta);
}
/**/
