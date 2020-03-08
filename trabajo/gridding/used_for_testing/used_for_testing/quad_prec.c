#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
int main()
{
  long double pi=4e00*(atanl(1e00));
  long double e=expl(1e00);
  fprintf(stdout,"\n\npi=%32.24Le;\n e=%32.24Le\n\n",pi,e); 
  return 0;
}
