
typedef struct  
{
  REAL val;
  INT  id; //for the permutation
} weights;

typedef struct  
{
  INT mask;
  INT val;
  INT  id; //for the permutation
} iweights;

typedef INT (*testit)(const void*, const void*);

