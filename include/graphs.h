/**
 * \struct stack
 * \brief  to be added -- Ludmil
 */
typedef struct stack {
  int top;
  int num_items;
  int null_item;
  int *items;
} stack;
//////////////////////////////////////////////////////////////
/**
 * \struct weights
 * \brief  to be added -- Xiaozhe
 */
typedef struct  
{
  REAL val;
  INT  id; //for the permutation
} weights;

/**
 * \struct iweights
 * \brief  to be added -- Xiaozhe
 */
typedef struct  
{
  INT mask;
  INT val;
  INT  id; //for the permutation
} iweights;

typedef INT (*testit)(const void*, const void*);

