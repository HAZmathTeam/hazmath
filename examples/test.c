//
//  test.c
//  
//
//  Created by Hu, Xiaozhe on 1/10/15.
//
//

#include "hazmat.h"

INT main (INT argc, const char * argv[])
{

    dCSRmat A = dcsr_create(10,10,20);
    
    dCSRmat AT;
    
    dcsr_trans(&A, &AT);
    
    printf("hello!\n");
    
    return 0;
    
}
