/*! \file LCData.h
 *  FOR HAZMATH
 *  Created by James Adler on 9/12/2017
 *  Copyright 2015_HAZMATH__. All rights reserved.
 *
 * 3D harmonic mapping solutions as used in
 * Error Estimators and Marking Strategies for Electrically Coupled Liquid Crystal Systems
 * 
 * The goal is to provide functions here that, for a given point, return the value of the true solution and the
 * gradients of the true solution at a provided point
*/

#include <math.h>

void solution_value(REAL *val, REAL *x) {

  REAL x_ = x[0] + 0.2;
  REAL y_ = x[1] + 0.1;
  REAL z_ = x[2];

  // unit length domain projection
  REAL p1 = x/(sqrt(x_*x_ + y_*y_ + z_*z_));
  REAL p2 = y/(sqrt(x_*x_ + y_*y_ + z_*z_));
  REAL p3 = z/(sqrt(x_*x_ + y_*y_ + z_*z_));  

  // forward stereographic projection
  REAL pi1 = p1/(1-p3);
  REAL pi2 = p2/(1-p3); 

  // Holomorphic Function
  REAL r1 = pi1/(pi1*pi1 + pi2*pi2) + pi1*pi1 - pi2*pi2;
  REAL r2 = 2*pi1*pi2 - pi2/(pi1*pi1 + pi2*pi2);  

  val[0] = (2*r1)/(1 + r1*r1 + r2*r2);
  val[1] = (2*r2)/(1 + r1*r1 + r2*r2);
  val[2] =  (r1*r1 + r2*r2 - 1)/(1 + r1*r1 + r2*r2); 
  return;
}

void n1_gradient(REAL *val, REAL *x) {

  REAL alpha1 = 0.1 + x[1];
  REAL beta1 = 0.2 + x[0];
  REAL z = x[2];

  REAL theta1 = z*z + alpha1*alpha1 + beta1*beta1;
  REAL zeta1 = 1 - z/sqrt(theta1);

  REAL zeta2 = alpha1*alpha1/(theta1*zeta1*zeta1) + beta1*beta1/(theta1*zeta1*zeta1);
  REAL zeta3 = (-2*z*alpha1*alpha1*alpha1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*alpha1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*alpha1/(theta1*theta1*zeta1*zeta1)
		  - 2*alpha1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*alpha1/(theta1*zeta1*zeta1));
  REAL zeta4 = (-2*z*alpha1*alpha1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*beta1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*beta1/(theta1*theta1*zeta1*zeta1)
		  - 2*beta1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*beta1/(theta1*zeta1*zeta1));
  REAL zeta5 = alpha1*alpha1 + beta1*beta1;
  REAL zeta6 = 5*alpha1*alpha1*alpha1*alpha1 - 10*alpha1*alpha1*beta1*beta1 + beta1*beta1*beta1*beta1;
  REAL zeta7 = alpha1-beta1;
  REAL zeta8 = alpha1+beta1;
  REAL zeta9 = -3*alpha1*alpha1 + beta1*beta1;
  REAL zeta10 = beta1*beta1 - alpha1*alpha1 + theta1;
  REAL zeta11 = alpha1*alpha1 - beta1*beta1;
  REAL zeta12 = beta1*beta1 + theta1;
  REAL zeta13 = beta1*beta1 - alpha1*alpha1;

  REAL denom1 = sqrt(theta1)*pow(theta1*zeta1*zeta1*zeta5
				   + zeta2*zeta2*(theta1*theta1*pow(zeta1, 4.0) + zeta5*zeta5)
				   + 2*beta1*sqrt(theta1)*zeta1*zeta2*zeta9, 2.0);

  val[0] = -(2*zeta1*(2*beta1*z*pow(zeta2, 4)*pow(zeta5, 2)*zeta7*zeta8 + 
			  pow(zeta1, 6)*pow(zeta2, 2)*pow(theta1, 3)*(zeta2*(pow(beta1, 2) - theta1) + 
								      beta1*zeta4*theta1) + 
			  beta1*pow(zeta1, 5)*pow(zeta2, 3)*pow(theta1, 5.0/2.0)
			  *(beta1*z - 2*zeta2*(zeta11 + theta1)) - pow(zeta1, 4)*pow(theta1, 2.0)
			  *(2*beta1*z*pow(zeta2, 4)*zeta7*zeta8 + beta1*zeta4*zeta5*theta1 + 
			    zeta2*(pow(beta1, 4) + pow(alpha1, 2)*zeta12 - pow(beta1, 2)*theta1))
			  - beta1*zeta1*pow(zeta2, 3)*sqrt(theta1)
			  *(beta1*(9*pow(alpha1, 4) - 2*pow(alpha1, 2)*pow(beta1, 2)
				   + 5*pow(beta1, 4))*z - 2*zeta2*zeta5*(pow(alpha1, 4) - pow(beta1, 4)
									 - 3*pow(alpha1, 2)*theta1
									 + pow(beta1, 2)*theta1))
			  - pow(zeta1, 3)*zeta2*pow(theta1, 3.0/2.0)
			  *(pow(beta1, 2.0)*z*zeta5+ 2*(-pow(alpha1, 4) + pow(beta1, 4))
			    *zeta4*theta1 + 4*beta1*zeta2*(pow(beta1, 4) - pow(beta1, 2)*theta1
							   + pow(alpha1, 2)*(-3*pow(beta1, 2) + theta1)))
			  + pow(zeta1, 2)*pow(zeta2, 2)*theta1
			  *(-4*pow(beta1, 3)*z*zeta9 - beta1*zeta4*zeta6*theta1
			    + zeta2*(-5*pow(beta1, 6) + 2*pow(alpha1, 2)*pow(beta1, 2)*zeta12
				     + 5*pow(beta1, 4)*theta1 + pow(alpha1, 4)
				     *(-9*pow(beta1, 2) + 5*theta1)))))/denom1;

  val[1] = -(2*zeta1*(2*alpha1*z*pow(zeta2, 4)*pow(zeta5, 2)*zeta7*zeta8 + 
			  alpha1*pow(zeta1, 5)*pow(zeta2, 3)*(beta1*z + 2*zeta10*zeta2)*pow(theta1, 5.0/2.0)
			  + beta1*pow(zeta1, 6)*pow(zeta2, 2)*pow(theta1, 3)*(alpha1*zeta2 + zeta3*theta1)
			  - pow(zeta1, 4)*pow(theta1, 2)*(2*alpha1*z*pow(zeta2, 4)*zeta7*zeta8
							  + alpha1*beta1*zeta2*(zeta5 - 2*theta1)
							  + beta1*zeta3*zeta5*theta1)
			  - pow(zeta1, 3)*zeta2*pow(theta1, 3.0/2.0)
			  *(alpha1*beta1*z*zeta5 + 2*(-pow(alpha1, 4) + pow(beta1, 4))*zeta3*theta1
			    + 4*alpha1*pow(beta1, 2)*zeta2*(zeta9 + 2*theta1))
			  + alpha1*zeta1*pow(zeta2, 3)*sqrt(theta1)
			  *(beta1*(-9*pow(alpha1, 4) + 2*pow(alpha1, 2)*pow(beta1, 2) - 5*pow(beta1, 4))*z
			    + 2*zeta2*zeta5*(pow(alpha1, 4) - pow(beta1, 4) - pow(alpha1, 2)*theta1
					     + 3*pow(beta1, 2)*theta1))
			  - beta1*pow(zeta1, 2)*pow(zeta2, 2)*theta1
			  *(4*alpha1*beta1*z*zeta9 + zeta3*zeta6*theta1
			    + alpha1*zeta2*(9*pow(alpha1, 4) + 5*pow(beta1, 4) + 4*pow(beta1, 2)*theta1
					    - 2*pow(alpha1, 2)*(pow(beta1, 2) + 2*theta1)))))/denom1;

  val[2] = -(2*(z*z + z*zeta1*sqrt(theta1) - theta1)
		 *(2*beta1*pow(zeta2, 2)*zeta5*zeta6 + 2*zeta1*zeta2*(-2 + pow(zeta2, 3))
		   *pow(zeta5, 2)*zeta7*zeta8*sqrt(theta1) + beta1*pow(zeta1, 2)
		   *((-9*pow(alpha1, 4) + 2*pow(alpha1, 2)*pow(beta1,  2)
		      - 5*pow(beta1, 4))*pow(zeta2, 3) + 2*pow(zeta5, 2))*theta1
		   - 4*pow(beta1, 2)*pow(zeta1, 3)*pow(zeta2,  2)*zeta9*pow(theta1, 3.0/2.0)
		   - beta1*pow(zeta1, 4)*zeta2*(1 + 2*zeta2)*zeta5*pow(theta1, 2)
		   + 2*pow(zeta1, 5)*zeta13*pow(zeta2, 4)*pow(theta1, 5.0/2.0)
		   + beta1*pow(zeta1, 6)*pow(zeta2, 3)*pow(theta1, 3)))
    /(theta1*pow((2*beta1*zeta1*zeta2*zeta9*sqrt(theta1) + pow(zeta1, 2)*zeta5*theta1
		  + pow(zeta2, 2)*(pow(zeta5, 2) + pow(zeta1, 4)*pow(theta1, 2))), 2));
  return;
}

void n2_gradient(REAL *val, REAL *x) {
  REAL alpha1 = 0.1 + x[1];
  REAL beta1 = 0.2 + x[0];
  REAL z = x[2];

  REAL theta1 = z*z + alpha1*alpha1 + beta1*beta1;
  REAL zeta1 = 1 - z/sqrt(theta1);

  REAL zeta2 = alpha1*alpha1/(theta1*zeta1*zeta1) + beta1*beta1/(theta1*zeta1*zeta1);
  REAL zeta3 = (-2*z*alpha1*alpha1*alpha1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*alpha1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*alpha1/(theta1*theta1*zeta1*zeta1)
		  - 2*alpha1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*alpha1/(theta1*zeta1*zeta1));
  REAL zeta4 = (-2*z*alpha1*alpha1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*beta1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*beta1/(theta1*theta1*zeta1*zeta1)
		  - 2*beta1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*beta1/(theta1*zeta1*zeta1));
  REAL zeta5 = alpha1*alpha1 + beta1*beta1;
  REAL zeta6 = 5*alpha1*alpha1*alpha1*alpha1 - 10*alpha1*alpha1*beta1*beta1 + beta1*beta1*beta1*beta1;
  REAL zeta7 = alpha1-beta1;
  REAL zeta8 = alpha1+beta1;
  REAL zeta9 = -3*alpha1*alpha1 + beta1*beta1;
  REAL zeta10 = beta1*beta1 - alpha1*alpha1 + theta1;
  REAL zeta11 = alpha1*alpha1 - beta1*beta1;
  REAL zeta12 = beta1*beta1 + theta1;
  REAL zeta13 = beta1*beta1 - alpha1*alpha1;

  REAL denom2 = sqrt(theta1)*pow((2*beta1*zeta1*zeta2*zeta9*sqrt(theta1) + pow(zeta1, 2)*zeta5*theta1
				    + pow(zeta2, 2)*(pow(zeta5, 2) + pow(zeta1, 4)*pow(theta1, 2))), 2);

  val[0] = (2*alpha1*zeta1*(4*pow(beta1, 2)*z*pow(zeta2, 4)*pow(zeta5, 2)
				+ pow(zeta1, 6)*pow(zeta2, 2)*pow(theta1, 3)*(beta1*zeta2 + zeta4*theta1)
				- pow(zeta1, 4)*pow(theta1, 2)*(beta1*zeta2*(4*beta1*z*pow(zeta2, 3) + zeta5
									     - 2*theta1) + zeta4*zeta5*theta1)
				+ pow(zeta1, 5)*pow(zeta2, 3)*pow(theta1, 5.0/2.0)
				*(beta1*z + 2*zeta2*(-2*pow(beta1, 2) + theta1))
				- pow(zeta1, 3)*zeta2*pow(theta1, 3.0/2.0)
				*(beta1*zeta5*(z - 4*zeta4*theta1)
				  +  4*zeta2*(pow(beta1, 4) - pow(beta1, 2)*theta1 + pow(alpha1, 2)
					      *(-3*pow(beta1, 2) + theta1)))
				+ zeta1*pow(zeta2, 3)*sqrt(theta1)
				*(beta1*(-3*pow(alpha1, 4) - 18*pow(alpha1, 2)*pow(beta1, 2) + pow(beta1, 4))*z
				  + 2*zeta2*zeta5*(2*pow(beta1, 4) - 3*pow(beta1, 2)*theta1
						   + pow(alpha1, 2)*(2*pow(beta1, 2) + theta1)))
				+ pow(zeta1, 2)*pow(zeta2, 2)*theta1
				*(-4*pow(beta1, 2)*z*zeta9 + (pow(alpha1, 4) - 10*pow(alpha1, 2)
							      *pow(beta1, 2) + 5*pow(beta1, 4))
				  *zeta4*theta1 + beta1*zeta2*(-3*pow(alpha1, 4) + pow(beta1, 4)
							       - 4*pow(beta1, 2)*theta1 + pow(alpha1, 2)
							       *(-18*pow(beta1, 2) + 4*theta1)))))/denom2;

  val[1] = (2*zeta1*(4*pow(alpha1, 2)*beta1*z*pow(zeta2, 4)*pow(zeta5, 2)
			 + pow(zeta1, 6)*pow(zeta2, 2)*pow(theta1, 3)*(zeta2*(pow(alpha1, 2) - theta1) + alpha1*zeta3*theta1)
			 + pow(zeta1, 5)*pow(zeta2, 3)*pow(theta1, 5.0/2.0)
			 *(pow(alpha1, 2)*z + 2*beta1*zeta2*(-2*pow(alpha1, 2) + theta1))
			 + zeta1*pow(zeta2, 3)*sqrt(theta1)
			 *(pow(alpha1, 2)*(-3*pow(alpha1, 4) - 18*pow(alpha1, 2)*pow(beta1, 2) + pow(beta1, 4))*z
			   + 2*beta1*zeta2*zeta5*(2*pow(alpha1, 4) + pow(alpha1, 2)
						  *(2*pow(beta1, 2) - 3*theta1) + pow(beta1, 2)*theta1))
			 - pow(zeta1, 4)*pow(theta1, 2)*(4*pow(alpha1, 2)*beta1*z*pow(zeta2, 4)
							 + alpha1*zeta3*zeta5*theta1
							 + zeta2*(pow(alpha1, 4) + pow(alpha1, 2)*(pow(beta1, 2) - theta1)
								  + pow(beta1, 2)*theta1))
			 - alpha1*pow(zeta1, 3)*zeta2*pow(theta1, 3.0/2.0)
			 *(4*alpha1*beta1*zeta2*(zeta9 + 2*theta1) + zeta5*(alpha1*z - 4*beta1*zeta3*theta1))
			 + pow(zeta1, 2)*pow(zeta2, 2)*theta1*(alpha1*(-4*alpha1*beta1*z*zeta9
								       + (pow(alpha1, 4) - 10*pow(alpha1, 2)*pow(beta1, 2)
									  + 5*pow(beta1, 4))*zeta3*theta1)
							       + zeta2*(-3*pow(alpha1, 6) + 3*pow(beta1, 4)*theta1
									+ 3*pow(alpha1, 4)*(-6*pow(beta1, 2) + theta1)
									+ pow(alpha1, 2)*(pow(beta1, 4)
											  + 14*pow(beta1, 2)*theta1)))))/denom2;

  val[2] = -((2*alpha1*(pow(z, 2) + z*zeta1*sqrt(theta1) - theta1)
		  *(2*(pow(alpha1, 4) - 10*pow(alpha1, 2)*pow(beta1, 2) + 5*pow(beta1, 4))*pow(zeta2, 2)*zeta5
		    - 4*beta1*zeta1*zeta2*(-2 + pow(zeta2, 3))*pow(zeta5, 2)*sqrt(theta1)
		    + pow(zeta1, 2)*((3*pow(alpha1, 4) + 18*pow(alpha1, 2)*pow(beta1, 2)
				      - pow(beta1, 4))*pow(zeta2, 3) - 2*pow(zeta5, 2))*theta1
		    + 4*beta1*pow(zeta1, 3)*pow(zeta2, 2)*zeta9*pow(theta1, 3.0/2.0)
		    + pow(zeta1, 4)*zeta2*(1 + 2*zeta2)*zeta5*pow(theta1, 2)
		    + 4*beta1*pow(zeta1, 5)*pow(zeta2, 4)*pow(theta1, 5.0/2.0)
		    - pow(zeta1, 6)*pow(zeta2, 3)*pow(theta1, 3)))
		 /(theta1*pow((2*beta1*zeta1*zeta2*zeta9*sqrt(theta1) + pow(zeta1, 2)*zeta5*theta1
			       + pow(zeta2, 2)*(pow(zeta5, 2) + pow(zeta1, 4)*pow(theta1, 2))), 2)));
  return;
}

void n3_gradient(REAL *val, REAL *x) {
  REAL alpha1 = 0.1 + x[1];
  REAL beta1 = 0.2 + x[0];
  REAL z = x[2];

  REAL theta1 = z*z + alpha1*alpha1 + beta1*beta1;
  REAL zeta1 = 1 - z/sqrt(theta1);

  REAL zeta2 = alpha1*alpha1/(theta1*zeta1*zeta1) + beta1*beta1/(theta1*zeta1*zeta1);
  REAL zeta3 = (-2*z*alpha1*alpha1*alpha1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*alpha1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*alpha1/(theta1*theta1*zeta1*zeta1)
		  - 2*alpha1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*alpha1/(theta1*zeta1*zeta1));
  REAL zeta4 = (-2*z*alpha1*alpha1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*z*beta1*beta1*beta1/(pow(theta1, 5.0/2.0)*zeta1*zeta1*zeta1)
		  - 2*alpha1*alpha1*beta1/(theta1*theta1*zeta1*zeta1)
		  - 2*beta1*beta1*beta1/(theta1*theta1*zeta1*zeta1)
		  + 2*beta1/(theta1*zeta1*zeta1));
  REAL zeta5 = alpha1*alpha1 + beta1*beta1;
  REAL zeta6 = 5*alpha1*alpha1*alpha1*alpha1 - 10*alpha1*alpha1*beta1*beta1 + beta1*beta1*beta1*beta1;
  REAL zeta7 = alpha1-beta1;
  REAL zeta8 = alpha1+beta1;
  REAL zeta9 = -3*alpha1*alpha1 + beta1*beta1;
  REAL zeta10 = beta1*beta1 - alpha1*alpha1 + theta1;
  REAL zeta11 = alpha1*alpha1 - beta1*beta1;
  REAL zeta12 = beta1*beta1 + theta1;
  REAL zeta13 = beta1*beta1 - alpha1*alpha1;

  REAL denom3 = pow((2*beta1*zeta1*zeta2*zeta9*sqrt(theta1) + pow(zeta1, 2)*zeta5*theta1
		       + pow(zeta2, 2)*(pow(zeta5, 2) + pow(zeta1, 4)*pow(theta1, 2))), 2);

  val[0] = -((4*pow(zeta1, 3)*zeta2*sqrt(theta1)
		  *(2*beta1*z*pow(zeta2, 3)*pow(zeta5, 2) + beta1*zeta1*pow(zeta2, 2)
		    *(3*beta1*z*zeta9 + 2*zeta2*zeta5*(zeta5 - theta1))*sqrt(theta1)
		    + pow(zeta1, 3)*pow(theta1, 3.0/2.0)*(beta1*zeta2*(zeta5 - theta1) + zeta4*zeta5*theta1)
		    + pow(zeta1, 2)*zeta2*theta1*(beta1*z*zeta5 + beta1*zeta4*zeta9*theta1
						  + 3*zeta2*(pow(beta1, 4) - pow(beta1, 2)*theta1
							     + pow(alpha1, 2)*(-3*pow(beta1, 2) + theta1)))))/denom3);

  val[1] = -((4*pow(zeta1, 3)*zeta2*sqrt(theta1)
		  *(2*alpha1*z*pow(zeta2, 3)*pow(zeta5, 2) + zeta1*pow(zeta2, 2)
		    *(3*alpha1*beta1*z*zeta9 + 2*alpha1*zeta2*zeta5*(zeta5 - theta1))*sqrt(theta1)
		    + pow(zeta1, 3)*pow(theta1, 3.0/2.0)*(alpha1*zeta2*(zeta5 - theta1) + zeta3*zeta5*theta1)
		    + pow(zeta1, 2)*zeta2*theta1*(alpha1*z*zeta5 + beta1*zeta3*zeta9*theta1
						  + 3*alpha1*beta1*zeta2*(zeta9 + 2*theta1))))/denom3);

  val[2] = (4*pow(zeta1, 2)*zeta2*(z*z + z*zeta1*sqrt(theta1) - theta1)
		*(2*beta1*zeta2*zeta5*zeta9 - 2*zeta1*(-1 + pow(zeta2, 3))*pow(zeta5, 2)*sqrt(theta1)
		  - 3*beta1*pow(zeta1, 2)*pow(zeta2, 2)*zeta9*theta1
		  - pow(zeta1, 3)*zeta2*zeta5*pow(theta1, 3.0/2.0)))/denom3;
  return;
}
