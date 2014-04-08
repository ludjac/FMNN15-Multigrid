#ifndef mg_2d
#define mg_2d

/* TYPEDEFS STRUCTS AND MACROS */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* In order to be able to change precision for the entire 
 * computation we used a typedef data_t.
 * 
 * pow_t is used instead of int to make sure there wont be a datatype 
 * overflow in any of the loops.   
 */

typedef double data_t;
typedef unsigned long pow_t; 

/* Struct containing the entire problem. */
typedef struct{
	pow_t k; 		//(2^k-1)^2 grid points for the finest grid
	pow_t c; 		// current grid
	pow_t n; 		//current grid as 2^c-1
	pow_t r; 		//roughest grid
	pow_t sum; 		// Size of allocated memory, used in debugging
	data_t* f; 		// Left-hand side 
	data_t* bc0; 	// beginning conditions, 
	data_t* bc1; 	// beginning conditions, 
	data_t* bc2; 	// beginning conditions, 
	data_t* bc3; 	// beginning conditions, 
	data_t* g;   	// the u grid
	data_t* rf;  	// the Residual grid
#if 0
       ----bc0----
       |	 	 |
       b	 	 b
       c	 	 c
       1	 	 2
       |	 	 |
       ----bc3----
#endif
} grid_t;

/* BASIC OUTLINE 
rf = Tdx*v - f; 		% Compute residual
v = v - omega*Ddx∖rf;	% Pre-smoothing Jacobi
rf = Tdx*v - f;			% Compute residual
rf = lowpass(rf);		% Optional LP ﬁlter
rc = FMGrestrict(rf);	% Restrict to coarse
ec =FMGV(0,T2dx,rc);	% Solve error equation
ef = FMGprolong(ec);	% Prolong to ﬁne grid
v = v - ef;				% Correct: remove error
rf = Tdx*v - f;			% Compute residual
v = v - omega*Ddx\rf;	% Post-smoothing Jacobi
*/

/* --------------------- MAIN-FUNCTIONS  --------------------- */

/* Recursive function */
void FMGVh(grid_t* g, data_t gamma);
/* Makes first pre-smooting and calls the recursive function. */
void FMGV(grid_t* g, data_t gamma);

/*  ------------------------- DRIVERS ------------------------ */

/* Allocates memory for the problem and sets pointers to correct 
 * initial positions. 
 */ 
grid_t* initiate_grid(long k, long r);

/* First guess of v 
 */
void guess_g(grid_t* g, data_t gamma);

/* Makes the calculations before the recursive call: 
	rf = Tdx*v - f; 		% Compute residual
	v = v - omega*Ddx∖rf;	% Pre-smoothing Jacobi
	rf = Tdx*v - f;			% Compute residual
 */
void prereq(grid_t* g, data_t gamma);

/* Makes the calculations after the recursive call 
	v = v - ef;				% Correct: remove error
	rf = Tdx*v - f;			% Compute residual
	v = v - omega*Ddx\rf;	% Post-smoothing Jacobi
 */
void postreq(grid_t* g, data_t gamma);

/*  ----------------------- ROUTINES  ----------------------- */

/* 2D convolution with the following kernel layout: 
		k3 k1 k3
		k1 k2 k1
		k3 k1 k3
 */
void conv2(grid_t* g, data_t* pt_org, data_t* pt_dest, data_t k1, data_t k2, data_t k3);

/* Prolong grid - used in inital guess of v. 
 */
void prolong_init(grid_t* g);

/* Iterative conjugate-gradient method for solving 
 * the linear system on the smallest grid. 
 */
void conjgrad(grid_t* g, data_t gamma);

/* Lowpass and restrict grid 
 */
void lp_rest2(grid_t* g);

/* Prolong grid 
 */
void prolong_grid_2d(grid_t* g);

/*  ---------------------- SUBROUTINES  ---------------------- */

/* Returns the index of an element on a finer grid. 
*/ 
pow_t lthmap(pow_t k, pow_t n, pow_t N);

/* Get element (x, y) from a.
 */
data_t* eget(data_t* a, pow_t c, pow_t x, pow_t y);

pow_t bcmap(pow_t i, pow_t d);

void vecsub2(data_t* a, data_t* b, pow_t c);

void vecadd(data_t* a, data_t* b, pow_t c);

/*  ----------------------- DEBUGGING ----------------------- */

void print_2d(grid_t* g);

void print_all(grid_t* g);
#endif
