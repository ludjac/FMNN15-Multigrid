#include "mg_2d.h"

/* TYPEDEFS STRUCTS AND MACROS */

#define omega (0.5)

// Se header file mg_2d.h

/* --------------------- MAIN-FUNCTIONS  --------------------- */

void FMGVh(grid_t* g, data_t gamma)
{
	if(g->c == g->r){
		conjgrad(g, gamma);
	}else{

		prereq(g, gamma);

		lp_rest2(g);

		FMGVh(g, gamma);

		prolong_grid_2d(g); 

		postreq(g, gamma);
	}
}

void FMGV(grid_t* g, data_t gamma)
{
	pow_t n = g->n;
	data_t dx2 = 1/(( (data_t) n+1)*( (data_t) n+1));

	data_t Tdx1 = 1/dx2;
	data_t Tdx2 =  -4/dx2 + gamma;
	data_t D_inv = 1/Tdx2;

	pow_t i;

	conv2(g, g->g, g->rf, Tdx1, Tdx2, 0);

	for(i = 0; i< n*n; ++i){
		g->rf[i] -= g->f[i];
		g->g[i] -= omega * D_inv * g->rf[i];
	}

	conv2(g, g->g, g->rf, Tdx1, Tdx2, 0);	

	for(i = 0; i< n*n; ++i){
		g->rf[i] -= g->f[i];
		g->g[i] = -g->g[i];
	}

	lp_rest2(g);

	FMGVh(g, gamma);

	prolong_grid_2d(g);

	for(i = 0; i < n*n; ++i){
		g->g[i] = g->rf[i];
	}	

	conv2(g, g->rf, g->rf, Tdx1, Tdx2, 0);

	for(i = 0; i < n*n; ++i){
		g->rf[i] = -g->rf[i];
		g->rf[i] -= g->f[i];
		g->g[i] += omega * D_inv * g->rf[i];
	}

}

/*  ------------------------- DRIVERS ------------------------ */

grid_t* initiate_grid(long k, long r)
{
	pow_t sum = 0;
	pow_t n = (1<<k)-1;
	sum += 4*n+4; //begin conds
	sum += 2*n*n;
	grid_t* g = malloc(sizeof(grid_t));
	g->sum = sum;
	g->k = k;
	g->c = k;
	g->n = n;
	g->r = r;
	g->bc0 = malloc(sizeof(data_t) * sum);
	if (g->bc0 == NULL){
		free(g);
		fprintf(stderr, "%s\n", "Malloc failed!");
		exit(0);
	}
	g->bc1 = g->bc0 + n + 2;
	g->bc2 = g->bc1 + n;
	g->bc3 = g->bc2 + n;
	g->g = g->bc3 + n + 2;
	g->rf = g->g+n*n;
	return g;
}

void guess_g(grid_t* g, data_t gamma)
{
	pow_t N = g->n;
	g->c = g->r;
	g->n = (1<<g->c)-1;
	pow_t i;
	for(i = 0; i < g->n * g->n; ++i){
		g->g[i] = g->f[lthmap(i, g->n, N)];
	}
	conjgrad(g, gamma);
	while(g->c < g->k) prolong_init(g);
}

void prereq(grid_t* g, data_t gamma)
{
	pow_t n = g->n;
	pow_t i;

	data_t dx2 = 1/(( (data_t) n+1)*( (data_t) n+1));
	data_t Tdx1 = 1/dx2;
	data_t Tdx2 =  -4/dx2 + gamma;
	data_t D_inv = 1/Tdx2;

	// ORG, DEST
	conv2(g, g->g, g->rf,  omega*D_inv*Tdx1, omega*D_inv*Tdx2, 0);

	for(i = 0; i < n*n; ++i){
		g->rf[i] -= g->g[i];
		g->g[i] *= omega*-D_inv; 		
	}
	

}

void postreq(grid_t* g, data_t gamma)
{
	pow_t n = g->n;
	pow_t i;
	data_t dx2 = 1/(( (data_t) n+1)*( (data_t) n+1));

	data_t Tdx1 = 1/dx2;
	data_t Tdx2 =  -4/dx2 + gamma;
	data_t D_inv = 1/Tdx2;

	// u = u -ef is in g now
	//org, dest

	conv2(g, g->rf, g->rf, -omega*D_inv*Tdx1, (1-omega*D_inv*Tdx2), 0);

	for(i = 0; i < n*n; ++i){
		g->rf[i] += omega*D_inv*g->g[i];
		g->g[i] += g->rf[i];

	}

}

/*  ----------------------- ROUTINES  ----------------------- */

void conv2(grid_t* g, data_t* pt_org, data_t* pt_dest, data_t k1, data_t k2, data_t k3)
{
	#if 0	
		k3 k1 k3
		k1 k2 k1
		k3 k1 k3
	#endif
	int n, m, i;
	
	/* Using Stack */ 
	data_t top_row_center[g->n+2];
	data_t top_row_left[g->n+2];
	data_t left; 	
	data_t center;	
	data_t temp4;
	data_t bottom_left;


#if 0
    ---bc0---
    |	 	|
    b	 	b
    c	 	c
    1	 	2
    |		|
    ---bc3---
#endif

	for (i = 0; i <= g->n+1; ++i)
	{
		top_row_center[i] = g->bc0[i];
		top_row_left[i] = g->bc0[i]; 
	}

	/* All rows but last */
	for (m = 0; m < g->n-1; ++m)
	{
		left = g->bc1[m];
		bottom_left = g->bc1[m+1];

		for (n = 0; n < g->n-1; n++)
		{
			center = pt_org[m*g->n + n];
			pt_dest[m*g->n + n] = 
					k3*top_row_left[n] +				// TOP-LEFT
					k1*top_row_center[n] +				// TOP
					k3*top_row_center[n+1] +			// TOP-RIGHT
					k1*left +							// LEFT
					k2*center +							// CENTER
					k1*pt_org[m*g->n + n + 1] +			// RIGHT
					k3*bottom_left + 					// BOTTOM-LEFT
					k1*pt_org[m*g->n + n + g->n] + 		// BOTTOM
					k3*pt_org[m*g->n + n + g->n + 1];	// BOTTOM-RIGHT

			top_row_left[n] = left;
			top_row_center[n] = center;
			left = center;
			bottom_left = pt_org[m*g->n + n + g->n];
		}

		left  = center;  
		center = pt_org[m*g->n + n];

		pt_dest[m*g->n + n] = 
				k3*top_row_left[n] + 				// TOP-LEFT
				k1*top_row_center[n] +				// TOP
				k3*g->bc2[m] + 						// TOP-RIGHT
				k1*left +							// LEFT
				k2*center +							// CENTER
				k1*g->bc2[m+1] +					// RIGHT
				k3*bottom_left + 					// BOTTOM-LEFT
				k1*pt_org[m*g->n + n + g->n] +		// BOTTOM
				k3*g->bc2[m+2];						// BOTTOM-RIGHT
		top_row_center[n] = center;
		top_row_left[n] = left;
	}

	/* Last row */
	m = g->n-1;
	left = g->bc1[m];
	bottom_left = g->bc1[m+1];

	for (n = 0; n < g->n-1; ++n)
	{
		center = pt_org[m*g->n + n];
		pt_dest[m*g->n + n] = 
				k3*top_row_left[n] + 		// TOP-LEFT
				k1*top_row_center[n] +		// TOP
				k3*top_row_center[n+1] + 	// TOP-RIGHT
				k1*left +					// LEFT
				k2*center +					// CENTER
				k1*pt_org[m*g->n + n + 1] +	// RIGHT
				k3*0 + 						// BOTTOM-LEFT
				k1*0;						// BOTTOM
				k3*0;						// BOTTOM-RIGHT
		top_row_left[n] = left;
		top_row_center[n] = center;
		left = center;
		bottom_left = pt_org[m*g->n + n + g->n];

	}

	center = pt_org[m*g->n + n];
	
	pt_dest[m*g->n + n] = 
			k3*top_row_left[n] + 				// TOP-LEFT
			k1*top_row_center[n] +				// TOP
			k3*g->bc2[m-1] + 					// TOP-RIGHT
			k1*left +							// LEFT
			k2*center +							// CENTER
			k1*g->bc2[m] +						// RIGHT
			k3*g->bc3[m-1] + 					// BOTTOM-LEFT
			k1*g->bc3[m] +						// BOTTOM
			k3*g->bc2[m+1];						// BOTTOM-RIGHT

}

void prolong_init(grid_t* g)
{
	pow_t i, k, N;
	N = (g->n << 1) + 1;

	// ok det här är kanske lite oläsligt,
	// eftersom i är unsigned kommer det bli underflow på --i när i är noll,
	// i kommer bli jättestort och loopen avbryts.
	for(i = (g->n)*(g->n)-1; i < (g->n)*(g->n); --i){ 
		g->g[lthmap(i, g->n, N)] = g->g[i];
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-x-x-x-      x = initated value
			-------
			-x-x-x-
			-------
			-x-x-x-
			-------
		If we go from a 3 to a 7 grid. 
	*/
	g->c++;
	g->n = (g->n << 1) + 1;
	data_t* t1;
	data_t* t2;
	for(k = 1; k<N; k += 2){
		for(i = 2; i<N-1; i += 2){
			t1 = eget(g->g, g->c, i, k);
			*t1 = t1[-1]/2 + t1[1]/2;
		}
	}
	for(k = 2; k<N-1; k += 2){
		for(i = 1; i<N; i += 2){
			t2 = eget(g->g, g->c, i, k);
			*t2 = t2[-N]/2 + t2[N]/2;
		}
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-xxxxx-      x = initated value
			-x-x-x-
			-xxxxx-
			-x-x-x-
			-xxxxx-
			-------
		If we go from a 3 to a 7 grid. 
	 */ 
	for(k = 2; k<N-2; k += 2){
		for(i = 2; i<N-2; i += 2){
			t1 = eget(g->g, g->c, i, k);
			*t1 = t1[1]/4 + t1[-1]/4 + t1[N]/4 + t1[-N]/4;
		}
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-xxxxx-      x = initated value
			-xxxxx-
			-xxxxx-
			-xxxxx-
			-xxxxx-
			-------
		If we go from a 3 to a 7 grid. 
	 */
	data_t t3;
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->g, g->c, i, 0);
		t3 = g->bc0[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->g, g->c, i, 1) /2;
	}
	
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->g, g->c, 0, i);
		t3 = g->bc1[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->g, g->c, 1, i) /2;
	}
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->g, g->c, N-1, i);
		t3 = g->bc2[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->g, g->c, N-2, i) /2;
	}
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->g, g->c, i, N-1);
		t3 = g->bc3[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->g, g->c, i, N-2) /2;
	}

	*eget(g->g, g->c, 0, 0) = 
		g->bc0[0]/4 + 
		g->bc1[0]/4 + 
		*eget(g->g, g->c, 1, 0)/4 + 
		*eget(g->g, g->c, 0, 1)/4;

	*eget(g->g, g->c, 0, N-1) = 
		g->bc0[0]/4 + 
		g->bc2[0]/4 + 
		*eget(g->g, g->c, N-2, 0)/4 + 
		*eget(g->g, g->c, N-1, 1)/4;
	
	*eget(g->g, g->c, N-1, 0) = 
		g->bc1[0]/4 + 
		g->bc3[0]/4 + 
		*eget(g->g, g->c, N-2, 0)/4 + 
		*eget(g->g, g->c, N-1, 1)/4;
	*eget(g->g, g->c, N-1, N-1) = 
		g->bc2[0]/4 + 
		g->bc3[0]/4 + 
		*eget(g->g, g->c, N-2, N-1)/4 + 
		*eget(g->g, g->c, N-1, N-2)/4;
}

void conjgrad(grid_t* g, data_t gamma)
{
	pow_t i, k, n;
	n = g->n;
	data_t p[n*n];
	data_t Ap[n*n];
	// r is ptr_b

	data_t rsold; 
	data_t rsnew;
	data_t pAp;
	data_t alpha;

	data_t dx2 = 1/(( (data_t) n+1)*( (data_t) n+1));

	data_t Tdx1 = 1/dx2;
	data_t Tdx2 =  -4/dx2 + gamma;
	data_t D_inv = 1/Tdx2;
	
	/* TOTAL STACK SIZE (not pow_t's): (2*(n*n) + 4) */

	// INITIATE
	rsold = 0;
	pAp = 0;
	alpha = 0;
	rsnew = 0;


	for(i = 0; i<n*n; ++i){
		g->rf[i] = 0.0;
		p[i] = g->g[i];
		rsold += g->g[i]*g->g[i];
	}

	// FOR-LOOP
	for (k = 0; k<1E6; ++k){
		conv2(g, p, Ap, Tdx1, Tdx2, 0);

		
		pAp = 0;
		for(i = 0; i<n*n; ++i){
			pAp += p[i]*Ap[i];
		}
		
		alpha = rsold/pAp;
		rsnew = 0;
		for(i = 0; i<n*n; ++i){
			g->rf[i] += alpha*p[i];
			g->g[i] -= alpha*Ap[i];
			rsnew += g->g[i]*g->g[i];
		}

		// BREAK IF rsnew is small enough
		if (sqrt(rsnew)<1E-2){
			for (i=0; i<n*n; ++i){
				g->g[i] = -g->rf[i];
			}
			return;
		}
		
		for(i = 0; i<n*n; ++i){
			p[i] = g->g[i] + (rsnew/rsold)*p[i];
		}
		rsold = rsnew;
	}
}

void lp_rest2(grid_t* g)
{
	pow_t N = g->n;
	conv2(g, g->rf, g->rf, 0.125, 0.25, 0.0625);

	g->c--;
	g->n = g->n>>1;
	pow_t n2 = (g->n)*(g->n);
	pow_t i;
	for(i = 0; i<n2; ++i){
		g->rf[i] = g->rf[lthmap(i, g->n, N)];
	}
	g->g = g->rf;
	g->rf += n2;
	
}

void prolong_grid_2d(grid_t* g)
{
	pow_t i, k, N;
	N = (g->n << 1) + 1;
	g->rf = g->g;
	g->g = g->g-(N)*(N);
	for(i = (g->n)*(g->n)-1; i < (g->n)*(g->n); --i){ // ok det här är kanske lite oläsligt, eftersom i är unsigned kommer det bli underflow på --i när i är noll, i kommer bli jättestort och loopen avbryts.
		g->rf[lthmap(i, g->n, N)] = g->rf[i];
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-x-x-x-      x = initated value
			-------
			-x-x-x-
			-------
			-x-x-x-
			-------
		If we go from a 3 to a 7 grid. 
	*/
	g->c++;
	g->n = (g->n << 1) + 1;
	data_t* t1;
	data_t* t2;
	for(k = 1; k<N; k += 2){
		for(i = 2; i<N-1; i += 2){
			t1 = eget(g->rf, g->c, i, k);
			*t1 = t1[-1]/2 + t1[1]/2;
		}
	}
	for(k = 2; k<N-1; k += 2){
		for(i = 1; i<N; i += 2){
			t2 = eget(g->rf, g->c, i, k);
			*t2 = t2[-N]/2 + t2[N]/2;
		}
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-xxxxx-      x = initated value
			-x-x-x-
			-xxxxx-
			-x-x-x-
			-xxxxx-
			-------
		If we go from a 3 to a 7 grid. 
	 */ 
	for(k = 2; k<N-2; k += 2){
		for(i = 2; i<N-2; i += 2){
			t1 = eget(g->rf, g->c, i, k);
			*t1 = t1[1]/4 + t1[-1]/4 + t1[N]/4 + t1[-N]/4;
		}
	}
	/*
		The grid should now look like this:
			-------      - = unitiated value
			-xxxxx-      x = initated value
			-xxxxx-
			-xxxxx-
			-xxxxx-
			-xxxxx-
			-------
		If we go from a 3 to a 7 grid. 
	 */
	data_t t3;
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->rf, g->c, i, 0);
		t3 = g->bc0[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->rf, g->c, i, 1) /2;
	}
	
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->rf, g->c, 0, i);
		t3 = g->bc1[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->rf, g->c, 1, i) /2;
	}
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->rf, g->c, N-1, i);
		t3 = g->bc2[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->rf, g->c, N-2, i) /2;
	}
	for(i = 1; i<N-1; ++i){
		t1 = eget(g->rf, g->c, i, N-1);
		t3 = g->bc3[bcmap(i, g->k - g->c)];
		*t1 = t3 / 2 + *eget(g->rf, g->c, i, N-2) /2;
	}

	*eget(g->rf, g->c, 0, 0) = g->bc0[0]/4 + g->bc1[0]/4 + *eget(g->rf, g->c, 1, 0)/4 + *eget(g->rf, g->c, 0, 1)/4;
	*eget(g->rf, g->c, 0, N-1)= g->bc0[0]/4 + g->bc2[0]/4 + *eget(g->rf, g->c, N-2, 0)/4 + *eget(g->rf, g->c, N-1, 1)/4;
	*eget(g->rf, g->c, N-1, 0)= g->bc1[0]/4 + g->bc3[0]/4 + *eget(g->rf, g->c, N-2, 0)/4 + *eget(g->rf, g->c, N-1, 1)/4;
	*eget(g->rf, g->c, N-1, N-1)= g->bc2[0]/4 + g->bc3[0]/4 + *eget(g->rf, g->c, N-2, N-1)/4 + *eget(g->rf, g->c, N-1, N-2)/4;

	vecsub2(g->rf, g->g, N*N);
}

/*  ---------------------- SUBROUTINES  ---------------------- */

pow_t lthmap(pow_t k, pow_t n, pow_t N)
{
	return N*((k/n)*2+1)+1+((k%n)*2);
}

data_t* eget(data_t* a, pow_t c, pow_t x, pow_t y)
{
	return (a + x + (y<<c) - y);
}

pow_t bcmap(pow_t i, pow_t d)
{
	return (i<<d)+(1<<d)-1;
}

void vecsub2(data_t* a, data_t* b, pow_t c)
{
	pow_t i;
	for(i = 0; i < c; ++i){
		a[i] = b[i] - a[i];
	}
}
void vecadd(data_t* a, data_t* b, pow_t c)
{
	pow_t i;
	for(i = 0; i < c; ++i){
		a[i] += b[i];
	}
}

/*  ----------------------- DEBUGGING ----------------------- */

void print_2d(grid_t* g)
{
	pow_t i, k;
	printf("N = %lu\n", g->n);
	printf("g = [\n");
	for(i = 0; i<g->n; ++i){
		for(k = 0; k<g->n; ++k){
			printf("%.8f, ", *eget(g->g, g->c, k, i));
		}
		printf(";\n");
	}
	printf("];\n");
#if 1
	printf("rf = [\n");
	for(i = 0; i<g->n; ++i){
		for(k = 0; k<g->n; ++k){
			printf("%.8f, ", *eget(g->rf, g->c, k, i));
		}
		printf(";\n");
	}
	printf("];\n");
#endif
}

void print_all(grid_t* g)
{
	unsigned long n = (1<<g->k)-1;
	pow_t i;
	printf("bc0\t\tbc1\t\tbc2\t\tbc3\n");
	for(i = 0; i<n; ++i){
		printf("%f\t%f\t%f\t%f\n", g->bc0[i], g->bc1[i], g->bc2[i], g->bc3[i]);
	}
	printf("\n");
	for(i = 0; i<1*n*n; ++i){
		printf("%f\t\t%f\t\t%f\n", g->bc0[4*n+i], g->bc0[4*n+i+n*n], g->bc0[4*n+i+2*n*n]);
	}
}
