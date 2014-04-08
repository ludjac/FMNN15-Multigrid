#include "mg_2d.h"

int main(int argc, char const *argv[])
{

	grid_t* g = initiate_grid(12, 5);
	pow_t i, k;

	// Set Boundary conditions
	for(i = 0; i<g->n; ++i){
		g->bc0[i] = 0;
		g->bc1[i] = 0;
		g->bc2[i] = 0;
		g->bc3[i] = 0;
	}

	g->bc0[g->n+1] = 0;
	g->bc0[g->n+2] = 0;
	g->bc3[g->n+1] = 0;
	g->bc3[g->n+2] = 0;

	g->f = malloc(sizeof(data_t) * g->n * g->n);
	for(i = 0; i<g->n; ++i){
		for(k = 0; k<g->n; ++k){
			*eget(g->f, g->c, k, i) = cos(((double)((i+k)/(g->n * g->n))));
		}
	}

	data_t* td = g->g;
	data_t* tr = g->rf;

	data_t gamma = 1;
	guess_g(g, gamma);
//	printf("går in i FMGV\n");
	//print_2d(g);
	printf("%.16f\n", g->g[(g->n * g->n) / 2]);
	for(i = 0; i<1; ++i){
		FMGV(g, 1);
		//for(data_t* tp = g->rf; tp < g->lastp; ++tp){
		//	*tp = 0;
		//}
	//	print_2d(g);
		printf("%.16f\n", g->g[(g->n * g->n) / 2]);
	}
	//print_2d(g);
	
	return !(td == g->g);
}