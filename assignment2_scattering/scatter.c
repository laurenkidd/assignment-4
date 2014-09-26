#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#ifndef N_SCATTER
#define N_SCATTER 1000000
#endif


int main(int argc, char **argv)
{
 //intialize pointers to all arrays
  srand(time(NULL));
  double *E = (double*) malloc(N_SCATTER * sizeof( double ));
  if (E == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *alpha = (double*) malloc(N_SCATTER * sizeof( double ));
  if (alpha == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *z1 = (double*) malloc(N_SCATTER * sizeof( double ));
  if (z1 == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *RND1 = (double*) malloc(N_SCATTER * sizeof( double ));
  if (RND1 == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *RND2 = (double*) malloc(N_SCATTER * sizeof( double ));
  if (RND2 == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
   double *phi = (double*) malloc(N_SCATTER * sizeof( double ));
  if (phi == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *psi = (double*) malloc(N_SCATTER * sizeof( double ));
  if (psi == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *sigma_e = (double*) malloc(N_SCATTER * sizeof( double ));
  if (sigma_e == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *lambda = (double*) malloc(N_SCATTER * sizeof( double ));
  if (lambda == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *step = (double*) malloc(N_SCATTER * sizeof( double ));
  if (step == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *dE_dS = (double*) malloc(N_SCATTER * sizeof( double ));
  if (dE_dS == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *delta_E = (double*) malloc(N_SCATTER * sizeof( double ));
  if (delta_E == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}
  double *e_target = (double*) malloc(N_SCATTER * sizeof( double ));
  if (e_target == NULL ) { puts("Cannot Allocate Result Array.");exit(1);}


  double pi = 3.14159;
  double Na = 6.022e23; //Avogadro's number
  double k = 0.1 ; //define something for beam spread
  double Z, A, rho ;
  int material;
  printf("Enter a material (gold = 1 or silicon = 2), or other for new material = 3)");
  scanf( "%d", &material ); 
  switch ( material ) {
    case 1:
      //gold
      Z = 79; //charge of nucleus  
      A = 196.9665; //atomic weight g/mol
      rho = 19.30; //density g/cm^3
      break;
    case 2:
      Z = 14; //charge of nucleus  
      A = 28.085; //atomic weight g/mol
      rho = 2.3290; //density g/cm^3
      break;
    case 3:
      printf("enter nucleus charge");
      scanf( "%le", &Z ); 
      getchar();
      printf("eneter atomic weight (g/mol)");
      scanf( "%le", &A ); 
      getchar();
      printf("enter density (g/cm^3)");
      scanf( "%le", &rho ); 
      getchar();
      break;
  }
//
  double J = (9.76 * Z + 58.5/pow(Z,0.19)) * pow(10.0,-3);
  double energy;
  printf("Enter a center energy (in kev)");
  scanf( "%le", &energy ); 
  getchar();
  double E_c = energy ; //center energy 
 
  for ( int i = 0; i < N_SCATTER; i++){
  // Random number between 0 and 1
    RND1[i] = rand()/(double)RAND_MAX ;
    RND2[i] = rand()/(double)RAND_MAX ;
// Gaussian deviate
    z1[i] = sqrt(-2 * log(RND1[i])) * cos(2 * pi * RND2[i]);
// Electron beam energy
    E[i] = E_c * ((double)1 + k * z1[i]);
// Screening parameter
    alpha[i] = 3.4e-3 * pow(Z,0.67)/E[i];
// phi
    phi[i] = acos(1e0 - (2e0 * alpha[i] * RND1[i]/(1e0+alpha[i] - RND1[i])));
// psi
    psi[i] = 2e0 * pi * RND1[i];
//sigma_e
    sigma_e[i] = 5.21e-21 * pow(Z,2)/pow(E[i],2) * 4 * pi/(alpha[i]*(1.0 + alpha[i])) * pow((E[i]+511),2)/pow((E[i]+1024),2); //cm
//mean free path
    lambda[i] = A/(Na * rho * sigma_e[i]); //cm
//step = scattering position
    step[i] = -lambda[i] * log(RND1[i]);
//dE/dS
    dE_dS[i] = -(double) 78500 * (rho * Z)/(A * E[i]) * log(1.166 * E[i]/J + 1.0); //keV
//delta E
    delta_E[i] = dE_dS[i]*step[i];
// energy as electron exits target
    e_target[i] = E[i] - delta_E[i];
    printf("%le %le %le %le %le\n", step[i], E[i], e_target[i], psi[i], phi[i]);
  }
   


 free(E);
 free(alpha);
 free(z1);
 free(RND1);
 free(RND2);
 free(phi);
 free(psi);
 free(sigma_e);
 free(lambda);
 free(dE_dS);
 free(delta_E);
 free(e_target);
  return 0;
}

