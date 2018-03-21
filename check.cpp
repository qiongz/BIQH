#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include<ctime>
#include<sstream>
#include<cstdlib>
#if __cplusplus > 199711L
#include<chrono>
#endif

using namespace std;

int main(int argc,char *argv[]) {
    int nsite,nel,nel_up,nel_down,lambda;
    int Sz;
    unsigned seed;
    double d,mu;
    double func_Lieb_Wu(double x, void *params);
    double Integrate_Lieb_Wu(double U);

    nsite=4;
    nel=4;
    nel_up=2;
    nel_down=2;
    d=1;
    lambda=200;

    init_argv(nsite,nel,nel_up,d,lambda,argc,argv);
    nel_down=nel-nel_up;

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    Sz=nel_up-nel_down;
    basis sector(nsite,nel_up,nel_down);
    sector.init();
    //sector.prlong();
    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    //lhamil lconfig(sector,d,lambda,seed);
    //config.coeff_explicit_update();
    //config.diag();
    //config.eigenstates_reconstruction();

    /* basis check */

    hamil config;
    //config.print_hamil();
    //for(int i=0;i<config.nHilbert;i++)
    //   cout<<config.eigenvalues[i]<<endl;
    //config.print_lhamil(6);
    //config.print_eigen(6);
    //config.save_to_file(filename.c_str());

    /* ground-state energy check */
    //for(int i=0;i<20;i++){
      //config.set_hamil(sector,i*0.1+0.01);
      config.set_hamil(sector,d);
      config.print_hamil();
      config.diag();
      //cout<<i*0.1+0.01<<" "<<config.ground_state_energy()/nsite<<endl;
      //cout<<d<<" "<<config.ground_state_energy()/nsite<<endl;
    //}

    return 0;
}
