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
    int nphi,nel,nel_up,nel_down,lambda;
    double lx,ly;
    int Sz;
    unsigned seed;
    double d,mu;
    double func_Lieb_Wu(double x, void *params);
    double Integrate_Lieb_Wu(double U);

    nphi=2;
    nel=2;
    nel_up=1;
    nel_down=1;
    d=1;
    lambda=200;

    init_argv(nphi,nel,nel_up,d,lambda,argc,argv);
    //ly=2.0;
    //lx=nphi/ly;
    lx=ly=sqrt(nphi*2.0*M_PI);
    nel_down=nel-nel_up;

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    Sz=nel_up-nel_down;
    basis sector(nphi,nel_up,nel_down);
    sector.init();
    //sector.prlong();
    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    //config.coeff_explicit_update();
    //config.diag();
    //config.eigenstates_reconstruction();

    /* basis check */

    hamil config;
    lhamil lconfig(sector, d, lambda, seed);
    //config.print_hamil();
    //for(int i=0;i<config.nHilbert;i++)
    //   cout<<config.eigenvalues[i]<<endl;
    //config.print_lhamil(6);
    //config.print_eigen(6);
    //config.save_to_file(filename.c_str());

    /* ground-state energy check */
    //for(int i=0;i<40;i++){
      lconfig.set_hamil(sector,lx,ly,nphi,d);
      //lconfig.set_hamil(sector,lx,ly,nphi,d);
      //lconfig.print_hamil(4);
      //lconfig.coeff_explicit_update();
      //lconfig.print_lhamil(4);
      //lconfig.diag();
      //lconfig.eigenstates_reconstruction();
      //config.set_hamil(sector,lx,ly,nphi,d);
      lconfig.print_hamil(lconfig.nHilbert);
      //config.diag();
      //cout<<i*0.1+0.01<<" "<<lconfig.ground_state_energy()/nsite<<endl;
     // cout<<i*0.1+0.01<<" "<<config.ground_state_energy()<<endl;
      //cout<<d<<" "<<config.ground_state_energy()<<endl;
      //cout<<d<<" "<<lconfig.ground_state_energy()<<endl;
    //}

    return 0;
}
