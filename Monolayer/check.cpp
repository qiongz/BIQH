#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<ctime>
#include<sstream>
#include<cstdlib>
#if __cplusplus > 199711L
#include<chrono>
#endif

using namespace std;

int main(int argc,char *argv[]) {
    int nphi,nel,lambda;
    double lx,ly;
    int Sz;
    unsigned seed;
    double d,mu;

    nphi=2;
    nel=2;
    lambda=200;

    init_argv(nphi,nel,lambda,argc,argv);
    lx=ly=sqrt(nphi*2.0*M_PI);

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    basis sector(nphi,nel);
    sector.init();

    /* basis check */
    //hamil config;
    //config.set_hamil(sector,lx,ly,nphi,d);
    //sector.prlong();
    //config.diag();
    //lhamil lconfig(sector, d, lambda, seed);
    //lconfig.set_hamil(sector,lx,ly,nphi,d);
    //lconfig.coeff_update();
    //lconfig.coeff_explicit_update();

    //cout<<"Hamiltonian in full diagonalization routine"<<endl;
    //config.print_hamil(9);

    //cout<<"Hamiltonian in Lanczos diagonalization routine"<<endl;
    //lconfig.print_hamil(lconfig.nHilbert);

    //lconfig.diag();
    //lconfig.eigenstates_reconstruction();
    //cout<<"# d  E_gs"<<endl;
    //cout<<d<<" "<<config.ground_state_energy()<<endl;
    //cout<<d<<" "<<lconfig.ground_state_energy()/nel<<endl;
    //for(int i=0;i<config.nHilbert;i++)
    //   cout<<config.eigenvalues[i]<<" "<<lconfig.eigenvalues[i]<<endl;
    //config.print_lhamil(6);
    //config.print_eigen(6);
    //config.save_to_file(filename.c_str());


    /*
    lconfig.set_hamil(sector,lx,ly,nphi,d);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<d<<" "<<lconfig.ground_state_energy()/nel<<endl;
    */

    /* ground-state energy check */



    lhamil lconfig(lambda,seed);
    for(nel=2;nel<nphi;nel++){
      basis sector(nphi,nel);
      sector.init();
      lconfig.set_hamil(sector,lx,ly,nphi);
      lconfig.coeff_explicit_update();
      lconfig.diag();
      lconfig.eigenstates_reconstruction();
      cout<<nel*1.0/nphi<<" "<<lconfig.ground_state_energy()/nel<<endl;
    }


    return 0;
}
