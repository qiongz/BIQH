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


    lhamil lconfig(lambda,seed);
    /*
    lconfig.set_hamil(sector,lx,ly,nphi,d);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<d<<" "<<lconfig.ground_state_energy()/nel<<endl;
    */
    /* ground-state energy check */

    for(int n=0;n<=20;n++){
      //config.set_hamil(sector,lx,ly,nphi,i*0.2+0.2);
      //basis sector(n,n/2,n/2);
      //sector.init();
      lconfig.set_hamil(sector,lx,ly,nphi,n*0.1+0.001);
      //lconfig.print_hamil(4);
      //lconfig.coeff_update();
      lconfig.coeff_explicit_update();
      lconfig.diag();
      lconfig.eigenstates_reconstruction();
      //config.set_hamil(sector,lx,ly,nphi,d);
      //config.print_hamil();
      //config.set_hamil(sector,lx,ly,nphi,n*0.1+0.05);
      //config.diag();
      //cout<<i*0.2+0.2<<" "<<config.ground_state_energy()<<endl;
      //cout<<d<<" "<<config.ground_state_energy()<<endl;
      // long nHilbert=pow(sector.factorial(nel,nel/2),2);
      // double mem_actual,mem_theory;
       //mem_theory=(2*nHilbert*nel+5*nHilbert)*8.0/1.0e9;
      cout<<n*0.1+0.001<<" "<<lconfig.ground_state_energy()/nel<<endl;
       //mem_actual=(lconfig.H.inner_indices.size()+lconfig.H.value.size()+lconfig.psir_0.size()+lconfig.Coulomb_matrix.size()*4)/1.0e9;
      //cout<<n*0.1+0.05<<" "<<mem_actual<<" "<<mem_theory<<endl;
    }



  /*
    cout<<"# nel    GB"<<endl;
    for(int n=8;n<=20;n++){
       long nHilbert=pow(sector.factorial(n,n/2),2);
       double memsize;
       memsize=(2*nHilbert*n+5*nHilbert)*8.0/1.0e9;
       cout<<n<<"  "<<memsize<<endl;
    }
    */



    return 0;
}
