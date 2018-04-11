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
    int nphi,nel,J,lambda;
    double lx,ly;
    unsigned seed;
    // aspect ratio lx/ly
    double gamma;

    nphi=4;
    nel=2;
    lambda=200;
    J=-1;
    gamma=1;

    init_argv(nphi,nel,J,gamma,lambda,argc,argv);
    // keep lx/ly=nel/4
    //gamma=nel/4.0;
    ly=sqrt(2.0*M_PI*nphi/gamma);
    lx=ly*gamma;

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    basis sector(nphi,nel,J);
    sector.init();

    /* basis check */

    //sector.prlong();
    hamil config;
    config.set_hamil(sector,lx,ly,nphi);
    //config.print_hamil(10);
    config.diag();
    config.print_eigen(10);
    cout<<"E_gs:= "<<config.ground_state_energy()/nel<<endl;
    for(int i=0;i<config.nHilbert;i++)
        if(abs(config.psi_0[i])>0.11){
          cout<<i<<" :   ";
          long c=sector.id[i];
          for(int n=0;n<nphi;n++)
             if((c>>n)%2==1)
                cout<<n+1<<" ";
          cout<<"    ";
          cout<<bitset<12>(sector.id[i]).to_string()<<" "<<config.psi_0[i]<<endl;
        }

    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi);
    cout<<"-----------Lanczos---------"<<endl;
    cout<<"lx: ="<<lx<<endl;
    cout<<"ly: ="<<ly<<endl;
    //lconfig.print_hamil(44);
    lconfig.coeff_explicit_update();
    //lconfig.print_lhamil(10);
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<"eigenvalues:= [";
    for(int i=0;i<10;i++)
       cout<<lconfig.eigenvalues[i]<<", ";
    cout<<", ...]"<<endl;
    cout<<"E_gs:="<<setprecision(5)<<lconfig.ground_state_energy()/nel<<endl;
    for(int i=0;i<lconfig.nHilbert;i++)
        if(abs(lconfig.psir_0[i])>0.11){
          cout<<i<<" :   ";
          long c=sector.id[i];
          for(int n=0;n<nphi;n++)
             if((c>>n)%2==1)
                cout<<n+1<<" ";
          cout<<"    ";
          cout<<bitset<12>(sector.id[i]).to_string()<<" "<<lconfig.psir_0[i]<<endl;
        }



    /* ground-state energy check */
    /*
    for(int nphi=nel*4-1;nphi>=nel*1.5;nphi-=1){
      ly=sqrt(8.0*M_PI*nphi/nel);
      lx=nel*ly/4.0;
      //double E_gs=0;
      //for(int j=0;j<nphi;j++){
        basis sector(nphi,nel);
        sector.init();
        seed=tmr.nanoseconds();
        lhamil lconfig(lambda,seed);
        lconfig.set_hamil(sector,lx,ly,nphi);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
       // if(lconfig.ground_state_energy()<E_gs)
        //  E_gs=lconfig.ground_state_energy();
      //}
      cout<<nel*1.0/nphi<<" "<<lconfig.ground_state_energy()/nel<<endl;
    }
    */




    return 0;
}
