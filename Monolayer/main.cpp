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
    int kx;
    unsigned seed;
    // aspect ratio lx/ly
    double gamma;

    nphi=4;
    nel=2;
    lambda=200;
    kx=-1;
    J=-1;
    gamma=1;

    init_argv(nphi,nel,J,kx,gamma,lambda,argc,argv);
    // keep lx/ly=nel/4
    gamma=nel/4.0;
    ly=sqrt(2.0*M_PI*nphi/gamma);
    lx=ly*gamma;

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    //basis sector(nphi,nel,J);
    //sector.init();

    /* basis check */
    //sector.prlong();
    //sector.generate(count,index,a);
    /*
    cout<<"nphi: = "<<nphi<<endl;
    cout<<"nel: = "<<nel<<endl;
    cout<<"lx: = "<<lx<<endl;
    cout<<"ly: = "<<ly<<endl;
    cout<<"J: = "<<J<<endl;
    cout<<"-----------Exact diag---------"<<endl;
*/
/*
    hamil config;
    config.set_hamil(sector,lx,ly,nphi);
    //config.print_hamil(10);
    config.diag();
    //config.print_eigen(10);
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
*/


    lhamil lconfig(lambda,seed);
  //  for(int j=0;j<nel;j++){
    basis sector(nphi,nel,J,kx);
    sector.init();
    //sector.prlong();

   // for(int k=0;k<nel;k++){
    lconfig.set_hamil(sector,lx,ly,nphi);
    //cout<<"-----------Lanczos-----------"<<endl;
    //lconfig.print_hamil(44);
    lconfig.coeff_explicit_update();
    //lconfig.print_lhamil(10);
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<"E_gs:="<<setprecision(5)<<lconfig.ground_state_energy()/nel<<endl;
    /*
    double K=sqrt(j*j+kx*kx*gamma*gamma)*2.0*M_PI/nphi/gamma;
    cout<<K<<" ";
    int i=0;
    int count=0;
    do{
       i++;
       if(lconfig.eigenvalues[i]-lconfig.eigenvalues[i-1]>1e-3) {
         cout<<lconfig.eigenvalues[i]<<" ";
         count++;
      }
    }while(count<2);
    cout<<endl;
    */
  // }
 //}


    /*
    cout<<"eigenvalues:= [";
    for(int i=0;i<10;i++)
       cout<<lconfig.eigenvalues[i]<<", ";
    cout<<", ...]"<<endl;
    */

    /*
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
        */







    /* ground-state energy check */

/*
    lhamil lconfig(lambda,seed);
    basis sector;
    for(int nphi=nel*4+1;nphi>=nel+1;nphi-=1){
      ly=sqrt(8.0*M_PI*nphi/nel);
      lx=nel*ly/4.0;
      //double E_gs=0;
      //for(int j=0;j<nphi;j++){
        sector.init(nphi,nel);
        seed=tmr.nanoseconds();
        lconfig.set_hamil(sector,lx,ly,nphi);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        //if(lconfig.ground_state_energy()<E_gs)
        //  E_gs=lconfig.ground_state_energy();
      //}
      cout<<nel*1.0/nphi<<" "<<lconfig.ground_state_energy()/nel<<endl;
      //cout<<nel*1.0/nphi<<" "<<E_gs/nel<<endl;
        sector.clear();
        lconfig.clear();
    }
*/




    return 0;
}
