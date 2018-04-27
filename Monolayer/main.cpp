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
    lambda=400;
    kx=-1;
    J=-1;
    gamma=1;

    init_argv(nphi,nel,J,kx,gamma,lambda,argc,argv);
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
    //basis sector(nphi,nel,J,kx);
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
    cout<<"K: = "<<kx<<endl;
    cout<<"nHilbert: ="<<sector.nbasis<<endl;
    cout<<"-----------Exact Diag---------"<<endl;



    hamil config;
    //basis sector(nphi,nel,J,kx);
    //sector.init();
    //sector.prlong();
    config.set_hamil(sector,lx,ly,nphi);
    //config.print_hamil(4);
    config.diag();
    config.print_eigen(10);
    cout<<"E_gs:= "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
    */

    /*
    for(int i=0;i<config.nHilbert;i++)
        if(abs(config.psi_0[i])>0.11){
          cout<<i<<" :   ";
          long c=sector.id[i];
          for(int n=0;n<nphi;n++)
             if((c>>n)%2==1)
                cout<<n+1<<" ";
          cout<<"    ";
          cout<<bitset<12>(sector.id[i]).to_string()<<" "<<abs(config.psi_0[i])<<endl;
        }

   */

   /*
   cout<<"-----------Lanczos---------"<<endl;
   //basis sector(nphi,nel,J,kx);
   //sector.init();
   lhamil lconfig(lambda,seed);
   lconfig.set_hamil(sector,lx,ly,nphi);
   //lconfig.print_hamil(4);
   lconfig.coeff_update();
   //lconfig.coeff_explicit_update();
   //lconfig.print_lhamil(4);
   lconfig.diag();
   lconfig.eigenstates_reconstruction();
   lconfig.print_eigen(10);
   cout<<"E_gs(LANCZOS):="<<setprecision(6)<<lconfig.ground_state_energy()/nel<<endl;
    
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


    double j0,k0;
    if(nel%2==0)
       j0=k0=nel/2;
    else
       j0=k0=0;

    for(int j=0;j<nphi;j++){
    for(int k=0;k<nel;k++){

    basis sector(nphi,nel,j,k);
    sector.init();
    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    double K=sqrt((j%sector.C-j0)*(j%sector.C-j0)+(k%sector.C-k0)*(k%sector.C-k0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);

    cout<<K<<" "<<lconfig.eigenvalues[0]<<" ";
    int i=0;
    int count=0;
    do{
       i++;
       if(lconfig.eigenvalues[i]-lconfig.eigenvalues[i-1]>1e-5) {
         cout<<lconfig.eigenvalues[i]<<" ";
         count++;
      }
    }while(count<10 && count<lconfig.norm.size());
    cout<<endl;
    }
    }

    /*
    for(int j=0;j<nphi;j++){
    for(int k=0;k<nel;k++){
    basis sector(nphi,nel,j,k);
    sector.init();
    hamil config;
    config.set_hamil(sector,lx,ly,nphi);
    config.diag();

    double K=sqrt((j%sector.C-j0)*(j%sector.C-j0)+(k%sector.C-k0)*(k%sector.C-k0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
   cout<<K<<" ";
   for(int s=0;s<10;s++)
     cout<<config.eigenvalues[s]<<" ";
    cout<<endl;
    }
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
