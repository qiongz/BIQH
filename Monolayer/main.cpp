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
    int nLL,nphi,nel,J,lambda;
    double lx,ly;
    int kx;
    unsigned seed;
    // aspect ratio lx/ly
    double gamma;

    nLL=0;
    nphi=4;
    nel=2;
    lambda=400;
    kx=-1;
    J=-1;
    gamma=1;

    init_argv(nLL,nphi,nel,J,kx,gamma,lambda,argc,argv);
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

    
        basis sector(nphi,nel,J,kx);
        sector.init();

        //sector.prlong();
        //sector.generate(count,index,a);

        cout<<"nLL: ="<<nLL<<endl;
        cout<<"nphi: = "<<nphi<<endl;
        cout<<"nel: = "<<nel<<endl;
        cout<<"lx: = "<<lx<<endl;
        cout<<"ly: = "<<ly<<endl;
        cout<<"J: = "<<J<<endl;
        cout<<"K: = "<<kx<<endl;
        cout<<"nHilbert: ="<<sector.nbasis<<endl;
        cout<<"-----------Exact Diag---------"<<endl;


/*
        hamil config;
        //basis sector(nphi,nel,J,kx);
        //sector.init();
        //sector.prlong();
        config.set_hamil(sector,lx,ly,nphi,nLL);
       // config.print_hamil(8);
        config.diag();
        //config.print_eigen(10);
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

    
     cout<<"-----------Lanczos---------"<<endl;
     //basis sector(nphi,nel,J,kx);
     //sector.init();
     lhamil lconfig(lambda,seed);
     lconfig.set_hamil(sector,lx,ly,nphi,nLL);
     //lconfig.print_hamil(4);
     lconfig.coeff_update();
     //lconfig.coeff_explicit_update();
     //lconfig.print_lhamil(4);
     lconfig.diag();
     lconfig.eigenstates_reconstruction();
     lconfig.print_eigen(10);
     cout<<"E_gs(LANCZOS):="<<setprecision(6)<<lconfig.ground_state_energy()/nel<<endl;
    



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

/*
    basis sector0(nphi,nel,0,0);
    double s0,t0;
    t0=0;
    int N=sector0.common_divisor(nphi,nel);
    s0=N/2*nLL;

    hamil config;
    //lhamil lconfig(lambda,seed);
    for(int t=0; t<=N/2; t++) {
        for(int s=0; s<N; s++) {
            basis sector(nphi,nel,t,s);
            sector.init();
            config.set_hamil(sector,lx,ly,nphi,nLL);
            config.diag();

            double K=sqrt(gamma*gamma*(t-t0)*(t-t0)+(s-s0)*(s-s0))*sqrt(2.0*M_PI/nphi/gamma);

            cout<<K<<" ";
            for(int n=0; n<4; n++)
                cout<<config.eigenvalues[n]<<" ";
            cout<<endl;
        }
    }
*/



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
