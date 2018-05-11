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
    int nLL,nphi,nel,nel_up,nel_down,J,kx,lambda;
    double lx,ly,gamma;
    int Sz;
    unsigned seed;
    double d,mu;

    nLL=0;
    nphi=8;
    nel=4;
    nel_up=2;
    nel_down=2;
    gamma=1.0;
    J=-1;
    kx=-1;

    d=100.0;
    lambda=100;
    init_argv(nLL,nphi,nel,nel_up,J,kx,d,gamma,lambda,argc,argv);
    //gamma=nel_up/4.0;
    ly=sqrt(nphi*2.0*M_PI/gamma);
    lx=ly*gamma;

    nel_down=nel-nel_up;

     #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    Sz=nel_up-nel_down;


    basis sector(nphi,nel_up,nel_down,J,kx);
    sector.init();

    
    cout<<"nphi: = "<<nphi<<endl;
    cout<<"nel_up: = "<<nel_up<<endl;
    cout<<"nel_down: = "<<nel_up<<endl;
    cout<<"d:="<<d<<endl;
    cout<<"lx: = "<<lx<<endl;
    cout<<"ly: = "<<ly<<endl;
    cout<<"J: = "<<J<<endl;
    cout<<"kx: = "<<kx<<endl;
    cout<<"C: ="<<sector.C<<endl;
    cout<<"nHilbert: ="<<sector.nbasis<<endl;
    //cout<<"-----------Ground state---------"<<endl;
    

    //sector.prlong();


  /* 
   // hamil config;
    lhamil lconfig(lambda,seed);
    for(int i=0;i<20;i++){
    //config.set_hamil(sector,lx,ly,nphi,nLL,i*0.1+0.1);
    //config.print_hamil(10);
    //config.diag();
    // cout<<"E_gs:= "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
    //cout<<i*0.1+0.1<<" "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi,nLL,0.1*i+0.1);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<i*0.1+0.1<<" "<<setprecision(6)<<lconfig.ground_state_energy()/nel<<endl;
    }
   */ 

    

    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi,nLL,d);
    cout<<"Stage-1: Hamiltonian matrix initialized !"<<endl;
    //lconfig.print_hamil(4);
    lconfig.coeff_explicit_update();
    cout<<"Stage-2: Lanczos update completed !"<<endl;
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<"Stage-3: Groundstate wavefunction reconstructed !"<<endl;
    //lconfig.print_eigen(10);
    cout<<"E_gs:= "<<setprecision(6)<<lconfig.ground_state_energy()/nel<<endl;
    cout<<"Ec:= "<<lconfig.Ec<<endl;
    cout<<"Ec_d:= "<<lconfig.Ec_d<<endl;

    for(int i=0;i<lconfig.nHilbert;i++)
        if(abs(lconfig.psir_0[i])>0.11){
          cout<<i<<" :   |";
          unsigned long long u=sector.id[i];
          for(int n=0;n<nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1<<" ";
          cout<<")|";
          for(int n=nphi;n<2*nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1-nphi<<" ";
          cout<<")   ";
          cout<<bitset<8>((sector.id[i])).to_string()<<": "<<bitset<8>((sector.id[i])>>nphi).to_string()<<"  "<<abs(lconfig.psir_0[i])<<endl;
        }
   

  /* 
        basis sector0(nphi,nel_up,nel_down,0,0);
        double t0,s0;
        t0=s0=0;

        int N=sector0.common_divisor(nphi,nel_up);
        hamil config;
        for(int t=0; t<N; t++)
            for(int s=0; s<N; s++){
                basis sector(nphi,nel_up,nel_down,t,s);
                sector.init();
                config.set_hamil(sector,lx,ly,nphi,nLL,d);
                config.diag();

                double K=sqrt((s-s0)*(s-s0)+(t-t0)*(t-t0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
                cout<<K<<" ";
                for(int n=0; n<10; n++)
                    cout<<config.eigenvalues[n]<<" ";
                cout<<endl;

             }
   */ 


    /*
     cout<<K<<" "<<lconfig.eigenvalues[0]<<" ";
     int i=0;
     int count=0;
     do{
        i++;
       if(lconfig.eigenvalues[i]-lconfig.eigenvalues[i-1]>1e-5) {
          cout<<lconfig.eigenvalues[i]<<" ";
          count++;
       }
     }while(count<3 && count<lconfig.norm.size());
     cout<<endl;
     */
//     }


    /*

        lhamil lconfig(lambda,seed);
        //basis sector(nphi,nel_up,nel_down,J_up,J_down);
        //sector.init();
        double Ec_gap=0;
        for(int n=0;n<=20;n++){
          basis sector_p(nphi,nel_up+1,nel_down+1,J_up,J_down);
          sector_p.init();
          lconfig.set_hamil(sector_p,lx,ly,nphi,n*0.15+0.01);
          lconfig.coeff_explicit_update();
          lconfig.diag();
          lconfig.eigenstates_reconstruction();
          Ec_gap=lconfig.ground_state_energy();

          basis sector_h(nphi,nel_up-1,nel_down-1,J_up,J_down);
          sector_h.init();
          lconfig.set_hamil(sector_h,lx,ly,nphi,n*0.15+0.01);
          lconfig.coeff_explicit_update();
          lconfig.diag();
          lconfig.eigenstates_reconstruction();
          Ec_gap+=lconfig.ground_state_energy();

          basis sector(nphi,nel_up,nel_down,J_up,J_down);
          sector.init();
          lconfig.set_hamil(sector,lx,ly,nphi,n*0.15+0.01);
          lconfig.coeff_explicit_update();
          lconfig.diag();
          lconfig.eigenstates_reconstruction();
          Ec_gap-=2*lconfig.ground_state_energy();
          cout<<n*0.15+0.01<<" "<<Ec_gap<<endl;
          //cout<<n*0.1+0.01<<" "<<E_gap<<endl;
        }

    */




    return 0;
}
