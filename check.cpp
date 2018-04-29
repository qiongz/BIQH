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
    int nphi,nel,nel_up,nel_down,J_up, J_down, kx_up, kx_down,lambda;
    double lx,ly,gamma;
    int Sz;
    unsigned seed;
    double d,mu;

    nphi=4;
    nel=4;
    nel_up=2;
    nel_down=2;
    gamma=1.0;
    J_up=0;
    J_down=0;
    kx_up=0;
    kx_down=0;

    d=10;
    lambda=200;
    init_argv(nphi,nel,nel_up,J_up,J_down,kx_up,kx_down,d,gamma,lambda,argc,argv);
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

    basis sector(nphi,nel_up,nel_down,J_up,J_down,kx_up,kx_down);
    sector.init();
    /*
    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;
    */

    cout<<"nphi: = "<<nphi<<endl;
    cout<<"nel_up: = "<<nel_up<<endl;
    cout<<"nel_down: = "<<nel_up<<endl;
    cout<<"d:="<<d<<endl;
    cout<<"lx: = "<<lx<<endl;
    cout<<"ly: = "<<ly<<endl;
    cout<<"J_u: = "<<J_up<<endl;
    cout<<"K_u: = "<<kx_up<<endl;
    cout<<"J_down: = "<<J_down<<endl;
    cout<<"K_down: = "<<kx_down<<endl;
    cout<<"nbasis_up:="<<sector.nbasis_up<<endl;
    cout<<"nbasis_down:="<<sector.nbasis_down<<endl;
    cout<<"nHilbert: ="<<sector.nbasis_up*sector.nbasis_down<<endl;
    cout<<"-----------Ground state---------"<<endl;
    //sector.prlong();
    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi,d);
    //config.print_hamil(10);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    lconfig.print_eigen(10);
    cout<<"E_gs:= "<<lconfig.ground_state_energy()/nel<<endl;

/*
    basis sector(nphi,nel_up,nel_down,J_up,J_down);
    sector.init();
    lhamil lconfig(lambda,seed);
    lconfig.set_hamil(sector,lx,ly,nphi,d);
    cout<<"-----------Lanczos---------"<<endl;
    //lconfig.print_hamil_CSR();
    //lconfig.print_hamil(10);
    lconfig.coeff_explicit_update();
    //lconfig.print_lhamil(10);
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<"eigenvalues:= [";
    for(int i=0;i<10;i++)
       cout<<lconfig.eigenvalues[i]<<",  ";
    cout<<", ...]"<<endl;
    cout<<"E_gs:="<<lconfig.ground_state_energy()/nel<<endl;
    */



    /*
    for(int i=0;i<sector.nbasis_up;i++)
     for(int j=0;j<sector.nbasis_down;j++)
        if(abs(lconfig.psir_0[i*sector.nbasis_down+j])>0.11){
          cout<<i*sector.nbasis_down+j<<" :   ";
          long u=sector.id_up[i];
          long d=sector.id_down[j];
          for(int n=0;n<nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1<<" ";
          cout<<":";
          for(int n=0;n<nphi;n++)
             if((d>>n)%2==1)
                cout<<n+1<<" ";
          cout<<"    ";
          cout<<bitset<12>(sector.id_up[i]).to_string()<<":"<<bitset<12>(sector.id_down[j]).to_string()<<" "<<lconfig.psir_0[i*sector.nbasis_down+j]<<endl;
        }
        */

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
