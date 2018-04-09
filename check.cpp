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
    int nphi,nel,nel_up,nel_down,J_up, J_down,lambda;
    double lx,ly,gamma;
    int Sz;
    unsigned seed;
    double d,mu;

    nphi=2;
    nel=2;
    nel_up=1;
    nel_down=1;
    gamma=1;
    J_up=-1;
    J_down=-1;
    d=1;
    lambda=200;


    init_argv(nphi,nel,nel_up,J_up,J_down,d,gamma,lambda,argc,argv);
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
    basis sector(nphi,nel_up,nel_down,J_up,J_down);
    sector.init();
    //sector.prlong();
    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    /* basis check */
    /*
    sector.prlong();
    hamil config;
    config.set_hamil(sector,lx,ly,nphi);
    config.print_hamil(10);
    config.diag();
    config.print_eigen(10);
    cout<<"E_gs:= "<<config.ground_state_energy()/nel<<endl;
    */



    /*
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
          cout<<bitset<12>(sector.id_up[i]).to_string()<<":"<<bitset<12>(sector.id_down[j]).to_string()<<lconfig.psir_0[i*sector.nbasis_down+j]<<endl;
        }
        */


    lhamil lconfig(lambda,seed);
    for(int n=0;n<=30;n++){
      lconfig.set_hamil(sector,lx,ly,nphi,n*0.1+0.01);
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
      cout<<n*0.1+0.01<<" "<<lconfig.ground_state_energy()/nel<<endl;
      /*
      double E0=lconfig.eigenvalues[0];
      double delta_E;
      for(int i=0;i<lambda;i++)
         if(lconfig.eigenvalues[i]-E0>0.001){
            delta_E=lconfig.eigenvalues[i]-E0;
            break;
          }
      cout<<n*0.1+0.01<<" "<<delta_E<<endl;
      */
       //mem_actual=(lconfig.H.inner_indices.size()+lconfig.H.value.size()+lconfig.psir_0.size()+lconfig.Coulomb_matrix.size()*4)/1.0e9;
      //cout<<n*0.1+0.05<<" "<<mem_actual<<" "<<mem_theory<<endl;
    }




/*
    cout<<"# nel    GB"<<endl;
    for(int n=8;n<=20;n++){
       long nHilbert=pow(sector.factorial(n,n/2),2)/n;
       double memsize;
       memsize=(2*nHilbert*n+5*nHilbert)*8.0/1.0e9;
       //memsize=nHilbert*nHilbert*8.0/1.0e9;
       cout<<n<<"  "<<memsize<<endl;
    }
    */


    return 0;
}
