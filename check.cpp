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
    int nLL,nphi,nel,nel_up,nel_down,J_up, J_down, kx_up, kx_down,lambda;
    double lx,ly,gamma;
    int Sz;
    unsigned seed;
    double d,mu;

    nLL=0;
    nphi=4;
    nel=4;
    nel_up=2;
    nel_down=2;
    gamma=1.25;
    J_up=0;
    J_down=0;
    kx_up=0;
    kx_down=0;

    d=10;
    lambda=200;
    init_argv(nLL,nphi,nel,nel_up,J_up,J_down,kx_up,kx_down,d,gamma,lambda,argc,argv);
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
    /*
        basis sector(nphi,nel_up,nel_down,J_up,J_down,kx_up,kx_down);
        sector.init();

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

    double j0,k0;
    //if(nel_up%2==0)
    //   j0=k0=nel_up/2;
    //else
    j0=k0=0;

    for(int J=0; J<nphi; J++)
        for(int kx=0; kx<nel_up; kx++){
            double Egs=0;
            int j_up_0=0;
            int kx_up_0=0;
            for(int j_up=0; j_up<nphi; j_up++) {
                for(int kx_up=0; kx_up<nel_up; kx_up++) {
                    int j_down=(J-j_up<0?J-j_up+nphi:J-j_up)%nphi;
                    int kx_down=(kx-kx_up<0?kx-kx_up+nel_up:kx-kx_up)%nel_up;
                    basis sector(nphi,nel_up,nel_down,j_up,j_down,kx_up,kx_down);
                    sector.init();
                    /*
                    lhamil lconfig(lambda,seed);
                    lconfig.set_hamil(sector,lx,ly,nphi,d);
                    lconfig.coeff_explicit_update();
                    lconfig.diag();
                    */
                    hamil config;
                    config.set_hamil(sector,lx,ly,nphi,nLL,d);
                    config.diag();
                    if(config.ground_state_energy()<Egs){
                       Egs=config.ground_state_energy();
                       j_up_0=j_up;
                       kx_up_0=kx_up;
                    }
                }
                int j_down=(J-j_up_0<0?J-j_up_0+nphi:J-j_up_0)%nphi;
                int kx_down=(kx-kx_up_0<0?kx-kx_up_0+nel_up:kx-kx_up_0)%nel_up;
                basis sector(nphi,nel_up,nel_down,j_up_0,j_down,kx_up_0,kx_down);
                sector.init();
                    /*
                    lhamil lconfig(lambda,seed);
                    lconfig.set_hamil(sector,lx,ly,nphi,d);
                    lconfig.coeff_explicit_update();
                    lconfig.diag();
                    */
                hamil config;
                config.set_hamil(sector,lx,ly,nphi,nLL,d);
                config.diag();
                double K=sqrt((J%sector.C_up-j0)*(J%sector.C_up-j0)+(kx%sector.C_up-k0)*(kx%sector.C_up-k0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
                cout<<K<<" ";
                for(int n=0; n<10; n++)
                    cout<<config.eigenvalues[n]<<" ";
                cout<<endl;
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
            }
     }
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
