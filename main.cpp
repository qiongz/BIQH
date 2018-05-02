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
    J=0;
    kx=0;

    d=10;
    lambda=200;
    init_argv(nLL,nphi,nel,nel_up,J,kx,d,gamma,lambda,argc,argv);
    gamma=nel_up/4.0;
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

      
        basis sector(nphi,nel_up,nel_down,J,2,kx,2);
        sector.init();

        cout<<"nphi: = "<<nphi<<endl;
        cout<<"nel_up: = "<<nel_up<<endl;
        cout<<"nel_down: = "<<nel_up<<endl;
        cout<<"d:="<<d<<endl;
        cout<<"lx: = "<<lx<<endl;
        cout<<"ly: = "<<ly<<endl;
        cout<<"J: = "<<J<<endl;
        cout<<"kx: = "<<kx<<endl;
        cout<<"Ju: ="<<sector.J_up<<endl;
        cout<<"Jd: ="<<sector.J_down<<endl;
        cout<<"Ku: ="<<sector.K_up<<endl;
        cout<<"Kd: ="<<sector.K_down<<endl;
        cout<<"Cu: ="<<sector.C_up<<endl;
        cout<<"Cd: ="<<sector.C_down<<endl;
        cout<<"nbasis_up:="<<sector.nbasis_up<<endl;
        cout<<"nbasis_down:="<<sector.nbasis_down<<endl;
        cout<<"nHilbert: ="<<sector.nbasis_up*sector.nbasis_down<<endl;
        cout<<"-----------Ground state---------"<<endl;
        sector.prlong();
        hamil config;
        config.set_hamil(sector,lx,ly,nphi,nLL,d); 
        config.diag();
        cout<<"E_gs:= "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
        
        

        lhamil lconfig(lambda,seed);
        lconfig.set_hamil(sector,lx,ly,nphi,nLL,d);
        lconfig.print_hamil(4);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        lconfig.print_eigen(10);
        cout<<"E_gs:= "<<setprecision(6)<<lconfig.ground_state_energy()/nel<<endl;
     
    

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
    basis sector0(nphi,nel_up,nel_down,0,0,0,0);
    double t0,s0;
    t0=s0=0;
     
    int N=sector0.common_divisor(nphi,nel_up);

    for(int t=0; t<N; t++)
        for(int s=0; s<N; s++){
            double Egs=0;
            int tu_0=0;
            int su_0=0;
            
            for(int tu=0; tu<N; tu++) 
                for(int su=0; su<N; su++) {
                    int td=(t-tu<0?t-tu+N:t-tu)%N;
                    int sd=(s-su<0?s-su+N:s-su)%N;
                    basis sector(nphi,nel_up,nel_down,tu,td,su,sd);
                    sector.init();
                    //lhamil lconfig(lambda,seed);
                    //lconfig.set_hamil(sector,lx,ly,nphi,d);
                    //lconfig.coeff_explicit_update();
                    //lconfig.diag();
                    hamil config;
                    config.set_hamil(sector,lx,ly,nphi,nLL,d);
                    config.diag();
                    if(config.ground_state_energy()<Egs){
                       Egs=config.ground_state_energy();
                       tu_0=tu;
                       su_0=su;
                    }
                }
                basis sector(nphi,nel_up,nel_down,t,tu_0,s,su_0);
                sector.init();
                    //lhamil lconfig(lambda,seed);
                    //lconfig.set_hamil(sector,lx,ly,nphi,d);
                    //lconfig.coeff_explicit_update();
                //    lconfig.diag();
                hamil config;
                config.set_hamil(sector,lx,ly,nphi,nLL,d);
                config.diag();
                double K=sqrt((s-s0)*(s-s0)+(t-t0)*(t-t0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
                cout<<K<<" ";
                for(int n=0; n<3; n++)
                    cout<<config.eigenvalues[n]<<" ";
                cout<<endl;
              
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
     //}


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
