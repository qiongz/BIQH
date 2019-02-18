#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<ctime>
#include<sstream>
#include<cstdlib>

using namespace std;

int main(int argc,char *argv[]) {
    int nLL,nphi,nel,nel_up,nel_su,J,kx,lambda,nthread;
    double lx,ly,gamma,Delta_SAS,Delta_V,Delta_Z,theta_B;
    unsigned long seed;
    double d;

    init_argv(nLL,nphi,nel,nel_up,nel_su,J,kx,d,Delta_SAS,Delta_V,Delta_Z,gamma,lambda,theta_B,nthread,seed,argc,argv);
    ly=sqrt(nphi*2.0*M_PI/gamma);
    lx=ly*gamma;


    hamil config;
    config.sector.init(nphi,nel,nel_up,nel_su,J,kx);
    cout<<"----------- ED results --------------"<<endl;
    cout<<"nHilbert: ="<<config.sector.nbasis<<endl;
    config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,Delta_Z,nthread);
    config.diag();
    cout<<"E_gs:= "<<setprecision(10)<<config.ground_state_energy()/nel<<endl;

    cout<<"----------- Lanczos results ---------"<<endl;
    lhamil lconfig(lambda,seed);
    lconfig.sector.init(nphi,nel,nel_up,nel_su,J,kx);
    lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,Delta_Z,0,0,0,nthread);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    double Egs=lconfig.ground_state_energy()/nel;
    cout<<"E_gs:= "<<setprecision(10)<<Egs<<endl;

    /*
    cout<<"eigenvalues:="<<lconfig.eigenvalues[0]<<" ";
    int i=0;
    int count=0;
     do{
        i++;
         if(lconfig.eigenvalues[i]-lconfig.eigenvalues[i-1]>1e-4) {
              cout<<lconfig.eigenvalues[i]<<" ";
              count++;
           }
     }while(count<8 && count<lconfig.norm.size());
    cout<<endl;
     */

    vector< pair<long,double> > wf_gs;
    for(long i=0; i<lconfig.nHilbert; i++)
        wf_gs.push_back(pair<long,double>(i,abs(lconfig.psir_0[i])));

    sort(wf_gs.begin(),wf_gs.end(),[](pair<long, double>& elem1,pair<long, double> &elem2) {
        return elem1.second>elem2.second;
    });

    int count=0;
    // print wave function in descending order

    for(auto const & p: wf_gs)
    {
        count++;
        cout<<setw(10)<<"("<<p.first<<"): "<<bitset<32>((lconfig.sector.id[p.first])).to_string()<<"  "<<p.second<<endl;
        cout<<setw(10)<<"upper layer:   |";
        unsigned long long u=lconfig.sector.id[p.first];
        for(int n=0; n<nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")|";
        for(int n=nphi; n<2*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n-nphi;
            else
                cout<<setw(3)<<"_";
        cout<<")|"<<endl;
        cout<<setw(10)<<" down layer:   |";
        for(int n=2*nphi; n<3*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")|";
        for(int n=3*nphi; n<4*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")   "<<endl;
        //cout<<bitset<32>((lconfig.sector.id[p.first])).to_string()<<"  "<<p.second<<endl;
        if(count>lconfig.nHilbert || p.second<1e-5 || count > 10)
            break;
    }

    /*
    cout<<"##########################         First Excited Wave Function       ############################"<<endl;

    vector< pair<long,double> > wf_e1;
    for(long i=0; i<lconfig.nHilbert; i++)
        wf_e1.push_back(pair<long,double>(i,abs(lconfig.psir_1[i])));

    sort(wf_e1.begin(),wf_e1.end(),[](pair<long, double>& elem1,pair<long, double> &elem2) {
        return elem1.second>elem2.second;
    });

    // print wave function in descending order
    count=0;

    for(auto const & p: wf_e1)
    {
        count++;
        cout<<setw(10)<<"("<<p.first<<"): "<<bitset<32>((lconfig.sector.id[p.first])).to_string()<<"  "<<p.second<<endl;
        cout<<setw(10)<<"upper layer:   |";
        unsigned long long u=lconfig.sector.id[p.first];
        for(int n=0; n<nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")|";
        for(int n=nphi; n<2*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n-nphi;
            else
                cout<<setw(3)<<"_";
        cout<<")|"<<endl;
        cout<<setw(10)<<" down layer:   |";
        for(int n=2*nphi; n<3*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")|";
        for(int n=3*nphi; n<4*nphi; n++)
            if((u>>n)%2==1)
                cout<<setw(3)<<n;
            else
                cout<<setw(3)<<"_";
        cout<<")   "<<endl;
        //cout<<bitset<32>((lconfig.sector.id[p.first])).to_string()<<"  "<<p.second<<endl;
        if(count>lconfig.nHilbert || p.second<1e-5)
            break;
    }
    */







    // calculate the excitation spectrum
    /*
             double t0,s0;
             t0=s0=nphi/2;

        int N=nphi/2;
             //lhamil config(lambda,seed);
             for(int t=0; t<=N; t++)
                 for(int s=0; s<=N; s++){
     	        //hamil config;
     	        lhamil config(lambda,seed);
                 config.sector.init(nphi,nel,nel_up,nel_su,t,s);
                 config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,Delta_Z,0,0,0,nthread);
     		config.coeff_explicit_update();
                 config.diag();
     		config.eigenstates_reconstruction();

                 double K=sqrt((s-s0)*(s-s0)+(t-t0)*(t-t0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
                 cout<<t<<" "<<s<<" "<<K<<" ";
     		cout<<config.ground_state_energy()/nel<<" "<<config.Sz()<<" ";
                     //for(int n=0; n<5; n++)
                      //   cout<<config.eigenvalues[n]<<" ";
     		//cout<<config.ground_state_energy()/nel;
                     //cout<<endl;
     		//config.sector.clear();
     	        cout<<config.eigenvalues[0]<<" ";

                     int i=0;
                     int count=0;
                     do{
                      i++;
                      if(config.eigenvalues[i]-config.eigenvalues[i-1]>1e-5) {
                        cout<<config.eigenvalues[i]<<" ";
                        count++;
                         }
                       }while(count<6 && count<config.norm.size());
                      cout<<endl;
    	    config.sector.clear();
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

    /*
        lhamil lconfig(lambda,seed);
        lconfig.sector.init(nphi,nel,J,kx);
        for(int i=0; i<201; i++) {
            Delta_V=i*0.008;
            for(int j=0; j<201; j++) {
                Delta_SAS=j*0.003;
    	    d=j*0.05;
                //theta_B=j/100.0*M_PI;
                lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,Delta_Z,theta_B,0,0,nthread);
                lconfig.coeff_explicit_update();
                lconfig.diag();
                lconfig.eigenstates_reconstruction();
    	    double Egs=lconfig.ground_state_energy()/nel;
                double Px=lconfig.pseudospin_Sx();
                double Pz=lconfig.pseudospin_Sz();
    	    double Sz=lconfig.Sz();
    	    double Sf=lconfig.spinflip_tunneling();
    */
    /*
    lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS+0.0001,Delta_V,Delta_Z,nthread);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    double Chi_sf=(lconfig.spinflip_tunneling()-Sf)/0.0001;
    */
    //cout<<d<<" "<<Egs<<" "<<Sz<<" "<<Px<<" "<<Pz<<" "<<Sf<<endl;
    //cout<<Delta_SAS<<" "<<Egs<<" "<<Sz<<" "<<Px<<" "<<Pz<<" "<<Sf<<endl;
    /*          cout<<Delta_SAS<<" "<<Delta_V<<" "<<Egs<<" "<<Sz<<" "<<Px<<" "<<Pz<<" "<<Sf<<endl;
         }
         cout<<endl;
       }

    */
    //Delta_SAS/=pow(nel,1.5);

    /*
    lhamil lconfig(lambda,seed);
    lconfig.sector.init(nphi,nel,nel_up,nel_su,J,kx);
    for(int j=0;j<21;j++){
    //Delta_V=j*1e-4;
    d=j*0.1;
    //Delta_SAS=0.002*j+1e-6;
    //theta_y=j*0.25/10*M_PI;
    lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,Delta_Z,0,0,0,nthread);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    cout<<d<<setprecision(6)<<" "<<lconfig.ground_state_energy()<<" "<<lconfig.Sz()<<endl;
    }
    */




    return 0;
}
