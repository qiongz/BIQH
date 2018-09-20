#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<ctime>
#include<sstream>
#include<cstdlib>

using namespace std;

int main(int argc,char *argv[]) {
    int nLL,nphi,nel,nel_up,nel_down,J,kx,lambda,nthread;
    double lx,ly,gamma,Delta_SAS,Delta_V;
    int Sz;
    unsigned long seed;
    double d,mu;

    init_argv(nLL,nphi,nel,nel_up,J,kx,d,Delta_SAS,Delta_V,gamma,lambda,nthread,seed,argc,argv);

    ly=sqrt(nphi*2.0*M_PI/gamma);
    lx=ly*gamma;

    nel_down=nel-nel_up;

    /*

        lhamil lconfig(lambda,seed);
        lconfig.sector.init(nphi,nel,nel_up,J,kx);
        for(int n=0;n<20;n++){
        d=n*0.1+2.5;
        for(int i=0;i<20;i++){
        Delta_V=i*0.01;
        for(int j=0;j<5;j++){
        Delta_SAS=j*0.0001;
        //d=i*0.1;
        lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,nthread);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        double Et=lconfig.ground_state_energy()/nel;
        double Sx=lconfig.pseudospin_Sx();
        double Sz=lconfig.pseudospin_Sz();

        lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS+0.00001,Delta_V,nthread);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        double Chi=(lconfig.pseudospin_Sx()-Sx)/(0.00001*nel);

        //cout<<d<<setprecision(6)<<" "<<lconfig.ground_state_energy()/nel<<" "<<lconfig.eigenvalues[0]<<endl;
        //cout<<d<<setprecision(6)<<" "<<lconfig.ground_state_energy()/nel<<" "<<lconfig.density_imbalance()<<" "<<lconfig.pseudospin_Sx()/nel<<endl;
        cout<<d<<" "<<Delta_V<<" "<<Delta_SAS<<" "<<Sx<<" "<<Chi<<" "<<Sz<<endl;
        //cout<<d<<" "<<lconfig.ground_state_energy()/nel<<" "<<lconfig.density_imbalance()<<endl;
        }
        }
        }
    */
    /*
        Delta_SAS/=pow(nel,1.5);

        lhamil lconfig(lambda,seed);
        lconfig.sector.init(nphi,nel,nel_up,J,kx);
        for(int j=0;j<50;j++){
        //Delta_V=j*1e-4;
        d=j*0.004+0.7;
        //Delta_SAS=0.002*j+1e-6;
        lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,nthread);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        double Sz=lconfig.pseudospin_Sz();
        double Et=lconfig.ground_state_energy();
        double Sx=lconfig.pseudospin_Sx();

        lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS+0.0001,Delta_V,nthread);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        double E0=lconfig.ground_state_energy();
        //double Sx=-(Et-E0)/Delta_SAS/nel;
        double Chi=(lconfig.pseudospin_Sx()-Sx)/0.0001/nel;



        //cout<<d<<setprecision(6)<<" "<<-(E0-Et)/0.0001/nel<<" "<<Sx<<" "<<Chi<<endl;
        cout<<d<<setprecision(6)<<" "<<Sx<<" "<<Chi<<endl;
        //cout<<Delta_V<<" "<<Sx<<" "<<Sz<<endl;
        }

    */










    
    hamil config;
    config.sector.init(nphi,nel,nel_up,J,kx);
    cout<<"nHilbert: ="<<config.sector.nbasis<<endl;
    cout<<"-----------Ground state---------"<<endl;

    auto _t1=std::chrono::high_resolution_clock::now();
    config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,nthread);
    auto _t2=std::chrono::high_resolution_clock::now();
    cout<<"Full Hamiltonian matrix initialized !"<<endl;
    cout<<"time cost: "<<chrono::duration_cast<chrono::microseconds>(_t2-_t1).count()/1.0e6<<" seconds."<<endl;

    //config.print_hamil(10);
    auto _t3=std::chrono::high_resolution_clock::now();
    config.diag();

    auto _t4=std::chrono::high_resolution_clock::now();
    cout<<"time cost: "<<chrono::duration_cast<chrono::microseconds>(_t4-_t3).count()/1.0e6<<" seconds."<<endl;
    cout<<"E_gs:= "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
    cout<<"E0:="<<config.eigenvalues[0]<<endl;
    cout<<"E1:="<<config.eigenvalues[1]<<endl;
    // cout<<"Pz:="<<config.pseudospin_Sz()<<endl;
    //    cout<<"Px:="<<config.pseudospin_Sx()<<endl;
    
    /*
    cout<<"# ground state wave function"<<endl;
    for(int i=0;i<config.nHilbert;i++)
        if(abs(config.psi_0[i])>0.05){
          cout<<i<<" :   |";
          unsigned long long u=config.sector.id[i];
          for(int n=0;n<nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1<<" ";
          cout<<")|";
          for(int n=nphi;n<2*nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1-nphi<<" ";
          cout<<")   ";
          cout<<bitset<8>((config.sector.id[i])).to_string()<<": "<<bitset<8>((config.sector.id[i])>>nphi).to_string()<<"  "<<abs(config.psi_0[i])<<endl;
    }
    */


    /*
    cout<<"# first excited state wave function"<<endl;
    for(int i=0;i<config.nHilbert;i++)
        if(abs(config.psi_1[i])>0.01){
          cout<<i<<" :   |";
          unsigned long long u=config.sector.id[i];
          for(int n=0;n<nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1<<" ";
          cout<<")|";
          for(int n=nphi;n<2*nphi;n++)
             if((u>>n)%2==1)
                cout<<n+1-nphi<<" ";
          cout<<")   ";
          cout<<bitset<8>((config.sector.id[i])).to_string()<<": "<<bitset<8>((config.sector.id[i])>>nphi).to_string()<<"  "<<abs(config.psi_1[i])<<endl;
    }
    */


    /*
       lhamil lconfig(lambda,seed);
       lconfig.sector.init(nphi,nel_up,nel_down,J,kx);
       //lconfig.sector.prlong();
       unsigned long basis_i,basis_j;
       int sign;
       basis_i=lconfig.sector.id[120];
       basis_j=lconfig.sector.translate(basis_i,kx,sign);

       cout<<"basis i:  ";
            unsigned long c=basis_i;
            cout<<" | ";
            for(int n=0; n<nphi; n++) {
                if((c>>n)%2==1)
                    cout<<n<<" ";
            }
            cout<<">| ";
            for(int n=nphi; n<2*nphi; n++) {
                if((c>>n)%2==1)
                    cout<<n-nphi<<" ";
            }
            cout<<"> ";
       cout<<"  "<<bitset<20>(basis_i).to_string()<<endl;
       cout<<"basis j:   ";
            c=basis_j;
            cout<<" | ";
            for(int n=0; n<nphi; n++) {
                if((c>>n)%2==1)
                    cout<<n<<" ";
            }
            cout<<">| ";
            for(int n=nphi; n<2*nphi; n++) {
                if((c>>n)%2==1)
                    cout<<n-nphi<<" ";
            }
            cout<<"> ";
       cout<<bitset<20>(basis_j).to_string()<<endl;
       cout<<"sign: "<<sign<<endl;
    */



    cout<<"-----------Lanczos results---------"<<endl;

    lhamil lconfig(lambda,seed);
    auto t1=std::chrono::high_resolution_clock::now();
    lconfig.sector.init(nphi,nel,nel_up,J,kx);
    //lconfig.sector.prlong();
    auto t2=std::chrono::high_resolution_clock::now();
    cout<<"Stage-1: Basis initialized !"<<endl;
    cout<<"nHilbert: ="<<lconfig.sector.nbasis<<endl;
    cout<<"time cost: "<<chrono::duration_cast<chrono::microseconds>(t2-t1).count()/1.0e6<<" seconds."<<endl;

    auto t3=std::chrono::high_resolution_clock::now();
    lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,nthread);
    auto t4=std::chrono::high_resolution_clock::now();
    cout<<"Stage-2: Hamiltonian matrix initialized !"<<endl;
    cout<<"time cost: "<<chrono::duration_cast<chrono::microseconds>(t4-t3).count()/1.0e6<<" seconds."<<endl;
    //lconfig.print_hamil(lconfig.nHilbert);
    //lconfig.print_hamil_CSR();
    t3=std::chrono::high_resolution_clock::now();
    lconfig.coeff_explicit_update();
     cout<<"Stage-2: Lanczos update completed !"<<endl;
    lconfig.diag();
    lconfig.eigenstates_reconstruction();


    double Egs=lconfig.ground_state_energy();
    t4=std::chrono::high_resolution_clock::now();
    cout<<"Stage-3: Groundstate wavefunction reconstructed !"<<endl;
    cout<<"E_gs:= "<<setprecision(10)<<Egs/nel<<endl;
    //cout<<"E0:="<<setprecision(6)<<lconfig.eigenvalues[0]<<endl;
    cout<<"time cost: "<<chrono::duration_cast<chrono::microseconds>(t4-t3).count()/1.0e6<<" seconds."<<endl;
    /*
    cout<<"Pz:="<<lconfig.pseudospin_Sz()<<endl;
    cout<<"Px:="<<lconfig.pseudospin_Sx()<<endl;

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
        cout<<setw(4)<<p.first<<" :   |";
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
        cout<<")   ";
        cout<<bitset<24>((lconfig.sector.id[p.first])).to_string()<<"  "<<p.second<<endl;
        if(count>10 || p.second<1e-3)
            break;
    }
    */


    /*
            basis sector0(nphi,nel_up,nel_down,0,0);
            double t0,s0;
            t0=s0=nphi/2;

            //int N=sector0.common_divisor(nphi,nel_up);
            //lhamil config(lambda,seed);
            for(int t=0; t<nphi; t++)
                for(int s=0; s<nphi; s++){
    	        //hamil config;
    	        lhamil config(lambda,seed);
    		config.sector.init(nphi,nel_up,nel_down,t,s);
                    config.set_hamil(lx,ly,nphi,nLL,d,nthread);
    		config.coeff_explicit_update();
                    config.diag();
    		config.eigenstates_reconstruction();

                    double K=sqrt((s-s0)*(s-s0)+(t-t0)*(t-t0)*gamma*gamma)*sqrt(2.0*M_PI/nphi/gamma);
                    cout<<t<<" "<<s<<" "<<K<<" ";
    		cout<<config.ground_state_energy()/nel<<" ";
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




    return 0;
}
