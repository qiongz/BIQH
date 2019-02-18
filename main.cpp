#include"init.h"
#include"basis.h"
#include"chern.h"
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

    // test comparing full ED with Lanczos
    hamil config;
    cout<<"----------- ED results --------------"<<endl;
    config.sector.init(nphi,nel,nel_up,J,kx);
    cout<<"nHilbert: ="<<config.sector.nbasis<<endl;
    config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,0,0,nthread);
    config.diag();
    cout<<"E_gs:= "<<setprecision(6)<<config.ground_state_energy()/nel<<endl;
    cout<<"----------- Lanczos results ---------"<<endl;
    lhamil lconfig(lambda,seed);
    lconfig.sector.init(nphi,nel,nel_up,J,kx);
    lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,0,0,nthread);
    lconfig.coeff_explicit_update();
    lconfig.diag();
    lconfig.eigenstates_reconstruction();
    double Egs=lconfig.ground_state_energy();
    cout<<"E_gs:= "<<setprecision(10)<<Egs/nel<<endl;
    vector< pair<long,double> > wf_gs;
    for(long i=0; i<lconfig.nHilbert; i++)
        wf_gs.push_back(pair<long,double>(i,abs(lconfig.psir_0[i])));
    sort(wf_gs.begin(),wf_gs.end(),[](pair<long, double>& elem1,pair<long, double> &elem2) {
        return elem1.second>elem2.second;
    });
    cout<<"# ground state wave function"<<endl;
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




    // calculate the excitation spectrum 
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

    // calculate the distance dependence propertities
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

     //calculate the spin stiffness    
     /*
        lhamil lconfig(lambda,seed);
        lconfig.sector.init(nphi,nel,nel_up,J,kx);
        double theta_x,theta_y;
        for(int i=0;i<40;i++){
        theta_y=i/20.0*M_PI; 
        lconfig.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,theta_x,theta_y,nthread);
        lconfig.coeff_explicit_update();
        lconfig.diag();
        lconfig.eigenstates_reconstruction();
        double E0=lconfig.ground_state_energy()/nel;
        cout<<theta_y<<setprecision(6)<<" "<<E0<<endl;
        }
    */


    // calculate the drag Hall conductance 
    /*
    int theta_1,theta_2,n_mesh=10;
    lhamil config(lambda,seed);
    //hamil config;
    config.sector.init(nphi,nel,nel_up,J,kx);
    long long nHilbert=config.sector.nbasis;
    cerr<<"nHilbert:="<<nHilbert<<endl;
    complex<double> *wfs_full = new complex<double>[(n_mesh+1)*2*nHilbert];
    double *chern_numbers_theta = new double [ n_mesh];
    double chern_number=0;
    double theta_x,theta_y;
    theta_1=0;
    theta_x=theta_1*2.0*M_PI/n_mesh;
    for(theta_2=0;theta_2<n_mesh;theta_2++){
	theta_y=theta_2*2.0*M_PI/n_mesh;
        config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,theta_x,theta_y,nthread);
        config.coeff_explicit_update();
        config.diag();
        config.eigenstates_reconstruction();
	for(long i=0;i<nHilbert;i++)
	   wfs_full[((n_mesh+1)*(theta_1%2)+theta_2)*nHilbert+i]=config.psir_0[i];
	cout<<theta_1*n_mesh+theta_2<<" "<<config.eigenvalues[0];
	for(long i=1;i<nHilbert;i++)
	  if(abs(config.eigenvalues[i]-config.eigenvalues[0])>1e-4*abs(config.eigenvalues[0])){
	    cout<<" "<<config.eigenvalues[i]<<endl;
		break;
	   }
	}
    for(long i=0;i<nHilbert;i++)
	wfs_full[((n_mesh+1)*(theta_1%2)+n_mesh)*nHilbert+i]=wfs_full[((n_mesh+1)*(theta_1%2))*nHilbert+i];

    for(theta_1=1;theta_1<=n_mesh;theta_1++){
	theta_x=theta_1*2.0*M_PI/n_mesh;
    for(theta_2=0;theta_2<n_mesh;theta_2++){
	theta_y=theta_2*2.0*M_PI/n_mesh;
        config.set_hamil(lx,ly,nphi,nLL,d,Delta_SAS,Delta_V,theta_x,theta_y,nthread);
        config.coeff_explicit_update();
        config.diag();
        config.eigenstates_reconstruction();
	cout<<theta_1*n_mesh+theta_2<<" "<<config.eigenvalues[0];
	for(long i=1;i<nHilbert;i++)
	  if(abs(config.eigenvalues[i]-config.eigenvalues[0])>1e-4*abs(config.eigenvalues[0])){
	    cout<<" "<<config.eigenvalues[i]<<endl;
	    break;
	}
	for(long i=0;i<nHilbert;i++)
	    wfs_full[((n_mesh+1)*(theta_1%2)+theta_2)*nHilbert+i]=config.psir_0[i];
	}
        for(long i=0;i<nHilbert;i++)
	   wfs_full[((n_mesh+1)*(theta_1%2)+n_mesh)*nHilbert+i]=wfs_full[((n_mesh+1)*(theta_1%2))*nHilbert+i];
        cal_Chern(wfs_full,chern_numbers_theta,nHilbert , n_mesh, theta_1);
        for(int i=0; i<n_mesh; i++)
           chern_number+=chern_numbers_theta[i];
    }
    
    cerr<<"chern_number:= "<<chern_number<<endl;
    */

    return 0;
}
