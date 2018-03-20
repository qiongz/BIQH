#include"hamiltonian.h"
hamil::hamil() {}

hamil::hamil(basis &sector,double _t, double _U) {
    long nsite,nbasis_up,nbasis_down,signu,signd;
    t=_t;
    U=_U;
    nsite=sector.nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector< complex<double> > matrix_elements;
    std::map<long,complex<double> >::iterator it;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    long b,mask,p,q,n,m,i,j,k,l,r,s,t,nsignu,nsignd;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            std::map<long,double> col_indices;
                // select two electrons in left-basis <m_1, m_2|
                // consider the upper-layer two electrons
                for(n=0;n<nsite;n++)
                  for(m=0;m<nsite;m++){
                     mask=(1<<n)+(1<<m);
                     // consider the upper-layer two electrons
                     // if there're two electrons on n and m;
                     if(i&mask==mask){
                        // b is the rest electon positions
                        b=i^mask;
                        long nt,mt,it,mask_t,occ_t;
                        nsignu=0;
                        // perform translation
                        for(t=0;t<nsite;t++){
                          // PBC, if one electron cross left boundary, sign change with -1
                          if(n-t<0){
                            nt=n-t+nsite;
                            nsignu++;
                          }
                          // PBC, if one electron cross right boundary, sign change with -1
                          if(m+t>=nsite){
                            mt=m+t-nsite;
                            nsignu++;
                          }
                          // the translated two electrons indices
                          mask_t=(1<<nt)+(1<<mt);
                          // occupation of electons on the translated position
                          occ_t=mask_t&b;
                          // if there're no electon on the translated position
                          // which is a valid translation, can be applied
                          if(occ_t==0){
                              // the translated indices
                              k=mask_t+b;
                              // calculate the Coulomb matrix contribution and
                              // add it to the hamiltonian matrix
                              if(k!=i) {
                                  // Coulomb matrix element
                                  long q_x,q_y;
                                  q_y=nt;
                                  complex<double> V_uu=0;
                                  for(q_x=0;q_x<L;q_x++)
                                     if(q_y!=0&&q_x!=0)
                                         V_uu+=Coulomb_interaction(0,0,q_x,q_y)*complex<double>(cos((n-m+t)*q_x*2.0*M_PI/N_phi),sin((n-m+t)*q_x*2.0*M_PI/N_phi)))/(4.0*M_PI*N_phi);
                                  it=col_indices.find(k*nbasis_down+j);
                                  if(it==col_indices.end())
                                    col_indices.insert(std::pair<long,double>(k*nbasis_down+j,V_uu*pow(-1,nsignu)));
                                  else
                                    it->second+=V_uu*pow(-1,nsignu);
                               }
                          // two electrons are occupied, and to be crossed next
                          else if(occ_t==mask_t)
                            signu+=2;
                          // one electron is occupied, and to be crossed next
                          else if(occ_t!=0 && occ_t!=mask_t)
                            signu++;
                           }
                         }
                      }

                     // consider the lower-layer two electrons
                     // if there're two electrons on n and m;
                     if(j&mask==mask){
                      {
                        // b is the rest electon positions
                        b=j^mask;
                        long nt,mt,jt,mask_t,occ_t;
                        nsignu=0;
                        // perform translation
                        for(t=0;t<nsite;t++){
                          // PBC, if one electron cross left boundary, sign change with -1
                          if(n-t<0){
                            nt=n-t+nsite;
                            nsignu++;
                          }
                          // PBC, if one electron cross right boundary, sign change with -1
                          if(m+t>=nsite){
                            mt=m+t-nsite;
                            nsignu++;
                          }
                          // the translated two electrons indices
                          mask_t=(1<<nt)+(1<<mt);
                          // occupation of electons on the translated position
                          occ_t=mask_t&b;
                          // if there're no electon on the translated position
                          // which is a valid translation, can be applied
                          if(occ_t==0){
                              // the translated indices
                              k=mask_t+b;
                              // calculate the Coulomb matrix contribution and
                              // add it to the hamiltonian matrix
                              if(k!=j) {
                                  // Coulomb matrix element
                                  long q_x,q_y;
                                  q_y=nt;
                                  complex<double> V_uu=0;
                                  for(q_x=0;q_x<L;q_x++)
                                     if(q_y!=0&&q_x!=0)
                                         V_uu+=Coulomb_interaction(1,1,q_x,q_y)*complex<double>(cos((n-m+t)*q_x*2.0*M_PI/N_phi),sin((n-m+t)*q_x*2.0*M_PI/N_phi)))/(4.0*M_PI*N_phi);
                                  it=col_indices.find(i*nbasis_down+k);
                                  if(it==col_indices.end())
                                    col_indices.insert(std::pair<long,double>(i*nbasis_down+k,V_uu*pow(-1,nsignu)));
                                  else
                                    it->second+=V_uu*pow(-1,nsignu);
                               }
                          // two electrons are occupied, and to be crossed next
                          else if(occ_t==mask_t)
                            signu+=2;
                          // one electron is occupied, and to be crossed next
                          else if(occ_t!=0 && occ_t!=mask_t)
                            signu++;
                           }
                         }
                      }
                  }

            for(it=col_indices.begin(); it!=col_indices.end(); it++) {
                inner_indices.push_back(it->first);
                matrix_elements.push_back(it->second);
            }
            row+=col_indices.size();
            outer_starts.push_back(row);
            col_indices.clear();
          }
    }
    H.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}

hamil::~hamil() {}

void hamil::init(basis &sector,double _t,double _U){
    long nsite,nbasis_up,nbasis_down,signu,signd;
    t=_t;
    U=_U;
    nsite=sector.nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector<double> matrix_elements;
    std::map<long,double>::iterator it;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    long n,m,i,j,k,l,s,nsignu,nsignd;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            std::map<long,double> col_indices;
            for(n=0; n<nsite; n++) {
                if(n+1>=nsite) {
                    signu=pow(-1,sector.nel_up-1);
                    signd=pow(-1,sector.nel_down-1);
                }
                else {
                    signu=1;
                    signd=1;
                }
                m=n+1;
                k=sector.hopping_up(i,n,m);
                if(k!=i) {
                    it=col_indices.find(k*nbasis_down+j);
                    if(it==col_indices.end())
                        col_indices.insert(std::pair<long,double>(k*nbasis_down+j,-t*signu));
                    else
                        it->second+=-t*signu;
                }
                if(sector.potential(i,j,n)) {
                    it=col_indices.find(i*nbasis_down+j);
                    if(it==col_indices.end())
                        col_indices.insert(std::pair<long,double>(i*nbasis_down+j,U));
                    else
                        it->second+=U;
                }
            }
            for(it=col_indices.begin(); it!=col_indices.end(); it++) {
                inner_indices.push_back(it->first);
                matrix_elements.push_back(it->second);
            }

            row+=col_indices.size();
            outer_starts.push_back(row);
            col_indices.clear();
        }
    }
    H.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this !=&_gs_hconfig) {
        seed=_gs_hconfig.seed;
        nHilbert=_gs_hconfig.nHilbert;
        H=_gs_hconfig.H;
        t=_gs_hconfig.t;
        U=_gs_hconfig.U;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(),_gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(),_gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(),_gs_hconfig.psi_n0.end());
    }
    return *this;
}

double hamil::Coulomb_interaction(int alpha, int beta, int q_x, int q_y){
  double q=sqrt(q_x*q_x+q_y*q_y)*2*M_PI/L;
  if(alpha==beta)
      return 2.0*M_PI/(q+1e-8)*exp(-q*q/2);
  else
      return 2.0*M_PI/(q+1e-8)*exp(-q*q/2-q*d);

}

double hamil::spectral_function(vector<double> &O_phi_0,double omega,double _E0, double eta, int annil) {
    complex<double> E;
    complex<double> G=0;
    for(int i=0; i<nHilbert; i++)
        // set annil==1, which gives hole-sector
        if(annil==1){
            E=complex<double>(omega,eta);
            G+=pow(psi_n0[i]*O_phi_0[i],2)/(E+eigenvalues[i]-_E0);
          }
        // else particle-sector
        else{
            E=complex<double>(omega,eta);
            G+=pow(psi_n0[i]*O_phi_0[i],2)/(E+_E0-eigenvalues[i]);
          }

    return -G.imag()/M_PI;
}

double hamil::ground_state_energy() {
    if(psi_0.size()==0) return 0;
    double E_gs=0;
    vector<double> psi_t;
    psi_t=H*psi_0;
    for(int i=0; i<nHilbert; i++)
        E_gs+=psi_t[i]*psi_0[i];
    return E_gs;
}

void hamil::diag() {
    int i,idx;
    double *hamiltonian=new double[nHilbert*nHilbert];
    double *en=new double[nHilbert];
    memset(hamiltonian,0,sizeof(double)*nHilbert*nHilbert);
    for(i=0; i<H.outer_starts.size()-1; i++)
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            hamiltonian[i*nHilbert+H.inner_indices[idx]]=H.value[idx];
    diag_dsyev(hamiltonian,en,nHilbert);
    psi_0.assign(nHilbert,0);
    psi_n0.assign(nHilbert,0);
    eigenvalues.assign(nHilbert,0);
    for(i=0; i<nHilbert; i++) {
        eigenvalues[i]=en[i];
        psi_0[i]=hamiltonian[i];
        psi_n0[i]=hamiltonian[i*nHilbert];
    }
    delete hamiltonian,en;
}


void hamil::print_hamil() {
    std::cout<<"hamiltonian in CSR format: "<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    H.print();
}

void hamil::print_eigen(){
  std::cout<<"Eigenvalues:=[ ";
   for(int i=0;i<nHilbert;i++)
      if(i!=nHilbert-1)
        std::cout<<eigenvalues[i]<<", ";
      else
         std::cout<<eigenvalues[i]<<" ]"<<std::endl;
}
