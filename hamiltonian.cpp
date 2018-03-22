#include"hamiltonian.h"
hamil::hamil() {}

void hamil::set_hamil(basis & sector, double _d) {
    long nbasis_up,nbasis_down;
    d=_d;
    nsite=sector.nsite;
    nphi=nsite*nsite/(2.0*M_PI);
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector< complex<double> > matrix_elements;
    std::map<long,complex<double> > ::iterator it;
    std::map<long,complex<double> > col_indices;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    long mask,mask_u,mask_d,b,p,n,m,i,j,k,l,t,nsignu,nsignd;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            // start of new row of nonzero elements
            // select two electrons in left-basis <m_1, m_2|
            for(n=0; n<nsite; n++)
                for(m=0; m<nsite; m++) {
                    mask=(1<<n)+(1<<m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((sector.id_up[i]&mask)==mask && m!=n) {
                        // b is the rest electon positions
                        b=sector.id_up[i]^mask;
                        long nt,mt,mask_ut,occ_ut;
                        nsignu=0;
                        // perform translation
                        for(t=1; t<nsite; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n-t<0) {
                                nt=n-t+nsite;
                                nsignu++;
                            }
                            else
                                nt=n-t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m+t>=nsite) {
                                mt=m+t-nsite;
                                nsignu++;
                            }
                            else
                                mt=m+t;
                            // the translated two electrons indices
                            mask_ut=(1<<nt)+(1<<mt);
                            // occupation of electons on the translated position
                            occ_ut=mask_ut&b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_ut==0) {
                                // looking up Lin's table, and find the corresponding index
                                if(sector.basis_up.find(mask_ut+b)!=sector.basis_up.end())
                                    k=sector.basis_up[mask_ut+b];
                                else
                                    k=i;
                                // calculate the Coulomb matrix contribution and
                                // add it to the hamiltonian matrix
                                // Coulomb matrix element, in Landau gauge
                                if(k!=i) {
                                    long q_x,q_y;
                                    q_y=t;
                                    complex<double> V_uu=0;
                                    for(q_x=0; q_x<nsite; q_x++)
                                        if(q_y!=0&&q_x!=0)
                                            V_uu+=Coulomb_interaction(0,0,q_x,q_y)*complex<double>(cos((n-m-t)*q_x*2.0*M_PI/nphi),sin((n-m-t)*q_x*2.0*M_PI/nphi))/(2.0*nphi);
                                    it=col_indices.find(k*nbasis_down+j);
                                    if(it==col_indices.end())
                                        col_indices.insert(std::pair<long,complex<double> >(k*nbasis_down+j,V_uu*pow(-1,nsignu)));
                                    else
                                        it->second+=V_uu*pow(-1,nsignu);
                                }
                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_ut==mask_ut)
                                nsignu+=2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_ut!=0 && occ_ut!=mask_ut)
                                nsignu++;
                        }
                    }

                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((sector.id_down[j]&mask)==mask && m!=n) {
                        // b is the rest electon positions
                        b=sector.id_down[j]^mask;
                        long nt,mt,mask_dt,occ_dt;
                        nsignd=0;
                        // perform translation
                        for(t=1; t<nsite; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n-t<0) {
                                nt=n-t+nsite;
                                nsignd++;
                            }
                            else
                                nt=n-t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m+t>=nsite) {
                                mt=m+t-nsite;
                                nsignd++;
                            }
                            else
                                mt=m+t;
                            // the translated two electrons indices
                            mask_dt=(1<<nt)+(1<<mt);
                            // occupation of electons on the translated position
                            occ_dt=mask_dt&b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_dt==0) {
                                // the translated indices
                                if(sector.basis_down.find(mask_dt+b)!=sector.basis_down.end())
                                    l=sector.basis_down[mask_dt+b];
                                else
                                    l=j;
                                // calculate the Coulomb matrix contribution and
                                // add it to the hamiltonian matrix
                                // Coulomb matrix element, in Landau gauge
                                if(l!=j) {
                                    long q_x,q_y;
                                    q_y=t;
                                    complex<double> V_dd=0;
                                    for(q_x=0; q_x<nsite; q_x++)
                                        if(q_y!=0&&q_x!=0)
                                            V_dd+=Coulomb_interaction(1,1,q_x,q_y)*complex<double>(cos((n-m-t)*q_x*2.0*M_PI/nphi),sin((n-m-t)*q_x*2.0*M_PI/nphi))/(2.0*nphi);
                                    it=col_indices.find(i*nbasis_down+k);
                                    if(it==col_indices.end())
                                        col_indices.insert(std::pair<long,complex<double> >(i*nbasis_down+l,V_dd*pow(-1,nsignd)));
                                    else
                                        it->second+=V_dd*pow(-1,nsignd);
                                }
                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_dt==mask_dt)
                                nsignd+=2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_dt!=0 && occ_dt!=mask_dt)
                                nsignd++;
                        }
                    }
                    // consider the one electron in the upper layer
                    // and one electron in the lower layer case
                    mask_u=(1<<n);
                    mask_d=(1<<m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((sector.id_up[i]&mask_u)==mask_u && (sector.id_down[j]&mask_d)==mask_d) {
                        // b is the rest electon positions for upper-layer electrons
                        b=sector.id_up[i]^mask_u;
                        p=sector.id_down[j]^mask_d;
                        long nt,mt,mask_ut,occ_ut,mask_dt,occ_dt;
                        nsignu=0;
                        nsignd=0;
                        // perform translation
                        for(t=1; t<nsite; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n-t<0) {
                                nt=n-t+nsite;
                                nsignu++;
                            }
                            else
                                nt=n-t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m+t>=nsite) {
                                mt=m+t-nsite;
                                nsignd++;
                            }
                            else
                                mt=m+t;
                            // the translated upper electron index
                            mask_ut=(1<<nt);
                            mask_dt=(1<<mt);
                            // occupation of electons on the translated position
                            occ_ut=mask_ut&b;
                            occ_dt=mask_dt&p;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_ut==0 && occ_dt==0) {
                                // the translated indices
                                if(sector.basis_up.find(mask_ut+b)!=sector.basis_up.end())
                                    k=sector.basis_up[mask_ut+b];
                                else
                                    k=i;
                                if(sector.basis_down.find(mask_d+p)!=sector.basis_down.end())
                                    l=sector.basis_down[mask_dt+p];
                                else
                                    l=j;
                                // calculate the Coulomb matrix contribution and
                                // add it to the hamiltonian matrix
                                // Coulomb matrix element, in Landau gauge
                                if(k!=i&& l!=j) {
                                    long q_x,q_y;
                                    q_y=t;
                                    complex<double> V_ud=0;
                                    for(q_x=0; q_x<nsite; q_x++)
                                        if(q_y!=0&&q_x!=0) {
                                            V_ud+=Coulomb_interaction(1,0,q_x,q_y)*complex<double>(cos((n-m-t)*q_x*2.0*M_PI/nphi),sin((n-m-t)*q_x*2.0*M_PI/nphi))/(2.0*nphi);
                                            //V_ud+=complex<double>(cos((n-m-t)*q_x*2.0*M_PI/nphi),sin((n-m-t)*q_x*2.0*M_PI/nphi))/(2.0*nphi);
                                            cout<<"q:="<<sqrt(q_x*q_x+q_y*q_y)*2.0*M_PI/nsite<<endl;
                                            cout<<"Coulomb interaction:="<<Coulomb_interaction(1,0,q_x,q_y)<<endl;
                                        }
                                    it=col_indices.find(k*nbasis_down+l);
                                    if(it==col_indices.end())
                                        col_indices.insert(std::pair<long,complex<double> >(k*nbasis_down+l,V_ud*pow(-1,nsignu+nsignd)));
                                    else
                                        it->second+=V_ud*pow(-1,nsignu+nsignd);
                                    cout<<"("<<sector.id_up[k]<<","<<sector.id_down[l]<<")="<<V_ud*pow(-1,nsignu+nsignd)<<endl;
                                }
                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_ut==mask_ut && occ_dt==mask_dt)
                                nsignu+=2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_ut==0 && occ_dt==mask_dt || occ_dt==0 && occ_ut==mask_ut)
                                nsignu++;
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

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this !=&_gs_hconfig) {
        seed=_gs_hconfig.seed;
        nHilbert=_gs_hconfig.nHilbert;
        H=_gs_hconfig.H;
        d=_gs_hconfig.d;
        nphi=_gs_hconfig.nphi;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(),_gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(),_gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(),_gs_hconfig.psi_n0.end());
    }
    return *this;
}

double hamil::Coulomb_interaction(int alpha, int beta, int q_x,int q_y) {
    double q=sqrt(q_x*q_x+q_y*q_y)*2*M_PI/nsite;
    if(alpha==beta)
        return 1.0/(q+1e-30)*exp(-q*q/2.0);
    else
        return 1.0/(q+1e-30)*exp(-q*q/2.0-q*d);
}

double hamil::spectral_function(vector<complex<double> > &O_phi_0,double omega,double _E0, double eta, int annil) {
    complex<double> E;
    complex<double> G=0;
    for(int i=0; i<nHilbert; i++)
        // set annil==1, which gives hole-sector
        if(annil==1) {
            E=complex<double>(omega,eta);
            G+=pow(conj(psi_n0[i])*O_phi_0[i],2)/(E+eigenvalues[i]-_E0);
        }
    // else particle-sector
        else {
            E=complex<double>(omega,eta);
            G+=pow(conj(psi_n0[i])*O_phi_0[i],2)/(E+_E0-eigenvalues[i]);
        }

    return -G.imag()/M_PI;
}

double hamil::ground_state_energy() {
    if(psi_0.size()==0) return 0;
    complex<double> E_gs=0;
    vector< complex<double> > psi_t;
    psi_t=H*psi_0;
    for(int i=0; i<nHilbert; i++)
        E_gs+=conj(psi_t[i])*psi_0[i];
    return E_gs.real();
}

void hamil::diag() {
    int i,idx;
    complex<double> *hamiltonian=new complex<double>[nHilbert*nHilbert];
    double *en=new double[nHilbert];
    memset(hamiltonian,0,sizeof(complex<double> )*nHilbert*nHilbert);
    for(i=0; i<H.outer_starts.size()-1; i++)
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            hamiltonian[i*nHilbert+H.inner_indices[idx]]=H.value[idx];
    diag_zheev(hamiltonian,en,nHilbert);
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

void hamil::print_eigen() {
    std::cout<<"Eigenvalues:=[ ";
    for(int i=0; i<nHilbert; i++)
        if(i!=nHilbert-1)
            std::cout<<eigenvalues[i]<<", ";
        else
            std::cout<<eigenvalues[i]<<" ]"<<std::endl;
}
