#include"hamiltonian.h"
hamil::hamil() {}

double hamil::Coulomb_interaction(int alpha,int q_x, int q_y) {
    double q = sqrt(q_x * q_x / (lx * lx) + q_y * q_y / (ly * ly)) * 2.0 * M_PI;
    if(alpha ==0)
        return 2.0*M_PI/q*exp(-q*q/2.0)*pow(1.0-exp(-q*q/2.0),nLL*2);
    else
        return 2.0*M_PI/q*exp(-q*q/2.0-q*d)*pow(1.0-exp(-q*q/2.0),nLL*2);
}

void hamil::init_Coulomb_matrix() {
    Coulomb_matrix.assign(2 * nphi*nphi, 0);
    for(int alpha = 0; alpha < 2; alpha++)
        // n=j_1, m=_j3
        for(int s = 0; s < nphi; s++)
            for(int q_y = 0; q_y < nphi; q_y++) {
                double V=0;
                for(int q_x = -nphi/2; q_x <=nphi/2; q_x++)
                    if(!(q_x==0 && q_y==0))
                        V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);
                // Coulomb matrix elements in Landau gauge
                Coulomb_matrix[alpha*nphi*nphi+s*nphi+q_y]=V;
            }
    // initialize classical Coulomb energy
    E_cl=-2.0;
    for(int i=0; i<nphi; i++)
        for(int j=0; j<nphi; j++)
            if(!(i==0 &&j==0))
                E_cl+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
    E_cl/=sqrt(lx*ly);
}

void hamil::set_hamil(basis & sector, double _lx, double _ly, long _nphi,long _nLL, double _d)
{
    d = _d;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    hamiltonian.assign(nHilbert*nHilbert,0);
    int kx=sector.K;
    unsigned long long lbasis,rbasis,mask,mask_t,occ_t,b;
    long i,j,k,l;
    int n,m,s,t,nt,mt,sign,signl,signr,kl,kr,Cl,Cr;

    for(i = 0; i < nHilbert; i++) {
            for(int C=1; C<=sector.C; C++)
            if(sector.translate(sector.id[i],C,sign)==sector.id[i]) {
                Cl=C;
                break;
            }
            for(kl=0; kl<Cl; kl++) {
                lbasis=sector.translate(sector.id[i],kl,signl);
                for(n = 0; n < nphi-1; n++)
                    for(m = n+1; m < nphi; m++) {
                        mask = (1 << n) + (1 << m);
                        // consider the upper-layer two electrons
                        // looking up the corresponding basis in id_up
                        // if there're two electrons on n and m;
                        if((lbasis &mask) == mask && n!=m) {
                            // b is the rest electon positions
                            b = lbasis ^ mask;
                            // mt=j3, nt=j4
                            // perform translation along x-direction (q_y), positive q_y
                            for(t = -nphi/2; t <= nphi/2; t++) {
                                if(n + t >=nphi)
                                    nt = n + t - nphi;
                                else if (n+t <0)
                                    nt = n + t +nphi;
                                else
                                    nt = n + t;
                                if(m - t <0)
                                    mt = m - t + nphi;
                                else if (m - t >=nphi)
                                    mt = m - t -nphi;
                                else
                                    mt = m - t;

                                s=fabs(mt-n);
                                // the translated two electrons indices
                                mask_t = (1 << nt) + (1 << mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                // looking up Lin's table, and find the corresponding index
                                if(occ_t == 0) {
                                    // determine the subbasis size of right side basis
                                    for(int C=1; C<=sector.C; C++)
                                        if(sector.translate(mask_t +b, C,sign)==(mask_t+b)) {
                                            Cr=C;
                                            break;
                                        }
                                    for(kr=0; kr<Cr; kr++) {
                                        rbasis=sector.inv_translate(mask_t+b,kr,signr);
                                        if(sector.basis_set.find(rbasis) != sector.basis_set.end())
                                        {
                                            j = sector.basis_set[rbasis];
                                            sign=sector.get_sign(lbasis,n,m,nt,mt);
                                            complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(kl-kr)/sector.C),sin(2.0*M_PI*kx*(kl-kr)/sector.C))/sqrt(Cl*Cr);
                                            hamiltonian[i*nHilbert+j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor*signl*signr;
                                        }
                                    }
                                }
                            }
                        }
                    }

            // down-layer
            for(n = nphi; n < 2*nphi-1; n++)
               for(m = n+1; m < 2*nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // consider the lower-layer two electrons
                        // if there're two electrons on n and m;
                        if((lbasis &mask) == mask && m!=n) {
                            // p is the rest electon positions
                            b = lbasis ^ mask;
                            // perform translation in x-direction, negative q_y
                            for(t = -nphi/2; t <= nphi/2; t++) {
                                if(n + t >=2*nphi)
                                    nt = n + t - nphi;
                                else if (n+t <nphi)
                                    nt = n + t +nphi;
                                else
                                    nt = n + t;
                                if(m - t <nphi)
                                    mt = m - t + nphi;
                                else if (m - t >=2*nphi)
                                    mt = m - t -nphi;
                                else
                                    mt = m - t;
                                s=fabs(mt-n);
                                // the translated two electrons indices
                                mask_t = (1 << nt) + (1 << mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                if(occ_t == 0) {
                                    // determine the subbasis size of the right side basis
                                    for(int C=1; C<=sector.C; C++)
                                        if(sector.translate(mask_t +b, C,sign)==(mask_t+b)) {
                                            Cr=C;
                                            break;
                                        }
                                    for(kr=0; kr<Cr; kr++) {
                                        rbasis=sector.inv_translate(mask_t+b,kr,signr);
                                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                            j = sector.basis_set[rbasis];
                                            sign=sector.get_sign(lbasis,n,m,nt,mt);
                                            complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(kl-kr)/sector.C),sin(2.0*M_PI*kx*(kl-kr)/sector.C))/sqrt(Cl*Cr);
                                            hamiltonian[i*nHilbert+j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor*signl*signr;
                                        }
                                    }

                                }
                            }
                        }
                    }
            
            // consider the one electron in the upper layer
            // and one electron in the lower layer case
            for(n = 0; n < nphi; n++)
               for(m = nphi; m < 2*nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // if there is one electron at site n in upper-layer
                        // and one electron at site m in lower-layer
                        if((lbasis &mask) == mask) {
                            // b is the rest electon positions for upper-layer electrons
                            b = lbasis ^ mask;
                            // perform translation along x-direction
                            for(t = -nphi/2; t <= nphi/2 ; t++) {
                                if(n + t>=nphi)
                                    nt = n + t - nphi;
                                else if (n+t <0)
                                    nt = n + t +nphi;
                                else
                                    nt = n + t;
                                if(m - t <nphi)
                                    mt = m - t + nphi;
                                else if (m - t >=2*nphi)
                                    mt = m - t -nphi;
                                else
                                    mt = m - t;
                                s=fabs(mt-nphi-n);
                                // the translated upper electron index
                                mask_t = (1 << nt)+(1<<mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                // the translated indices
                                if(occ_t == 0) {
                                    // determine the subbasis size of the right side up-basis
                                    for(int C=1; C<=sector.C; C++)
                                        if(sector.translate(mask_t +b, C,sign)==(mask_t+b)) {
                                            Cr=C;
                                            break;
                                        }
                                    for(kr=0; kr<Cr; kr++) {
                                        rbasis=sector.inv_translate(mask_t+b,kr,signr);
                                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                            j = sector.basis_set[rbasis];
                                            sign=sector.get_sign(lbasis,n,m,nt,mt);
                                            complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(kl-kr)/sector.C),sin(2.0*M_PI*kx*(kl-kr)/sector.C))/sqrt(Cl*Cr);
                                            hamiltonian[i*nHilbert+j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT_factor*signl*signr;
                                        }
                                    }
                                }
                            }
                        }
                  }
           }
           // diagonal Coulomb classical energy term
           hamiltonian[i*nHilbert+i]+=E_cl*(sector.nel_up+sector.nel_down);
    }
}

hamil::~hamil() {}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this != &_gs_hconfig) {
        nHilbert = _gs_hconfig.nHilbert;
        d = _gs_hconfig.d;
        nphi = _gs_hconfig.nphi;
        hamiltonian.assign(_gs_hconfig.hamiltonian.begin(),_gs_hconfig.hamiltonian.end());
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(), _gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(), _gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(), _gs_hconfig.psi_n0.end());
    }
    return *this;
}

double hamil::ground_state_energy() {
    if(psi_0.size() == 0) return 0;
    complex<double> E_gs = 0;
    vector< complex<double> > psi_t;
    psi_t.assign(nHilbert,0);
    for(int i=0; i<nHilbert; i++)
        for(int j=0; j<nHilbert; j++)
            psi_t[i]+=hamiltonian[i*nHilbert+j]*psi_0[j];
    for(int i = 0; i < nHilbert; i++)
        E_gs += conj(psi_t[i]) * psi_0[i];
    psi_t.clear();
    return E_gs.real();
}

void hamil::diag() {
    int i;
    complex<double> *h = new complex<double>[nHilbert * nHilbert];
    double *en = new double[nHilbert];
    memset(h, 0, sizeof( complex<double>)*nHilbert * nHilbert);
    for(i = 0; i <nHilbert*nHilbert; i++)
        h[i]=hamiltonian[i];
    diag_zheev(h, en, nHilbert);
    psi_0.assign(nHilbert, 0);
    psi_n0.assign(nHilbert, 0);
    eigenvalues.assign(nHilbert, 0);
    for(i = 0; i < nHilbert; i++) {
        eigenvalues[i] = en[i];
        psi_0[i] = h[i];
        psi_n0[i] = h[i * nHilbert];
    }
    delete h, en;
}

void hamil::print_hamil(int range) {
    int i, j, count;
    if(range>nHilbert)
        range=nHilbert;
    for(i = 0; i < range; i++) {
        if(i == 0)
            cout <<setw(2)<< "[[";
        else cout <<setw(2)<< " [";
        // count is the No. of nonzero elements in the row
        for(j=0; j<range; j++)
            cout<<setw(5)<<setprecision(2)<<hamiltonian[i*nHilbert+j]<<", ";
        if(i == range - 1)
            cout << ",...]]" << endl;
        else cout << ",...]" << endl;
    }
}
void hamil::print_eigen(int range) {
    if(range>=nHilbert)
        range=nHilbert;
    std::cout << "Eigenvalues:=[ ";
    for(int i = 0; i < range; i++)
        if(i != range - 1)
            std::cout << eigenvalues[i] << ", ";
        else
            std::cout << eigenvalues[i] << " , ...]" << std::endl;
}
