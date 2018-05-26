#include"hamiltonian.h"
static std::mutex mutex_update;
hamil::hamil() {}

double hamil::Coulomb_interaction(int alpha,int q_x, int q_y) {
    double q = sqrt(q_x * q_x / (lx * lx) + q_y * q_y / (ly * ly)) * 2.0 * M_PI;
    if(alpha ==0)
        return 2.0*M_PI/q*exp(-q*q/2.0)*pow(1.0-exp(-q*q/2.0),nLL*2);
    else
        return 2.0*M_PI/q*exp(-q*q/2.0-q*d)*pow(1.0-exp(-q*q/2.0),nLL*2);
}

void hamil::init_Coulomb_matrix() {
    Ec=-2.0;
    for(int i=0; i<nphi; i++)
        for(int j=0; j<nphi; j++)
            if(!(i==0 &&j==0))
                //Ec+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
                Ec+=erfc(sqrt(M_PI*(i*i*lx/ly+j*j*ly/lx)));
    Ec/=sqrt(lx*ly);
    // classical Coulomb energy for interwell interaction
    Coulomb_matrix.assign(2 * nphi*nphi, 0);
    for(int alpha = 0; alpha < 2; alpha++)
        // n=j_1, m=_j3
        for(int s = 0; s < nphi; s++)
            for(int q_y = 0; q_y < nphi; q_y++) {
                double V=0;
                for(int q_x = -nphi/2; q_x <=nphi/2; q_x++)
                    if(!(q_x==0 && q_y==0))
                        V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);

                if(alpha==1) {
                    V=0;
                    for(int q_x = -50*nphi/d; q_x <50*nphi/d; q_x++)
                        if(!(q_x==0 &&q_y==0))
                            V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);
                }
                // Coulomb matrix elements in Landau gauge
                Coulomb_matrix[alpha*nphi*nphi+s*nphi+q_y]=V;
            }
    // store FT coefficients
    int kx=sector.K;
    FT.assign(sector.C*sector.C,0);
    for(int kl=0; kl<sector.C; kl++)
        for(int kr=0; kr<sector.C; kr++)
            FT[kl*sector.C+kr]=complex<double>(cos(2.0*M_PI*(kl-kr)*kx/sector.C),sin(2.0*M_PI*(kl-kr)*kx/sector.C));
}

inline void hamil::peer_set_hamil(int id, long nbatch,long nrange) {
    int kx=sector.K;
    unsigned long lbasis,rbasis,rbasis_0,mask,mask_t,occ_t,b;
    int n,m,s,t,nt,mt,sign,signl,signr,kl,kr,Cl,Cr,C;
    int ql,qr,Dl,Dr,D,_signl,_signr;
    unsigned long _lbasis,_rbasis;


    long i,j,k,l;
    vector<complex<double> > matrix_elements;
    for(int _i = 0; _i < nrange; _i++) {
        // i for each thread
        i=_i+id*nbatch;
        matrix_elements.assign(nHilbert,0);

        // determine the size of relative translate of this basis
        Dl=sector.C;
        for(D=1; D<sector.C; D++)
            if(sector.translate(sector.id[i],D,_signl)==sector.id[i]) {
                Dl=D;
                break;
            }
        // if parameter kx<0, do not perform basis translation
        Dl=(kx<0?1:Dl);
        for(ql=0; ql<Dl; ql++) {
            _signl=1;
            _lbasis=(ql==0?sector.id[i]:sector.translate(sector.id[i],ql,_signl));
            // determine the size of translate of this basis
            Cl=sector.C;
            for(C=1; C<sector.C; C++)
                if(sector.translate(_lbasis,C,signl)==_lbasis) {
                    Cl=C;
                    break;
                }
            // if parameter kx<0, do not perform basis translation
            Cl=(kx<0?1:Cl);
	    Cl=1;
            for(kl=0; kl<Cl; kl++) {
                signl=1;
                // translation with 0 bits is not necessary
                lbasis=(kl==0?_lbasis:sector.translate(_lbasis,kl,signl));
                //down-layer

                for(n=0; n<nphi-1; n++)
                    for(m = n+1; m < nphi; m++) {
                        mask = (1 << n) + (1 << m);
                        // consider the upper-layer two electrons
                        // looking up the corresponding basis in id_up
                        // if there're two electrons on n and m;
                        if((lbasis &mask) == mask ) {
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

                                s=abs(mt-n);
                                // the translated two electrons indices
                                mask_t = (1 << nt) + (1 << mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                rbasis_0=mask_t+b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                if(occ_t == 0 ) {
                                    // determine the right side size of the translation
                                    Cr=sector.C;
                                    for(C=1; C<sector.C; C++)
                                        if(sector.translate(rbasis_0,C,signr)==rbasis_0) {
                                            Cr=C;
                                            break;
                                        }
                                    // if parameter kx<0, do not perform basis translation
                                    Cr=(kx<0?1:Cr);
				    Cr=1;
                                    for(kr=0; kr<Cr; kr++) {
                                        signr=1;
                                        // 0 bits shifting is not performed
                                        _rbasis=(kr==0?rbasis_0:sector.inv_translate(rbasis_0,kr,signr));
                                        // determine the right side size of the relative_translation

                                        Dr=sector.C;
                                        for(D=1; D<sector.C; D++)
                                            if(sector.inv_translate(_rbasis,D,_signr)==_rbasis) {
                                                Dr=D;
                                                break;
                                            }

                                        // if parameter kx<0, do not perform basis translation
                                        Dr=(kx<0?1:Dr);
                                        for(qr=0; qr<Dr; qr++){
                                            _signr=1;
                                            rbasis=(qr==0?_rbasis:sector.inv_translate(_rbasis,qr,_signr));
                                            //rbasis=sector.inv_translate(_rbasis,ql,_signr);
                                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                                j = sector.basis_set[rbasis];
                                                sign=sector.get_sign(lbasis,n,m,nt,mt)*signl*signr*_signl*_signr;
                                                matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*sector.C+qr]/sqrt(Cl*Cr*Dl*Dr);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                // upper-layer
                for( n=nphi; n<2*nphi-1; n++)
                    for( m = n+1; m < 2*nphi; m++) {
                        mask = (1 << n) + (1 << m);
                        // consider the lower-layer two electrons
                        // if there're two electrons on n and m;
                        if((lbasis &mask) == mask && n!=m) {
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
                                s=abs(mt-n);
                                // the translated two electrons indices
                                mask_t = (1 << nt) + (1 << mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                rbasis_0=mask_t+b;
                                if(occ_t == 0 && nt!=mt) {
                                    // determine the right side size of the translation
                                    Cr=sector.C;
                                    for(C=1; C<sector.C; C++)
                                        if(sector.translate(rbasis_0,C,signr)==rbasis_0) {
                                            Cr=C;
                                            break;
                                        }
                                    // if parameter kx<0, do not perform basis translation
                                    Cr=(kx<0?1:Cr);
				    Cr=1;
                                    for(kr=0; kr<Cr; kr++) {
                                        signr=1;
                                        // 0 bits shifting is not performed
                                        _rbasis=(kr==0?rbasis_0:sector.inv_translate(rbasis_0,kr,signr));
                                        // determine the right side size of the relative_translation

                                        Dr=sector.C;
                                        for(D=1; D<sector.C; D++)
                                            if(sector.translate(_rbasis,D,_signr)==_rbasis) {
                                                Dr=D;
                                                break;
                                            }

                                        // if parameter kx<0, do not perform basis translation
                                        Dr=(kx<0?1:Dr);
                                        for(qr=0; qr<Dr; qr++) {
                                            _signr=1;
                                            rbasis=(qr==0?_rbasis:sector.inv_translate(_rbasis,qr,_signr));
                                            //rbasis=sector.inv_translate(_rbasis,ql,_signr);
                                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                                j = sector.basis_set[rbasis];
                                                sign=sector.get_sign(lbasis,n,m,nt,mt)*signl*signr*_signl*_signr;
                                                matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*sector.C+qr]/sqrt(Cl*Cr*Dl*Dr);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                for(n=0; n<nphi; n++)
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
                                s=abs(mt-nphi-n);
                                // the translated electron index
                                mask_t = (1 << nt)+(1<<mt);
                                // occupation of electons on the translated position
                                occ_t = mask_t & b;
                                // if there're no electon on the translated position
                                // which is a valid translation, can be applied
                                // the translated indices
                                rbasis_0=mask_t+b;
                                if(occ_t == 0 ) {
                                    // determine the right side size of the translation
                                    Cr=sector.C;
                                    for(C=1; C<sector.C; C++)
                                        if(sector.translate(rbasis_0,C,signr)==rbasis_0) {
                                            Cr=C;
                                            break;
                                        }
                                    // if parameter kx<0, do not perform basis translation
                                    Cr=(kx<0?1:Cr);
				    Cr=1;
                                    for(kr=0; kr<Cr; kr++) {
                                        signr=1;
                                        // 0 bits shifting is not performed
                                        _rbasis=(kr==0?rbasis_0:sector.inv_translate(rbasis_0,kr,signr));
                                        // determine the right side size of the translation
                                        Dr=sector.C;
                                        for(D=1; D<sector.C; D++)
                                            if(sector.translate(_rbasis,D,_signr)==_rbasis) {
                                                Dr=D;
                                                break;
                                            }
                                        // if parameter kx<0, do not perform basis translation
                                        Dr=(kx<0?1:Dr);
                                        for(qr=0; qr<Dr; qr++) {
                                            _signr=1;
                                            rbasis=(qr==0?_rbasis:sector.inv_translate(_rbasis,qr,_signr));
                                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                                j = sector.basis_set[rbasis];
                                                sign=sector.get_sign(lbasis,n,m,nt,mt)*signl*signr*_signl*_signr;
                                                matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*sector.C+qr]/sqrt(Cl*Cr*Dl*Dr);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
               }
        }

        matrix_elements[i]+=Ec*(sector.nel_up+sector.nel_down);
        long count=0;
        mutex_update.lock();
        for(k=0; k<nHilbert; k++)
	    hamiltonian[i*nHilbert+k]=matrix_elements[k]; 
        mutex_update.unlock();
    }
    matrix_elements.clear();
}

void hamil::set_hamil(double _lx, double _ly, long _nphi,long _nLL, double _d,int nthread){
    d = _d;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    hamiltonian.assign(nHilbert*nHilbert,0);
    std::vector<std::thread> threads;
    long  nbatch=nHilbert/nthread;
    long nresidual=nHilbert%nthread;
    for(int id = 0; id < nthread; id++)
        threads.push_back(std::thread(&hamil::peer_set_hamil,this,id,nbatch,nbatch));
    for(auto &th:threads)
        if(th.joinable())
            th.join();
    if(nresidual!=0)
        peer_set_hamil(nthread,nbatch,nresidual);
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
