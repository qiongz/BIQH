#include"lanczos_hamiltonian.h"
static std::mutex mutex_update;

inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,const long &_nHilbert,const long &_lambda, const unsigned &_seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(const long &_lambda,const unsigned &_seed) {
    lambda=_lambda;
    seed=_seed;
}

lhamil::~lhamil() {
}

const lhamil & lhamil::operator =(const lhamil & _config) {
    if(this !=&_config) {
        nHilbert=_config.nHilbert;
        H=_config.H;
        seed=_config.seed;
        lambda=_config.lambda;
        E0=_config.E0;
        norm.assign(_config.norm.begin(),_config.norm.end());
        overlap.assign(_config.overlap.begin(),_config.overlap.end());
        psir_0.assign(_config.psir_0.begin(),_config.psir_0.end());
        psi_n0.assign(_config.psi_n0.begin(),_config.psi_n0.end());
        eigenvalues.assign(_config.eigenvalues.begin(),_config.eigenvalues.end());
    }
    return *this;
}

double lhamil::Coulomb_interaction(const int &alpha,const int &q_x, const int &q_y) {
    double q=sqrt(q_x*q_x/(lx*lx)+q_y*q_y/(ly*ly))*2.0*M_PI;
    if(alpha==1)
        return 2.0*M_PI/q*exp(-q*q/2.0-q*d)*pow(1.0-exp(-q*q/2.0),nLL*2);
    else
        return 2.0*M_PI/q*exp(-q*q/2.0)*pow(1.0-exp(-q*q/2.0),nLL*2);
}

void lhamil::init_Coulomb_matrix(const double &theta_x) {
    double theta_1,theta_2;
    Ec=-2.0;
    for(int i=0; i<nphi; i++)
        for(int j=0; j<nphi; j++)
            if(!(i==0 &&j==0))
                //Ec+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
                Ec+=erfc(sqrt(M_PI*(i*i*lx/ly+j*j*ly/lx)));
    Ec/=sqrt(lx*ly);
    // classical Coulomb energy for interwell interaction
    Coulomb_matrix.assign(2 * nphi*nphi, 0);
    for(int alpha = 0; alpha < 2; alpha++) {
        // n=j_1, m=_j3
        // s=j_1-j_3
        for(int s = 0; s < nphi; s++)
            for(int q_y = 0; q_y <nphi; q_y++) {
                double V=0;
                for(int q_x = -nphi; q_x <=nphi; q_x++)
                    if(!(q_x==0 && q_y==0))
                        V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);

                if(alpha==1) {
                    V=0;
                    for(int q_x = -10*nphi/(d+0.01); q_x <10*nphi/(d+0.01); q_x++)
                        //for(int q_x = -nphi; q_x <=nphi; q_x++)
                        if(!(q_x==0 &&q_y==0))
			// theta_x=theta_x^L-theta_x^R, twisted phase to calculate spin stiffness
                            V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos((2.0*M_PI*s-theta_x)*q_x/nphi)/(2.0*lx*ly);
                }
                // Coulomb matrix elements in Landau gauge
                Coulomb_matrix[alpha*nphi*nphi+s*nphi+q_y]=V;
            }
    }
    // store FT coefficients
    int kx=sector.K;
    FT.assign(nphi*nphi,0);
    for(int kl=0; kl<nphi; kl++)
        for(int kr=0; kr<nphi; kr++)
            FT[kl*nphi+kr]=complex<double>(cos(2.0*M_PI*(kl-kr)*kx/nphi),sin(2.0*M_PI*(kl-kr)*kx/nphi));
}

inline void lhamil::peer_set_hamil(const double &Delta_SAS,const double &Delta_V,const double &Delta_Z,const double &_theta_B,const double &theta_x, const double &theta_y,const int &id, const long &nbatch,const long &nrange) {
    int kx=sector.K;
    unsigned long lbasis,rbasis,rbasis_0,mask,mask_t,occ_t,b;
    int n,m,s,t,nt,mt,sign,signl,signr;
    int ql,qr,Dl,Dr,D;

    long i,j,k,l;
    vector<complex<double> > matrix_elements;
    vector<long> inner_indices;
    vector<complex<double> > value;
    inner_indices.reserve(nphi);
    value.reserve(nphi);

    for(int _i = 0; _i < nrange; _i++) {
        // i for each thread
        i=_i+id*nbatch;
        matrix_elements.assign(nHilbert,0);

        Dl=(kx<0?1:sector.basis_C[i]);
        for(ql=0; ql<Dl; ql++) {
            signl=1;
            lbasis=(ql==0?sector.id[i]:sector.translate(sector.id[i],ql,signl));

            // upper-layer spin-up interaction
            for(n=0; n<nphi-1; n++)
                for(m = n+1; m < nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // b is the rest electon positions
                        b = lbasis ^ mask;
                        // mt=j3, nt=j4
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <0) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <0) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;

                            s=abs(mt-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            rbasis_0=mask_t+b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-layer spin-down interaction
            for(n=nphi; n<2*nphi-1; n++)
                for(m = n+1; m < 2*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // b is the rest electon positions
                        b = lbasis ^ mask;
                        // mt=j3, nt=j4
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=2*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=2*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;

                            s=abs(mt-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            rbasis_0=mask_t+b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-layer spin-up and spin-down interspin interaction
            for(n=0; n<nphi; n++)
                for(m = nphi; m < 2*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // b is the rest electon positions
                        b = lbasis ^ mask;
                        // mt=j3, nt=j4
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <0) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=2*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;

                            s=abs(mt-nphi-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            rbasis_0=mask_t+b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // down-layer spin-up interaction
            for( n=2*nphi; n<3*nphi-1; n++)
                for( m = n+1; m < 3*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // p is the rest electon positions
                        b = lbasis ^ mask;
                        // perform translation in x-direction, negative q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=3*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <2*nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <2*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=3*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }

                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }
            //down-layer spin-down interaction
            for( n=3*nphi; n<4*nphi-1; n++)
                for( m = n+1; m < 4*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // p is the rest electon positions
                        b = lbasis ^ mask;
                        // perform translation in x-direction, negative q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=4*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <3*nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <3*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=4*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }

                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            //down-layer spin-up and spin-down interspin interaction
            for( n=2*nphi; n<3*nphi; n++)
                for( m = 3*nphi; m < 4*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask) {
                        // p is the rest electon positions
                        b = lbasis ^ mask;
                        // perform translation in x-direction, negative q_y
                        for(t = -nphi+1; t < nphi; t++) {
                            if(n + t >=3*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <2*nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <3*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=4*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-nphi-n);
                            // the translated two electrons indices
                            mask_t = (1UL << nt) + (1UL << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 && nt!=mt) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }

                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-down layer spin-up electrons interaction
            for(n=0; n<nphi; n++)
                for(m = 2*nphi; m < 3*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((lbasis &mask) == mask ) {
                        // b is the rest electon positions for upper-layer electrons
                        b = lbasis ^ mask;
                        // perform translation along x-direction
                        for(t = -nphi+1; t < nphi ; t++) {
                            if(n + t>=nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <0) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <2*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=3*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-2*nphi-n);
                            // the translated electron index
                            mask_t = (1UL << nt)+(1UL<<mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // the translated indices
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 ) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-down layer spin-down electrons interaction
            for(n=nphi; n<2*nphi; n++)
                for(m = 3*nphi; m < 4*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((lbasis &mask) == mask ) {
                        // b is the rest electon positions for upper-layer electrons
                        b = lbasis ^ mask;
                        // perform translation along x-direction
                        for(t = -nphi+1; t < nphi ; t++) {
                            if(n + t>=2*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <3*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=4*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-2*nphi-n);
                            // the translated electron index
                            mask_t = (1UL << nt)+(1UL<<mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // the translated indices
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 ) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-layer spin up to down-layer spin-down interaction
            for(n=0; n<nphi; n++)
                for(m = 3*nphi; m < 4*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((lbasis &mask) == mask ) {
                        // b is the rest electon positions for upper-layer electrons
                        b = lbasis ^ mask;
                        // perform translation along x-direction
                        for(t = -nphi+1; t < nphi ; t++) {
                            if(n + t>=nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <0) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <3*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=4*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-3*nphi-n);
                            // the translated electron index
                            mask_t = (1UL << nt)+(1UL<<mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // the translated indices
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 ) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
					complex<double> FT_twisted=FT[ql*nphi+qr]*complex<double>(cos(theta_y*t/nphi),sin(theta_y*t/nphi));
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT_twisted/sqrt(Dl*Dr);
                                        //matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            // upper-layer spin down to down-layer spin-up interaction
            for(n=nphi; n<2*nphi; n++)
                for(m = 2*nphi; m < 3*nphi; m++) {
                    mask = (1UL << n) + (1UL << m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((lbasis &mask) == mask ) {
                        // b is the rest electon positions for upper-layer electrons
                        b = lbasis ^ mask;
                        // perform translation along x-direction
                        for(t = -nphi+1; t < nphi ; t++) {
                            if(n + t>=2*nphi) {
                                nt = n + t - nphi;
                            }
                            else if (n+t <nphi) {
                                nt = n + t +nphi;
                            }
                            else
                                nt = n + t;
                            if(m - t <2*nphi) {
                                mt = m - t + nphi;
                            }
                            else if (m - t >=3*nphi) {
                                mt = m - t -nphi;
                            }
                            else
                                mt = m - t;
                            s=abs(mt-nphi-n);
                            // the translated electron index
                            mask_t = (1UL << nt)+(1UL<<mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // the translated indices
                            rbasis_0=mask_t+b;
                            if(occ_t == 0 ) {
                                // determine the right side size of the translation
                                Dr=nphi;
                                for(D=1; D<nphi; D++)
                                    if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                        Dr=D;
                                        break;
                                    }
                                // if parameter kx<0, do not perform basis translation
                                Dr=(kx<0?1:Dr);
                                for(qr=0; qr<Dr; qr++) {
                                    signr=1;
                                    rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        j = sector.basis_set.at(rbasis);
                                        sign=sector.get_sign(lbasis,n,m,nt,mt,t)*signl*signr;
					complex<double> FT_twisted=FT[ql*nphi+qr]*complex<double>(cos(theta_y*t/nphi),sin(theta_y*t/nphi));
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT_twisted/sqrt(Dl*Dr);
                                        //matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                                    }
                                }
                            }
                        }
                    }
                }

            if(Delta_SAS>0) {
                //spin-up up-down layer tunneling
                for(n=0; n<nphi; n++) {
                    nt=n+2*nphi;
		    complex<double> Q_phase(cos(_theta_B*n*d),sin(_theta_B*n*d));
                    mask = (1UL << n)+(1UL<<nt);
                    unsigned long Kn=mask & lbasis;
                    if(Kn!=mask && Kn!=0) {
                        unsigned long Ln=Kn ^ mask;
                        rbasis_0=lbasis-Kn+Ln;
                        // determine the right side size of the translation
                        Dr=nphi;
                        for(D=1; D<nphi; D++)
                            if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                Dr=D;
                                break;
                            }
                        // if parameter kx<0, do not perform basis translation
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]*Q_phase/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
                //spin-up down-up layer tunneling
                for(n=2*nphi; n<3*nphi; n++) {
                    nt=n-2*nphi;
		    // phase induced by parallel magnetic field
		    complex<double> Q_phase(cos(_theta_B*nt*d),-sin(_theta_B*nt*d));
                    mask = (1UL << n)+(1UL<<nt);
                    unsigned long Kn=mask & lbasis;
                    if(Kn!=mask && Kn!=0) {
                        unsigned long Ln=Kn ^ mask;
                        rbasis_0=lbasis-Kn+Ln;
                        // determine the right side size of the translation
                        Dr=nphi;
                        for(D=1; D<nphi; D++)
                            if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                Dr=D;
                                break;
                            }
                        // if parameter kx<0, do not perform basis translation
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]*Q_phase/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
		
                //spin-down up-down layer tunneling
                for(n=nphi; n<2*nphi; n++) {
                    nt=n+2*nphi;
		    // phase induced by parallel magnetic field
		    complex<double> Q_phase(cos(_theta_B*(n%nphi)*d),sin(_theta_B*(n%nphi)*d));
                    mask = (1UL << n)+(1UL<<nt);
                    unsigned long Kn=mask & lbasis;
                    if(Kn!=mask && Kn!=0) {
                        unsigned long Ln=Kn ^ mask;
                        rbasis_0=lbasis-Kn+Ln;
                        // determine the right side size of the translation
                        Dr=nphi;
                        for(D=1; D<nphi; D++)
                            if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                Dr=D;
                                break;
                            }
                        // if parameter kx<0, do not perform basis translation
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]*Q_phase/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
                //spin-down down-up layer tunneling
                for(n=3*nphi; n<4*nphi; n++) {
                    nt=n-2*nphi;
		    complex<double> Q_phase(cos(_theta_B*(n%nphi)*d),-sin(_theta_B*(n%nphi)*d));
                    mask = (1UL << n)+(1UL<<nt);
                    unsigned long Kn=mask & lbasis;
                    if(Kn!=mask && Kn!=0) {
                        unsigned long Ln=Kn ^ mask;
                        rbasis_0=lbasis-Kn+Ln;
                        // determine the right side size of the translation
                        Dr=nphi;
                        for(D=1; D<nphi; D++)
                            if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                                Dr=D;
                                break;
                            }
                        // if parameter kx<0, do not perform basis translation
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]*Q_phase/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
		
                // spin-up spin-down interlayer tunneling
                for(n=0; n<nphi; n++) {
		nt=n+3*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set[rbasis];
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            matrix_elements[j]+=-0.5*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                        }
                    }
                }
            }
            // spin-down spin-up interlayer tunneling
            for(n=nphi; n<2*nphi; n++) {
		nt=n+nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            matrix_elements[j]+=-0.5*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                        }
                    }
                }
            }
            
            }

        }
        // background charge energy
        matrix_elements[i]+=Ec*sector.nel;
        // charging energy
        matrix_elements[i]+=-d*(sector.get_nel(0,i)+sector.get_nel(1,i))*(sector.get_nel(2,i)+sector.get_nel(3,i))/sector.nphi;
        // bias-voltage (pseudospin Zeemann energy)
        matrix_elements[i]+=-0.5*Delta_V*(sector.get_nel(0,i)+sector.get_nel(1,i)-sector.get_nel(2,i)-sector.get_nel(3,i));
	// Zeemann energy
        matrix_elements[i]+=-0.5*Delta_Z*(sector.get_nel(0,i)+sector.get_nel(2,i)-sector.get_nel(1,i)-sector.get_nel(3,i));

        long count=0;
        for(k=0; k<nHilbert; k++)
            if(abs(matrix_elements[k])>1e-10) {
                inner_indices.push_back(k);
                value.push_back(matrix_elements[k]);
                count++;
            }
        mutex_update.lock();
        H.outer_starts[i]=H.inner_indices.size();
        H.inner_indices.insert(H.inner_indices.end(),inner_indices.begin(),inner_indices.end());
        H.value.insert(H.value.end(),value.begin(),value.end());
        H.outer_size[i]=count;
        mutex_update.unlock();
        inner_indices.clear();
        value.clear();
    }
    matrix_elements.clear();
}


void lhamil::set_hamil(const double &_lx,const double &_ly,const long &_nphi,const long &_nLL,const double &_d,const double &Delta_SAS,const double &Delta_V,const double &Delta_Z,const double &theta_B,const double &theta_x, const double &theta_y,const int &nthread) {
    d = _d;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix(theta_x);
    nHilbert = sector.nbasis;
    H.clear();
    H.inner_indices.reserve(nHilbert * nphi);
    H.value.reserve(nHilbert * nphi);
    H.outer_starts.assign(nHilbert,0);
    H.outer_size.assign(nHilbert,0);

    std::vector<std::thread> threads;
    long  nbatch=nHilbert/nthread;
    long nresidual=nHilbert%nthread;
    for(int id = 0; id < nthread; id++)
        threads.push_back(std::thread(&lhamil::peer_set_hamil,this,Delta_SAS,Delta_V,Delta_Z,theta_B,theta_x,theta_y,id,nbatch,nbatch));
    for(auto &th:threads)
        if(th.joinable())
            th.join();
    if(nresidual!=0)
        peer_set_hamil(Delta_SAS,Delta_V,Delta_Z,theta_B,theta_x,theta_y,nthread,nbatch,nresidual);
}


void lhamil::coeff_update() {
    Vec phi_0,phi_1,phi_2;

    phi_0.init_random(nHilbert,seed);
    norm.push_back(1);
    phi_1 = H*phi_0;
    overlap.push_back( (phi_0 * phi_1).real() );

    phi_1 -= phi_0 * overlap[0];
    norm.push_back( phi_1.normed());

    /* reorthogonalize basis */
    //double q= (phi_0*phi_1).real();
    //phi_1 = (phi_1 - phi_0*q)/(1-q*q);


    for(int i = 1; i < lambda; i++) {
        phi_2 = H * phi_1;
        overlap.push_back( (phi_1 * phi_2).real());
        #pragma ivdep
        phi_2 -= phi_1 * overlap[i] + phi_0 * norm[i];
        norm.push_back( phi_2.normed());

        /* reorthogonalize basis */
        //Gram_Schmidt_orthogonalization( phi_2, i);

        if(norm[i]<1e-15)
            break;

        swap(&phi_0,&phi_1,&phi_2);
    }

    phi_0.clear();
    phi_1.clear();
    phi_2.clear();
}

void lhamil::coeff_explicit_update()
{
    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);
    int i,j,idx;
    double norm_factor;
    complex<double> overlap_factor;
    complex<double> *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new complex<double>[nHilbert];
    phi_1=new complex<double>[nHilbert];
    phi_2=new complex<double>[nHilbert];
    phi_t=new complex<double>[nHilbert];

    //phi_0.init_random(nHilbert,seed);
    #if __cplusplus > 199711L
    std::mt19937 rng(seed);
    for(i=0; i<nHilbert; i++)
        phi_0[i]=complex<double>(rng()*1.0/rng.max()-0.5,rng()*1.0/rng.max()-0.5);
    #else
    init_genrand64(seed);
    for(i=0; i<nHilbert; i++)
        phi_0[i]=complex<double>(genrand64_real3()-0.5,genrand64_real3()-0.5);
    #endif
    norm_factor=0;
    for(i=0; i<nHilbert; i++)
        norm_factor+=std::norm(phi_0[i]);
    norm_factor=sqrt(norm_factor);
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(complex<double>)*nHilbert);
    for(i=0; i<H.outer_starts.size(); i++) {
        //for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
        for(idx=H.outer_starts[i]; idx< H.outer_starts[i]+H.outer_size[i]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    for(i=0; i<nHilbert; i++)
        overlap_factor+=conj(phi_1[i])*phi_0[i];

    //overlap.push_back(overlap_factor.real());
    overlap[0]=overlap_factor.real();

    //phi_1 -= phi_0 * overlap[0];
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    for(i=0; i<nHilbert; i++)
        norm_factor+=std::norm(phi_1[i]);
    norm_factor=sqrt(norm_factor);

    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor;


    //norm.push_back(1);
    //norm.push_back(norm_factor);
    norm[0]=1;
    norm[1]=norm_factor;

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(complex<double>)*nHilbert);
        for(i=0; i<H.outer_starts.size(); i++) {
            //for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            for(idx=H.outer_starts[i]; idx< H.outer_starts[i]+H.outer_size[i]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }
        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        for(i=0; i<nHilbert; i++)
            overlap_factor+=conj(phi_1[i])*phi_2[i];

        //overlap.push_back(overlap_factor.real());
        overlap[j]=overlap_factor.real();

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normed();
        norm_factor=0;
        for(i=0; i<nHilbert; i++)
            norm_factor+=std::norm(phi_2[i]);
        norm_factor=sqrt(norm_factor);

        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor;


        //norm.push_back(norm_factor);
        norm[j+1]=norm_factor;
        //if(norm.back()<1e-8)
        //  break;

        phi_s=phi_0;
        phi_0=phi_1;
        phi_1=phi_2;
        phi_2=phi_s;

    }
    delete phi_0,phi_1,phi_2,phi_t;
}

void lhamil::diag()
{
    int l=overlap.size();
    if(l==0) coeff_update();
    double *e = new double[l];
    double *h = new double [l * l];
    memset(h, 0, sizeof(double)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_dsyevd(h,e,l);

    psi_0.assign(l,0);
    psi_1.assign(l,0);
    psi_n0.assign(l,0);
    eigenvalues.assign(l,0);
    int E1_idx=1;
    for(int i=1;i<l;i++)
	  if(abs(eigenvalues[i]-eigenvalues[0])>1e-2){
		E1_idx=i;
	  	break;
		}

    for(int i=0; i<l; i++) {
        psi_0[i]=h[i];
        psi_1[i]=h[i+l*E1_idx];
        psi_n0[i]=h[i*l];
        eigenvalues[i]=e[i];
    }
    delete h,e;
}


// given the seed of phi_0, generate the Lanczos basis again, and transform the eigenstate to the desired eigenstate representation
void lhamil::eigenstates_reconstruction() {
    int l,n;
    if(overlap.size()==0) {
        coeff_update();
        diag();
    }
    Vec phi_0,phi_1,phi_2;
    E0=eigenvalues[0];
    // repeat the iteration, but with normalization and overlap values at hand
    phi_0.init_random(nHilbert,seed);
    phi_1 = H * phi_0;
    phi_1 -= phi_0*overlap[0];
    phi_1 /= norm[1];
    psir_0.assign(nHilbert,0);
    psir_1.assign(nHilbert,0);
    for(n=0; n<nHilbert; n++){
        psir_0[n] += psi_0[0]*conj(phi_0.value[n])+psi_0[1]*conj(phi_1.value[n]);
        psir_1[n] += psi_1[0]*conj(phi_0.value[n])+psi_1[1]*conj(phi_1.value[n]);
    }

    for(l=2; l<overlap.size(); l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++){
            psir_0[n]+= psi_0[l]*conj(phi_2.value[n]);
            psir_1[n]+= psi_1[l]*conj(phi_2.value[n]);
	}
        swap(&phi_0,&phi_1,&phi_2);
    }

    double Norm=0;
    for(n=0; n<nHilbert; n++)
        Norm+=std::norm(psir_0[n]);
    Norm=sqrt(Norm);
    for(n=0; n<nHilbert; n++)
        psir_0[n]/=Norm;
    
    Norm=0;
    for(n=0; n<nHilbert; n++)
        Norm+=std::norm(psir_1[n]);
    Norm=sqrt(Norm);
    for(n=0; n<nHilbert; n++)
        psir_1[n]/=Norm;

    phi_0.clear();
    phi_1.clear();
    phi_2.clear();
}

double lhamil::ground_state_energy() const{
    vector<complex<double> > H_psir0;
    complex<double> overlap=0;
    if(psir_0.size()!=0) {
        H_psir0=H*psir_0;
        for(int i=0; i<nHilbert; i++)
            overlap+=conj(psir_0[i])*H_psir0[i];
    }
    H_psir0.clear();
    return overlap.real();
}

double lhamil::first_excited_state_energy()const {
    vector<complex<double> > H_psir1;
    complex<double> overlap=0;
    if(psir_1.size()!=0) {
        H_psir1=H*psir_1;
        for(int i=0; i<nHilbert; i++)
            overlap+=conj(psir_1[i])*H_psir1[i];
    }
    H_psir1.clear();
    return overlap.real();
}

double lhamil::occupatation_number(const int &alpha,const int &j)const {
    complex<double> occ=0;
    unsigned long mask,lbasis;
    int sign,q,D;
    // '0' for upper-layer, '1' for down-layer
    if(alpha==0)
        mask = (1UL << j);
    else
        mask = (1UL << (j+nphi));

    for(int i=0; i<nHilbert; i++) {
        D=(sector.K<0?1:sector.basis_C[i]);
        for(q=0; q<D; q++) {
            sign=1;
            lbasis=(q==0?sector.id[i]:sector.translate(sector.id[i],q,sign));
            if((mask& lbasis)==mask)
                occ+=conj(psir_0[i])*psir_0[i]*FT[q*nphi+q]/(1.0*D);
        }
    }
    return occ.real();
}

double lhamil::pseudospin_Sz()const {
    double nel_up=0;
    for(long i=0; i<nHilbert; i++)
	nel_up+=(sector.get_nel(0,i)+sector.get_nel(1,i))*std::norm(psir_0[i]);
    return (2.0*nel_up-sector.nel)/(2.0*sector.nel);
}

double lhamil::Sz() const{
    double nel_up=0;
    for(long i=0; i<nHilbert; i++)
	nel_up+=(sector.get_nel(0,i)+sector.get_nel(2,i))*std::norm(psir_0[i]);
    return (2.0*nel_up-sector.nel)/(2.0*sector.nel);
}

double lhamil::upper_Sz()const
{
    double Sz=0,nel_up;
    for(long i=0; i<nHilbert; i++){
	nel_up=sector.get_nel(0,i)+sector.get_nel(1,i);
	Sz+=(2.0*sector.get_nel(0,i)-nel_up)/(2.0*nel_up+1e-8)*std::norm(psir_0[i]);
    }
    return Sz;
}

double lhamil::down_Sz()const
{
    double Sz=0,nel_up;
    for(long i=0; i<nHilbert; i++){
	nel_up=sector.get_nel(2,i)+sector.get_nel(3,i);
	Sz+=(2.0*sector.get_nel(2,i)-nel_up)/(2.0*nel_up+1e-8)*std::norm(psir_0[i]);
    }
    return Sz;
}

double lhamil::pseudospin_Sx() const{
    unsigned long mask,mask_t,b,occ_t,lbasis,rbasis,rbasis_0;
    long i,j;
    int n,nt,sign,signl,signr;
    int ql,qr,Dl,Dr,D;
    complex<double> Sx_mean=0;
    for(i=0; i<nHilbert; i++) {
        Dl=(sector.K<0?1:sector.basis_C[i]);
        for(ql=0; ql<Dl; ql++) {
            signl=1;
            lbasis=(ql==0?sector.id[i]:sector.translate(sector.id[i],ql,signl));
            //lbasis=sector.id[i];
            // spin-up interlayer tunneling 
            for(n=0; n<nphi; n++) {
		nt=n+2*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            for(n=2*nphi; n<3*nphi; n++) {
		nt=n-2*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            // spin-down interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+2*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            for(n=3*nphi; n<4*nphi; n++) {
		nt=n-2*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            // spin-up spin-down interlayer tunneling 
            for(n=0; n<nphi; n++) {
		nt=n+3*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            for(n=3*nphi; n<4*nphi; n++) {
		nt=n-3*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            
            // spin-down spin-up interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            for(n=2*nphi; n<3*nphi; n++) {
		nt=n-nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                            Sx_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.25/sqrt(Dl*Dr));
                        }
                    }
                }
            }
        }
    }
    return abs(Sx_mean)/sector.nel;
}

double lhamil::spinflip_tunneling() const{
    unsigned long mask,mask_t,b,occ_t,lbasis,rbasis,rbasis_0;
    long i,j;
    int n,nt,sign,signl,signr;
    int ql,qr,Dl,Dr,D;
    complex<double> Sf_mean=0;
    for(i=0; i<nHilbert; i++) {
        Dl=(sector.K<0?1:sector.basis_C[i]);
        for(ql=0; ql<Dl; ql++) {
            signl=1;
            lbasis=(ql==0?sector.id[i]:sector.translate(sector.id[i],ql,signl));
            // spin-up spin-down interlayer tunneling 
            for(n=0; n<nphi; n++) {
		nt=n+3*nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sf_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            
            // spin-down spin-up interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+nphi;
                mask = (1UL << n)+(1UL<<nt);
                unsigned long Kn=mask & lbasis;
                if(Kn!=mask && Kn!=0) {
                    unsigned long Ln=Kn ^ mask;
                    rbasis_0=lbasis-Kn+Ln;
                    // determine the right side size of the translation

                    Dr=nphi;
                    for(D=1; D<nphi; D++)
                        if(sector.translate(rbasis_0,D,signr)==rbasis_0) {
                            Dr=D;
                            break;
                        }
                    // if parameter kx<0, do not perform basis translation
                    Dr=(sector.K<0?1:Dr);
                    for(qr=0; qr<Dr; qr++) {
                        signr=1;
                        rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                        if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                            j = sector.basis_set.at(rbasis);
                            sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                            Sf_mean+=FT[ql*nphi+qr]*conj(psir_0[i])*psir_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
        }
    }
    return abs(Sf_mean)/sector.nel;
}


void lhamil::save_to_file(const char* filename) {
    if(eigenvalues.size()==0) return;
    ofstream odf;
    odf.open(filename,ofstream::out);
    odf<<"nHilbert:="<<setw(10)<<" "<<nHilbert<<endl;
    odf<<"lambda:="<<setw(10)<<" "<<lambda<<endl;
    odf<<"seed:="<<setw(10)<<" "<<seed<<endl;
    odf<<"eigenvalues:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<overlap.size(); i++)
        odf<<eigenvalues[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"norm:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<=overlap.size(); i++)
        odf<<norm[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"overlap:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<overlap.size(); i++)
        odf<<overlap[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<overlap.size(); i++)
        odf<<psi_0[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_n0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<overlap.size(); i++)
        odf<<psi_n0[i]<<" , ";
    odf<<" ] "<<endl;
    odf.close();
}

void lhamil::read_from_file(const char* filename) {
    ifstream idf;
    idf.open(filename,ifstream::in);
    string buffer;
    idf>>buffer;
    idf>>nHilbert;
    idf>>buffer;
    idf>>lambda;
    idf>>buffer;
    idf>>seed;
    double value;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        eigenvalues.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<=lambda; i++) {
        idf>>value;
        idf>>buffer;
        norm.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        overlap.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        psi_0.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        psi_n0.push_back(value);
    }
    idf.close();
}

void lhamil::print_hamil_CSR()const {
    std::cout << "hamiltonian in CSR format: " << std::endl;
    std::cout << "------------------------------" << std::endl;
    H.print();
}

void lhamil::print_hamil(int range)const {
    int i,k,row_starts,col_index;
    if(range>=nHilbert)
        range=nHilbert;
    for(i=0; i<range; i++) {
        if(i==0)
            cout<<"[[";
        else cout<<" [";
        row_starts=H.outer_starts[i];
        for(k=0; k<range; k++) {
            col_index=H.inner_indices[row_starts];
            //if(k==col_index && row_starts<H.outer_starts[i+1]) {
            if(k==col_index && row_starts<H.outer_size[i]) {
                cout<<setw(4)<<setprecision(2)<<H.value[row_starts]<<",";
                row_starts++;
            }
            else
                cout<<setw(4)<<setprecision(2)<<0<<",";
        }
        if(i==range-1)
            cout<<"]]"<<endl;
        else cout<<"]"<<endl;
    }
}

void lhamil::print_lhamil(int range)const {
    if(range>=overlap.size()) range=overlap.size();
    for(int i=0; i<range; i++) {
        if(i==0)
            cout<<"[[";
        else cout<<" [";
        for(int j=0; j<range; j++) {
            if(j==i+1)
                cout<<norm[j]<<",";
            else if(j==i-1)
                cout<<norm[i]<<",";
            else if(j==i)
                cout<<overlap[i]<<",";
            else
                cout<<0<<",";
        }
        if(i==range-1)
            cout<<"]]"<<endl;
        else cout<<"]"<<endl;
    }
}

void lhamil::print_eigen( int range) const{
    if(range>=overlap.size()) range=overlap.size();
    cout<<"Eigenvalues:= [";
    for(int i=0; i<range; i++)
        cout<<eigenvalues[i]<<", ";
    cout<<",...]"<<endl;
}
