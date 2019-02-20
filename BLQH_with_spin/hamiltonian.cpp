#include"hamiltonian.h"
static std::mutex mutex_update;
hamil::hamil() {}

double hamil::Coulomb_interaction(const int &alpha,const int &q_x,const  int &q_y) {
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
                for(int q_x = -nphi; q_x <=nphi; q_x++)
                    if(!(q_x==0 && q_y==0))
                        V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);

                if(alpha==1) {
                    V=0;
                    for(int q_x = -10*nphi/(d+0.01); q_x <10*nphi/(d+0.01); q_x++)
                        if(!(q_x==0 &&q_y==0))
                            V+=2.0*Coulomb_interaction(alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);
                }
                // Coulomb matrix elements in Landau gauge
                Coulomb_matrix[alpha*nphi*nphi+s*nphi+q_y]=V;
            }
    // store FT coefficients
    int kx=sector.K;
    FT.assign(nphi*nphi,0);
    for(int kl=0; kl<nphi; kl++)
        for(int kr=0; kr<nphi; kr++)
            FT[kl*nphi+kr]=complex<double>(cos(2.0*M_PI*(kl-kr)*kx/nphi),sin(2.0*M_PI*(kl-kr)*kx/nphi));
}

inline void hamil::peer_set_hamil(const double &Delta_SAS,const double &Delta_V,const double &Delta_Z,const int &id, const long &nbatch,const long &nrange) {
    int kx=sector.K;
    unsigned long lbasis,rbasis,rbasis_0,mask,mask_t,occ_t,b;
    int n,m,s,t,nt,mt,sign,signl,signr;
    int ql,qr,Dl,Dr,D;

    long i,j,k,l;
    vector<complex<double> > matrix_elements;
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
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
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
                                        matrix_elements[j]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
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
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
                //spin-up down-up layer tunneling
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
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
		
                //spin-down up-down layer tunneling
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
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,n,nt)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
                            }
                        }
                    }
                }
                //spin-down down-up layer tunneling
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
                        Dr=(kx<0?1:Dr);
                        for(qr=0; qr<Dr; qr++) {
                            signr=1;
                            rbasis=(qr==0?rbasis_0:sector.inv_translate(rbasis_0,qr,signr));
                            if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                j = sector.basis_set.at(rbasis);
                                sign=sector.get_sign(lbasis,nt,n)*signl*signr;
                                matrix_elements[j]+=-0.25*Delta_SAS*sign*FT[ql*nphi+qr]/sqrt(Dl*Dr);
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

        mutex_update.lock();
        for(k=0; k<nHilbert; k++)
	    hamiltonian[i*nHilbert+k]=matrix_elements[k]; 
        mutex_update.unlock();
    }
    matrix_elements.clear();
}

void hamil::set_hamil(const double &_lx,const double & _ly, const long &_nphi, const long &_nLL,const double& _d, const double &Delta_SAS,const double &Delta_V,const double &Delta_Z,const int &nthread) {
    d = _d;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    eigenvalues=new double[nHilbert]; 
    hamiltonian=new complex<double>[nHilbert*nHilbert];
    memset(hamiltonian,0,nHilbert*nHilbert*sizeof(complex<double>));
    std::vector<std::thread> threads;
    long  nbatch=nHilbert/nthread;
    long nresidual=nHilbert%nthread;
    for(int id = 0; id < nthread; id++)
        threads.push_back(std::thread(&hamil::peer_set_hamil,this,Delta_SAS,Delta_V,Delta_Z,id,nbatch,nbatch));
    for(auto &th:threads)
        if(th.joinable())
            th.join();
    if(nresidual!=0)
        peer_set_hamil(Delta_SAS,Delta_V,Delta_Z,nthread,nbatch,nresidual);
}

hamil::~hamil() {}


double hamil::ground_state_energy() const {
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

double hamil::pseudospin_Sz() const{
    double nel_up=0;
    for(long i=0; i<nHilbert; i++)
	nel_up+=(sector.get_nel(0,i)+sector.get_nel(1,i))*std::norm(psi_0[i]);
    return (2.0*nel_up-sector.nel)/(2.0*sector.nel);
}

double hamil::Sz()const{
    double nel_up=0;
    for(long i=0; i<nHilbert; i++)
	nel_up+=(sector.get_nel(0,i)+sector.get_nel(2,i))*std::norm(psi_0[i]);
    return (2.0*nel_up-sector.nel)/(2.0*sector.nel);
}

double hamil::pseudospin_Sx()const{
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
                //nt=(n<nphi?n+nphi:n-nphi);
		nt=n+2*nphi;
                mask = (1UL << n)+(1UL <<nt);
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
                            Sx_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            // spin-down interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+2*nphi;
                mask = (1UL << n)+(1UL <<nt);
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
                            Sx_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }

            // spin-up spin-down interlayer tunneling 
            for(n=0; n<nphi; n++) {
		nt=n+3*nphi;
                mask = (1UL << n)+(1UL <<nt);
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
                            Sx_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            
            // spin-down spin-up interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+nphi;
                mask = (1UL << n)+(1UL <<nt);
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
                            Sx_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
        }
    }
    return abs(Sx_mean)/sector.nel;
}

double hamil::spinflip_tunneling() const{
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
                mask = (1UL << n)+(1UL <<nt);
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
                            Sf_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
            
            // spin-down spin-up interlayer tunneling 
            for(n=nphi; n<2*nphi; n++) {
		nt=n+nphi;
                mask = (1UL << n)+(1UL <<nt);
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
                            Sf_mean+=FT[ql*nphi+qr]*conj(psi_0[i])*psi_0[j]*double(sign*0.5/sqrt(Dl*Dr));
                        }
                    }
                }
            }
        }
    }
    return abs(Sf_mean)/sector.nel;
}

void hamil::diag() {
    complex<double>* H=new complex<double>[nHilbert*nHilbert];
    memcpy(H,hamiltonian,nHilbert*nHilbert*sizeof(complex<double>));
    diag_zheevd(H,eigenvalues,nHilbert);
    psi_0.assign(nHilbert, 0);
    psi_1.assign(nHilbert, 0);
    psi_n0.assign(nHilbert, 0);
    for(int i = 0; i < nHilbert; i++) {
        psi_0[i] = H[i];
        psi_1[i] = H[i+nHilbert];
        psi_n0[i] = H[i * nHilbert];
    }
    delete H;
}

void hamil::print_hamil(int &range) const{
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

void hamil::print_eigen(int &range) const{
    if(range>=nHilbert)
        range=nHilbert;
    std::cout << "Eigenvalues:=[ ";
    for(int i = 0; i < range; i++)
        if(i != range - 1)
            std::cout << eigenvalues[i] << ", ";
        else
            std::cout << eigenvalues[i] << " , ...]" << std::endl;
}
