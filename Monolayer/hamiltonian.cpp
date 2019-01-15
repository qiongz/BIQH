#include"hamiltonian.h"

hamil::hamil() {}

double hamil::Coulomb_interaction(int k_x, int k_y) {
    double q=sqrt(k_x*k_x/(lx*lx)+k_y*k_y/(ly*ly))*2.0*M_PI;
    return 2.0*M_PI/(q+1e-30)*exp(-q*q/2.0)*pow(1.0-exp(-q*q/2.0),nLL*2);
}

void hamil::init_Coulomb_matrix() {
    Coulomb_matrix.assign(nphi*nphi, 0);
    for(int s = 0; s < nphi; s++)
        for(int k_y = 0; k_y <nphi; k_y++) {
            double V=0;
            for(int k_x=-nphi/2; k_x<=nphi/2; k_x++)
                if(!(k_y==0 && k_x==0))
                    V+=Coulomb_interaction(k_x,k_y)*cos((2.0*M_PI*s)*k_x/nphi)/(2.0*lx*ly);
            Coulomb_matrix[s*nphi+k_y]=V;
        }
    // initialize classical Coulomb energy
    E_cl=-2.0;
    for(int i=0; i<nphi; i++)
        for(int j=0; j<nphi; j++)
            if(!(i==0 &&j==0))
                E_cl+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
    E_cl/=sqrt(lx*ly);
}

void hamil::set_hamil(basis sector, double _lx, double _ly, int _nphi, int _nLL) {
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    hamiltonian = new complex<double>[nHilbert * nHilbert];
    memset(hamiltonian, 0, sizeof(complex<double>)*nHilbert * nHilbert);
    int sign,signl,signr,ncross;
    long mask, b, n, m, i, k, s,t;
    long k1,k2,lbasis,rbasis;
    long kx=sector.K;
    long Cl,Cr;
          
    for(i = 0; i < nHilbert; i++) {
        // select two electrons in left-basis <m_1, m_2|
        // n=j1, m=j2
        Cl=(kx<0?1:sector.basis_C[i]);
        for(k1=0; k1<Cl; k1++) {
            signl=1;
            lbasis=(k1==0?sector.id[i]:sector.translate(sector.id[i],k1,signl,ncross));
            for(n = 0; n < nphi-1; n++)
                for(m = n+1; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // looking up the corresponding basis in id
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask && n!=m) {
                        // b is the rest electon positions
                        b = lbasis ^ mask;
                        // mt=j3, nt=j4
                        long nt, mt, mask_t, occ_t;
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = -nphi/2; t<=nphi/2; t++) {
                            if(n + t >=nphi)
                                nt = n + t - nphi;
                            else if(n+t<0)
                                nt = n + t + nphi;
                            else
                                nt = n + t;
                            if(m - t <0)
                                mt = m - t + nphi;
                            else if (m-t>=nphi)
                                mt = m - t - nphi;
                            else
                                mt = m - t;
                            s=abs(n-mt);
                            // the translated two electrons indices
                            mask_t = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_t == 0){
                                // determine the subbasis size of right side basis
                                for(k2=1; k2<=sector.C; k2++)
                                  if(sector.translate(mask_t+b,k2,signr,ncross)==(mask_t+b)) {
                                    Cr=k2;
                                    break;
                                  }
                                Cr=(kx<0?1:Cr);
                                for(k2=0; k2<Cr; k2++) {
                                    signr=1;
				    rbasis=(k2==0?mask_t+b:sector.inv_translate(mask_t+b,k2,signr,ncross));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        k = sector.basis_set[rbasis];
                                        sign=sector.get_sign(lbasis,n,m,nt,mt)*signl*signr;
					// if there're crossing between two electrons and no BC-cross
        				if(nt>mt && n<0)
        				  sign*=-1;
                                        complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(k1-k2)/sector.C),sin(2.0*M_PI*kx*(k1-k2)/sector.C))/sqrt(Cl*Cr);
                                        hamiltonian[i*nHilbert+k]+=2.0*Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor;
                                    }
                                }
                            }
                        }
                    }
                }

        }
        hamiltonian[i*nHilbert+i]+=E_cl*sector.nel;
    }
}

void hamil::set_hamil(basis sector, double _lx, double _ly, int _nphi, int _nLL, double theta_x, double theta_y) {
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    hamiltonian = new complex<double>[nHilbert * nHilbert];
    memset(hamiltonian, 0, sizeof(complex<double>)*nHilbert * nHilbert);
    int sign,signl,signr,ncrossl,ncrossr,ncross;
    long mask, b, n, m, i, k, s,t;
    long k1,k2,lbasis,rbasis;
    long kx=sector.K;
    long Cl,Cr,J,c;
    complex<double> FT_twisted,FT_twisted_x,FT_twisted_y; 
    for(i = 0; i < nHilbert; i++) {
        // select two electrons in left-basis <m_1, m_2|
        c=sector.id[i];
    	J=0;
	for(int Ji=0;Ji<nphi;Ji++){
	   mask=(1<<Ji);
	   if((mask&c)==mask)
	   J+=Ji;
	}
        // n=j1, m=j2
        Cl=(kx<0?1:sector.basis_C[i]);
        for(k1=0; k1<Cl; k1++) {
            //signl=1;
            lbasis=(k1==0?sector.id[i]:sector.translate(sector.id[i],k1,signl,ncrossl));
            for(n = 0; n < nphi-1; n++)
                for(m = n+1; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // looking up the corresponding basis in id
                    // if there're two electrons on n and m;
                    if((lbasis &mask) == mask && n!=m) {
                        // b is the rest electon positions
                        b = lbasis ^ mask;
                        // mt=j3, nt=j4
                        long nt, mt, mask_t, occ_t;
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = -nphi/2; t<=nphi/2; t++) {
			    FT_twisted=1;
			    // j2 fixed, and if j4 across the boundary, apply the twisted BC along x-direction
                            if(n + t >=nphi){
                                nt = n + t - nphi;
				FT_twisted*=complex<double>(cos(theta_x),sin(theta_x));
			    }
                            else if(n+t<0){
                                nt = n + t + nphi;
				FT_twisted*=complex<double>(cos(theta_x),sin(-theta_x));
			    }
                            else
                                nt = n + t;
			    // j1 fixed, and if j3 across the boundary, apply the twisted BC along y-direction
                            if(m - t <0){
                                mt = m - t + nphi;
				FT_twisted*=complex<double>(cos(theta_x),sin(-theta_x));
			    }				
                            else if (m-t>=nphi){
                                mt = m - t - nphi;
				FT_twisted*=complex<double>(cos(theta_x),sin(theta_x));
			    }
                            else
                                mt = m - t;
			    // s=j1-j3
                            s=abs(n-mt);
                            // the translated two electrons indices
                            mask_t = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_t == 0){
                                // determine the subbasis size of right side basis
                                for(k2=1; k2<=sector.C; k2++)
                                  if(sector.translate(mask_t+b,k2,signr,ncross)==(mask_t+b)) {
                                    Cr=k2;
                                    break;
                                  }
                                Cr=(kx<0?1:Cr);
                                for(k2=0; k2<Cr; k2++) {
                                    signr=1;
				    rbasis=(k2==0?mask_t+b:sector.inv_translate(mask_t+b,k2,signr,ncrossr));
				    FT_twisted_x=complex<double>(cos(theta_x*(ncrossl-ncrossr)),sin(theta_x*(ncrossl-ncrossr)));
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()) {
                                        k = sector.basis_set[rbasis];
                                        sign=sector.get_sign(lbasis,n,m,nt,mt);
					// if there're crossing between two electrons and no BC-cross
				        // add twisted phase by x-direction translation 
                                        complex<double> FT_factor=FT_twisted*FT_twisted_x*complex<double>(cos(((2.0*M_PI*kx)*(k1-k2)/sector.C)),sin(((2.0*M_PI*kx)*(k1-k2)/sector.C)))/sqrt(Cl*Cr);
                                        hamiltonian[i*nHilbert+k]+=2.0*Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor;
                                    }
                                }
                            }
                        }
                    }
                }

        }
        hamiltonian[i*nHilbert+i]+=E_cl*sector.nel;
    }
}
hamil::~hamil() {}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this != &_gs_hconfig) {
        nHilbert = _gs_hconfig.nHilbert;
        nphi = _gs_hconfig.nphi;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(), _gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(), _gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(), _gs_hconfig.psi_n0.end());
    }
    return *this;
}


double hamil::spectral_function(vector< complex<double> > &O_phi_0, double omega, double _E0, double eta, int annil) {
    complex<double> E;
    complex<double> G = 0;
    for(int i = 0; i < nHilbert; i++)
        // set annil==1, which gives hole-sector
        if(annil == 1) {
            E = complex<double>(omega, eta);
            G += pow(conj(psi_n0[i]) * O_phi_0[i], 2) / (E + eigenvalues[i] - _E0);
        }
    // else particle-sector
        else {
            E = complex<double>(omega, eta);
            G += pow( conj(psi_n0[i]) * O_phi_0[i], 2) / (E + _E0 - eigenvalues[i]);
        }

    return -G.imag() / M_PI;
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
    long i, idx;
    wfs = new complex<double>[nHilbert * nHilbert];
    double *en = new double[nHilbert];
    memset(wfs, 0, sizeof(complex<double>)*nHilbert * nHilbert);
    for(i = 0; i <nHilbert*nHilbert; i++)
        wfs[i]=hamiltonian[i];
    diag_zheev(wfs, en, nHilbert);
    psi_0.assign(nHilbert, 0);
    psi_n0.assign(nHilbert, 0);
    eigenvalues.assign(nHilbert, 0);
    for(i = 0; i < nHilbert; i++) {
        eigenvalues[i] = en[i];
        psi_0[i] = wfs[i];
        psi_n0[i] = wfs[i * nHilbert];
    }
    delete en;
}

void hamil::twist(basis sector,double theta_y){ 
    int J;
    unsigned long mask,c;
    long i,j,k; 
    for(i=0;i<nHilbert;i++){
      // get the total momentum of the basis
      c=sector.id[i];
      J=0;
      for(k=0;k<nphi;k++){
	mask=(1<<k);
	if((mask&c)==mask)
	   J+=k;
      }
      J=J%nphi;
      for(j=0;j<nHilbert;j++)
         wfs[j*nHilbert+i]*=complex<double>(cos(theta_y*J*sector.nel/nphi),sin(theta_y*J*sector.nel/nphi));
    }	    
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
            cout<<setw(5)<<setprecision(2)<<hamiltonian[i*nHilbert+j].real()<<", ";
        if(i == range -1)
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

void hamil::clear(){
    delete hamiltonian,wfs;
    eigenvalues.clear();
    psi_0.clear();
    psi_n0.clear();
}
