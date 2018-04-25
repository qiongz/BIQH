#include"hamiltonian.h"

hamil::hamil() {}

double hamil::Coulomb_interaction(int k_x, int k_y) {
    double q=sqrt(k_x*k_x/(lx*lx)+k_y*k_y/(ly*ly))*2.0*M_PI;
    return 2.0*M_PI/(q+1e-30)*exp(-q*q/2.0);
}

void hamil::init_Coulomb_matrix() {
    Coulomb_matrix.assign(nphi*nphi, 0);
    for(int s = 0; s < nphi; s++)
      for(int k_y = 0; k_y <nphi; k_y++){
        double V=0;
        for(int k_x=-nphi/2;k_x<=nphi/2;k_x++)
          if(!(k_y==0 && k_x==0))
            V+=Coulomb_interaction(k_x,k_y)*cos(2.0*M_PI*s*k_x/nphi)/(2.0*lx*ly);
        Coulomb_matrix[s*nphi+k_y]=V;
      }
    // initialize classical Coulomb energy
    E_cl=-2.0;
    for(int i=0;i<nphi;i++)
      for(int j=0;j<nphi;j++)
        if(!(i==0 &&j==0))
        E_cl+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
    E_cl/=sqrt(lx*ly);
}

void hamil::set_hamil(basis & sector, double _lx, double _ly, int _nphi) {
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    long kx=sector.K;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    hamiltonian = new complex<double>[nHilbert * nHilbert];
    memset(hamiltonian, 0, sizeof(complex<double>)*nHilbert * nHilbert);
    long mask, b, n, m, i, k, s,t, sign;
    for(i = 0; i < nHilbert; i++){
            // k-space basis is a linear combination of translated basis subset
            // sum over left-basis after translation (k1 index) and right-basis after translation (k2 index)
            long k1,k2;
            for(k1=0;k1<sector.C;k1++){
            long left_basis=sector.translate(sector.id[i],k1);
            // select two electrons in left-basis <m_1, m_2|
            // n=j1, m=j2
            for(n = 0; n < nphi-1; n++)
                for(m = n+1; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // looking up the corresponding basis in id
                    // if there're two electrons on n and m;
                    if((left_basis &mask) == mask && n!=m) {
                        // b is the rest electon positions
                        b = left_basis ^ mask;
                        // mt=j3, nt=j4
                        long nt, mt, mask_t, occ_t;
                        // perform translation along x-direction (q_y), positive q_y
                        for(t =-nphi/2 ; t <= nphi/2; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n + t >=nphi)
                                nt = n + t - nphi;
                            else if(n+t<0)
                                nt = n +t+nphi;
                            else
                                nt = n + t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m - t <0)
                                mt = m - t + nphi;
                            else if(m-t>=nphi)
                                mt =m-t-nphi;
                            else
                                mt = m - t;
                            s=abs(mt-n);
                            // the translated two electrons indices
                            mask_t = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_t == 0)
                             for(k2=0;k2<sector.C;k2++){
                                long right_basis=sector.inv_translate(mask_t+b,k2);
                                if(sector.basis_set.find(right_basis) != sector.basis_set.end())
                                {
                                k = sector.basis_set[right_basis];
                                sign=sector.get_sign(i,n,m,nt,mt);
                                complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(k1-k2)/sector.C),sin(2.0*M_PI*kx*(k1-k2)/sector.C))/sector.C;
                                hamiltonian[i*nHilbert+k]+=2.0*Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor;
                               }
                           }
                        }
                    }
                }

           }
            // diagonal Coulomb classical energy term
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

double hamil::ground_state_energy(){
    if(psi_0.size() == 0) return 0;
    complex<double> E_gs = 0;
    vector< complex<double> > psi_t;
    psi_t.assign(nHilbert,0);
    for(int i=0;i<nHilbert;i++)
       for(int j=0;j<nHilbert;j++)
        psi_t[i]+=hamiltonian[i*nHilbert+j]*psi_0[j];
    for(int i = 0; i < nHilbert; i++)
        E_gs += conj(psi_t[i]) * psi_0[i];
    psi_t.clear();
    return E_gs.real();
}

void hamil::diag() {
    int i, idx;
    complex<double> *h = new complex<double>[nHilbert * nHilbert];
    double *en = new double[nHilbert];
    memset(h, 0, sizeof(complex<double>)*nHilbert * nHilbert);
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
        for(j=0;j<range;j++)
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
