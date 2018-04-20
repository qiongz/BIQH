#include"hamiltonian.h"
hamil::hamil() {}

double hamil::Coulomb_interaction(int alpha, int beta, int q_x, int q_y) {
    double q = sqrt(q_x * q_x / (lx * lx) + q_y * q_y / (ly * ly)) * 2.0 * M_PI;
    if(alpha == beta)
        return 2.0*M_PI/q*exp(-q*q/2.0);
    else
        return 2.0*M_PI/q*exp(-q*q/2.0-q*d);
}

void hamil::init_Coulomb_matrix() {
    Coulomb_matrix.assign(2 * nphi*nphi, 0);
    for(int alpha = 0; alpha < 2; alpha++)
      // n=j_1, m=_j3
      for(int s = 0; s < nphi; s++)
          for(int q_y = 0; q_y < nphi; q_y++){
            double V=0;
            for(int q_x = -nphi/2; q_x <=nphi/2; q_x++)
                if(!(q_x==0 && q_y==0))
                  V+=Coulomb_interaction(0,alpha,q_x,q_y)*cos(2.0*M_PI*s*q_x/nphi)/(2.0*lx*ly);
              // Coulomb matrix elements in Landau gauge
              Coulomb_matrix[alpha*nphi*nphi+s*nphi+q_y]=V;
            }
    // initialize classical Coulomb energy
    Ec=-2.0;
    for(int i=0;i<nphi;i++)
      for(int j=0;j<nphi;j++)
        if(!(i==0 &&j==0))
        Ec+=Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
    Ec/=sqrt(lx*ly);
}

void hamil::set_hamil(basis & sector, double _lx, double _ly, int _nphi, double _d)
{
    long nbasis_up, nbasis_down;
    d = _d;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    init_Coulomb_matrix();
    nbasis_up = sector.nbasis_up;
    nbasis_down = sector.nbasis_down;
    nHilbert = nbasis_up * nbasis_down;
    hamiltonian.assign(nHilbert*nHilbert,0);
    long mask, mask_u, mask_d, b, p, n, m, i, j, k, l, s,t,sign;
    for(i = 0; i < nbasis_up; i++)
        for(j = 0; j < nbasis_down; j++) {
            // select two electrons in left-basis <m_1, m_2|
            // n=j1, m=j2
            for(n = 0; n < nphi; n++)
                for(m = 0; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((sector.id_up[i]&mask) == mask && n!=m) {
                        // b is the rest electon positions
                        b = sector.id_up[i] ^ mask;
                        // mt=j3, nt=j4
                        long nt, mt, mask_ut, occ_ut;
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
                            mask_ut = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_ut = mask_ut & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_ut == 0 && sector.basis_up.find(mask_ut + b) != sector.basis_up.end())
                            {
                                k = sector.basis_up[mask_ut + b];
                                sign=sector.get_signu(i,n,m,nt,mt);
                                hamiltonian[(i*nbasis_down+j)*nHilbert+k*nbasis_down+j]+=Coulomb_matrix[s*nphi+abs(t)]*sign;
                            }
                        }
                    }
                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((sector.id_down[j]&mask) == mask && m!=n) {
                        // b is the rest electon positions
                        b = sector.id_down[j] ^ mask;
                        long nt, mt, mask_dt, occ_dt;
                        // perform translation in x-direction, negative q_y
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
                            mask_dt = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_dt = mask_dt & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_dt == 0 && sector.basis_down.find(mask_dt + b) != sector.basis_down.end()) {
                                l = sector.basis_down[mask_dt + b];
                                sign=sector.get_signd(j,n,m,nt,mt);
                                hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+l]+=Coulomb_matrix[s*nphi+abs(t)]*sign;
                            }
                        }
                    }
                    // consider the one electron in the upper layer
                    // and one electron in the lower layer case
                    mask_u = (1 << n);
                    mask_d = (1 << m);
                    // if there is one electron at site n in upper-layer
                    // and one electron at site m in lower-layer
                    if((sector.id_up[i]&mask_u) == mask_u && (sector.id_down[j]&mask_d) == mask_d) {
                        // b is the rest electon positions for upper-layer electrons
                        b = sector.id_up[i] ^ mask_u;
                        p = sector.id_down[j] ^ mask_d;
                        long nt, mt, mask_ut, occ_ut, mask_dt, occ_dt;
                        // perform translation along x-direction
                        for(t = -nphi/2; t <= nphi/2 ; t++) {
                            if(n + t>=nphi)
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
                            // the translated upper electron index
                            mask_ut = (1 << nt);
                            mask_dt = (1 << mt);
                            // occupation of electons on the translated position
                            occ_ut = mask_ut & b;
                            occ_dt = mask_dt & p;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // the translated indices
                            if(occ_ut == 0 && occ_dt == 0 && sector.basis_up.find(mask_ut + b) != sector.basis_up.end() && sector.basis_down.find(mask_dt + p) != sector.basis_down.end()) {
                                k = sector.basis_up[mask_ut + b];
                                l = sector.basis_down[mask_dt + p];
                                sign=sector.get_signud(i,j,n,m,nt,mt);
                                hamiltonian[(i*nbasis_down+j)*nHilbert+k*nbasis_down+l]+=Coulomb_matrix[nphi*nphi+s*nphi+abs(t)]*sign;
                            }
                        }
                    }
          }
          // diagonal Coulomb classical energy term
          hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=Ec*(sector.nel_up+sector.nel_down);
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
    double E_gs = 0;
    vector<double> psi_t;
    psi_t.assign(nHilbert,0);
    for(int i=0;i<nHilbert;i++)
       for(int j=0;j<nHilbert;j++)
        psi_t[i]+=hamiltonian[i*nHilbert+j]*psi_0[j];
    for(int i = 0; i < nHilbert; i++)
        E_gs += psi_t[i] * psi_0[i];
    psi_t.clear();
    return E_gs;
}

void hamil::diag() {
    int i;
    double *h = new double[nHilbert * nHilbert];
    double *en = new double[nHilbert];
    memset(h, 0, sizeof(double)*nHilbert * nHilbert);
    for(i = 0; i <nHilbert*nHilbert; i++)
        h[i]=hamiltonian[i];
    diag_dsyev(h, en, nHilbert);
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
