#include"hamiltonian.h"

hamil::hamil() {}

double hamil::Coulomb_interaction(int q_x, int q_y) {
    double q = sqrt(q_x * q_x / (lx * lx) + q_y * q_y / (ly * ly)) * 2.0 * M_PI;
    return 2.0*M_PI/q*exp(-q*q/2.0);
}

void hamil::init_Coulomb_matrix() {
    off_head=nphi/2+nphi%2;
    long dim1 = off_head;
    long dim2 = off_head* dim1;
    long dim3 = dim2* nphi;
    Coulomb_matrix.assign(dim3*nphi, 0);
    for(int n = 0; n < nphi; n++)
        for(int m = 0; m < nphi; m++)
          for(int q_y = 0; q_y < off_head; q_y++)
            for(int q_x = 0; q_x < off_head; q_x++)
              //Coulomb_matrix[alpha * dim4 + n * dim3 + m * dim2 + q_y * dim1 + q_x] = Coulomb_interaction(0, alpha, q_x, q_y) * cos(-2.0 * M_PI * q_x * q_y / nphi + 2.0 * M_PI * (m - n) * q_x / nphi)/(2.0*lx*ly);
              Coulomb_matrix[ n * dim3 + m* dim2 + q_y * dim1 + q_x] = Coulomb_interaction(q_x, q_y) * cos(2.0 * M_PI * (n-m) * q_x / nphi)/(2.0*lx*ly);

    // initialize classical Coulomb energy
    E_cl=-2.0/sqrt(2.0*M_PI*nphi);
    for(int i=0;i<nphi;i++)
      for(int j=0;j<nphi;j++)
        if(!(i==0 &&j==0))
        E_cl+=1.0/sqrt(2.0*M_PI*nphi)*Integrate_ExpInt((i*i*lx/ly+j*j*ly/lx)*M_PI);
}

void hamil::set_hamil(basis & sector, double _lx, double _ly, int _nphi) {
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    off_head=nphi/2+nphi%2+1;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    vector<double> matrix_elements;
    H.clear();
    H.inner_indices.reserve(nHilbert * nphi);
    H.value.reserve(nHilbert * nphi);
    H.outer_starts.reserve(nHilbert + 1);
    long mask, b, p, n, m, i, j, k, l, t, nsignu, nsign;
    long row = 0;
    H.outer_starts.push_back(0);
    for(i = 0; i < nHilbert; i++){
            // start of new row of nonzero elements
            matrix_elements.assign(nHilbert,0);
            // select two electrons in left-basis <m_1, m_2|
            long dim1= off_head;
            long dim2 = dim1* off_head;
            long dim3 = dim2* nphi;
            // n=j1, m=j2
            for(n = 0; n < nphi; n++)
                for(m = 0; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((sector.id[i]&mask) == mask && m != n) {
                        // b is the rest electon positions
                        b = sector.id[i] ^ mask;
                        // mt=j3, nt=j4
                        long nt, mt, mask_t, occ_t;
                        nsignu = 0;
                        bool ncross=false;
                        bool mcross=false;
                        // perform translation along x-direction (q_y), positive q_y
                        for(t = 0; t < off_head; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n - t < 0) {
                                nt = n - t + nphi;
                                // crossing boundary= crossing all rest electrons
                                ncross=true;
                            }
                            else
                                nt = n - t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m + t >= nphi) {
                                mt = m + t - nphi;
                                // crossing boundary= crossing all rest electrons
                                mcross=true;
                            }
                            else
                                mt = m + t;
                            // the translated two electrons indices
                            mask_t = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_t == 0 && sector.basis_set.find(mask_t + b) != sector.basis_set.end())
                            {
                                k = sector.basis_set[mask_t + b];
                                long q_x, q_y;
                                // t=-off_head corresponds to q_y= -Pi
                                // t= off_head corresponds to q_y= +Pi
                                q_y = t;
                                double V_uu = 0;
                                // only positive part of q_x is added, q_x-> -q_x reflection gives cosine term
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_uu += Coulomb_matrix[n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignu;
                                if(ncross&& mcross)
                                   nsign+=sector.nel+sector.nel-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel-1;

                                matrix_elements[k]+=V_uu*pow(-1,nsign);
                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_t == mask_t)
                                nsignu += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_t != 0 && occ_t != mask_t)
                                nsignu++;
                        }

                        nsignu = 0;
                        ncross=false;
                        mcross=false;
                        // perform translation along x-direction (q_y), negative q_y
                        for(t = 1; t < off_head; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if (n + t >= nphi){
                                nt = n +t -nphi;
                                // crossing boundary= crossing all rest electrons
                                ncross=true;
                            }
                            else
                                nt = n + t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if (m-t <0 ){
                                mt =m -t +nphi;
                                // crossing boundary= crossing all rest electrons
                                mcross=true;
                            }
                            else
                                mt = m - t;
                            // the translated two electrons indices
                            mask_t = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_t = mask_t & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_t == 0 && sector.basis_set.find(mask_t + b) != sector.basis_set.end())
                            {
                                k = sector.basis_set[mask_t + b];
                                long q_x, q_y;
                                // t=-off_head corresponds to q_y= -Pi
                                // t= off_head corresponds to q_y= +Pi
                                q_y = t;
                                double V_uu = 0;
                                // only positive part of q_x is added, q_x-> -q_x reflection gives cosine term
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_uu += Coulomb_matrix[n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignu;
                                if(ncross&& mcross)
                                  nsign+=sector.nel+sector.nel-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel-1;
                                matrix_elements[k]+=V_uu*pow(-1,nsign);

                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_t == mask_t)
                                nsignu += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_t != 0 && occ_t != mask_t)
                                nsignu++;
                        }
                    }
                }
            // diagonal Coulomb classical energy term
           matrix_elements[i]+=E_cl*(sector.nel);

            long count=0;
            for(k=0;k<nHilbert;k++)
               if(abs(matrix_elements[k])>1e-6){
                 H.inner_indices.push_back(k);
                 H.value.push_back(matrix_elements[k]);
                 count++;
               }
            row += count;
            H.outer_starts.push_back(row);
          }
    matrix_elements.clear();
}

hamil::~hamil() {}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this != &_gs_hconfig) {
        seed = _gs_hconfig.seed;
        nHilbert = _gs_hconfig.nHilbert;
        H = _gs_hconfig.H;
        nphi = _gs_hconfig.nphi;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(), _gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(), _gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(), _gs_hconfig.psi_n0.end());
    }
    return *this;
}


double hamil::spectral_function(vector<double > &O_phi_0, double omega, double _E0, double eta, int annil) {
    complex<double> E;
    complex<double> G = 0;
    for(int i = 0; i < nHilbert; i++)
        // set annil==1, which gives hole-sector
        if(annil == 1) {
            E = complex<double>(omega, eta);
            G += pow(psi_n0[i] * O_phi_0[i], 2) / (E + eigenvalues[i] - _E0);
        }
    // else particle-sector
        else {
            E = complex<double>(omega, eta);
            G += pow(psi_n0[i] * O_phi_0[i], 2) / (E + _E0 - eigenvalues[i]);
        }

    return -G.imag() / M_PI;
}

double hamil::ground_state_energy() {
    if(psi_0.size() == 0) return 0;
    double E_gs = 0;
    vector<double> psi_t;
    psi_t = H * psi_0;
    for(int i = 0; i < nHilbert; i++)
        E_gs += psi_t[i] * psi_0[i];
    return E_gs;
}

void hamil::diag() {
    int i, idx;
    double *hamiltonian = new double[nHilbert * nHilbert];
    double *en = new double[nHilbert];
    memset(hamiltonian, 0, sizeof(double)*nHilbert * nHilbert);
    for(i = 0; i < H.outer_starts.size() - 1; i++)
        for(idx = H.outer_starts[i]; idx < H.outer_starts[i + 1]; idx++)
            hamiltonian[i * nHilbert + H.inner_indices[idx]] = H.value[idx];
    diag_dsyev(hamiltonian, en, nHilbert);
    psi_0.assign(nHilbert, 0);
    psi_n0.assign(nHilbert, 0);
    eigenvalues.assign(nHilbert, 0);
    for(i = 0; i < nHilbert; i++) {
        eigenvalues[i] = en[i];
        psi_0[i] = hamiltonian[i];
        psi_n0[i] = hamiltonian[i * nHilbert];
    }
    delete hamiltonian, en;
}


void hamil::print_hamil_CSR() {
    std::cout << "hamiltonian in CSR format: " << std::endl;
    std::cout << "------------------------------" << std::endl;
    H.print();
}

void hamil::print_hamil() {
    int i, j, count;
    for(i = 0; i < nHilbert; i++) {
        if(i == 0)
            cout <<setw(2)<< "[[";
        else cout <<setw(2)<< " [";
        count = 0;
        for(j = 0; j < nHilbert; j++) {
            if(j == H.inner_indices[H.outer_starts[i] + count])
                cout << setprecision(1) <<setw(7)<< H.value[H.outer_starts[i] + count++] << ",";
            else
                cout << setw(7)<<0 << ",";
        }
        if(i == nHilbert - 1)
            cout << "]]" << endl;
        else cout << "]" << endl;
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
        count = 0;
        for(j = 0; j < range; j++) {
            if(j == H.inner_indices[H.outer_starts[i] + count])
                cout << setprecision(1) <<setw(7)<< H.value[H.outer_starts[i] + count++] << ",";
            else
                cout << setw(7)<<0 << ",";
        }
        if(i == range - 1)
            cout << "]]" << endl;
        else cout << "]" << endl;
    }
}

void hamil::print_eigen() {
    std::cout << "Eigenvalues:=[ ";
    for(int i = 0; i < nHilbert; i++)
        if(i != nHilbert - 1)
            std::cout << eigenvalues[i] << ", ";
        else
            std::cout << eigenvalues[i] << " ]" << std::endl;
}
