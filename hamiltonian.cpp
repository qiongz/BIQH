#include"hamiltonian.h"
hamil::hamil() {}

double hamil::Coulomb_interaction(int alpha, int beta, int q_x, int q_y) {
    double q = sqrt(q_x * q_x / (lx * lx) + q_y * q_y / (ly * ly)) * 2.0 * M_PI;
    if(alpha == beta)
        // symmetric gauge
        //return 2.0*M_PI/ (q + 1e-30) * exp(-M_PI * (q_x * q_x * ly * 1.0 / lx + q_y * q_y * lx * 1.0 / ly) / nphi);
        return 2.0*M_PI/q*exp(-q*q/2.0);
    else
        //return 2.0*M_PI / (q + 1e-30) * exp(-M_PI * (q_x * q_x * ly * 1.0 / lx + q_y * q_y * lx * 1.0 / ly) / nphi - q * d);
        return 2.0*M_PI/q*exp(-q*q/2.0-q*d);
}

void hamil::init_Coulomb_matrix() {
    off_head=nphi/2+nphi%2;
    long dim1 = off_head;
    long dim2 = off_head* dim1;
    long dim3 = dim2* nphi;
    long dim4 = dim3* nphi;
    Coulomb_matrix.assign(2 * dim4, 0);
    for(int alpha = 0; alpha < 2; alpha++)
      for(int n = 0; n < nphi; n++)
        for(int m = 0; m < nphi; m++)
          for(int q_y = 0; q_y < off_head; q_y++)
            for(int q_x = 0; q_x < off_head; q_x++)
              //Coulomb_matrix[alpha * dim4 + n * dim3 + m * dim2 + q_y * dim1 + q_x] = Coulomb_interaction(0, alpha, q_x, q_y) * cos(-2.0 * M_PI * q_x * q_y / nphi + 2.0 * M_PI * (m - n) * q_x / nphi)/(2.0*lx*ly);
              Coulomb_matrix[alpha * dim4 + n * dim3 + m* dim2 + q_y * dim1 + q_x] = Coulomb_interaction(0, alpha, q_x, q_y) * cos(2.0 * M_PI * (n-m) * q_x / nphi)/(2.0*lx*ly);
}

void hamil::set_hamil(basis & sector, double _lx, double _ly, int _nphi, double _d) {
    long nbasis_up, nbasis_down;
    d = _d;
    nsite = sector.nsite;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    off_head=nphi/2+nphi%2+1;
    init_Coulomb_matrix();
    nbasis_up = sector.nbasis_up;
    nbasis_down = sector.nbasis_down;
    nHilbert = nbasis_up * nbasis_down;
    vector<double> matrix_elements;
    H.clear();
    H.inner_indices.reserve(nHilbert * nsite);
    H.value.reserve(nHilbert * nsite);
    H.outer_starts.reserve(nHilbert + 1);
    long mask, mask_u, mask_d, b, p, n, m, i, j, k, l, t, nsignu, nsignd, nsign;
    long row = 0;
    H.outer_starts.push_back(0);
    for(i = 0; i < nbasis_up; i++) {
        for(j = 0; j < nbasis_down; j++) {
            // start of new row of nonzero elements
            matrix_elements.assign(nHilbert,0);
            // select two electrons in left-basis <m_1, m_2|
            long dim1= off_head;
            long dim2 = dim1* off_head;
            long dim3 = dim2* nphi;
            long dim4 = dim3* nphi;
            for(n = 0; n < nphi; n++)
                for(m = 0; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // consider the upper-layer two electrons
                    // looking up the corresponding basis in id_up
                    // if there're two electrons on n and m;
                    if((sector.id_up[i]&mask) == mask && m != n) {
                        // b is the rest electon positions
                        b = sector.id_up[i] ^ mask;
                        //cout<<"#i n m mask b=i^mask b+mask"<<endl;
                        //cout<<bitset<4>(sector.id_up[i]).to_string()<<" "<<bitset<4>(1<<n).to_string()<<" "<<bitset<4>(1<<m).to_string()<<" "<<bitset<4>(mask).to_string()<<" "<<bitset<4>(b).to_string()<<" "<<bitset<4>(mask+b)<<" "<<endl;
                        long nt, mt, mask_ut, occ_ut;
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
                            mask_ut = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_ut = mask_ut & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_ut == 0 && sector.basis_up.find(mask_ut + b) != sector.basis_up.end())
                            {
                                k = sector.basis_up[mask_ut + b];
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
                                   nsign+=sector.nel_up+sector.nel_up-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel_up-1;

                                matrix_elements[k*nbasis_down+j]+=V_uu*pow(-1,nsign);
                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_ut == mask_ut)
                                nsignu += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_ut != 0 && occ_ut != mask_ut)
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
                            mask_ut = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_ut = mask_ut & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            // looking up Lin's table, and find the corresponding index
                            if(occ_ut == 0 && sector.basis_up.find(mask_ut + b) != sector.basis_up.end())
                            {
                                k = sector.basis_up[mask_ut + b];
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
                                  nsign+=sector.nel_up+sector.nel_up-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel_up-1;
                                matrix_elements[k*nbasis_down+j]+=V_uu*pow(-1,nsign);

                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_ut == mask_ut)
                                nsignu += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_ut != 0 && occ_ut != mask_ut)
                                nsignu++;
                        }
                    }

                    // consider the lower-layer two electrons
                    // if there're two electrons on n and m;
                    if((sector.id_down[j]&mask) == mask && m != n) {
                        // b is the rest electon positions
                        b = sector.id_down[j] ^ mask;
                        long nt, mt, mask_dt, occ_dt;
                        nsignd = 0;
                        bool ncross=false;
                        bool mcross=false;
                        // perform translation in x-direction, negative q_y
                        for(t = 0; t < off_head; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n - t < 0) {
                                nt = n - t + nphi;
                                ncross=true;
                            }
                            else
                                nt = n - t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m + t >= nphi) {
                                mt = m + t - nphi;
                                mcross=true;
                            }
                            else
                                mt = m + t;
                            // the translated two electrons indices
                            mask_dt = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_dt = mask_dt & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_dt == 0 && sector.basis_down.find(mask_dt + b) != sector.basis_down.end()) {
                                l = sector.basis_down[mask_dt + b];
                                long q_x, q_y;
                                // t=-off_head corresponds to q_y= -Pi
                                // t= off_head corresponds to q_y= +Pi
                                q_y = t;
                                double V_dd = 0;
                                // only positive part of q_x is added, q_x-> -q_x reflection gives cosine term
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_dd += Coulomb_matrix[n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_down+sector.nel_down-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel_down-1;
                                matrix_elements[i*nbasis_down+l]+=V_dd*pow(-1,nsign);

                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_dt == mask_dt)
                                nsignd += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_dt != 0 && occ_dt != mask_dt)
                                nsignd++;
                        }
                        nsignd=0;
                        ncross=false;
                        mcross=false;
                        // positive q_y
                        for(t = 1; t < off_head; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n+t>=nphi){
                                nt = n+t -nphi;
                                ncross=true;
                            }
                            else
                                nt = n + t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m - t <0) {
                                mt = m - t + nphi;
                                mcross=true;
                            }
                            else
                                mt = m - t;
                            // the translated two electrons indices
                            mask_dt = (1 << nt) + (1 << mt);
                            // occupation of electons on the translated position
                            occ_dt = mask_dt & b;
                            // if there're no electon on the translated position
                            // which is a valid translation, can be applied
                            if(occ_dt == 0 && sector.basis_down.find(mask_dt + b) != sector.basis_down.end()) {
                                l = sector.basis_down[mask_dt + b];
                                long q_x, q_y;
                                // t=-off_head corresponds to q_y= -Pi
                                // t= off_head corresponds to q_y= +Pi
                                q_y = t;
                                double V_dd = 0;
                                // only positive part of q_x is added, q_x-> -q_x reflection gives cosine term
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_dd += Coulomb_matrix[n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_down+sector.nel_down-2;
                                else if(ncross && !mcross || !ncross && mcross)
                                   nsign+=sector.nel_down-1;

                                matrix_elements[i*nbasis_down+l]+=V_dd*pow(-1,nsign);

                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_dt == mask_dt)
                                nsignd += 2;
                            // one electron is occupied, and to be crossed next
                            else if(occ_dt != 0 && occ_dt != mask_dt)
                                nsignd++;
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
                        nsignu = 0;
                        nsignd = 0;
                        bool ncross=false;
                        bool mcross=false;
                        // perform translation along x-direction
                        for(t = 0; t < off_head ; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if(n - t < 0) {
                                ncross= true;
                                nt = n - t + nphi;
                            }
                            else
                                nt = n - t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if(m + t >= nphi) {
                                mt = m + t - nphi;
                                mcross=true;
                            }
                            else
                                mt = m + t;
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
                                long q_x, q_y;
                                q_y = t;
                                double V_ud = 0;
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_ud += Coulomb_matrix[dim4 + n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignu+nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_up+sector.nel_down-2;
                                else if(ncross && !mcross)
                                   nsign+=sector.nel_up-1;
                                else if(!ncross && mcross)
                                   nsign+=sector.nel_down-1;

                                matrix_elements[k*nbasis_down+l]+=V_ud*pow(-1,nsign);

                            }
                            // two electrons are occupied, and to be crossed next
                            else if(occ_ut == mask_ut && occ_dt == mask_dt){
                                nsignu ++;
                                nsignd ++;
                            }
                            else if(occ_ut == 0 && occ_dt == mask_dt && t!=0)
                                nsignd++;
                            // one electron is occupied, and to be crossed next
                            else if(occ_dt == 0 && occ_ut == mask_ut && t!=0)
                                nsignu++;

                        }

                        nsignu = 0;
                        nsignd = 0;
                        ncross=false;
                        mcross=false;
                        // perform translation along x-direction, positive q_y
                        for(t = 1; t < off_head ; t++) {
                            // PBC, if one electron cross left boundary, sign change with -1
                            if (n+t>=nphi){
                                nt = n+t -nphi;
                                ncross=true;
                            }
                            else
                                nt = n + t;
                            // PBC, if one electron cross right boundary, sign change with -1
                            if (m-t <0){
                                mt = m-t +nphi;
                                mcross=true;
                            }
                            else
                                mt = m - t;
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
                                long q_x, q_y;
                                q_y = t;
                                double V_ud = 0;
                                for(q_x = 0; q_x < off_head; q_x++)
                                    if(!(q_x==0 && q_y==0))
                                        V_ud += Coulomb_matrix[dim4 + n * dim3 + mt* dim2 + q_y * dim1 + q_x];

                                nsign=nsignu+nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_up+sector.nel_down-2;
                                else if(ncross && !mcross)
                                   nsign+=sector.nel_up-1;
                                else if(!ncross && mcross)
                                   nsign+=sector.nel_down-1;

                                matrix_elements[k*nbasis_down+l]+=V_ud*pow(-1,nsign);
                            }
                            // both electrons are occupied, and to be crossed next
                            else if(occ_ut == mask_ut && occ_dt == mask_dt){
                                nsignu ++;
                                nsignd ++;
                            }
                            else if(occ_ut == 0 && occ_dt == mask_dt ){
                                nsignd++;
                            }
                            // one electron is occupied, and to be crossed next
                            else if(occ_dt == 0 && occ_ut == mask_ut){
                                nsignu++;
                            }
                        }
                    }
                }
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
    }
    matrix_elements.clear();
}

hamil::~hamil() {}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this != &_gs_hconfig) {
        seed = _gs_hconfig.seed;
        nHilbert = _gs_hconfig.nHilbert;
        H = _gs_hconfig.H;
        d = _gs_hconfig.d;
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
    return E_gs / nsite;
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
