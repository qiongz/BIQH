#include"lanczos_hamiltonian.h"
inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,long _nHilbert,long _lambda, unsigned _seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(basis &_sector,double _d, long _lambda,unsigned _seed) {
    sector=_sector;
    lambda=_lambda;
    d=_d;
    seed=_seed;
}

lhamil::~lhamil() {
}

void lhamil::init(basis &_sector,double _d, long _lambda,unsigned _seed) {
    sector=_sector;
    lambda=_lambda;
    seed=_seed;
    d=_d;
}
const lhamil & lhamil::operator =(const lhamil & _config) {
    if(this !=&_config) {
        nHilbert=_config.nHilbert;
        sector=_config.sector;
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

double lhamil::Coulomb_interaction(int alpha,int beta, int q_x, int q_y){
    double q=sqrt(q_x*q_x/(lx*lx)+q_y*q_y/(ly*ly))*2.0*M_PI;
    if(alpha==beta)
        // symmetric gauge
        return 2.0*M_PI/(q+1e-30)*exp(-M_PI*(q_x*q_x*ly*1.0/lx+q_y*q_y*lx*1.0/ly)/nphi);
    else
        return 2.0*M_PI/(q+1e-30)*exp(-M_PI*(q_x*q_x*ly*1.0/lx+q_y*q_y*lx*1.0/ly)/nphi-q*d);
}

void lhamil::init_Coulomb_matrix(){
    off_head=nphi/2+nphi%2;
    long dim2 = pow(nphi, 2);
    long dim3 = pow(nphi,3);
    long dim4=dim3*(off_head+1);
    Coulomb_matrix.assign(2 * dim4, 0);
    for(int alpha = 0; alpha < 2; alpha++)
        for(int q_y = 0; q_y <= off_head; q_y++)
            for(int q_x = 0; q_x < nphi; q_x++)
                for(int n = 0; n < nphi; n++)
                    for(int m = 0; m < nphi; m++)
                        Coulomb_matrix[alpha * dim4 + q_y * dim3 + q_x * dim2 + n * nphi + m] = Coulomb_interaction(0, alpha, q_x, q_y) * cos(-2.0 * M_PI * q_x * q_y / nphi + 2.0 * M_PI * (m - n) * q_x / nphi)/2.0;
}

void lhamil::set_hamil(basis & _sector ,double _lx, double _ly, long _nphi,double _d){
    long nbasis_up, nbasis_down;
    d = _d;
    sector=_sector;
    nsite = sector.nsite;
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    off_head=nphi/2+nphi%2;
    init_Coulomb_matrix();
    nbasis_up = sector.nbasis_up;
    nbasis_down = sector.nbasis_down;
    nHilbert = nbasis_up * nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector<double> matrix_elements;
    std::map<long, double> ::iterator it;
    std::map<long, double> col_indices;
    inner_indices.reserve(nHilbert * nsite);
    matrix_elements.reserve(nHilbert * nsite);
    outer_starts.reserve(nHilbert + 1);
    long mask, mask_u, mask_d, b, p, n, m, i, j, k, l, t, nsignu, nsignd, nsign;
    long row = 0;
    outer_starts.push_back(0);
    for(i = 0; i < nbasis_up; i++) {
        for(j = 0; j < nbasis_down; j++) {
            // start of new row of nonzero elements
            // select two electrons in left-basis <m_1, m_2|
            long dim3 = pow(nphi, 3);
            long dim2 = pow(nphi, 2);
            long dim4 = dim3 * (off_head+1);
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
                        for(t = 0; t <= off_head; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    // q=0 is the uniform background charge, which is canceled out
                                    if( q_y!= 0 && q_x != 0) {
                                        // Coulomb matrix element in symmetric gauge
                                        V_uu += Coulomb_matrix[q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                if(abs(V_uu) > 1e-6) {
                                    nsign=nsignu;
                                    if(ncross&& mcross)
                                       nsign+=sector.nel_up+sector.nel_up-2;
                                    else if(ncross && !mcross || !ncross && mcross)
                                       nsign+=sector.nel_up-1;

                                    it = col_indices.find(k * nbasis_down + j);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(k * nbasis_down + j, V_uu * pow(-1, nsign)));
                                    else
                                        it->second += V_uu * pow(-1, nsign);
                                }
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
                        for(t = 1; t <= off_head; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    // q=0 is the uniform background charge, which is canceled out
                                    if( q_y!= 0 && q_x != 0) {
                                        // Coulomb matrix element in symmetric gauge
                                        V_uu += Coulomb_matrix[q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                if(abs(V_uu) > 1e-6) {
                                    nsign=nsignu;
                                    if(ncross&& mcross)
                                       nsign+=sector.nel_up+sector.nel_up-2;
                                    else if(ncross && !mcross || !ncross && mcross)
                                       nsign+=sector.nel_up-1;
                                    it = col_indices.find(k * nbasis_down + j);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(k * nbasis_down + j, V_uu * pow(-1, nsign)));
                                    else
                                        it->second += V_uu * pow(-1, nsign);
                                }
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
                        for(t = 0; t <= off_head; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    // q=0 is the uniform background charge, which is canceled out
                                    if(q_y != 0 && q_x != 0) {
                                        V_dd += Coulomb_matrix[q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                if(abs(V_dd) > 1e-6) {
                                    nsign=nsignd;
                                    if(ncross&& mcross)
                                       nsign+=sector.nel_down+sector.nel_down-2;
                                    else if(ncross && !mcross || !ncross && mcross)
                                       nsign+=sector.nel_down-1;
                                    it = col_indices.find(i * nbasis_down + l);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(i * nbasis_down + l, V_dd * pow(-1, nsign)));
                                    else
                                        it->second += V_dd * pow(-1, nsign);
                                }
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
                        for(t = 1; t <= off_head; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    // q=0 is the uniform background charge, which is canceled out
                                    if(q_y != 0 && q_x != 0) {
                                        V_dd += Coulomb_matrix[q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                if(abs(V_dd) > 1e-6) {
                                    nsign=nsignd;
                                    if(ncross&& mcross)
                                       nsign+=sector.nel_down+sector.nel_down-2;
                                    else if(ncross && !mcross || !ncross && mcross)
                                       nsign+=sector.nel_down-1;
                                    it = col_indices.find(i * nbasis_down + l);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(i * nbasis_down + l, V_dd * pow(-1, nsign)));
                                    else
                                        it->second += V_dd * pow(-1, nsign);
                                }
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
                        for(t = 0; t <=off_head ; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    if(q_y != 0 && q_x != 0) {
                                        // Coulomb matrix element, in symmetric gauge
                                        V_ud += Coulomb_matrix[dim4 + q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                nsign=nsignu+nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_up+sector.nel_down-2;
                                else if(ncross && !mcross)
                                   nsign+=sector.nel_up-1;
                                else if(!ncross && mcross)
                                   nsign+=sector.nel_down-1;

                                if(abs(V_ud) > 1e-6) {
                                    it = col_indices.find(k * nbasis_down + l);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(k * nbasis_down + l, V_ud * pow(-1, nsign)));
                                    else
                                        it->second += V_ud * pow(-1, nsign);
                                }
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
                        for(t = 1; t <=off_head ; t++) {
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
                                for(q_x = 0; q_x < nphi; q_x++)
                                    if(q_y != 0 && q_x != 0) {
                                        // Coulomb matrix element, in symmetric gauge
                                        V_ud += Coulomb_matrix[dim4 + q_y * dim3 + q_x * dim2 + n * nphi + mt];
                                    }
                                nsign=nsignu+nsignd;
                                if(ncross&& mcross)
                                   nsign+=sector.nel_up+sector.nel_down-2;
                                else if(ncross && !mcross)
                                   nsign+=sector.nel_up-1;
                                else if(!ncross && mcross)
                                   nsign+=sector.nel_down-1;
                                if(abs(V_ud) > 1e-6) {
                                    it = col_indices.find(k * nbasis_down + l);
                                    if(it == col_indices.end())
                                        col_indices.insert(std::pair<long, double>(k * nbasis_down + l, V_ud * pow(-1, nsign)));
                                    else
                                        it->second += V_ud * pow(-1, nsign);
                                }
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
            for(it = col_indices.begin(); it != col_indices.end(); it++) {
                inner_indices.push_back(it->first);
                matrix_elements.push_back(it->second);
            }
            row += col_indices.size();
            outer_starts.push_back(row);
            col_indices.clear();
        }
    }
    H.init(outer_starts, inner_indices, matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}



void lhamil::coeff_update() {
    double eigenvalues_0=1;
    double epsilon=1e-6;

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    Vec phi_0,phi_1,phi_2;
    phi_0.init_random(nHilbert,seed);
    norm[0]=1;
    phi_1 = H*phi_0;
    overlap[0] = phi_0 * phi_1;
    phi_1 -= phi_0 * overlap[0];
    norm[1] =phi_1.normalize();

    for(int i = 1; i < lambda; i++) {
        phi_2= H * phi_1;
        overlap[i] = phi_1 * phi_2;
        #pragma ivdep
        phi_2 -= phi_1 * overlap[i] + phi_0 * norm[i];
        norm[i+1] = phi_2.normalize();
        swap(&phi_0,&phi_1,&phi_2);

        // checking if the iteration is converged
        /*
        if(i>10 and i%5==0) {
            diag(i);
            if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon) {
                lambda=i+1;
                break;
            }
            else
                eigenvalues_0=eigenvalues[0];
                //cout<<i<<" "<<eigenvalues_0<<endl;
        }
        */

    }
}

void lhamil::coeff_explicit_update()
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-6;

    double norm_factor,overlap_factor;
    double *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new double[nHilbert];
    phi_1=new double[nHilbert];
    phi_2=new double[nHilbert];
    phi_t=new double[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    //phi_0.init_random(nHilbert,seed);
    #if __cplusplus > 199711L
    std::mt19937 rng(seed);
    for(i=0; i<nHilbert; i++)
        phi_0[i]=rng()*1.0/rng.max()-0.5;
    #else
    init_genrand64(seed);
    for(i=0; i<nHilbert; i++)
        phi_0[i]=genrand64_real3()-0.5;
    #endif

    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_0[i]*phi_0[i];
    norm_factor=sqrt(norm_factor);

    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor;
    norm[0]=1;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(double)*nHilbert);
    #pragma omp parallel for schedule(static)
    for(i=0; i<H.outer_starts.size()-1; i++) {
        #pragma ivdep
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    #pragma omp parallel for reduction(+:overlap_factor)
    for(i=0; i<nHilbert; i++)
        overlap_factor+=phi_1[i]*phi_0[i];
    overlap[0]=overlap_factor;

    //phi_1 -= phi_0 * overlap[0];
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_1[i]*phi_1[i];
    norm_factor=sqrt(norm_factor);

    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor;
    norm[1]=norm_factor;

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(double)*nHilbert);
        #pragma omp parallel for schedule(static)
        for(i=0; i<H.outer_starts.size()-1; i++) {
            #pragma ivdep
            for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }

        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        #pragma omp parallel for reduction(+:overlap_factor)
        for(i=0; i<nHilbert; i++)
            overlap_factor+=phi_1[i]*phi_2[i];
        overlap[j]=overlap_factor;

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normalize();
        norm_factor=0;
        #pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=phi_2[i]*phi_2[i];
        norm_factor=sqrt(norm_factor);

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor;
        norm[j+1]=norm_factor;

        phi_s=phi_0;
        phi_0=phi_1;
        phi_1=phi_2;
        phi_2=phi_s;

        // checking if the iteration is converged
        // the overhead of diagonalization is small
         /*
        if(j>10 and j%5==0){
          diag(j);
          if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon){
              lambda=j+1;
              break;
          }
          else
             eigenvalues_0=eigenvalues[0];
             //cout<<j<<" "<<eigenvalues_0<<endl;
        }
        */

    }
    delete phi_0,phi_1,phi_2,phi_t;
}

//Lanczos update version for spectral function calculation
void lhamil::coeff_update_wopt(vector<double> O_phi_0)
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-8;

    double norm_factor,overlap_factor;
    double *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new double[nHilbert];
    phi_1=new double[nHilbert];
    phi_2=new double[nHilbert];
    phi_t=new double[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    for(i=0; i<nHilbert; i++)
        phi_0[i]=O_phi_0[i];

    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_0[i]*phi_0[i];
    norm_factor=sqrt(norm_factor)+1e-24;
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor;

    norm[0]=1;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(double)*nHilbert);
    #pragma omp parallel for schedule(static)
    for(i=0; i<H.outer_starts.size()-1; i++) {
        #pragma ivdep
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    #pragma omp parallel for reduction(+:overlap_factor)
    for(i=0; i<nHilbert; i++)
        overlap_factor+=phi_1[i]*phi_0[i];
    overlap[0]=overlap_factor;

    //phi_1 -= phi_0 * overlap[0];
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_1[i]*phi_1[i];
    norm_factor=sqrt(norm_factor);

    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor;
    norm[1]=norm_factor;

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(double)*nHilbert);
        #pragma omp parallel for schedule(static)
        for(i=0; i<H.outer_starts.size()-1; i++) {
            #pragma ivdep
            for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }

        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        #pragma omp parallel for reduction(+:overlap_factor)
        for(i=0; i<nHilbert; i++)
            overlap_factor+=phi_1[i]*phi_2[i];
        overlap[j]=overlap_factor;

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normalize();
        norm_factor=0;
        #pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=phi_2[i]*phi_2[i];
        norm_factor=sqrt(norm_factor);

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor;
        norm[j+1]=norm_factor;

        phi_s=phi_0;
        phi_0=phi_1;
        phi_1=phi_2;
        phi_2=phi_s;

        // checking if the iteration is converged
        // the overhead of diagonalization is small
        /*
        if(j>10 and j%5==0){
          diag(j);
          if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon){
              lambda=j+1;
              break;
          }
          else
             eigenvalues_0=eigenvalues[0];
        }
        */
    }
    delete phi_0,phi_1,phi_2,phi_t;
}

void lhamil::diag()
{
    if(norm.size()==0) coeff_update();
    int l=lambda;
    double *e = new double[l];
    double *h = new double [l * l];
    memset(h, 0, sizeof(double)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_dsyev(h,e,l);

    psi_0.assign(l,0);
    psi_n0.assign(l,0);
    eigenvalues.assign(l,0);

    for(int i=0; i<l; i++) {
        psi_0[i]=h[i];
        psi_n0[i]=h[i*l];
        eigenvalues[i]=e[i];
    }
    delete h,e;
}

void lhamil::diag(int l)
{
    if(norm.size()==0) coeff_update();
    double *e = new double[l];
    double *h = new double [l * l];
    memset(h, 0, sizeof(double)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_dsyev(h,e,l);
    eigenvalues.assign(l,0);
    psi_0.assign(l,0);
    psi_n0.assign(l,0);
    for(int i=0; i<l; i++) {
        eigenvalues[i]=e[i];
        psi_0[i]=h[i];
        psi_n0[i]=h[i*l];
    }
    delete h,e;
}

// given the seed of phi_0, generate the Lanczos basis again, and transform the eigenstate to the desired eigenstate representation
void lhamil::eigenstates_reconstruction() {
    int l,n;
    double Norm=0;
    Vec phi_0,phi_1,phi_2;
    E0=eigenvalues[0];
    // repeat the iteration, but with normalization and overlap values at hand
    phi_0.init_random(nHilbert,seed);
    phi_1 = H * phi_0;
    phi_1 -= phi_0*overlap[0];
    phi_1 /= norm[1];
    psir_0.assign(nHilbert,0);
    for(n=0; n<nHilbert; n++)
        psir_0[n]+=psi_0[0]*phi_0.value[n]+psi_0[1]*phi_1.value[n];

    for(l=2; l<overlap.size(); l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++)
            psir_0[n]+=psi_0[l]*phi_2.value[n];
        swap(&phi_0,&phi_1,&phi_2);
    }
    for(n=0; n<nHilbert; n++)
        Norm+=psir_0[n]*psir_0[n];
    for(n=0; n<nHilbert; n++)
        psir_0[n]/=sqrt(Norm);
}

double lhamil::ground_state_energy() {
    vector<double> H_psir0;
    double overlap=0;
    if(psir_0.size()!=0) {
        H_psir0=H*psir_0;
        for(int i=0; i<nHilbert; i++)
            overlap+=psir_0[i]*H_psir0[i];
    }
    return overlap;
}


// spectral function continued fraction version
double lhamil::spectral_function(double omega, double eta) {
    // calculation continued fraction using modified Lentz method
    complex<double> E(omega,eta);
    vector< complex<double> >  f,c,d,delta;
    complex<double> a,b,G;
    double I;
    b=E-overlap[0];
    if(abs(b)<1e-25)
        f.push_back(1e-25);
    else
        f.push_back(b);
    c.push_back(f[0]);
    d.push_back(0);
    delta.push_back(0);
    for(int n=1; n<lambda; n++) {
        b=E-overlap[n];
        a=-norm[n]*norm[n];
        d.push_back(b+a*d[n-1]);
        if(abs(d[n])<1e-25)
            d[n]=1e-25;
        c.push_back(b+a/c[n-1]);
        if(abs(c[n])<1e-25)
            c[n]=1e-25;
        d[n]=1.0/d[n];
        delta.push_back(c[n]*d[n]);
        f.push_back(f[n-1]*delta[n]);
        if(abs(delta.back()-1)<1e-15)
            break;
    }
    G=1.0/f.back();

    I=-G.imag()/M_PI;
    c.clear();
    d.clear();
    f.clear();
    delta.clear();
    return I;
}


void lhamil::save_to_file(const char* filename) {
    if(eigenvalues.size()==0) return;
    ofstream odf;
    odf.open(filename,ofstream::out);
    odf<<"nHilbert:="<<setw(10)<<" "<<nHilbert<<endl;
    odf<<"lambda:="<<setw(10)<<" "<<lambda<<endl;
    odf<<"seed:="<<setw(10)<<" "<<seed<<endl;
    odf<<"eigenvalues:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<eigenvalues[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"norm:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<=lambda; i++)
        odf<<norm[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"overlap:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<overlap[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<psi_0[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_n0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
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

void lhamil::print_hamil(int range) {
    int i,j,count;
    for(i=0; i<range; i++) {
        if(i==0)
            cout<<"[[";
        else cout<<" [";
        count=0;
        for(j=0; j<range; j++) {
            if(j==H.inner_indices[H.outer_starts[i]+count])
                cout<<H.value[H.outer_starts[i]+count++]<<",";
            else
                cout<<0<<",";
        }
        if(i==range-1)
            cout<<"]]"<<endl;
        else cout<<"]"<<endl;
    }
}

void lhamil::print_lhamil(int range) {
    for(int i=0; i<range; i++) {
        if(i==0)
            cout<<"[[";
        else cout<<" [";
        for(int j=0; j<nHilbert; j++) {
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

void lhamil::print_eigen( int range) {
    if(range>=lambda) range=lambda;
    cout<<"Eigenvalues:= [";
    for(int i=0; i<range; i++)
        cout<<eigenvalues[i]<<", ";
    cout<<",...]"<<endl;
}
