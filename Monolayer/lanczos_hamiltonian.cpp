#include"lanczos_hamiltonian.h"
inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,long _nHilbert,long _lambda, unsigned _seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(long _lambda,unsigned _seed) {
    lambda=_lambda;
    seed=_seed;
}

lhamil::~lhamil() {
}

void lhamil::clear() {
    H.clear();
    if(psir_0.size()!=0)
        psir_0.clear();
    if(psi_n0.size()!=0)
        psi_n0.clear();
    if(psi_0.size()!=0)
        psi_0.clear();
    if(eigenvalues.size()!=0)
        eigenvalues.clear();
    if(Coulomb_matrix.size()!=0)
        Coulomb_matrix.clear();
    if(norm.size()!=0)
        norm.clear();
    if(overlap.size()!=0)
        overlap.clear();
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

double lhamil::Coulomb_interaction(int k_x, int k_y) {
    double q=sqrt(k_x*k_x/(lx*lx)+k_y*k_y/(ly*ly))*2.0*M_PI;
    return 2.0*M_PI/(q+1e-30)*exp(-q*q/2.0)*pow(1.0-exp(-q*q/2.0),nLL*2);
}

void lhamil::init_Coulomb_matrix() {
    Coulomb_matrix.assign(nphi*nphi, 0);
    for(int s = 0; s < nphi; s++)
        for(int k_y = 0; k_y <nphi; k_y++) {
            double V=0;
            for(int k_x=-nphi/2; k_x<=nphi/2; k_x++)
                if(!(k_y==0 && k_x==0))
                    V+=Coulomb_interaction(k_x,k_y)*cos(2.0*M_PI*s*k_x/nphi)/(2.0*lx*ly);
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


void lhamil::set_hamil(basis & sector ,double _lx, double _ly,long _nphi, long _nLL) {
    lx = _lx;
    ly = _ly;
    nphi = _nphi;
    nLL = _nLL;
    init_Coulomb_matrix();
    nHilbert = sector.nbasis;
    vector<complex<double> > matrix_elements;
    H.clear();
    H.inner_indices.reserve(nHilbert * nphi);
    H.value.reserve(nHilbert * nphi);
    H.outer_starts.reserve(nHilbert + 1);
    long mask, b, n, m, i, k, s,t, sign;
    long k1,k2,lbasis,rbasis,signl,signr;
    long kx=sector.K;
    long Cl,Cr;
    long row = 0;
    H.outer_starts.push_back(0);
    for(i = 0; i < nHilbert; i++) {
        // start of new row of nonzero elements
        matrix_elements.assign(nHilbert,0);
        
        // determine the subbasis size of basis i
        for(k1=1; k1<=sector.C; k1++)
            if(sector.translate(sector.id[i],k1,signl)==sector.id[i]) {
                Cl=k1;
                break;
            }
        // select two electrons in left-basis <m_1, m_2|
        // n=j1, m=j2
        for(k1=0; k1<Cl; k1++) {
            lbasis=sector.translate(sector.id[i],k1,signl);
            for(n = 0; n < nphi-1; n++)
                for(m = n+1; m < nphi; m++) {
                    mask = (1 << n) + (1 << m);
                    // looking up the corresponding basis in id
                    // if there're two electrons on n and m;
                    if((lbasis&mask) == mask && n!=m) {
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
                                  if(sector.translate(mask_t+b,k2,signr)==(mask_t+b)) {
                                    Cr=k2;
                                    break;
                                  }
                                for(k2=0; k2<Cr; k2++) {
                                    rbasis=sector.inv_translate(mask_t+b,k2,signr);
                                    if(sector.basis_set.find(rbasis) != sector.basis_set.end()){
                                        k = sector.basis_set[rbasis];
                                        sign=sector.get_sign(lbasis,n,m,nt,mt);
                                        complex<double> FT_factor=complex<double>(cos(2.0*M_PI*kx*(k1-k2)/sector.C),sin(2.0*M_PI*kx*(k1-k2)/sector.C))/sqrt(Cl*Cr);
                                        matrix_elements[k]+=2.0*Coulomb_matrix[s*nphi+abs(t)]*sign*FT_factor*signl*signr;
                                        
                                    }
                                }
                             }
                        }
                    }
                }
        }
        matrix_elements[i]+=E_cl*sector.nel;
        long count=0;
        for(k=0; k<nHilbert; k++)
            if(abs(matrix_elements[k])>1e-10) {
                H.inner_indices.push_back(k);
                H.value.push_back(matrix_elements[k]);
                count++;
            }
        row += count;
        H.outer_starts.push_back(row);
    }
    matrix_elements.clear();
}

void lhamil::Gram_Schmidt_orthogonalization(Vec &phi, int n){
    Vec phi_0,phi_1,phi_2;
    phi_0.init_random(nHilbert,seed);
    phi_1 = H* phi_0;
    phi_1 -= phi_0 *overlap[0];
    phi_1 /= norm[1]; 

    double q= (phi_0*phi).real();
    phi = (phi-phi_0*q)/(1-q*q);
    q=(phi_1*phi).real();
    phi = (phi-phi_1*q)/(1-q*q);

    for(int i = 1; i < n; i++) {
        phi_2= H * phi_1;
        #pragma ivdep
        phi_2 -= phi_1 * overlap[i] + phi_0 * norm[i];
        phi_2 /= norm[i+1];

        q=(phi_2*phi).real();
        phi = (phi-phi_2*q)/(1-q*q);
        
        swap(&phi_0,&phi_1,&phi_2);
        }
    phi_0.clear();
    phi_1.clear();
    phi_2.clear();
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
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=std::norm(phi_0[i]);
    norm_factor=sqrt(norm_factor);
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(complex<double>)*nHilbert);
    #pragma omp parallel for schedule(static)
    for(i=0; i<H.outer_starts.size()-1; i++) {
        #pragma ivdep
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    //#pragma omp parallel for reduction(+:overlap_factor)
    for(i=0; i<nHilbert; i++)
        overlap_factor+=conj(phi_1[i])*phi_0[i];

    //overlap.push_back(overlap_factor.real());
    overlap[0]=overlap_factor.real();

    //phi_1 -= phi_0 * overlap[0];
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=std::norm(phi_1[i]);
    norm_factor=sqrt(norm_factor);

    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor;


    //norm.push_back(1);
    //norm.push_back(norm_factor);
     norm[0]=1;
    norm[1]=norm_factor;

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(complex<double>)*nHilbert);
        #pragma omp parallel for schedule(guided,4)
        for(i=0; i<H.outer_starts.size()-1; i++) {
            #pragma ivdep
            for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }
        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        //#pragma omp parallel for reduction(+:overlap_factor)
        for(i=0; i<nHilbert; i++)
            overlap_factor+=conj(phi_1[i])*phi_2[i];

        //overlap.push_back(overlap_factor.real());
        overlap[j]=overlap_factor.real();

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

        //norm[j+1] = phi_2.normed();
        norm_factor=0;
        #pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=std::norm(phi_2[i]);
        norm_factor=sqrt(norm_factor);

        #pragma omp parallel for schedule(static)
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

//Lanczos update version for spectral function calculation
/*
void lhamil::coeff_update_wopt(vector<complex<double> > O_phi_0)
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
    }
    delete phi_0,phi_1,phi_2,phi_t;
}
*/

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
    for(n=0; n<nHilbert; n++)
        psir_0[n] += psi_0[0]*conj(phi_0.value[n])+psi_0[1]*conj(phi_1.value[n]);

    for(l=2; l<overlap.size(); l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++)
            psir_0[n]+= psi_0[l]*conj(phi_2.value[n]);
        swap(&phi_0,&phi_1,&phi_2);
    }
    double Norm=0;
    for(n=0; n<nHilbert; n++)
        Norm+=std::norm(psir_0[n]);
    Norm=sqrt(Norm);
    for(n=0; n<nHilbert; n++)
        psir_0[n]/=Norm;
    phi_0.clear();
    phi_1.clear();
    phi_2.clear();
}

double lhamil::ground_state_energy() {
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
    for(int n=1; n<overlap.size(); n++) {
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

void lhamil::print_hamil_CSR() {
    std::cout << "hamiltonian in CSR format: " << std::endl;
    std::cout << "------------------------------" << std::endl;
    H.print();
}

void lhamil::print_hamil(int range) {
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
            if(k==col_index && row_starts<H.outer_starts[i+1]) {
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

void lhamil::print_lhamil(int range) {
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

void lhamil::print_eigen( int range) {
    if(range>=overlap.size()) range=overlap.size();
    cout<<"Eigenvalues:= [";
    for(int i=0; i<range; i++)
        cout<<eigenvalues[i]<<", ";
    cout<<",...]"<<endl;
}
