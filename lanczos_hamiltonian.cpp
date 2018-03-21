#include"lanczos_hamiltonian.h"
inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,long _nHilbert,long _lambda, unsigned _seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(basis &_sector,double _d,long _lambda,unsigned _seed) {
    sector=_sector;
    lambda=_lambda;
    seed=_seed;
    set_hamil(sector,_d);
}

lhamil::~lhamil() {
}

void lhamil::init(basis &_sector,double _d,long _lambda,unsigned _seed) {
    sector=_sector;
    lambda=_lambda;
    seed=_seed;
    set_hamil(sector,_d);
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

double lhamil::Coulomb_interaction(int alpha, int beta, int q_x, int q_y){
  double q=sqrt(q_x*q_x+q_y*q_y)*2*M_PI/nsite;
  if(alpha==beta)
      return 2.0*M_PI/(q+1e-8)*exp(-q*q/2);
  else
      return 2.0*M_PI/(q+1e-8)*exp(-q*q/2-q*d);
}

void lhamil::set_hamil(basis &_sector, double _d)
{
    long nbasis_up,nbasis_down;
    d=_d;
    sector=_sector;
    nsite=sector.nsite;
    nphi=nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector< complex<double> > matrix_elements;
    std::map<long,complex<double> > ::iterator it;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    long mask,mask_u,mask_d,b,p,n,m,i,j,k,l,t,nsignu,nsignd;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
      for(j=0; j<nbasis_down; j++) {
        std::map<long,complex<double> > col_indices;
        // select two electrons in left-basis <m_1, m_2|
        for(n=0;n<nsite;n++)
          for(m=n;m<nsite;m++){
            mask=(1<<n)+(1<<m);
            // consider the upper-layer two electrons
            // if there're two electrons on n and m;
            if(i&mask==mask && m!=n){
              // b is the rest electon positions
              b=i^mask;
              long nt,mt,mask_ut,occ_ut;
              nsignu=0;
              // perform translation
              for(t=0;t<nsite;t++){
                // PBC, if one electron cross left boundary, sign change with -1
                if(n-t<0){
                  nt=n-t+nsite;
                  nsignu++;
                }
                // PBC, if one electron cross right boundary, sign change with -1
                if(m+t>=nsite){
                  mt=m+t-nsite;
                  nsignu++;
                }
                // the translated two electrons indices
                mask_ut=(1<<nt)+(1<<mt);
                // occupation of electons on the translated position
                occ_ut=mask_ut&b;
                // if there're no electon on the translated position
                // which is a valid translation, can be applied
                if(occ_ut==0){
                  // the translated indices
                  k=mask_ut+b;
                  // calculate the Coulomb matrix contribution and
                  // add it to the hamiltonian matrix
                  if(k!=i) {
                    // Coulomb matrix element, in Landau gauge
                    long q_x,q_y;
                    q_y=nt;
                    complex<double> V_uu=0;
                    for(q_x=0;q_x<nsite;q_x++)
                      if(q_y!=0&&q_x!=0)
                        V_uu+=Coulomb_interaction(0,0,q_x,q_y)*complex<double>(cos((n-m+t)*q_x*2.0*M_PI/nphi),sin((n-m+t)*q_x*2.0*M_PI/nphi))/(4.0*M_PI*nphi);
                    it=col_indices.find(k*nbasis_down+j);
                    if(it==col_indices.end())
                      col_indices.insert(std::pair<long,complex<double> >(k*nbasis_down+j,V_uu*pow(-1,nsignu)));
                    else
                      it->second+=V_uu*pow(-1,nsignu);
                    }
                  }
                // two electrons are occupied, and to be crossed next
                else if(occ_ut==mask_ut)
                  nsignu+=2;
                 // one electron is occupied, and to be crossed next
                else if(occ_ut!=0 && occ_ut!=mask_ut)
                  nsignu++;
                }
               }

               // consider the lower-layer two electrons
               // if there're two electrons on n and m;
               if(j&mask==mask && m!=n){
                 // b is the rest electon positions
                 b=j^mask;
                 long nt,mt,mask_dt,occ_dt;
                 nsignd=0;
                 // perform translation
                 for(t=0;t<nsite;t++){
                   // PBC, if one electron cross left boundary, sign change with -1
                   if(n-t<0){
                     nt=n-t+nsite;
                     nsignd++;
                    }
                   // PBC, if one electron cross right boundary, sign change with -1
                   if(m+t>=nsite){
                     mt=m+t-nsite;
                     nsignd++;
                     }
                   // the translated two electrons indices
                   mask_dt=(1<<nt)+(1<<mt);
                   // occupation of electons on the translated position
                   occ_dt=mask_dt&b;
                   // if there're no electon on the translated position
                   // which is a valid translation, can be applied
                   if(occ_dt==0){
                     // the translated indices
                     k=mask_dt+b;
                     // calculate the Coulomb matrix contribution and
                     // add it to the hamiltonian matrix
                     if(k!=j) {
                       // Coulomb matrix element, in Landau gauge
                       long q_x,q_y;
                       q_y=nt;
                       complex<double> V_dd=0;
                       for(q_x=0;q_x<nsite;q_x++)
                         if(q_y!=0&&q_x!=0)
                           V_dd+=Coulomb_interaction(1,1,q_x,q_y)*complex<double>(cos((n-m+t)*q_x*2.0*M_PI/nphi),sin((n-m+t)*q_x*2.0*M_PI/nphi))/(4.0*M_PI*nphi);
                       it=col_indices.find(i*nbasis_down+k);
                       if(it==col_indices.end())
                         col_indices.insert(std::pair<long,complex<double> >(i*nbasis_down+k,V_dd*pow(-1,nsignd)));
                       else
                         it->second+=V_dd*pow(-1,nsignd);
                       }
                     }
                     // two electrons are occupied, and to be crossed next
                   else if(occ_dt==mask_dt)
                     nsignd+=2;
                    // one electron is occupied, and to be crossed next
                   else if(occ_dt!=0 && occ_dt!=mask_dt)
                     nsignd++;
                  }
                }
                // consider the one electron in the upper layer
                // and one electron in the lower layer case
                mask_u=(1<<n);
                mask_d=(1<<m);
                if(i&mask_u==mask_u && j&mask_d==mask_d){
                  // b is the rest electon positions for upper-layer electrons
                  b=i^mask_u;
                  p=j^mask_d;
                  long nt,mt,mask_ut,occ_ut,mask_dt,occ_dt;
                  nsignu=0;
                  nsignd=0;
                  // perform translation
                  for(t=0;t<nsite;t++){
                    // PBC, if one electron cross left boundary, sign change with -1
                    if(n-t<0){
                      nt=n-t+nsite;
                      nsignu++;
                      }
                    // PBC, if one electron cross right boundary, sign change with -1
                    if(m+t>=nsite){
                      mt=m+t-nsite;
                      nsignd++;
                    }
                    // the translated upper electron index
                    mask_ut=(1<<nt);
                    mask_dt=(1<<mt);
                    // occupation of electons on the translated position
                    occ_ut=mask_ut&b;
                    occ_dt=mask_dt&p;
                    // if there're no electon on the translated position
                    // which is a valid translation, can be applied
                    if(occ_ut==0 && occ_dt==0){
                      // the translated indices
                      k=mask_ut+b;
                      l=mask_dt+p;
                      // calculate the Coulomb matrix contribution and
                      // add it to the hamiltonian matrix
                      if(k!=i && l!=j) {
                        // Coulomb matrix element, in Landau gauge
                        long q_x,q_y;
                        q_y=nt;
                        complex<double> V_ud=0;
                        for(q_x=0;q_x<nsite;q_x++)
                          if(q_y!=0&&q_x!=0)
                            V_ud+=Coulomb_interaction(1,0,q_x,q_y)*complex<double>(cos((n-m+t)*q_x*2.0*M_PI/nphi),sin((n-m+t)*q_x*2.0*M_PI/nphi))/(4.0*M_PI*nphi);
                        it=col_indices.find(k*nbasis_down+l);
                        if(it==col_indices.end())
                          col_indices.insert(std::pair<long,complex<double> >(k*nbasis_down+l,V_ud*pow(-1,nsignu+nsignd)));
                        else
                          it->second+=V_ud*pow(-1,nsignu+nsignd);
                        }
                      // two electrons are occupied, and to be crossed next
                      else if(occ_ut==mask_ut && occ_dt==mask_dt)
                        nsignu+=2;
                      // one electron is occupied, and to be crossed next
                      else if(occ_ut==0 && occ_dt==mask_dt || occ_dt==0 && occ_ut==mask_ut)
                        nsignu++;
                      }
                    }
                  }
               }
              for(it=col_indices.begin(); it!=col_indices.end(); it++) {
                inner_indices.push_back(it->first);
                matrix_elements.push_back(it->second);
              }
             row+=col_indices.size();
             outer_starts.push_back(row);
             col_indices.clear();
           }
    }
    H.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}

void lhamil::coeff_update() {
    double eigenvalues_0=1;
    double epsilon=1e-8;

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
        if(i>10 and i%5==0) {
            diag(i);
            if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon) {
                lambda=i+1;
                break;
            }
            else
                eigenvalues_0=eigenvalues[0];
        }
    }
}

void lhamil::coeff_explicit_update()
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-8;

    complex<double> norm_factor,overlap_factor;
    complex<double> *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new complex<double>[nHilbert];
    phi_1=new complex<double>[nHilbert];
    phi_2=new complex<double>[nHilbert];
    phi_t=new complex<double>[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

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
    //#pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=conj(phi_0[i])*phi_0[i];
    norm_factor=abs(norm_factor);

    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor.real();
    norm[0]=1;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(complex<double> )*nHilbert);
    //#pragma omp parallel for schedule(static)
    for(i=0; i<H.outer_starts.size()-1; i++) {
        //#pragma ivdep
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    //#pragma omp parallel for reduction(+:overlap_factor)
    for(i=0; i<nHilbert; i++)
        overlap_factor+=conj(phi_1[i])*phi_0[i];
    overlap[0]=overlap_factor;

    //phi_1 -= phi_0 * overlap[0];
    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    //#pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=conj(phi_1[i])*phi_1[i];
    norm_factor=abs(norm_factor);

    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor.real();
    norm[1]=norm_factor.real();

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(complex<double>)*nHilbert);
        //#pragma omp parallel for schedule(static)
        for(i=0; i<H.outer_starts.size()-1; i++) {
            //#pragma ivdep
            for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }

        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        //#pragma omp parallel for reduction(+:overlap_factor)
        for(i=0; i<nHilbert; i++)
            overlap_factor+=conj(phi_1[i])*phi_2[i];
        overlap[j]=overlap_factor;

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normalize();
        norm_factor=0;
        //#pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=conj(phi_2[i])*phi_2[i];
        norm_factor=abs(norm_factor);

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor.real();
        norm[j+1]=norm_factor.real();

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

//Lanczos update version for spectral function calculation
void lhamil::coeff_update_wopt(vector<complex<double> > O_phi_0)
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-8;

    complex<double> norm_factor,overlap_factor;
    complex<double> *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new complex<double>[nHilbert];
    phi_1=new complex<double>[nHilbert];
    phi_2=new complex<double>[nHilbert];
    phi_t=new complex<double>[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    for(i=0; i<nHilbert; i++)
        phi_0[i]=O_phi_0[i];

    norm_factor=0;
    //#pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=conj(phi_0[i])*phi_0[i];
    norm_factor=abs(norm_factor)+1e-24;
    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor.real();

    norm[0]=1;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(complex<double> )*nHilbert);
    //#pragma omp parallel for schedule(static)
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
    overlap[0]=overlap_factor;

    //phi_1 -= phi_0 * overlap[0];
    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    //#pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=conj(phi_1[i])*phi_1[i];
    norm_factor=abs(norm_factor);

    //#pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor.real();
    norm[1]=norm_factor.real();

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(double)*nHilbert);
        //#pragma omp parallel for schedule(static)
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
        overlap[j]=overlap_factor;

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normalize();
        norm_factor=0;
        //#pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=conj(phi_2[i])*phi_2[i];
        norm_factor=abs(norm_factor);

        //#pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor.real();
        norm[j+1]=norm_factor.real();

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
    complex<double> *h = new complex<double> [l * l];
    memset(h, 0, sizeof(complex<double>)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_zheev(h,e,l);

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
    complex<double> *h = new complex<double> [l * l];
    memset(h, 0, sizeof(complex<double>)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_zheev(h,e,l);
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
    complex<double> Norm=0;
    Vec phi_0,phi_1,phi_2;
    E0=eigenvalues[0];
    // repeat the iteration, but with normalization and overlap values at hand
    phi_0.init_random(nHilbert,seed);
    phi_1 = H * phi_0;
    phi_1 -= phi_0*overlap[0];
    phi_1 /= norm[1];
    psir_0.assign(nHilbert,0);
    for(n=0; n<nHilbert; n++)
        psir_0[n]+=conj(psi_0[0])*phi_0.value[n]+conj(psi_0[1])*phi_1.value[n];

    for(l=2; l<overlap.size(); l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++)
            psir_0[n]+=conj(psi_0[l])*phi_2.value[n];
        swap(&phi_0,&phi_1,&phi_2);
    }
    for(n=0; n<nHilbert; n++)
        Norm+=conj(psir_0[n])*psir_0[n];
    for(n=0; n<nHilbert; n++)
        psir_0[n]/=abs(Norm);
}

double lhamil::ground_state_energy() {
    vector<complex<double> > H_psir0;
    complex<double> overlap=0;
    if(psir_0.size()!=0) {
        H_psir0=H*psir_0;
        for(int i=0; i<nHilbert; i++)
            overlap+=psir_0[i]*H_psir0[i];
    }
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

void lhamil::print_hamil() {
    int i,j,count;
    for(i=0; i<nHilbert; i++) {
        if(i==0)
            cout<<"[[";
        else cout<<" [";
        count=0;
        for(j=0; j<nHilbert; j++) {
            if(j==H.inner_indices[H.outer_starts[i]+count])
                cout<<H.value[H.outer_starts[i]+count++]<<",";
            else
                cout<<0<<",";
        }
        if(i==nHilbert-1)
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
