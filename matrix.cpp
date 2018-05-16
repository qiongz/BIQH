#include"matrix.h"
Vec::Vec() {}

Vec::Vec(long _size) {
    size = _size;
    value.reserve(size);
}

Vec::Vec(long _size,const complex<double> _init) {
    size = _size;
    value.assign(size, _init);
}

Vec::Vec(const Vec &rhs) {
    size = rhs.size;
    for(int i = 0; i < size; i++)
        value.assign((rhs.value).begin(), (rhs.value).end());
}

Vec::~Vec() {
    if(size!=0) {
        value.clear();
        size=0;
    }
}

void Vec::assign(long  _size,const complex<double> _init) {
    size=_size;
    value.assign(size,_init);
}

void Vec::init_random(unsigned seed) {
    double norm=0;
    #if __cplusplus > 199711L
    std::mt19937 rng(seed);
    for(int i = 0; i < size; i++) {
        value[i] = complex<double>(rng() * 1.0 / rng.max()-0.5,rng()*1.0/rng.max()-0.5 );
        norm += std::norm(value[i]);
    }
    #else
    init_genrand64(seed);
    for(int i = 0; i < size; i++) {
        value[i]=complex<double>(genrand64_real3()-0.5,genrand64_real3()-0.5);
        //value[i]=complex<double>(genrand64_real3()-0.5,genrand64_real3()-0.5);
        norm += std::norm(value[i]);
    }
    #endif
    norm = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= norm;
}

void Vec::init_random(long _size,unsigned seed) {
    size=_size;
    value.resize(size);
    double norm=0;
    #if __cplusplus > 199711L
    std::mt19937 rng(seed);
    for(int i = 0; i < size; i++) {
        value[i] = complex<double>(rng() * 1.0 / rng.max()-0.5,rng()*1.0/rng.max()-0.5 );
        norm += std::norm(value[i]);
    }
    #else
    init_genrand64(seed);
    for(int i = 0; i < size; i++) {
        value[i]=complex<double>(genrand64_real3()-0.5,genrand64_real3()-0.5);
        norm += std::norm(value[i]);
    }
    #endif
    norm = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= norm;
}

void Vec::clear() {
    if(size!=0) {
        value.clear();
        size=0;
    }
}

double Vec::normed() {
    int i;
    double norm = 0;
    #pragma omp parallel for reduction(+:norm)
    for(i = 0; i < size; i++)
        norm += std::norm(value[i]);
    double normsq = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(i = 0; i < size; i++)
        value[i] /= normsq;
    return normsq;
}

void Vec::normalize() {
    int i;
    double norm = 0;
    #pragma omp parallel for reduction(+:norm)
    for(i = 0; i < size; i++)
        norm += std::norm(value[i]);
    double normsq = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(i = 0; i < size; i++)
        value[i] /= normsq;
}

Vec & Vec::operator=(const Vec & rhs) {
    if(this==&rhs) return *this;
    size = rhs.size;
    value.assign((rhs.value).begin(), (rhs.value).end());
    return *this;
}

Vec & Vec::operator-=(const Vec & rhs) {
    if(this==&rhs)
        return *this;
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] -= rhs.value[i];
    return *this;
}

Vec & Vec::operator+=(const Vec & rhs) {
    if(this==&rhs)
        return *this;
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] += rhs.value[i];
    return *this;
}

Vec & Vec::operator*=(const double & rhs) {
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] *=rhs;
    return *this;
}

Vec & Vec::operator*=(const complex<double> & rhs) {
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] *=rhs;
    return *this;
}

Vec & Vec::operator/=(const double & rhs) {
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= rhs;
    return *this;
}

Vec Vec::operator+(const Vec & rhs) {
    Vec rt(rhs.size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] + (rhs.value)[i];
    return rt;
}

Vec Vec::operator-(const Vec & rhs) {
    Vec rt(rhs.size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] - (rhs.value)[i];
    return rt;
}

Vec Vec::operator/(const double &rhs) {
    Vec rt(size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] / rhs;
    return rt;
}

Vec Vec::operator*(const double &rhs) {
    Vec rt(size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i]* rhs;
    return rt;
}

Vec Vec::operator*(const complex<double> &rhs) {
    Vec rt(size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i]* rhs;
    return rt;
}

complex<double> Vec::operator*(const Vec &rhs) {
    if(size != rhs.size) return 0 ;
    complex<double> overlap = 0;
    //#pragma omp parallel for reduction(+:overlap)
    for(int i = 0; i < size; i++)
        overlap += conj(value[i]) * (rhs.value)[i];
    return overlap;
}

ostream & operator<<(ostream & os, const Vec & _b) {
    os << "[ ";
    for(int i = 0; i < _b.size - 1; i++)
        os << setprecision(6) << setw(10) << (_b.value)[i] << ", ";
    os << setprecision(6) << setw(10) << (_b.value)[_b.size - 1] << " ]";
    return os;
}

Mat::Mat() {}

Mat::Mat(const Mat &rhs) {
    inner_indices.assign(rhs.inner_indices.begin(),rhs.inner_indices.end());
    value.assign(rhs.value.begin(),rhs.value.end());
    outer_starts.assign(rhs.outer_starts.begin(),rhs.outer_starts.end());
}

Mat::~Mat() {
    if(value.size()!=0) {
        value.clear();
        inner_indices.clear();
        outer_starts.clear();
    }
}

void Mat::clear() {
    if(value.size()!=0) {
        value.clear();
        inner_indices.clear();
        outer_starts.clear();
    }
}

Mat & Mat::operator=(const Mat & rhs) {
    if(this==&rhs) return *this;
    value.assign((rhs.value).begin(), (rhs.value).end());
    inner_indices.assign(rhs.inner_indices.begin(),rhs.inner_indices.end());
    outer_starts.assign(rhs.outer_starts.begin(),rhs.outer_starts.end());
    return *this;
}

Vec Mat::operator*(const Vec &rhs)const {
    Vec phi(rhs.size,0);
    if(rhs.size!=outer_starts.size()) return phi;
    #pragma omp parallel for schedule(guided,4)
    for(int i=0; i<outer_starts.size(); i++) {
        #pragma ivdep
       // for(int idx=outer_starts[i]; idx<outer_starts[i+1]; idx++)
       for(int idx=outer_starts[i];idx<outer_starts[i]+outer_size[i];idx++)
            phi.value[i]+=value[idx]*rhs.value[inner_indices[idx]];
    }
    return phi;
}

vector<complex<double> > Mat::operator*(const vector< complex<double> > &rhs)const {
    vector< complex<double> > phi;
    phi.assign(rhs.size(),0);
    if(rhs.size()!=outer_starts.size()) return phi;
    #pragma omp parallel for schedule(guided,4)
    for(int i=0; i<outer_starts.size(); i++) {
        #pragma ivdep
       // for(int idx=outer_starts[i]; idx<outer_starts[i+1]; idx++)
       for(int idx=outer_starts[i];idx<outer_starts[i]+outer_size[i];idx++)
            phi[i]+=value[idx]*rhs[inner_indices[idx]];
    }
    return phi;
}

void Mat::init(const vector<long> & _outer,const vector<long> & _inner, const vector< complex<double> > &_value) {
    outer_starts.assign(_outer.begin(),_outer.end());
    inner_indices.assign(_inner.begin(),_inner.end());
    value.assign(_value.begin(),_value.end());
}

void Mat::print() {
    std::cout<<"value:=         [";
    for(int i=0; i<value.size(); i++) std::cout<<setw(4)<<setprecision(2)<<value[i]<<" ";
    std::cout<<" ]"<<std::endl;
    std::cout<<"inner_indices:= [";
    for(int i=0; i<inner_indices.size(); i++) std::cout<<setw(4)<<setprecision(2)<<inner_indices[i]<<" ";
    std::cout<<" ]"<<std::endl;
    std::cout<<"outer_starts:=  [";
    for(int i=0; i<outer_starts.size(); i++) std::cout<<setw(4)<<setprecision(2)<<outer_starts[i]<<" ";
    std::cout<<" ]"<<std::endl;
}

/*
void diag_dsyev(double *h, double *e, int l) {
    char jobz, uplo;
    int info;
    jobz = 'V';
    uplo = 'U';
    int lda = l;
    int lwork = 3 * l - 1;
    double *work = new double[lwork];
    dsyev_(&jobz, &uplo, &l, h, &lda, e, work, &lwork, &info);
    delete [] work;
}
*/

void diag_dsyev(double *hamiltonian, double *energy, int l){
    char jobz,uplo;
    int info,lda;
    jobz = 'V';
    uplo = 'U';
    lda=l;
    info=LAPACKE_dsyev(LAPACK_COL_MAJOR,jobz, uplo, l,hamiltonian,lda, energy);
}

void diag_zheev(complex<double> *hamiltonian, double *energy, int l){
    char jobz,uplo;
    int info,lda;
    jobz = 'V';
    uplo = 'U';
    lda=l;
    info=LAPACKE_zheevd(LAPACK_COL_MAJOR,jobz, uplo, l,hamiltonian,lda, energy);
}


double func_ExpInt(double t, void *params) {
    double z = *(double *)params;
    double ret_value = 1.0/sqrt(t)*exp(-z*t);
    return ret_value;
}

struct my_params{
    double r;
    double d;
};
   
double func_BesselInt(double k, void *p){
    struct my_params *params =(struct my_params *)p;
     
    double f = gsl_sf_bessel_J0(k*params->r)*exp(-k*params->d)/(k+1e-10);
    return f;
}

double Integrate_ExpInt(double x) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(50000);
    double result,error;

    gsl_function F;
    F.function = &func_ExpInt;
    F.params= & x;

    gsl_integration_qagiu(& F,1.0,1.0e-6, 1.0e-6, 50000,w, &result, & error);
    return result;
}

double Integrate_BesselInt(double r,double d) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(50000);
    double result,error;

    my_params params;
    params.r=r;
    params.d=d;

    gsl_function F;
    F.function = &func_BesselInt;
    F.params= & params;

    gsl_integration_qagiu(& F,0.0,1.0e-6, 1.0e-6, 50000,w, &result, & error);
    return result;
}
