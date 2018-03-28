#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(long _nphi,long _nel):nphi(_nphi),nel(_nel){
  K=-1;
}

basis::basis(long _nphi,long _nel, long _K):nphi(_nphi),nel(_nel),K(_K) {
}

const basis & basis::operator =(const basis & _basis) {
    if(this !=&_basis) {
        nphi=_basis.nphi;
        nel=_basis.nel;
        nbasis=_basis.nbasis;
        basis_set=_basis.basis_set;
        id=_basis.id;
    }
    return *this;
}

basis::~basis() {}

long basis::factorial(long N, long m) {
    unsigned long num,denum;
    long i;
    num=1;
    for(i=N-m+1; i<=N; i++)
        num*=i;
    denum=1;
    for(i=1; i<=m; i++)
        denum*=i;
    return num/denum;
}

void basis::init() {
    long i,config_init;
    nbasis=factorial(nphi,nel);
    config_init=0;
    for(i=0; i<nel; i++)
        config_init+=(1<<i);
    generate(config_init);

    std::map<long,long>::iterator it;
    for(it=basis_set.begin(); it!=basis_set.end(); it++)
        id.push_back(it->first);

    sort(id.begin(),id.end());
    basis_set.clear();
    for(i=0; i<nbasis; i++)
        basis_set[id[i]]=i;
}

void basis::init(long _nphi, long _nel){
   nphi=_nphi;
   nel=_nel;
   init();
}


long basis::onsite_potential(long i,long n) {
    long mask,b;
    mask=(1<<n);
    b=id[i]&mask;
    if(b==mask)
        return 1;
    else
        return 0;
}


long basis::creation(long s,long n)
{
    long mask,b;
    mask=(1<<n);
    b=s&mask;
    // which means s_n=0
    if(b!=mask)
        return s+mask;
    // there's already electron on site n
    else
        return s;
}

long basis::annihilation(long s,long n)
{
    long mask,b;
    mask=(1<<n);
    b=s&mask;
    if(b==mask)
        return s-mask;
    else
        return s;
}


void basis::generate(long a) {
    long mask,K,L,b,j;
    #if __cplusplus > 199711L
    basis_set.emplace(a,a);
    #else
    basis_set[a]=a;
    #endif
    for(long i=0; i<nphi; i++) {
        j=(i+1>=nphi)?i+1-nphi:i+1;
        mask=(1<<i)+(1<<j);
        K=mask&a;
        L=K^mask;
        if(L!=0 && L!=mask) {
            b=a-K+L;
            if(basis_set.find(b)==basis_set.end())
                generate(b);
            if(basis_set.size()==nbasis)
                return;
        }
    }
    return;
}


void basis::prlong() {
    std::map<long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_set.begin(); it!=basis_set.end(); it++)
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis: "<<setw(6)<<nbasis<<endl;
    /*
    cout<<"---------------------------------------"<<endl;
    cout<<"Lin's Table:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: id_up)
       cout<<x<<endl;
    */ 
}
