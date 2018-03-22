#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(long _nsite,long _nel_up, long _nel_down):nsite(_nsite),nel_up(_nel_up),nel_down(_nel_down) {
}

const basis & basis::operator =(const basis & _basis) {
    if(this !=&_basis) {
        nsite=_basis.nsite;
        nel_up=_basis.nel_up;
        nel_down=_basis.nel_down;
        nbasis_up=_basis.nbasis_up;
        nbasis_down=_basis.nbasis_down;
        basis_up=_basis.basis_up;
        basis_down=_basis.basis_down;
        id_up=_basis.id_up;
        id_down=_basis.id_down;
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
    nbasis_up=factorial(nsite,nel_up);
    nbasis_down=factorial(nsite,nel_down);
    config_init=0;
    for(i=0; i<nel_up; i++)
        config_init+=(1<<i);
    generate_up(config_init);

    config_init=0;
    for(i=0; i<nel_down; i++)
        config_init+=(1<<i);
    generate_down(config_init);

    std::map<long,long>::iterator it;
    for(it=basis_up.begin(); it!=basis_up.end(); it++)
        id_up.push_back(it->first);
    for(it=basis_down.begin(); it!=basis_down.end(); it++)
        id_down.push_back(it->first);

    sort(id_up.begin(),id_up.end());
    sort(id_down.begin(),id_down.end());
    basis_up.clear();
    basis_down.clear();
    for(i=0; i<nbasis_up; i++)
        basis_up[id_up[i]]=i;
    for(i=0; i<nbasis_down; i++)
        basis_down[id_down[i]]=i;
}

void basis::init(long _nsite, long _nel_up, long _nel_down){
   nsite=_nsite;
   nel_up=_nel_up;
   nel_down=_nel_down;
   init();
}

long basis::interlayer_hopping(long i,long n,long m) {
    long mask,K,L,b;
    if(m<0) m+=nsite;
    else if (m>=nsite) m-=nsite;

    mask=(1<<n)+(1<<m);
    K=mask&id_up[i];
    L=K^mask;

    if(L!=0 && L!=mask) {
        b=id_up[i]-K+L;
        if(basis_up.find(b)!=basis_up.end())
            return basis_up[b];
        else
            return basis_up[id_up[i]];
    }
    else
        return basis_up[id_up[i]];
}


long basis::onsite_potential(long i,long j,long n) {
    long mask,bu,bd;
    mask=(1<<n);
    bu=id_up[i]&mask;
    bd=id_down[j]&mask;
    if(bu==mask && bd==mask)
        return 1;
    else
        return 0;
}


long basis::creation(long s,long n)
{
    long mask,bu;
    mask=(1<<n);
    bu=s&mask;
    // there's no electron on site n
    // which means s_n=0
    if(bu!=mask)
        return s+mask;
    // there's already electron on site n
    else
        return s;
}

long basis::annihilation(long s,long n)
{
    long mask,bu;
    mask=(1<<n);
    bu=s&mask;
    if(bu==mask)
        return s-mask;
    else
        return s;
}


void basis::generate_up(long a) {
    long mask,K,L,b,j;
#if __cplusplus > 199711L
    basis_up.emplace(a,a);
#else
    basis_up[a]=a;
#endif
    for(long i=0; i<nsite; i++) {
        j=(i+1>=nsite)?i+1-nsite:i+1;
        mask=(1<<i)+(1<<j);
        K=mask&a;
        L=K^mask;
        if(L!=0 && L!=mask) {
            b=a-K+L;
            if(basis_up.find(b)==basis_up.end())
                generate_up(b);
            if(basis_up.size()==nbasis_up)
                return;
        }
    }
    return;
}

void basis::generate_down(long a) {
    long mask,K,L,b,j;
#if __cplusplus > 199711L
    basis_down.emplace(a,a);
#else
    basis_down[a]=a;
#endif
    for(long i=0; i<nsite-1; i++) {
        j=(i+1>=nsite)?i+1-nsite:i+1;
        mask=(1<<i)+(1<<j);
        K=mask&a;
        L=K^mask;
        if(L!=0 && L!=mask) {
            b=a-K+L;
            if(basis_down.find(b)==basis_down.end())
                generate_down(b);
            if(basis_down.size()==nbasis_down)
                return;
        }
    }
}

void basis::prlong() {
    std::map<long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_up.begin(); it!=basis_up.end(); it++)
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_down.begin(); it!=basis_down.end(); it++)
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis for spin-up electrons: "<<setw(6)<<nbasis_up<<endl;
    cout<<"No. basis for spin-down electrons: "<<setw(6)<<nbasis_down<<endl;
    /*
    cout<<"---------------------------------------"<<endl;
    cout<<"Lin's Table:"<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: id_up)
       cout<<x<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: id_down)
       cout<<x<<endl;
    cout<<"---------------------------------------"<<endl;
    */
}
