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

    //cout<<nbasis*4/1.0e9<<endl;
    std::map<long,long>::iterator it;
    for(it=basis_set.begin(); it!=basis_set.end(); it++){
      if(K<0)
        id.push_back(it->first);
      else{
        long sum_j=0;
        long c=it->first;
        for(i=0;i<nphi;i++)
           if((c>>i)%2==1)
             sum_j+=i;
         if(sum_j%nphi==K)
            id.push_back(it->first);
      }
    }

    sort(id.begin(),id.end());
    basis_set.clear();
    for(i=0; i<id.size(); i++)
        basis_set[id[i]]=i;
    nbasis=id.size();
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

int basis::get_sign(long i,long n, long m, long nt, long mt){
     long b,k,kl,kr,mask,mask_k,nsign;
     mask=(1<<n)+(1<<m);
      // get the rest electrons
     b=id[i]^mask;
     // if there're no crossing between two electrons
     nsign=0;
        kl=nt<n?nt:n;
        kr=nt<n?n:nt;
        for(k=kl+1;k<kr;k++){
          mask_k=(1<<k);
          if((b&mask_k)==mask_k)
             nsign++;
        }
        kl=mt<m?mt:m;
        kr=mt<m?m:mt;
        for(k=kl+1;k<kr;k++){
          mask_k=(1<<k);
          if((b&mask_k)==mask_k)
             nsign++;
        }
        // if there're crossings between two electrons
        if(nt>mt && m>n || mt>nt && m<n)
          nsign++;

     return pow(-1,nsign);
}

void basis::prlong() {
    std::map<long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_set.begin(); it!=basis_set.end(); it++){
        long c=it->first;
        for(int n=0;n<nphi;n++){
            if((c>>n)%2==1)
               cout<<n<<" ";
        }
        cout<<"   ";
        cout<<bitset<12>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
      }
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
