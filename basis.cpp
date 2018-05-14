#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(int _nphi,int _nel_up, int _nel_down):nphi(_nphi),nel_up(_nel_up),nel_down(_nel_down) {
    K=-1;
    J=-1;
    nel=nel_up+nel_down;
}


basis::basis(int _nphi,int _nel_up, int _nel_down,int _J,int _K):nphi(_nphi),nel_up(_nel_up),nel_down(_nel_down),J(_J),K(_K) {
    C=common_divisor(nphi,nel_up);
    nel=nel_up+nel_down;
}

const basis & basis::operator =(const basis & _basis) {
    if(this !=&_basis) {
        nphi=_basis.nphi;
        nel_up=_basis.nel_up;
        nel_down=_basis.nel_down;
        nbasis=_basis.nbasis;
        basis_set=_basis.basis_set;
        id=_basis.id;
        K=_basis.K;
        J=_basis.J;
    }
    return *this;
}

basis::~basis() {}

void basis::clear() {
    basis_set.clear();
    id.clear();
}

unsigned long basis::factorial(int N, int m) {
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

long basis::common_divisor(int n, int m) {
    long N=1;
    for(int i=1; i<=n; i++)
        if(n%i == 0 && m%i==0)
            N=i;
    return N;
}

void basis::init() {
    long n,i,j,Ki,count,q,C;
    unsigned long long config;
    std::map<unsigned long long,unsigned long long>::iterator it;
    count=j=config=Ki=0;
    generate(count,j,Ki,config);

    // shrink the up-layer basis to an unique subset L
    // and with nel_up/nel_down in each layer
    C=common_divisor(nphi,nel_up);
    q=nphi/C;

    if(K>=0) {
        vector<long> cv_up,cv_down;
        for(it=basis_set.begin(); it!=basis_set.end();) {
            unsigned long long c=it->first;
            for(n=0; n<nphi; n++)
                if((c>>n)%2==1)
                    cv_up.push_back(n);
            for(n=nphi; n<2*nphi; n++)
                if((c>>n)%2==1)
                    cv_down.push_back(n);
            // if translation could generate this configuration, delete it
            bool delete_flag=false;
            for(n=1; n<C; n++) {
                config=0;
                // translate in the upper-layer
                for(i=0; i<cv_up.size(); i++) {
                    j=(cv_up[i]+q*n>=nphi?cv_up[i]+q*n-nphi:cv_up[i]+q*n);
                    config+=(1<<j);
                }
                // inverse translate in the down-layer
                for(i=0; i<cv_down.size(); i++) {
                    j=(cv_down[i]-q*n<nphi?cv_down[i]-q*n+nphi:cv_down[i]-q*n);
                    config+=(1<<j);
                }
                if(basis_set.find(config)!=basis_set.end()&& config!=c) {
                    delete_flag=true;
                    break;
                }
            }
            if(delete_flag)
                basis_set.erase(it++);
            else
                ++it;
            cv_up.clear();
            cv_down.clear();
        }
    }

    for(it=basis_set.begin(); it!=basis_set.end(); it++)
        id.push_back(it->second);
    sort(id.begin(),id.end());
    basis_set.clear();
    for(i=0; i<id.size(); i++)
        basis_set[id[i]]=i;
    nbasis=id.size();
}

void basis::init(int _nphi, int _nel_up, int _nel_down) {
    nphi=_nphi;
    nel_up=_nel_up;
    nel_down=_nel_down;
    init();
}

void basis::init(int _nphi,int _nel_up, int _nel_down,int _J,int _K){
    nphi=_nphi;
    nel_up=_nel_up;
    nel_down=_nel_down;
    J=_J;
    K=_K;
    C=common_divisor(nphi,nel_up);
    nel=nel_up+nel_down;
    init();
}

/*
long basis::interlayer_hopping(long i,long n,long m) {
    long mask,K,L,b;
    if(m<0) m+=nphi;
    else if (m>=nphi) m-=nphi;

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
*/

/*
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
*/


int basis::get_sign(unsigned long long c,int n, int m, int nt, int mt){
     int k,kl,kr,nsign;
     unsigned long long b,mask, mask_k;
     mask=(1<<n)+(1<<m);
      // get the rest electrons
     b=c^mask;
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

/*
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
*/


void basis::generate(long count,long j, long Ji, unsigned long long config) {
    long i,k,range;
    unsigned long long c;

    count++;
    config=(count==1?0:config);
    // start of upper-layer electron momentum
    if(count==1 && nel_up>0)
        j=0;
    // start of down-layer electron momentum
    else if(count==nel_up+1)
        j=nphi;
    else
        j=j+1;

    // upper-layer electron momentum range
    if(j<nphi && nel_up>0)
        range=nphi-(nel_up-count);
    // down-layer electron momentum range
    else
        range=2*nphi-(nel_down-(count-nel_up));

    for(i=j; i<range; i++) {
        // total sum of k
        k=Ji+i%nphi;
        c=config+(1<<i);
        // if the No. of upper-layer electron != nel_up
        if(count<nel) {
            generate(count,i,k,c);
        }
        else {
            if(J<0)
                basis_set[c]=c;
            else if(k%nphi==J)
                basis_set[c]=c;
        }
    }
}


unsigned long long basis::translate(unsigned long long c, int k, int &sign){
     long j,n;
     unsigned long long config;
     vector<long> cv_up,cv_down;
     int q=nphi/C;
     for(n=0;n<nphi;n++)
        if((c>>n)%2==1)
           cv_up.push_back(n);
     for(n=nphi;n<2*nphi;n++)
        if((c>>n)%2==1)
           cv_down.push_back(n);

     config=0;
     sign=0;
     for(n=0;n<cv_up.size();n++){
      if(cv_up[n]+q*k>=nphi){
        j=cv_up[n]+q*k-nphi;
        sign+=nel_up-1;
        }
      else
         j=cv_up[n]+q*k;
      config+=(1<<j);
     }

     for(n=0;n<cv_down.size();n++){
      if(cv_down[n]-q*k<nphi){
        j=cv_down[n]-q*k+nphi;
        sign+=nel_down-1;
        }
      else
         j=cv_down[n]-q*k;
      config+=(1<<j);
     }

     sign=pow(-1,sign);
     cv_up.clear();
     cv_down.clear();
     return config;
}

unsigned long long basis::inv_translate(unsigned long long c, int k , int &sign){
     long j,n;
     unsigned long long config;
     vector<long> cv_up,cv_down;
     int q=nphi/C;
     for(n=0;n<nphi;n++)
        if((c>>n)%2==1)
           cv_up.push_back(n);
     for(n=nphi;n<2*nphi;n++)
        if((c>>n)%2==1)
           cv_down.push_back(n);

     config=0;
     sign=0;
     for(n=0;n<cv_up.size();n++){
      if(cv_up[n]-q*k<0){
        j=cv_up[n]-q*k+nphi;
        sign+=nel_up-1;
        }
      else
         j=cv_up[n]-q*k;
      config+=(1<<j);
     }

     for(n=0;n<cv_down.size();n++){
      if(cv_down[n]+q*k>=2*nphi){
        j=cv_down[n]+q*k-nphi;
        sign+=nel_down-1;
        }
      else
         j=cv_down[n]+q*k;
      config+=(1<<j);
     }

     sign=pow(-1,sign);
     cv_up.clear();
     cv_down.clear();
     return config;
}


void basis::prlong() {
    std::map<unsigned long long,unsigned long long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_set.begin(); it!=basis_set.end(); it++) {
        unsigned long long c=it->first;
        cout<<" | ";
        for(int n=0; n<nphi; n++) {
            if((c>>n)%2==1)
                cout<<n<<" ";
        }
        cout<<">| ";
        for(int n=nphi; n<2*nphi; n++) {
            if((c>>n)%2==1)
                cout<<n<<" ";
        }
        cout<<"> ";
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    }
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis: "<<setw(6)<<nbasis<<endl;
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
