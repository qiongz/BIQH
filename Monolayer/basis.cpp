#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(long _nphi,long _nel):nphi(_nphi),nel(_nel){
  J=-1;
  K=-1;
}

basis::basis(long _nphi,long _nel, long _J,long _K):nphi(_nphi),nel(_nel),J(_J),K(_K) {
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

void basis::clear(){
   if(nbasis!=0){
     basis_set.clear();
     id.clear();
   }
}

unsigned long basis::factorial(long N, long m) {
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

long basis::common_divisor(long n, long m){
    long N=1;
    for(int i=1;i<=n;i++)
      if(n%i==0 && m%i==0)
         N=i;
    return N;
}

void basis::init() {
    long i,j,n,Ji,count;
    unsigned long config;
    C=common_divisor(nel,nphi);
    long q=nphi/C;
    count=j=config=Ji=0;
    generate(count,j,Ji,config);
    std::unordered_map<unsigned long,long>::iterator it;
    // shrink to an unique subset L, all other subset could be generated with translation
    if(K>=0)
    {
        unsigned long mask;
        for(it=basis_set.begin(); it!=basis_set.end();) {
            unsigned long c=it->first;
            // generate the read mask with field length nphi, and shift to the position
            mask=(1<<nphi)-1;
            c=c & mask;
            bool delete_flag=false;
            for(n=1; n<C; n++) {
                config=((c>>(q*n))|(c<<(nphi-q*n)))&mask ;
                // right rotate the bits in the down-layer: c_down
                // if translation could generate this configuration, delete it
                if(basis_set.find(config)!=basis_set.end()&& config!=c) {
                    delete_flag=true;
                    break;
                }
            }
            if(delete_flag)
                basis_set.erase(it++);
            else
                ++it;
        }
    }

    for(it=basis_set.begin(); it!=basis_set.end(); it++)
        id.push_back(it->second);
    sort(id.begin(),id.end());
    basis_set.clear();
    for(i=0; i<id.size(); i++)
        basis_set.insert(pair<unsigned long, long>(id[i],i));
    nbasis=id.size();
    

    // initialize the bit sets count table
    for(i=0; i<pow(2,nphi); i++) {
        int count=0;
        for(n=0; n<nphi; n++)
            if((i>>n)%2==1)
                count++;
        popcount_table.push_back(count);
    }
    // initialize the sector translation size for each basis
   int sign;
    basis_C.assign(nbasis,C);
    for(i=0;i<id.size();i++)
       for(n=1; n<C; n++)
            if(translate(id[i],n,sign)==id[i]) {
                basis_C[i]=n;
                break;
    }
}

void basis::init(long _nphi, long _nel){
   nphi=_nphi;
   nel=_nel;
   J=-1;
   K=-1;
   init();
}

void basis::init(long _nphi, long _nel, long _J,long _K){
   nphi=_nphi;
   nel=_nel;
   J=_J;
   K=_K;
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

int basis::get_sign(long c,long n, long m, long nt, long mt){
     long b,k,kl,kr,mask,mask_k,nsign;
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

// for loop with range and number of iterations as an argument
void basis::for_sum(long &a,int count,int index,int range,int n_iter){
  int i,j,id;
  count++;
  j=(count==1?index%range:index%range+1);
  for(i=j;i<range-(n_iter-count);i++){
    id=i+index*range;
    if(count<n_iter){
       for_sum(a,count,id,range,n_iter);
     }
    else{
      a+=id;
      //cout<<"sum: "<<i<<" "<<j<<" "<<id<<endl;
    }
  }
}


void  basis::generate(long count,long j,long Ji,unsigned long config){
  int i,id,k,c;
  count++;
  config=(count==1?0:config);
  j=(count==1?0:j+1);
  for(i=j;i<nphi-(nel-count);i++){
    // total sum of k
    k=Ji+i;
    c=config+(1<<i);
    if(count<nel){
       generate(count,i,k,c);
     }
    else{
      //cout<<count<<" "<<i<<" "<<j<<" "<<c<<" "<<bitset<12>(c).to_string()<<endl;
      if(J<0)
        basis_set[c]=c;
      else if(k%nphi==J)
        basis_set[c]=c;
    }
  }
}


void basis::prlong() {
    std::unordered_map<unsigned long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_set.begin(); it!=basis_set.end(); it++){
        long c=it->first;
        for(int n=0;n<nphi;n++){
            if((c>>n)%2==1)
               cout<<n<<" ";
        }
        cout<<"   ";
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
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

unsigned long basis::translate(unsigned long c, int k, int &sign) {
        unsigned long config,mask,mask_sign;
        int nsign;
        int bits=k*nphi/C;
        int reverse_bits=nphi-bits;
        mask_sign=(1<<bits)-1;
        mask=(1<<nphi)-1;
        c= c & mask;
        nsign=popcount_table[(mask_sign<<reverse_bits) & c]*(nel-1);
        config=((c<<bits)|(c>>reverse_bits))&mask;
        sign=(nsign%2==0?1:-1);
        return config;
    }

unsigned long basis::inv_translate(unsigned long c, int k,int &sign) {
        unsigned long config,mask,mask_sign;
        int nsign;
        int bits=k*nphi/C;
        int reverse_bits=nphi-bits;
        mask_sign=(1<<bits)-1;
        mask=(1<<nphi)-1;
        c=c&mask;
        nsign=popcount_table[mask_sign & c]*(nel-1);
        config=((c>>bits)|(c<<reverse_bits))&mask;
        sign=(nsign%2==0?1:-1);
        return config;
    }
