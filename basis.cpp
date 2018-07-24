#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(int _nphi,int _nel):nphi(_nphi),nel(_nel){
    nel_up=-1;
    nel_down=-1;
    K=-1;
    J=-1;
}

basis::basis(int _nphi,int _nel, int _nel_up):nphi(_nphi),nel(_nel),nel_up(_nel_up) {
    K=-1;
    J=-1;
    nel_down=nel-nel_up;
}


basis::basis(int _nphi,int _nel, int _nel_up,int _J,int _K):nphi(_nphi),nel(_nel),nel_up(_nel_up),J(_J),K(_K) {
    nel_down=nel-nel_up;
}

const basis & basis::operator =(const basis & _basis) {
    if(this !=&_basis) {
        nphi=_basis.nphi;
        nel=_basis.nel;
        nel_up=_basis.nel_up;
        nbasis=_basis.nbasis;
        basis_set=_basis.basis_set;
        id=_basis.id;
        K=_basis.K;
        J=_basis.J;
    }
    return *this;
}

basis::~basis() {}

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
    int sign;
    long n,i,j,Ki,count;
    unsigned long config;
    //C=common_divisor(nphi,nel_up);
    std::unordered_map<unsigned long,long>::iterator it;
    count=j=config=Ki=0;
    if(nel_up<0)
      generate_all_density(count,j,Ki,config);
    else 
      generate(count,j,Ki,config);
    // initialize the bit sets count table
    for(i=0; i<pow(2,nphi); i++) {
        int count=0;
        for(n=0; n<nphi; n++)
            if((i>>n)%2==1)
                count++;
        popcount_table.push_back(count);
    }
    // shrink the up-layer basis to an unique subset L
    // and with nel_up/nel_down in each layer
    if(K>=0) {
        for(it=basis_set.begin(); it!=basis_set.end();) {
            unsigned long c=it->first;
            bool delete_flag=false;
            for(n=1; n<nphi; n++) {
		config=inv_translate(c,n,sign);
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

    for(it=basis_set.begin(); it!=basis_set.end(); it++){
        id.push_back(it->second);
    }

    sort(id.begin(),id.end());
    basis_set.clear();
    for(i=0; i<id.size(); i++){
        basis_set.insert(pair<unsigned long, long>(id[i],i));

	}
    nbasis=id.size();
    basis_C.assign(nbasis,nphi);
    for(i=0;i<id.size();i++)
       for(n=1; n<nphi; n++)
            if(inv_translate(id[i],n,sign)==id[i]) {
                basis_C[i]=n;
                break;
    }

}

void basis::init(int _nphi,int _nel){
    nphi=_nphi;
    nel=_nel; 
    nel_up=-1;
    nel_down=-1;
    K=-1;
    J=-1;
    init();
}

void basis::init(int _nphi, int _nel, int _nel_up) {
    nphi=_nphi;
    nel_up=_nel_up;
    nel=_nel;
    K=-1;
    J=-1;
    init();
}

void basis::init(int _nphi,int _nel, int _nel_up,int _J,int _K) {
    nphi=_nphi;
    nel_up=_nel_up;
    nel=_nel;
    J=_J;
    K=_K;
    nel_down=nel-nel_up;
    init();
}

void basis::clear(){
    basis_set.clear();
    id.clear();
    popcount_table.clear(); 
}


int basis::get_sign(unsigned long c,int n, int m, int nt, int mt,int t) {
    int k,kl,kr,nsign,sign;
    unsigned long b,mask, mask_k;
    mask=(1<<n)+(1<<m);
    // get the rest electrons
    b=c^mask;
    // if there're no crossing between two electrons
    nsign=0;
    kl=nt<n?nt:n;
    kr=nt<n?n:nt;
    for(k=kl+1; k<kr; k++) {
      if((b>>k)%2==0)
        nsign++;
    }
    kl=mt<m?mt:m;
    kr=mt<m?m:mt;
    for(k=kl+1; k<kr; k++) {
      if((b>>k)%2==0)
        nsign++;
    }
    // if there're crossings between two electrons
    if(nt>mt && m>n || mt>nt && m<n)
       nsign++;
     
    sign=(nsign%2==0?1:-1);
    return sign;
}

// sign change function for electron hopping from one layer to another
int basis::get_sign(unsigned long basis_i, int n,int nt){
    int nsign,sign;
    nsign=0;
    for(int i=n+1;i<nt;i++)
	if((basis_i>>i)%2==0)
	 nsign++;

    sign=(nsign%2==0?1:-1);
    return sign;
}


void basis::generate(long count,long j, long Ji, unsigned long config) {
    long i,k,range;
    unsigned long c;

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

void basis::generate_all_density(long count,long j, long Ji, unsigned long config) {
  int i,id,k,c;
  count++;
  config=(count==1?0:config);
  j=(count==1?0:j+1);
  for(i=j;i<2*nphi-(nel-count);i++){
    // total sum of k
    k=Ji+i%nphi;
    c=config+(1<<i);
    if(count<nel){
       generate_all_density(count,i,k,c);
     }
    else{
      if(J<0)
        basis_set[c]=c;
      else if(k%nphi==J)
        basis_set[c]=c;
    }
  }
}

int basis::get_nel_upper(long i){
    unsigned long mask_u,c_u;
    mask_u =(1<<nphi)-1;
    c_u=id[i] & mask_u;
    return popcount_table[c_u & mask_u]; 
}

void basis::prlong() {
    std::unordered_map<unsigned long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    int count=0;
    for(it=basis_set.begin(); it!=basis_set.end(); it++) {
        unsigned long c=it->first;
        cout<<" | ";
        for(int n=0; n<nphi; n++) {
            if((c>>n)%2==1)
                cout<<n<<" ";
        }
        cout<<">| ";
        for(int n=nphi; n<2*nphi; n++) {
            if((c>>n)%2==1)
                cout<<n-nphi<<" ";
        }
        cout<<"> ";
        cout<<bitset<12>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<" "<<basis_C[count]<<endl;
	count++;
    }
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis: "<<setw(6)<<nbasis<<endl;
}

unsigned long basis::translate(unsigned long c, int k, int &sign) {
        unsigned long config,mask_d,mask_u,c_d,c_u,mask_sign;
        int nsign;
        int bits=k;
        int inv_bits=nphi-bits;

        mask_u=(1<<nphi)-1;
        mask_d=mask_u<<nphi;
        c_u= c & mask_u;
        c_d= (c & mask_d)>>nphi;

        mask_sign=(1<<bits)-1;
        int ncross_up=popcount_table[(mask_sign<<inv_bits) & c_u];
        int ncross_down=popcount_table[(mask_sign<<inv_bits) & c_d];
        int _nel_up=popcount_table[c_u];
        int _nel_down=nel-_nel_up;
         
        nsign=(_nel_up-ncross_up)*ncross_up+(_nel_down-ncross_down)*ncross_down;
        c_u=((c_u<<bits)|(c_u>>inv_bits))&mask_u;
        c_d=((c_d<<bits)|(c_d>>inv_bits))&mask_u;

        config=((c_d<<nphi)|c_u);
        sign=(nsign%2==0?1:-1);
        return config;
    }

unsigned long basis::inv_translate(unsigned long c, int k, int &sign) {
        unsigned long config,mask_d,mask_u,c_d,c_u,mask_sign;
        int nsign;
        int bits=k;
        int inv_bits=nphi-bits;

        mask_u=(1<<nphi)-1;
        mask_d=mask_u<<nphi;
        c_u= c & mask_u;
        c_d=(c & mask_d)>>nphi;

        mask_sign=(1<<bits)-1;
        int ncross_up=popcount_table[mask_sign & c_u];
        int ncross_down=popcount_table[mask_sign & c_d];
        int _nel_up=popcount_table[c_u];
        int _nel_down=nel-_nel_up;
        nsign=(_nel_up-ncross_up)*ncross_up+(_nel_down-ncross_down)*ncross_down;
        c_u=((c_u>>bits)|(c_u<<inv_bits))&mask_u;
        c_d=((c_d>>bits)|(c_d<<inv_bits))&mask_u;
        config=((c_d<<nphi)|c_u);
        sign=(nsign%2==0?1:-1);
        return config;
    }

