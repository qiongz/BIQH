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
		config=translate(c,n,sign);
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
    for(i=0; i<id.size(); i++)
        basis_set.insert(pair<unsigned long, long>(id[i],i));
    nbasis=id.size();

    basis_C.assign(nbasis,nphi);
    for(i=0;i<id.size();i++)
       for(n=1; n<nphi; n++)
            if(translate(id[i],n,sign)==id[i]) {
                basis_C[i]=n;
                break;
    }
}

void basis::init(int _nphi, int _nel_up, int _nel_down) {
    nphi=_nphi;
    nel_up=_nel_up;
    nel_down=_nel_down;
    init();
}

void basis::init(int _nphi,int _nel_up, int _nel_down,int _J,int _K) {
    nphi=_nphi;
    nel_up=_nel_up;
    nel_down=_nel_down;
    J=_J;
    K=_K;
    nel=nel_up+nel_down;
    init();
}

void basis::clear(){
    basis_set.clear();
    id.clear();
    popcount_table.clear(); 
}


int basis::get_sign(unsigned long c,int n, int m, int nt, int mt) {
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
        mask_k=(1<<k);
        if((b&mask_k)==mask_k)
            nsign++;
    }
    kl=mt<m?mt:m;
    kr=mt<m?m:mt;
    for(k=kl+1; k<kr; k++) {
        mask_k=(1<<k);
        if((b&mask_k)==mask_k)
            nsign++;
    }
    // if there're crossings between two electrons
    if(nt>mt && m>n || mt>nt && m<n)
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
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<" "<<basis_C[count]<<endl;
	count++;
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
unsigned long basis::translate(unsigned long c, int k_u, int &sign) {
        unsigned long config,mask_d,mask_u,c_d,c_u,mask_sign_u,mask_sign_d;
        int nsign;
	//int k_d=nel_down*k_u/nel_up;
	int k_d=k_u;
        int bits_u=k_u;
        int bits_d=k_d;

        mask_d=(1<<nphi)-1;
        mask_u=mask_d<<nphi;
        c_d= c & mask_d;
        c_u=(c & mask_u)>>nphi;

        mask_sign_u=(1<<bits_u)-1;
        mask_sign_d=(1<<bits_d)-1;
        nsign=popcount_table[(mask_sign_u<<(nphi-bits_u)) & c_u]*(nel_up-1)+popcount_table[(mask_sign_d<<(nphi-bits_d)) & c_d]*(nel_down-1);

        c_u=((c_u<<bits_u)|(c_u>>(nphi-bits_u)))&mask_d;
        c_d=((c_d<<bits_d)|(c_d>>(nphi-bits_d)))&mask_d;

        config=((c_u<<nphi)|c_d);
        sign=(nsign%2==0?1:-1);
        return config;
    }

unsigned long basis::inv_translate(unsigned long c, int k_u, int &sign) {
        unsigned long config,mask_d,mask_u,c_d,c_u,mask_sign_u,mask_sign_d;
        int nsign;
	//int k_d=nel_down*k_u/nel_up;
	int k_d=k_u;
        int bits_u=k_u;
        int bits_d=k_d;

        mask_d=(1<<nphi)-1;
        mask_u=mask_d<<nphi;
        c_d= c & mask_d;
        c_u=(c & mask_u)>>nphi;

        mask_sign_u=(1<<bits_u)-1;
        mask_sign_d=(1<<bits_d)-1;
        nsign=popcount_table[mask_sign_u & c_u]*(nel_up-1)+popcount_table[mask_sign_d & c_d]*(nel_down-1);

        c_u=((c_u>>bits_u)|(c_u<<(nphi-bits_u)))&mask_d;
        c_d=((c_d>>bits_d)|(c_d<<(nphi-bits_d)))&mask_d;
        config=((c_u<<nphi)|c_d);
        sign=(nsign%2==0?1:-1);
        return config;
    }

