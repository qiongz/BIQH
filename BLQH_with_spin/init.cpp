#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -L                       nLL\n";
    std::cout<<"  -n                       nphi\n";
    std::cout<<"  -e                       No. of electrons\n";
    std::cout<<"  -U                       No. of upper-layer electrons\n";
    std::cout<<"  -u                       No. of spin-up electrons\n";
    std::cout<<"  -S                       Delta_SAS: tunnelling amplitude\n";
    std::cout<<"  -V                       Delta_V: bias voltage\n";
    std::cout<<"  -Z                       Delta_Z: Zeemann energy\n";
    std::cout<<"  -J                       total J (K_y)\n";
    std::cout<<"  -K                       Kx \n";
    std::cout<<"  -g                       gamma=lx/ly  aspect ratio\n";
    std::cout<<"  -d                       interlayer distance\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"  -s                       seed\n";
    std::cout<<"  -t                       No. threads\n";
}

void init_argv(int &nLL,int &nphi, int& nel, int &nel_up,int &nel_su,int &J, int &kx, double &d,double &Delta_SAS, double &Delta_V, double &Delta_Z,double &gamma ,int &lambda,double & theta,int &nthread,unsigned long & seed,int argc,char *argv[])
{
    nLL=0;
    nphi=4;
    nel=8;
    nel_up=-1;
    nel_su=-1;
    gamma=1.0;
    theta=0;
    nthread=8;
    J=-1;
    kx=-1;
    Delta_SAS=0;
    Delta_Z=0;
    Delta_V=0;
    d=1.0;
    lambda=200;

    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif

    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"L:e:n:U:u:d:J:K:g:Z:m:t:S:s:V:h:"))!=-1) {
        switch(ch) {
        case 'L':
            nLL=atoi(optarg);
            break;
        case 'n':
            nphi=atoi(optarg);
            break;
        case 'e':
            nel=atoi(optarg);
            break;
        case 'U':
            nel_up=atoi(optarg);
            break;
        case 'u':
            nel_su=atoi(optarg);
            break;
        case 'K':
            kx=atoi(optarg);
            break;
        case 'J':
            J=atoi(optarg);
            break;
        case 'g':
            gamma=atof(optarg);
            break;
        case 'd':
            d=atof(optarg);
            break;
        case 'S':
	    Delta_SAS=atof(optarg);
	    break;
        case 's':
	    seed=atof(optarg);
	    break;
        case 'V':
	    Delta_V=atof(optarg);
	    break;
        case 'Z':
	    Delta_Z=atof(optarg);
	    break;
        case 'm':
            lambda=atoi(optarg);
            break;
	case 'T':
	    theta=atof(optarg);
	    break;
        case 't':
            nthread=atoi(optarg);
            break;
        case 'h':
            errFlag++;
            break;
        default:
            errFlag++;
            break;
        }
    }
    try {
        if(nphi<0 )
            throw std::logic_error("-n: positive value required !");
    } 
    catch(std::logic_error &n) {
        std::cout<<n.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    try {
	if(nel<=0 || nel>=4*nphi)
            throw std::logic_error("-e: Total No. of electrons should be in (0,4*nphi)");
    }
    catch(std::logic_error &e) {
        std::cout<<e.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    try {
	if(nel_up<0 || nel_up>=nel)
            throw std::logic_error("-U: No. of upper-layer electrons should be in [0,nel_up] !");
    }
    catch(std::logic_error &u) {
        std::cout<<u.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    try {
	if(nel_su<0 || nel_su>=nel)
            throw std::logic_error("-u: No. of spin-up electrons should be in [0,nel_up] !");
    }
    catch(std::logic_error &s) {
        std::cout<<s.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    try{
	if(nphi>21 && nel>2*nphi/3 && nel<nphi*4/3)
            throw std::overflow_error("Hilbert space too large !");

    }
    catch(std::overflow_error &h) {
        std::cout<<h.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }

    if(errFlag) {
        usage(argv[0]);
        exit(0);
    }
}
