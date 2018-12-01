#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -l                       nLL\n";
    std::cout<<"  -n                       nphi\n";
    std::cout<<"  -e                       Total No. of electrons in upper-layer\n";
    std::cout<<"  -u                       No. of electrons in upper-layer\n";
    std::cout<<"  -k                       kx in upper-layer\n";
    std::cout<<"  -S                       Delta_SAS: tunnelling amplitude\n";
    std::cout<<"  -v                       Delta_V: bias voltage\n";
    std::cout<<"  -j                       total J in upper-layer or down-layer\n";
    std::cout<<"  -g                       gamma=lx/ly  aspect ratio\n";
    std::cout<<"  -d                       interlayer distance\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"  -t                       nthread\n";
    std::cout<<"  -s                       seed\n";
    std::cout<<"Default: (l,n,u,d,lambda) = (4,4,2,1,200)\n";
}

void init_argv(int &nLL,int &nphi, int& nel, int &nel_up, int &J, int &kx, double &d,double &Delta_SAS, double &Delta_V,double &gamma ,int &lambda,int &nthread,unsigned long seed,int argc,char *argv[])
{
    nLL=0;
    nphi=8;
    nel=8;
    nel_up=4;
    gamma=1.0;
    nthread=32;
    J=-1;
    kx=-1;
    Delta_SAS=0;
    Delta_V=0;
    d=100.0;
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
    while((ch=getopt(argc,argv,"L:e:n:u:d:j:k:g:m:t:s:v:h:S:"))!=-1) {
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
        case 'u':
            nel_up=atoi(optarg);
            break;
        case 'k':
            kx=atoi(optarg);
            break;
        case 'j':
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
        case 'v':
	    Delta_V=atof(optarg);
	    break;
        case 'm':
            lambda=atoi(optarg);
            break;
        case 't':
            nthread=atoi(optarg);
            break;
        case 's':
            seed=atoi(optarg);
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
        if(nel<0)
            throw std::logic_error("-n: positive value required !");
        //if(nel_up<0)
        //    throw std::logic_error("-u: positive value required !");
        if(nel_up>nel)
            throw std::logic_error("-u: nel_up < nel !");
        //if(fabs(d)<1e-8)
        //    throw std::logic_error("-d: at least one finite coupling constant required !");
    } catch(std::logic_error &e) {
        std::cout<<e.what()<<std::endl;
        usage(argv[0]);
        exit(2);
    }
    catch(std::overflow_error &e) {
        std::cout<<e.what()<<std::endl;
        exit(2);
    }

    if(errFlag) {
        usage(argv[0]);
        exit(0);
    }
}
