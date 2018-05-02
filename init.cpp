#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -l                       nLL\n";
    std::cout<<"  -p                       nphi\n";
    std::cout<<"  -n                       Total No. of electrons in upper-layer\n";
    std::cout<<"  -u                       No. of electrons in upper-layer\n";
    std::cout<<"  -s                       kx in upper-layer\n";
    std::cout<<"  -t                       total J in upper-layer or down-layer\n";
    std::cout<<"  -g                       gamma=lx/ly  aspect ratio\n";
    std::cout<<"  -d                       interlayer distance\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"Default: (l,n,u,d,lambda) = (4,4,2,1,200)\n";
}

void init_argv(int &nLL,int &nphi, int& nel, int &nel_up, int &J, int &kx, double &d,double &gamma ,int &lambda,int argc,char *argv[])
{
    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:p:n:u:d:s:t:g:m:h:"))!=-1) {
        switch(ch) {
        case 'l':
            nLL=atoi(optarg);
            break;
        case 'p':
            nphi=atoi(optarg);
            break;
        case 'n':
            nel=atoi(optarg);
            break;
        case 'u':
            nel_up=atoi(optarg);
            break;
        case 's':
            kx=atoi(optarg);
            break;
        case 't':
            J=atoi(optarg);
            break;
        case 'g':
            gamma=atof(optarg);
            break;
        case 'd':
            d=atof(optarg);
            break;
        case 'm':
            lambda=atoi(optarg);
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
        if(nel_up<0)
            throw std::logic_error("-u: positive value required !");
        if(nel_up>nel)
            throw std::logic_error("-u: nel_up < nel !");
        if(fabs(d)<1e-8)
            throw std::logic_error("-d: at least one finite coupling constant required !");
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
