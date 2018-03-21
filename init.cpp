#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -l                       No. of sites\n";
    std::cout<<"  -n                       Total No. of electrons in upper-layer\n";
    std::cout<<"  -u                       No. of electrons in upper-layer\n";
    std::cout<<"  -d                       interlayer distance\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"Default: (l,n,u,d,lambda) = (4,4,2,1,200)\n";
}

void init_argv(int& nsite,int& nel, int &nel_up, double &d, int &lambda,int argc,char *argv[])
{
    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:n:u:d:m:h:"))!=-1) {
        switch(ch) {
        case 'l':
            nsite=atoi(optarg);
            break;
        case 'n':
            nel=atoi(optarg);
            break;
        case 'u':
            nel_up=atoi(optarg);
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
        if(nsite>=16 && nel==nsite/2)
            throw std::overflow_error("-l -n:Hilbert space larger than the memory space !");
        if(nsite<=0)
            throw std::logic_error("-l: positive value required !");
        if(nel<=0)
            throw std::logic_error("-n: positive value required !");
        if(nel_up<=0)
            throw std::logic_error("-u: positive value required !");
        if(nel_up>nel)
            throw std::logic_error("-u: nel_up < nel !");
        if(nsite*2<nel)
            throw std::logic_error("-l -n: electron filling nu_t<=1 required !");
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
