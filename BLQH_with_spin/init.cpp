#include"init.h"
void usage(char *target) {
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -L                       nLL\n";
    std::cout<<"  -n                       nphi\n";
    std::cout<<"  -e                       No. of electrons\n";
    std::cout<<"  -s                       Delta_SAS: tunnelling amplitude\n";
    std::cout<<"  -v                       Delta_V: bias voltage\n";
    std::cout<<"  -z                       Delta_Z: Zeemann energy\n";
    std::cout<<"  -j                       total J (k_y)\n";
    std::cout<<"  -k                       kx \n";
    std::cout<<"  -g                       gamma=lx/ly  aspect ratio\n";
    std::cout<<"  -d                       interlayer distance\n";
    std::cout<<"  -m                       Lambda\n";
    std::cout<<"  -t                       nthread\n";
}

void init_argv(int &nLL,int &nphi, int& nel, int &J, int &kx, double &d,double &Delta_SAS, double &Delta_V, double &Delta_Z,double &gamma ,int &lambda,double & theta,int &nthread,int argc,char *argv[])
{
    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"L:e:n:u:d:j:k:g:z:m:t:s:v:h:"))!=-1) {
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
        case 's':
	    Delta_SAS=atof(optarg);
	    break;
        case 'v':
	    Delta_V=atof(optarg);
	    break;
        case 'z':
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
        if(nel<0)
            throw std::logic_error("-n: positive value required !");
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
