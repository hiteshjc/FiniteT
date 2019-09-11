#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER

#include"global.h"
using namespace std;

class Ham{

    public:
    double spin;
    int num_sites;
    std::vector< std::vector<int> > pairs_list;
    std::vector< std::vector<int> > pairs1;
    std::vector< std::vector<int> > pairs2;
    std::vector< std::vector<int> > model_pairs_list;
    std::vector< std::vector<int> > model_triangles;
    std::vector< std::vector<int> > model_strong_triangles;
    std::vector< std::vector<int> > model_weak_triangles;
    virtual void operator()
                       (std::vector<int> const &config, 
                        std::vector< std::vector<int> > &touched_sites_list, 
                        std::vector< std::vector<int> > &vals_on_touched_list, 
                        std::vector< complex<double> > &hints_list) {};
    
    virtual void init(){};
    
    virtual Ham* clone() const=0;     
            
};

#endif
