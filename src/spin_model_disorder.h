#ifndef SPIN_MODEL_DISORDER_HEADER
#define SPIN_MODEL_DISORDER_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 SPIN MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"global.h"
#include"hamiltonian.h"
#include"number_functions.h"
#include"math_utilities.h"
using namespace std;

class Spin_Model_Disorder: public Ham
{
    public:
    double j_x,j_z,jbq_x,jbq_z,h,d;
    std::vector<double> jx_bonds,jz_bonds,jbq_bonds,jbq_x_bonds,jbq_z_bonds;
    //double spin;
    int max_ham_spin_change;
    std::vector< complex<double> > omega,omegadag; 
    public:
    void init(double spin, double j_x, double j_z, 
	      double jbq_x, double jbq_z, double delta, 
	      double h, double d, 
	      std::vector< std::vector<int> > list_of_pairs,
	      std::vector< std::vector<int> > model_pairs)
    {
	this->omega.resize(3);
	this->omegadag.resize(3);
        int max=0;
        this->spin=spin;
        this->j_x=j_x;this->j_z=j_z;this->h=h;
	this->d=d;this->jbq_x=jbq_x;this->jbq_z=jbq_z;
        this->pairs_list=list_of_pairs;
        this->model_pairs_list=model_pairs;
 
        for (int i=0;i<list_of_pairs.size();i++)
        {
            for (int j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
	    double rj=delta*(2.0*uniform_rnd()-1)*0.5;
	    this->jx_bonds.push_back(this->j_x + rj);
	    this->jz_bonds.push_back(this->j_z + rj);
	    this->jbq_bonds.push_back(0.0);
	    cout<<list_of_pairs[i][0]<<","<<list_of_pairs[i][1]<<" -- "<<this->jx_bonds[i]<<","<<this->jz_bonds[i]<<endl;
        }
        this->num_sites=max+1;
	if (j_z!=0.0)   this->max_ham_spin_change=0;
	if (j_x!=0.0)   this->max_ham_spin_change=1;
	if (jbq_x!=0.0) this->max_ham_spin_change=2;
   	cout<<"Jx       = "<<j_x<<endl; 
   	cout<<"Jz       = "<<j_z<<endl; 
   	cout<<"Delta    = "<<delta<<endl; 
   	cout<<"Jbq      = "<<jbq_x<<endl; 
   	cout<<"Nsites   = "<<this->num_sites<<endl; 
   	cout<<"Max /ham = "<<this->max_ham_spin_change<<endl;
	this->omega[0]=complex<double>(-0.5,-0.86602540378);
	this->omega[1]=complex<double>(1.0,0.00);
	this->omega[2]=complex<double>(-0.5,0.86602540378); 
	this->omegadag[0]=complex<double>(-0.5,0.86602540378);
	this->omegadag[1]=complex<double>(1.0,0.00);
	this->omegadag[2]=complex<double>(-0.5,-0.86602540378); 
   }
    
   void operator() (std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list, 
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);
   
    Ham* clone() const
    {return new Spin_Model_Disorder(*this);}    
};

void spin_model_disorder_setup(std::string filename, 
               Spin_Model_Disorder &model);

#endif
