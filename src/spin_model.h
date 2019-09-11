#ifndef SPIN_MODEL_HEADER
#define SPIN_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 SPIN MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"global.h"
#include"hamiltonian.h"
using namespace std;

class Spin_Model: public Ham
{
    public:
    //double j_ring,j_x,j_z,jbq_x,jbq_z,h,d,q_x,q_y;
    double j_strong, j_weak, jbq_strong, jbq_weak;
    //double spin;
    int max_ham_spin_change;
    std::vector< complex<double> > omega,omegadag; 
    public:
    void init(double spin, double j_strong, double jbq_strong, double j_weak, double jbq_weak, 
	      std::vector< std::vector<int> > model_strong_triangles,
	      std::vector< std::vector<int> > model_weak_triangles)
    {
        int max=0;
	this->spin=spin;
	this->j_strong=j_strong;
	this->j_weak=j_weak;
	this->jbq_strong=jbq_strong;
	this->jbq_weak=jbq_weak;
        this->model_strong_triangles=model_strong_triangles;
        this->model_weak_triangles=model_weak_triangles;
 
        for (int i=0;i<model_strong_triangles.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                if (model_strong_triangles[i][j]>max) {max=model_strong_triangles[i][j];}
            }
        }
        for (int i=0;i<model_weak_triangles.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                if (model_weak_triangles[i][j]>max) {max=model_weak_triangles[i][j];}
            }
        }
        this->num_sites=max+1;
	this->max_ham_spin_change=2;
   	cout<<"Jstrong   ="<<this->j_strong<<endl; 
   	cout<<"Jweak     ="<<this->j_weak<<endl; 
   	cout<<"Jbq_strong="<<this->jbq_strong<<endl; 
   	cout<<"Jbq_weak  ="<<this->jbq_weak<<endl; 
   	cout<<"Nsites    ="<<this->num_sites<<endl; 
   }
    /*void init(double spin, double j_x, double j_z, 
	      double jbq_x, double jbq_z, double h, double d, 
	      std::vector< std::vector<int> > list_of_pairs,
	      std::vector< std::vector<int> > model_pairs)
    {
	this->omega.resize(3);
	this->omegadag.resize(3);
        int max=0;
	this->q_x=0.0;
        this->q_y=0.0;
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
        }
        this->num_sites=max+1;
	if (j_z!=0.0) this->max_ham_spin_change=0;
	if (j_x!=0.0) this->max_ham_spin_change=1;
	if (jbq_x!=0.0) this->max_ham_spin_change=2;
   	cout<<"Jx="<<j_x<<endl; 
   	cout<<"Jz="<<j_z<<endl; 
   	cout<<"Jbq="<<jbq_x<<endl; 
   	cout<<"Qx="<<q_x<<endl; 
   	cout<<"Qy="<<q_y<<endl; 
   	cout<<"Nsites="<<this->num_sites<<endl; 
   	cout<<"Max /ham"<<this->max_ham_spin_change<<endl;
	this->omega[0]=complex<double>(-0.5,-0.86602540378);
	this->omega[1]=complex<double>(1.0,0.00);
	this->omega[2]=complex<double>(-0.5,0.86602540378); 
	this->omegadag[0]=complex<double>(-0.5,0.86602540378);
	this->omegadag[1]=complex<double>(1.0,0.00);
	this->omegadag[2]=complex<double>(-0.5,-0.86602540378); 
   }
    
   void init(double spin, double j_x, double j_z, 
	      double jbq_x, double jbq_z, double h, double d, 
	      double q_x, double q_y,
	      std::vector< std::vector<int> > list_of_pairs,
	      std::vector< std::vector<int> > model_pairs)
    {
	this->omega.resize(3);
	this->omegadag.resize(3);
        int max=0;
	this->q_x=q_x;
        this->q_y=q_y;
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
        }
        this->num_sites=max+1;
	if (j_z!=0.0) this->max_ham_spin_change=0;
	if (j_x!=0.0) this->max_ham_spin_change=1;
	if (jbq_x!=0.0) this->max_ham_spin_change=2;
   	cout<<"Jx="<<j_x<<endl; 
   	cout<<"Jz="<<j_z<<endl; 
   	cout<<"Jbq="<<jbq_x<<endl; 
   	cout<<"Qx="<<q_x<<endl; 
   	cout<<"Qy="<<q_y<<endl; 
   	cout<<"Nsites="<<this->num_sites<<endl; 
   	cout<<"Max /ham"<<this->max_ham_spin_change<<endl; 
	this->omega[0]=complex<double>(-0.5,-0.86602540378);
	this->omega[1]=complex<double>(1.0,0.00);
	this->omega[2]=complex<double>(-0.5,0.86602540378); 
	this->omegadag[0]=complex<double>(-0.5,0.86602540378);
	this->omegadag[1]=complex<double>(1.0,0.00);
	this->omegadag[2]=complex<double>(-0.5,-0.86602540378); 
   }
   
   void init(double spin, double j_ring, double j_x, double j_z, 
	      double jbq_x, double jbq_z, double h, double d, 
	      double q_x, double q_y,
	      std::vector< std::vector<int> > list_of_pairs,
	      std::vector< std::vector<int> > model_pairs,
	      std::vector< std::vector<int> > model_triangles)
    {
	this->omega.resize(3);
	this->omegadag.resize(3);
        int max=0;
	this->j_ring=j_ring;
	this->q_x=q_x;
        this->q_y=q_y;
        this->spin=spin;
        this->j_x=j_x;this->j_z=j_z;this->h=h;
	this->d=d;this->jbq_x=jbq_x;this->jbq_z=jbq_z;
        this->pairs_list=list_of_pairs;
        this->model_pairs_list=model_pairs;
        this->model_triangles=model_triangles;
 
        for (int i=0;i<list_of_pairs.size();i++)
        {
            for (int j=0;j<2;j++)
            {
                if (list_of_pairs[i][j]>max) {max=list_of_pairs[i][j];}
            }
        }
        this->num_sites=max+1;
	if (j_z!=0.0) this->max_ham_spin_change=0;
	if (j_x!=0.0) this->max_ham_spin_change=1;
	if (jbq_x!=0.0 or j_ring!=0.0) this->max_ham_spin_change=2;
   	cout<<"Jx="<<this->j_x<<endl; 
   	cout<<"Jz="<<this->j_z<<endl; 
   	cout<<"Jring="<<this->j_ring<<endl; 
   	cout<<"Jbq="<<jbq_x<<endl; 
   	cout<<"Qx="<<q_x<<endl; 
   	cout<<"Qy="<<q_y<<endl; 
   	cout<<"Nsites="<<this->num_sites<<endl; 
   	cout<<"Max /ham"<<this->max_ham_spin_change<<endl; 
	this->omega[0]=complex<double>(-0.5,-0.86602540378);
	this->omega[1]=complex<double>(1.0,0.00);
	this->omega[2]=complex<double>(-0.5,0.86602540378); 
	this->omegadag[0]=complex<double>(-0.5,0.86602540378);
	this->omegadag[1]=complex<double>(1.0,0.00);
	this->omegadag[2]=complex<double>(-0.5,-0.86602540378); 
   }
   */
   void operator() (std::vector<int> const &config,
                    std::vector< std::vector<int> > &touched_sites_list, 
                    std::vector< std::vector<int> > &vals_on_touched_list,
                    std::vector< complex<double> > &hints_list);
   
    Ham* clone() const
    {return new Spin_Model(*this);}    
};

void spin_model_setup(std::string filename, 
               Spin_Model &model);

void spin_model_six_spin_setup(Spin_Model &model);

#endif
