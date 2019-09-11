#include"global.h"
#include"spin_model.h"
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"search_for.h"
#include"printing_functions.h"
#include"spin_functions.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 SPIN MODEL
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Spin_Model::operator()
                    (std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list, 
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)
{
    std::vector< std::vector<int> >:: iterator p;
    double spin=this->spin;
    double j_strong=this->j_strong; 
    double jbq_strong=this->jbq_strong; 
    //double q_x=this->q_x;
    double j_weak=this->j_weak; 
    double jbq_weak=this->jbq_weak; 
    complex<double> diaghint=0.0;

    //cout<<"Spin = "<<spin<<endl;
    //cout<<"J strong = "<<j_strong<<endl; 
    //cout<<"J weak   = "<<j_weak<<endl; 
    for (int p=0;p<this->model_strong_triangles.size();p++)
    {
	for (int n=0;n<3;n++)
	{
		int first=this->model_strong_triangles[p][(0+n)%3];int second=this->model_strong_triangles[p][(1+n)%3];
		//cout<<"Strong first,second="<<first<<","<<second<<endl;
		calc_hints_any_spin_all_sip_sjm(spin,j_strong,jbq_strong,
				     first,second,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sip_sjm(spin,j_strong,jbq_strong,
				     second,first,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sipp_sjmm(spin,0.25*jbq_strong,
				     first,second,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sipp_sjmm(spin,0.25*jbq_strong,
				     second,first,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		int a=config[first];
        	int b=config[second];
		double sz1=double(a)-spin;
		double sz2=double(b)-spin;
		diaghint+=(j_strong*sz1*sz2)+(jbq_strong*sz1*sz1*sz2*sz2);
		diaghint+=(0.5*jbq_strong*sminus_splus_fn(spin,a)*sminus_splus_fn(spin,b));
		diaghint+=(0.5*jbq_strong*sminus_splus_fn(spin,a)*sz2);
		diaghint+=(0.5*jbq_strong*sminus_splus_fn(spin,b)*sz1);
	}
    }
    
    for (int p=0;p<this->model_weak_triangles.size();p++)
    {
	for (int n=0;n<3;n++)
	{
		int first=this->model_weak_triangles[p][(0+n)%3];int second=this->model_weak_triangles[p][(1+n)%3];
		//cout<<"Weak first,second="<<first<<","<<second<<endl;
		calc_hints_any_spin_all_sip_sjm(spin,j_weak,jbq_weak,
				     first,second,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sip_sjm(spin,j_weak,jbq_weak,
				     second,first,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sipp_sjmm(spin,0.25*jbq_weak,
				     first,second,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		calc_hints_any_spin_all_sipp_sjmm(spin,0.25*jbq_weak,
				     second,first,config,
				     touched_sites_list,vals_on_touched_list,
				     hints_list);
		int a=config[first];
        	int b=config[second];
		double sz1=double(a)-spin;
		double sz2=double(b)-spin;
		diaghint+=(j_weak*sz1*sz2)+(jbq_weak*sz1*sz1*sz2*sz2);
		diaghint+=(0.5*jbq_weak*sminus_splus_fn(spin,a)*sminus_splus_fn(spin,b));
		diaghint+=(0.5*jbq_weak*sminus_splus_fn(spin,a)*sz2);
		diaghint+=(0.5*jbq_weak*sminus_splus_fn(spin,b)*sz1);
	}
    }

   
    /*if (abs(this->j_ring)>1.0e-8)
    { 
	    for (p=this->model_triangles.begin();p!=this->model_triangles.end();p++)
	    {
		int first=(*p)[0];int second=(*p)[1];int third=(*p)[2];
		//cout<<"first = "<<first<<" second = "<<second<<" third = "<<third<<endl;

		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,first,second,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		
		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,first,third,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		
		
		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,second,third,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		
		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,second,first,third,config,touched_sites_list,vals_on_touched_list,hints_list);
		
		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,third,first,second,config,touched_sites_list,vals_on_touched_list,hints_list);
		 
		calc_hints_any_spin_all_s1pp_s2m_s3m(spin,0.5*0.25*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pm_s2m_s3p(spin,0.5*0.25*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1pz_s2m_s3z(spin,0.5*0.50*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mp_s2p_s3m(spin,0.5*0.25*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mm_s2p_s3p(spin,0.5*0.25*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1mz_s2p_s3z(spin,0.5*0.50*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zp_s2z_s3m(spin,0.5*0.50*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
		calc_hints_any_spin_all_s1zm_s2z_s3p(spin,0.5*0.50*this->j_ring,third,second,first,config,touched_sites_list,vals_on_touched_list,hints_list);
	   }
   }*/
    // Diagonal terms
    //complex<double> hint=0.0;
    /*for (int i=0;i<this->num_sites;i++)
    {
	double sz=config[i]-spin;
	diaghint+=(this->h*sz)+(this->d*sz*sz);
    }*/
    
    std::vector<int> touched_sites, vals_on_touched;
    hints_list.push_back(diaghint);
    touched_sites_list.push_back(touched_sites);
    vals_on_touched_list.push_back(vals_on_touched);

    /*
    //Diagonal terms continued
    for (p=this->model_pairs_list.begin();p!=this->model_pairs_list.end();p++)
    {
	first=(*p)[0];second=(*p)[1];
	int a=config[first];
        int b=config[second];
	double sz1=double(a)-spin;
	double sz2=double(b)-spin;
	//cout<<"first = "<<first<<endl;
	//cout<<"second = "<<second<<endl;
        //cout<<"sz1="<<sz1<<endl;
	//cout<<"sz2="<<sz2<<endl;
	hint+=(this->j_z*sz1*sz2)+(j_bq*sz1*sz1*sz2*sz2);
	hint+=(0.5*j_bq*sminus_splus_fn(spin,a)*sminus_splus_fn(spin,b));
	hint+=(0.5*j_bq*sminus_splus_fn(spin,a)*sz2);
	hint+=(0.5*j_bq*sminus_splus_fn(spin,b)*sz1);
   	if (abs(spin-1.0)<1.0e-6) 
	{
	 //cout<<endl;
	 //cout<<"a,b = "<<a<<" "<<b<<endl;
	 //cout<<"sz1,sz2 = "<<sz1<<" "<<sz2<<endl;
	 //complex<double> thint2=complex<double>(0.0,0.0);
	 //thint2=(complex<double>(q_x,0.0)*this->omega[a]*this->omegadag[b]);
	 //thint2=thint2+(complex<double>(q_x,0.0)*this->omegadag[a]*this->omega[b]);
	 double thint=0.0;
	 thint=thint+((+1.5*q_x)*sz1*sz2)+( (4.5*q_x)*sz1*sz1*sz2*sz2);
	 thint=thint+((-3.0*q_x)*sz1*sz1)+( (-3.0*q_x)*sz2*sz2);
	 thint=thint+(2.0*q_x);
	 //cout<<"thint = "<<thint<<endl;
	 //cout<<"thint2 = "<<thint2<<endl;
	 //hint=hint+(this->q_x*this->omegadag[a]*this->omega[b]);
	 hint=hint+thint;
	}
    
    }

    //Diagonal terms continued AGAIN for the triangular piece
    for (p=this->model_triangles.begin();p!=this->model_triangles.end();p++)
    {
	first=(*p)[0];second=(*p)[1];third=(*p)[2];
	int a=config[first];
        int b=config[second];
        int c=config[third];
	double sz1=double(a)-spin;
	double sz2=double(b)-spin;
	double sz3=double(c)-spin;
	hint+=(1.0*this->j_ring*sz1*sz2*sz1*sz3);
	hint+=(1.0*this->j_ring*sz2*sz3*sz2*sz1);
	hint+=(1.0*this->j_ring*sz3*sz1*sz3*sz2);
    }
    std::vector<int> touched_sites, vals_on_touched;
    //cout<<"hint="<<hint<<endl;
    hints_list.push_back(hint);
    touched_sites_list.push_back(touched_sites);
    vals_on_touched_list.push_back(vals_on_touched);


    if (abs(spin-1.0)<1.0e-6 and abs(this->q_y)>1.0e-6)
    {
	// Add the Qy terms
    	for (p=this->model_pairs_list.begin();p!=this->model_pairs_list.end();p++)
    	{
		touched_sites.clear();
		vals_on_touched.clear();
		first=(*p)[0];second=(*p)[1];
		int a=config[first];
        	int b=config[second];
		int new_a=(a-1+3)%3;
		int new_b=(b+1)%3;
		touched_sites.push_back(first);
		touched_sites.push_back(second);
		vals_on_touched.push_back(new_a);
		vals_on_touched.push_back(new_b);
    		touched_sites_list.push_back(touched_sites);
    		vals_on_touched_list.push_back(vals_on_touched);
		hints_list.push_back(complex<double>(this->q_y));
		
		touched_sites.clear();
		vals_on_touched.clear();
                new_a=(a+1)%3;
		new_b=(b-1+3)%3;
		touched_sites.push_back(first);
		touched_sites.push_back(second);
		vals_on_touched.push_back(new_a);
		vals_on_touched.push_back(new_b);
    		touched_sites_list.push_back(touched_sites);
    		vals_on_touched_list.push_back(vals_on_touched);
		hints_list.push_back(complex<double>(this->q_y));
	}
    }
    */
}

////////////////////////////////////////////////////////////////////////////////
void spin_model_setup(string filename, 
               Spin_Model &model)

{    
    // Set the Spin Hamiltonian
    bool found=true;
    string str,str_ret;
    //double spin,j_ring,j_x,j_z,jbq_x,jbq_z,h,d,q_x,q_y;
    double spin,jbq_strong,jbq_weak,j_strong,j_weak;
    int L;
    std::vector<int> pair;
    //std::vector< std::vector<int> > pairs,model_pairs,model_triangles;
    std::vector< std::vector<int> > model_strong_triangles,model_weak_triangles;

    // search for coupling parameters
    search_for("spin",filename,str_ret,found);
    if (found){spin=str_to_d(str_ret);} else {spin=0.5;}
    
    //search_for("j_ring",filename,str_ret,found);
    //if (found){j_ring=str_to_d(str_ret);} else {j_ring=0.0;}
    
    search_for("j_strong",filename,str_ret,found);
    if (found){j_strong=str_to_d(str_ret);} else {j_strong=0.0;}

    search_for("jbq_strong",filename,str_ret,found);
    if (found){jbq_strong=str_to_d(str_ret);} else {jbq_strong=0.0;}
    
    search_for("j_weak",filename,str_ret,found);
    if (found){j_weak=str_to_d(str_ret);} else {j_weak=0.0;}

    search_for("jbq_weak",filename,str_ret,found);
    if (found){jbq_weak=str_to_d(str_ret);} else {jbq_weak=0.0;}
    
    /*search_for("h",filename,str_ret,found);
    if (found){h=str_to_d(str_ret);} else {h=0.0;}
    
    search_for("d",filename,str_ret,found);
    if (found){d=str_to_d(str_ret);} else {d=0.0;}
   
    bool found1,found2;
    search_for("q_x",filename,str_ret,found1);
    if (found1){q_x=str_to_d(str_ret);} else {q_x=0.0;}
    
    search_for("q_y",filename,str_ret,found2);
    if (found2){q_y=str_to_d(str_ret);} else {q_y=0.0;}
*/
    /*if (found1==false and found2==false)
    {
	double magq,phi; 
    	search_for("phi",filename,str_ret,found);
    	if (found){phi=str_to_d(str_ret);} else {phi=0.0;}
    
    	search_for("magq",filename,str_ret,found);
    	if (found){magq=str_to_d(str_ret);} else {magq=0.0;}
    
	q_x=magq*cos(phi);
	q_y=magq*sin(phi);
    }

    search_for("oned",filename,str_ret,found);
    if (str_to_bool(str_ret) and found)
    {
    	search_for("L",filename,str_ret,found);
    	if (found){L=str_to_int(str_ret);} else {L=4;}
	for (int i=0;i<L-1;i++) {
	pair.push_back(i); pair.push_back(i+1);
	pairs.push_back(pair);
	pair.clear();}
	model_pairs=pairs;
	bool periodic;
	search_for("periodic",filename,str_ret,found);
    	if (found){periodic=str_to_bool(str_ret);} else {periodic=false;}
	if (periodic and L>2)
	{
		pair.clear();
		pair.push_back(0);
		pair.push_back(L-1);
		model_pairs.push_back(pair);
	}
	cout<<"MODEL PAIRS"<<endl;
	print_mat_int(model_pairs);
    }    
    else
    {    
    	search_for("pairs",filename,str_ret,found);
    	if (found)
    	{if (str_ret.substr(0,1)==string("[")) pairs=convert_string_to_vec_of_vec(str_ret);}  
    	
	search_for("model_pairs",filename,str_ret,found);
    	if (found)
    	{if (str_ret.substr(0,1)==string("[")) model_pairs=convert_string_to_vec_of_vec(str_ret);} 
	else {model_pairs=pairs;} 
	cout<<"MODEL PAIRS"<<endl;
	print_mat_int(model_pairs); 
*/
	search_for("model_strong_triangles",filename,str_ret,found);
    	if (found)
    	{if (str_ret.substr(0,1)==string("[")) model_strong_triangles=convert_string_to_vec_of_vec(str_ret);}
	cout<<"MODEL STRONG TRIANGLES"<<endl;
	print_mat_int(model_strong_triangles); 
	
	search_for("model_weak_triangles",filename,str_ret,found);
    	if (found)
    	{if (str_ret.substr(0,1)==string("[")) model_weak_triangles=convert_string_to_vec_of_vec(str_ret);}
	cout<<"MODEL WEAK TRIANGLES"<<endl;
	print_mat_int(model_weak_triangles); 
//    }  
    //if (abs(q_x)<1.0e-6 and abs(q_y)<1.0e-6) {model.init(spin,j_x,j_z,jbq_x,jbq_z,h,d,pairs,model_pairs);}
    //else                                     {model.init(spin,j_x,j_z,jbq_x,jbq_z,h,d,q_x,q_y,pairs,model_pairs);}
    //if (abs(q_x)<1.0e-6 and abs(q_y)<1.0e-6) {model.init(spin,j_x,j_z,jbq_x,jbq_z,h,d,pairs,model_pairs);}
    model.init(spin,j_strong, jbq_strong, j_weak, jbq_weak, model_strong_triangles,model_weak_triangles);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void spin_model_six_spin_setup(Spin_Model &model)
{    
    // Set the Spin Hamiltonian
    bool found=true;
    string str,str_ret;
    double spin,j_x,j_z,jbq_x,jbq_z,h,d,q_x,q_y;
    int L;
    std::vector<int> pair;
    std::vector< std::vector<int> > pairs,model_pairs;

    spin=0.5;
    j_x=1.0;
    j_z=1.0;
    jbq_x=0.0;
    jbq_z=0.0;
    h=0.0;
    d=0.0;
    q_x=0.0;
    q_y=0.0;
    L=6;
    bool periodic=true;
    for (int i=0;i<L-1;i++) {
    pair.push_back(i); pair.push_back(i+1);
    pairs.push_back(pair);
    pair.clear();}
    model_pairs=pairs;
    if (periodic and L>2)
    {
		pair.clear();
		pair.push_back(0);
		pair.push_back(L-1);
		model_pairs.push_back(pair);
    }
    cout<<"MODEL PAIRS (six spin)"<<endl;
    print_mat_int(model_pairs);
    if (abs(q_x)<1.0e-6 and abs(q_y)<1.0e-6) {model.init(spin,j_x,j_z,jbq_x,jbq_z,h,d,pairs,model_pairs);}
}*/

