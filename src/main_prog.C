#include"global.h"
#include"tmp_info.h"

// Test
#include"math_utilities.h"
#include"mtrand.h"

// Search
#include"search_for.h"

// Hamiltonian
#include"hamiltonian.h"

// Heisenberg
#include"spin_model.h"
#include"spin_model_disorder.h"

// ED
#include"ed.h"
#include"lanczos_more.h"
#include"qddm.h"

using namespace std;

void do_model(std::string model,std::string read_file,std::string dump_file);

int main(int argc, char *argv[])
{
    time_t start,end;
    double dif;
    bool found;
    //std::vector<int> dets;
    string str_ret;
    string filename;
    int seed;
    MTRand irand;

    if (argc <= 1)
    {
        cout << "Usage: " << argv[0] << " <Filename>" << endl;
        exit(1);
    }
 
    filename=argv[1];
    search_for(string("seed"),filename,str_ret,found);
    if (found) 
    {
		if (found) {seed=str_to_int(str_ret);}
		else {seed=1;}
    }	
    irand.seed(seed);

    //Tsr_Network tn;
    //read_tensor_network("tn.txt",tn);
    //cout<<"Done with TTN read"<<endl;
    //for (int i=0;i<100;i++) cout<<uniform_rnd()<<endl;
    // Read the TN, now evaluate wavefunction    
    std::vector<int64_t> dets;

    // Too big to evaluate exactly in a short time.
    // Do Monte Carlo evaluation instead  
    //constrained_dets(36,18,dets);
    //int64_t zero=0;
    //int64_t sv=0;
    //cout<<"dets. size() ="<<dets.size()<<endl;
    //# pragma omp parallel for reduction (+:sv,zero)
    /*std::vector<int> config;
    int sign_neel=+1;
    for (int j=0;j<36;j++) 
    {
	if (j%2==0) config.push_back(1);
	if (j%2!=0) config.push_back(0);
    }
    double neel=wavefunction_value(config,tn);
    cout<< "Neel ="<<neel<<endl;*/

    /*int n_b=0;
    for (int j=1;j<36;j+=2) n_b+=btest64(dets[0],j);
    if (n_b%2!=0 and neel/abs(neel)>0) sign_neel=-1;
    if (n_b%2!=0 and neel/abs(neel)<0) sign_neel=+1;
    if (n_b%2==0 and neel/abs(neel)<0) sign_neel=-1;
    if (n_b%2==0 and neel/abs(neel)>0) sign_neel=+1;
    int64_t mco=0;
    double norm=0;*/
    // -1 rule check and check which negative
    /*for (int i=0;i<dets.size();i++)
    {
	std::vector<int> config;
	int n_b=0;
	int marshall_sign=+1;
	for (int j=0;j<36;j++) config.push_back(btest64(dets[i],j));
	for (int j=1;j<36;j+=2) n_b+=btest64(dets[i],j);
	int sign=+1;
	double wf=wavefunction_value(config,tn);
	//if (abs(wf)>1.0e-36) cout<<"state, wf ="<<i<<"  "<<wf<<endl;
	cout<<"state, wf ="<<i<<"  "<<wf<<endl;
	norm+=(wf*wf);
	if (abs(wf)/abs(neel)<1.0e-10) {zero=zero+1; cout<<"state, N zeros ="<<i<<"  "<<zero<<endl;}
	if ( wf/abs(wf) < 0 ) sign=-1;
	if (n_b%2!=0) marshall_sign=-1;
	if (sign != marshall_sign*sign_neel) {sv=sv+1; cout<<"state, N sv ="<<i<<"  "<<sv<<endl;}
	if (sign != sign_neel) {mco=mco+1; cout<<"state, N minus ="<<i<<"  "<<mco<<endl;}
    }*/
   
    /* 
    // -1 rule check and check which negative
    int mc=10;
    std::vector<int> old_config=config;
    double oldwf=neel;
    int prev_sign=+1;
    double rnd,avg_sign=0.0;
    //#pragma omp parallel for 
    for (int i=0;i<mc;i++)
    {
	int n_b=0;
	std::vector<int> config=old_config;
    	// Flip 2 spins
	int rnd1=0;
	int rnd2=0;
  	do
	{
		rnd1=int(uniform_rnd()*36.0);
		rnd2=int(uniform_rnd()*36.0);
		//cout<<"rnd1,rnd2 ="<<rnd1<<" "<<rnd2<<endl;
		//cout<<"config[rnd1],config[rnd2] ="<<config[rnd1]<<" "<<config[rnd2]<<endl;
	} while (rnd1==rnd2);
	//cout<<"Out of loop"<<endl;
        int tmp=config[rnd1];
        config[rnd1]=config[rnd2];
        config[rnd2]=tmp;

	//for (int j=1;j<36;j+=2) n_b+=config[j]; // incorrect way 
	for (int j=0;j<36;j++) 
	{
		int x=j%6;
		int y=j/6;
		if (((x+y)%2) !=0 ) n_b=n_b+(1-config[j]);
	}
	int sign=+1;
	double wf=wavefunction_value(config,tn);
        double rnd=uniform_rnd();
	if (abs(wf*wf)/abs(oldwf*oldwf)>rnd) 
	{
		//Accept
		int marshall_sign=+1;
		int wf_sign=+1;
		if (n_b%2!=0) marshall_sign=-1;
		if (wf<0) wf_sign=-1; 
		avg_sign=avg_sign+double(marshall_sign*wf_sign);
		cout <<"Wf value, wf sign , marshall sign = "<<wf<<"  "<<wf_sign<<"  "<<marshall_sign<<endl;
 		old_config=config;
		oldwf=wf;
		prev_sign=marshall_sign*wf_sign;
	}
	else    // Reject
	{
		avg_sign=avg_sign+double(prev_sign);
	}
    }
    avg_sign=avg_sign/double(mc);
    cout<<"Average sign = "<<avg_sign<<endl;    
    //return 0;

    cout<<sizeof(long)<<endl;
    cout<<sizeof(int64_t)<<endl;
  
    std::vector<int> dets2; 
    cout<<n_choose_k(16,8)<<endl;
    constrained_dets(30,15,dets2);
    //print_vec(dets2,true);*/
///////////////////////////////////////////////////////////////////////////////
////                              ENUMERATE CLUSTERS
///////////////////////////////////////////////////////////////////////////////
//    
//    std::vector< std::vector<Bethe_Cluster> > distinct_clusters; 
//    std::string enum_clus,max_size_str;
//    int max_size;
// 
//    search_for(string("enumerate_clusters"),filename,enum_clus,found);
//    if (found) 
//    {
//        if (str_to_bool(enum_clus))
//        {    
//    		search_for(string("max_size"),filename,max_size_str,found);
//		if (found) {max_size=str_to_int(max_size_str);}
//		else {max_size=10;}
//		enumerate_distinct_clusters(max_size,distinct_clusters); 
//        }
//    }	
//
 /////////////////////////////////////////////////////////////////////////////
 //                    HAMILTONIAN READ AND SETUP calls
 /////////////////////////////////////////////////////////////////////////////
    string hamiltonian;
    Ham *ham=NULL;
    Spin_Model model;

    bool ham_found;

    /*cout<<endl;
    cout<<"========================================="<<endl;
    cout<<"I am reading the Hamiltonian information "<<endl;
    cout<<"========================================="<<endl;*/

    search_for(string("hamiltonian"),filename,hamiltonian,ham_found);
    if (ham_found) 
    {
        if (hamiltonian.compare("spin_model")==0)
        {
            spin_model_setup(filename,model);cout<<endl;
            cout<<"Done with Spin model setup"<<endl;
            ham=model.clone();
            cout<<"Done with cloning Spin model"<<endl;
        }
        else
        {
            cout<<endl;cout<<"I could not find the requested hamiltonian"<<endl;
            //return 0;
        }
    }
 
/////////////////////////////////////////////////////////////////////////////
//                              EXACT DIAGONALIZATION
/////////////////////////////////////////////////////////////////////////////
    std::vector<double> eigs,szs,spins;

    if (ham_found)
    {
	    //print_mathematica_pairs((*ham).pairs_list);
    	    //cout<<"FINISHED Printing pairs list"<<endl;
	    
            search_for(string("diagonalize_all"),filename,str_ret,found);
	    if (str_to_bool(str_ret))
	    {		
		cout<<"Ham sites"<<(*ham).num_sites<<endl;
	    	if ((*ham).num_sites<=13)
		{ ed_get_eigs(*ham,eigs);
		  cout<<"ED energies (full)"<<endl;
		  print_vec(eigs,true);
		}
	    }   

	    search_for(string("diagonalize"),filename,str_ret,found);
	    if (str_to_bool(str_ret))
	    {		
		if (hamiltonian.compare("spin_model")==0)
		{
		    cout<<"STARTING Lanczos spin1, conserved quantum numbers"<<endl;
		    if ((*ham).spin==1) {lanczos_requested_sz_spin_1(*ham,100,eigs,0.0,spins,false);}
		    cout<<"ED energies (Sz=0)"<<endl;
		    print_vec_acc(eigs,true);
		}
	    }
     }
		
/////////////////////////////////////////////////////////////////////////////
//                              Other Analyses
/////////////////////////////////////////////////////////////////////////////
    std::vector< std::vector<double> > 	evecs;
    std::vector<double> 		three_cis;
    std::vector<Matrix> 		Si_z_matrices;
    std::vector<Matrix> 		Si_plus_matrices;
      std::string 			other_analyses;	     
//    int 				arrows;
//    std::vector<double>                 a_i_10,c_i_10;
//    std::string                         file_corrs,c_i_10string;
//
      search_for(string("other_analyses"),filename,other_analyses,found);
//
      if (found)
      {
//	    if (other_analyses.compare("singlet_triplet")==0) singlet_triplet(*ham,128,eigs,szs,spins,evecs,three_cis);
//	    else if (other_analyses.compare("two_d")==0) two_dangling_analysis(*ham,256,eigs,szs,spins,evecs,Si_z_matrices,Si_plus_matrices);
//	    else if (other_analyses.compare("qddm")==0) calc_qddm_fermions(*ham,256,eigs);
//	    else if (other_analyses.compare("qddm_spins")==0) calc_qddm_spins(*ham,256,eigs);
//	    else if (other_analyses.compare("qddm_spin1b2")==0) 
//	    {
//		std::string nlowstr;
//		int nlow=1;
//		search_for(string("nlow"),filename,nlowstr,found);
//		if (found) nlow=str_to_int(nlowstr);
//		
//		int sz=0;
//
//		search_for(string("sz"),filename,nlowstr,found);
//		if (found) sz=str_to_int(nlowstr);
//
//	    	std::vector<int> region_qddm;
//    		search_for(string("region_qddm"),filename,str_ret,found);
//    	        if (found){if (str_ret.substr(0,1)==string("[")) region_qddm=convert_string_to_vec(str_ret);}
//	    	calc_qddm_spin1b2(sz,*ham,256,region_qddm,eigs,nlow);
//	    }
	    if (other_analyses.compare("qddm_spin1")==0) 
	    {
		std::string nlowstr;
		int nlow=1;
		search_for(string("nlow"),filename,nlowstr,found);
		if (found) nlow=str_to_int(nlowstr);
		
		int sz=0;

		search_for(string("sz"),filename,nlowstr,found);
		if (found) sz=str_to_int(nlowstr);

	    	std::vector<int> region_qddm;
    		search_for(string("region_qddm"),filename,str_ret,found);
    	        if (found){if (str_ret.substr(0,1)==string("[")) region_qddm=convert_string_to_vec(str_ret);}
	    	calc_qddm_spin1(sz,*ham,256,region_qddm,eigs,nlow);
	   }
//	    else if (other_analyses.compare("four_d")==0) four_dangling_analysis(*ham,256,eigs,szs,spins,evecs,Si_z_matrices,Si_plus_matrices);
//	    else if (other_analyses.compare("sumi_ct")==0) sumiran_ct(*ham,256,eigs,szs,spins,evecs,Si_z_matrices,Si_plus_matrices);
//	    else if (other_analyses.compare("arrow_count")==0) arrow_count((*ham).pairs_list,arrows);
//	    else if (other_analyses.compare("arrow_count_ens")==0) make_and_analyze_balanced_clusters();
//	    //else if (other_analyses.compare("hida")==0) test_hida();
//	    else if (other_analyses.compare("test_ferro")==0) test_henley_algo();
//	    else if (other_analyses.compare("get_ai")==0) 
//	    {
//    		search_for(string("c_i"),filename,file_corrs,found);
//                if (not found)
//		{
//			cout<<"I can not obtain the a_i since c_i is not specified"<<endl;
//			return 0;
//		}
//                else
//                {
//            	        c_i_10=convert_string_to_vec_double(c_i_10string);
//                }
//
//    		search_for(string("file_corrs"),filename,file_corrs,found);
//                if (not found)
//		{
//			cout<<"I can not obtain the a_i since the file of correlations is not specified"<<endl;
//			return 0;
//		}
//
//		get_ai(file_corrs,c_i_10,a_i_10);
//	    }
//           
//	    else cout<<"Analysis requested not found in this code"<<endl;
     }
/////////////////////////////////////////////////////////////////////////////
////                              Numerical RG (Wilson)
///////////////////////////////////////////////////////////////////////////////
//    int max_states;
//    int size;
//    Bethe_Cluster cluster;
//    std::vector< std::vector<int> > o_adj;
//    std::vector< std::vector<int> > in_out;
//    std::vector< std::vector<int> > out_in;
//    std::vector< std::vector< Matrix > > tr;
//    std::vector< std::vector< Matrix > > s_p;
//    std::vector< std::vector< Matrix > > s_m;
//    std::vector< std::vector< Matrix > > s_z;
//    std::vector< std::vector< Matrix > > hams;
//
//    std::vector< std::vector< std::vector< Matrix > > > sp_tr;
//    std::vector< std::vector< std::vector< Matrix > > > sp_s_p;
//    std::vector< std::vector< std::vector< Matrix > > > sp_s_m;
//    std::vector< std::vector< std::vector< Matrix > > > sp_s_z;
//    std::vector< std::vector< std::vector< Matrix > > > sp_hams;
//    std::vector< std::vector< std::vector<int> > >      inverse_subspace_map;
//    std::vector< std::vector< std::vector<int> > >      inverse_sz_map;
//    std::vector< std::vector< std::vector<double> > >   szs_compact_on_block;
//    std::vector< std::vector< std::vector<double> > >   szs_full_on_block;
//    std::vector< std::vector<double> > 		   	eigs_by_spin;
//    std::vector<Block>  				blocks;
//    std::vector<Universe>  				universes;
//    std::vector<Matrix> l_trans,r_trans;
//    std::vector< std::vector< std::vector<Matrix> > >   left_A_matrices;
//    bool						one_d_rg=false;
//    	    
//    search_for(string("one_d_rg"),filename,str_ret,found);
//    if (found) {one_d_rg=str_to_bool(str_ret);}
//	
//    if (ham_found)
//    {
//    	    search_for(string("nrg"),filename,str_ret,found);
//	    if (str_to_bool(str_ret))
//	    {
//    	    	search_for(string("max_states"),filename,str_ret,found);
//		if (found) {max_states=str_to_int(str_ret);}
//		else {max_states=4;}
//		cout<<"max_states = "<<max_states<<endl;
//		cout<<"Making cluster from pairs"<<endl;
//		make_cluster_from_pairs((*ham).pairs_list,cluster);
//		cout<<"FINISHED making cluster from pairs"<<endl;
//		time(&start);
//		//nrg_refined(cluster,hamiltonian,eigs,o_adj,in_out,out_in,tr,s_p,s_m,s_z,hams,max_states);		
//		if (hamiltonian.compare("xxz")==0) nrg_spin(cluster,xxz.j_z,xxz.j_x,eigs_by_spin,spins,o_adj,in_out,out_in,szs_compact_on_block,szs_full_on_block,inverse_subspace_map,inverse_sz_map,sp_tr,sp_s_p,sp_s_m,sp_s_z,sp_hams,max_states,true);
//		if (hamiltonian.compare("spin_model")==0) 
//		{
//		Simulation_Params sp;
//		sp.ipr=true;
//		sp.max_trunc_states=max_states;
//		if (one_d_rg) {generic_nrg_spin(cluster,model, sp, eigs_by_spin, spins, universes,l_trans,r_trans, left_A_matrices);}
//		else          {tree_nrg_spin(cluster,model, sp, eigs_by_spin, spins, universes,l_trans,r_trans);}
//		}	
//		time(&end);
//		dif=difftime(end,start);
//
//		cout<<"==================================================================="<<endl;
//        	cout<<"Total time to do NRG was "<<dif<<" seconds"<<endl;
//        	cout<<"==================================================================="<<endl;
//	  }   
//     }
//
//////////////////////////////////////////////////////////////////////////////
////                           Calculate Gaps
//////////////////////////////////////////////////////////////////////////////
//    int nsamples;	    
//    bool read_clusters;
//
//    search_for(string("calc_gaps"),filename,str_ret,found);
//    if (found and str_to_bool(str_ret))
//    {	
//    	search_for(string("size"),filename,str_ret,found);
//	if (found)
//	{size=str_to_int(str_ret);}
//	else
//	{size=16;}
//    	
//	search_for(string("nsamples"),filename,str_ret,found);
//	if (found)
//	{nsamples=str_to_int(str_ret);}
//	else
//	{nsamples=1000;}
//
//        search_for(string("read_clusters"),filename,str_ret,found);
//	if (found)
//	{read_clusters=str_to_bool(str_ret);}
//	else
//	{read_clusters=false;}
//
//	calc_exact_or_dmrg_gaps_from_file(size,nsamples,filename,read_clusters);
//
//    }
//
///////////////////////////////////////////////////////////////////////////////
////                              DENSITY MATRIX RG (SR White)
///////////////////////////////////////////////////////////////////////////////
//    std::vector<double> 		target_szs;
//    std::vector<int>  			target_numbers;
//    int 				nsweeps;
//
//    if (ham_found)
//    {
//    	    search_for(string("dmrg"),filename,str_ret,found);
//	    if (str_to_bool(str_ret))
//	    {
//    	    	search_for(string("max_states"),filename,str_ret,found);
//		if (found) {max_states=str_to_int(str_ret);}
//		else {max_states=4;}
//		//cout<<"max_states = "<<max_states<<endl;
//    	    	
//		search_for(string("nsweeps"),filename,str_ret,found);
//		if (found) {nsweeps=str_to_int(str_ret);}
//		else {nsweeps=5;}
//		//cout<<"nsweeps = "<<nsweeps<<endl;
//
//    		search_for(string("target_szs"),filename,str_ret,found);
//    		if (found)
//    		{
//          		if (str_ret.substr(0,1)==string("["))
//          		{
//            			target_szs=convert_string_to_vec_double(str_ret);
//          		}
//    		}
//		else
//		{
//			cout<<"Could not find target s_z for DMRG. I will chose only Ground state from NRG warmup"<<endl;
//		}    
//    		search_for(string("target_numbers"),filename,str_ret,found);
//    		if (found)
//    		{
//          		if (str_ret.substr(0,1)==string("[")) target_numbers=convert_string_to_vec(str_ret);
//    		}
//		else
//		{
//			cout<<"Could not find target numbers for DMRG. I will chose only Ground state from NRG warmup"<<endl;
//		}    
//		//cout<<"Making cluster from pairs"<<endl;
//		make_cluster_from_pairs((*ham).pairs_list,cluster);
//		//cout<<"FINISHED making cluster from pairs"<<endl;
//		//cout<<"Entering DMRG calculation"<<endl;	
//		time(&start);
//		//dmrg_spin(cluster,xxz.j_z,xxz.j_x,target_szs,target_numbers,eigs_by_spin,szs_full_on_block,inverse_subspace_map,inverse_sz_map,sp_tr,sp_s_p,sp_s_m,sp_s_z,sp_hams,max_states,nsweeps,true);
//		if (hamiltonian.compare("xxz")==0)
//		{ 
//			dmrg_spin(cluster,xxz.j_z,xxz.j_x,target_szs,target_numbers,eigs_by_spin,szs_full_on_block,inverse_subspace_map,inverse_sz_map,sp_tr,sp_s_p,sp_s_m,sp_s_z,sp_hams,max_states,nsweeps,true,true,true);
//		}
//
//		if (hamiltonian.compare("spin_model")==0) 
//		{
//			Simulation_Params sp;
//			sp.ipr=true;
//			sp.max_trunc_states=max_states;
//			sp.target_szs=target_szs;
//			sp.target_states=target_numbers;
//    			search_for(string("sweep_history"),filename,str_ret,found);
//    			if (found){if (str_ret.substr(0,1)==string("[")) sp.sweep_history=convert_string_to_vec(str_ret);}
//			sp.n_sweeps=sp.sweep_history.size();
//    			search_for(string("region_1"),filename,str_ret,found);
//    			if (found){if (str_ret.substr(0,1)==string("[")) sp.region_1=convert_string_to_vec(str_ret);}
//    			search_for(string("region_2"),filename,str_ret,found);
//    			if (found){if (str_ret.substr(0,1)==string("[")) sp.region_2=convert_string_to_vec(str_ret);}
//
//			if (one_d_rg) {generic_dmrg_spin(cluster,model, sp, eigs_by_spin, spins, universes,l_trans,r_trans);}
//			else {cout<<"Carrying out Tree DMRG"<<endl; tree_dmrg_spin(cluster,model, sp, eigs_by_spin, spins, universes,l_trans,r_trans);}
//	        }	
//		//cout<<"Out of DMRG calculation"<<endl;	
//		//dmrg(cluster,hamiltonian,eigs,max_states,3);		
//	    	time(&end);	
//		dif=difftime(end,start);
//        
//		cout<<"==================================================================="<<endl;
//        	cout<<"Total time to do DMRG was "<<dif<<" seconds"<<endl;
//        	cout<<"==================================================================="<<endl;
//	     	//return 0;
//	 }   
//     }
//    	    
//
///////////////////////////////////////////////////////////////////////////////
////                              SPIN WAVE ANALYSIS
///////////////////////////////////////////////////////////////////////////////
//    std::vector<double> omegas;
//    int generation;
//    bool spin_wave_found;
//
//    if (ham_found)
//    {
//    	    search_for(string("spin_wave"),filename,str_ret,spin_wave_found);
//	    if (spin_wave_found and str_to_bool(str_ret)) 
//	    {spin_wave_neel((*ham).pairs_list,0.5,omegas);}
//    } 	
//	
//    search_for(string("spin_wave_ensemble"),filename,str_ret,found);
//    if (found and str_to_bool(str_ret)) {make_and_analyze_spin_wave_clusters();}
//
//    
//    if (spin_wave_found)
//    {
//    	    search_for(string("undiluted_Bethe"),filename,str_ret,found);
//	    if (found and str_to_bool(str_ret)) 
//	    {
//		search_for(string("gen"),filename,str_ret,found);
//		if (found) {generation=str_to_int(str_ret);}
//		else {generation=8;}
//		analyze_spin_wave_on_bethe(generation);
//	    }
//    }
//
///////////////////////////////////////////////////////////////////////////////
////                              LAPLACIAN
///////////////////////////////////////////////////////////////////////////////
//    bool laplacian_found;
//
//    if (ham_found)
//    {
//    	    search_for(string("laplacian"),filename,str_ret,laplacian_found);
//	    if (laplacian_found and str_to_bool(str_ret)) {laplacian((*ham).pairs_list,omegas);}
//    } 	
//
//    if (laplacian_found)
//    {
//	    search_for(string("undiluted_Bethe"),filename,str_ret,found);
//	    if (found and str_to_bool(str_ret)) 
//	    {
//		search_for(string("gen"),filename,str_ret,found);
//		if (found) {generation=str_to_int(str_ret);}
//		else {generation=8;}
//		analyze_spin_wave_on_bethe(generation,true,false);
//	    }
//    }	
//   	
///////////////////////////////////////////////////////////////////////////////
////                  DO MODEL calculation on many clusters
///////////////////////////////////////////////////////////////////////////////
//     	    
//    //do_model(string("tfim"),string("study_fm_clusters.txt"),string("./data_tfim_8/data_tfim"));
//    //do_model(string("xxz"),string("study_fm_clusters.txt"),string("./data_xxz/data_xxz_8/data_xxz"));
//
///////////////////////////////////////////////////////////////////////////////
////                          CONVERSION to g6 format
///////////////////////////////////////////////////////////////////////////////
//    std::string txt_file_clusters;
//    search_for(string("g6_convert"),filename,str_ret,found);
//    { 
//	if (found)
//	{
//    		if (str_to_bool(str_ret))
//		{
//			search_for(string("clusters_filename"),filename,txt_file_clusters,found);		
//    			if (found){convert_all_to_g6(txt_file_clusters);}
//		}
//	}
//     }	
//
///////////////////////////////////////////////////////////////////////////////
////                              CLEAN UP 
///////////////////////////////////////////////////////////////////////////////
//
    if (ham!=NULL) {delete ham;ham=NULL;}
    return 0;

}

////////////////////////////////////////////////////////////////////////////////
//
//void do_model(std::string model,
//	      std::string read_file,
//              std::string eigs_dump_file)
//{
//    Spin_Half_TFIM tfim;
//    Spin_Half_XXZ xxz;
//
//    Ham *ham=NULL;
//    bool found;
//    int i,j,ctr;
//    std::vector< std::vector<double> > eigs_all;
//    std::vector<double> hs,eigs,gaps,spins,szs;
//    std::string str,str_ret,tmp_string,new_file,gap_file;
//    std::vector< std::vector<int> > pairs;
//    bool measure_s;
//
//    ctr=0;found=true;
//    
//    hs.push_back(0.05);hs.push_back(0.10);hs.push_back(0.15);hs.push_back(0.20);hs.push_back(0.25);
//    hs.push_back(0.30);hs.push_back(0.50);hs.push_back(0.70);hs.push_back(1.00);
//
//    while (found)
//    {
//        str="pairs_cluster_";
//        str+=to_string(ctr);
//        search_for(str,read_file,str_ret,found);
//        
//        if (found)
//        {
//              if (str_ret.substr(0,1)==string("[")){pairs=convert_string_to_vec_of_vec(str_ret);}
//              eigs_all.clear();
//              gaps.clear();
//
//              for (i=0;i<hs.size();i++)
//              {
//                    if (model.compare("tfim")==0)
//                    {   tfim.init(-1.0,hs[i],pairs);
//                        ham=tfim.clone();
//                    }
//                    else if (model.compare("xxz")==0)
//                    {   
//			measure_s=false;
//			xxz.init(hs[i],1.0,pairs);
//                        ham=xxz.clone();
//			if (hs[i]==1.0) {measure_s=true;}
//                    }
//		    
//                    if ((*ham).num_sites<=10){ed_get_eigs(*ham,eigs);}
//                    else if (model.compare("tfim")==0){lanczos_no_sym_get_eigs(*ham,512,eigs);}
//		    else if (model.compare("xxz")==0) {lanczos_spin_sym(*ham,512,eigs,szs,spins,measure_s);}	
//                    eigs_all.push_back(eigs);
//                    gaps.push_back(eigs[1]-eigs[0]);
//                    if (ham!=NULL) {delete ham;ham=NULL;}
//              }
//              
//              new_file=eigs_dump_file+string(".txt.")+to_string(ctr);
//              dump_multiple_functions(hs,eigs_all,new_file);
//
//              new_file=eigs_dump_file+string(".gaps.txt.")+to_string(ctr);
//              dump_function(hs,gaps,new_file);
//              
//              ctr+=1;
//        }
//    }
//    if (ham!=NULL) {delete ham;ham=NULL;}
//}
//
////////////////////////////////////////////////////////////////////////////////
