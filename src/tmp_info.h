#ifndef TMP_INFO_HEADER
#define TMP_INFO_HEADER

#include"global.h"
#include"matrix.h"

using namespace std;

class TI
{
	public:
	Matrix								couplings;
	std::vector< std::vector<int> > 				list_of_converted_vecs;
	std::vector<int> 		      				connected_site_locations;
	std::vector< std::vector<int> >      				connected_site_systems;
	std::vector< std::vector< std::vector<int> > >     		connected_system_system;
	std::vector< std::vector<int> >      				connected_site_sys_env;
	std::vector< std::vector<int> >      				system_pairs;
        std::vector<int>      						block_num_states;
	std::vector<double>   						tmp_eigs;
	std::vector< std::vector<double> >   				eigs;
	std::vector< std::vector<int> > 				map_for_states;
	std::vector< std::vector< std::vector<int> > > 			map_for_hints;
	std::vector< std::vector< std::vector<double> > >		hints;
	int								max_spin_index;
	int 								max_ham_spin_change;
	int 								max_trunc_states;
	int 								lanc_dav_it;
	int 								num_spaces;
	bool								diag;
	std::vector<Matrix>                                             tmp_opt_evecs;
	bool								subspace_check;
	std::vector<int>						target_states;
	std::vector<int>						list;
	std::vector<int>						inverse_sz_map,inverse_subspace_map;
	std::vector<double>						distinct_szs;
	std::vector<double>						all_szs;
	std::vector<int>						sites;

};

class Simulation_Params
{
	public:
	int 								lanc_dav_it;
	int 								max_trunc_states;
	bool		                                                ipr;
 	std::vector<double>                                             target_szs;
 	std::vector<int>	                                        target_states;
 	std::vector<int>	                                        sweep_history;
	int 								n_sweeps;
	int 								direction;
        std::vector< std::vector<int> > 		   	        contraction_order;
        std::vector< std::vector<int> > 		   	        in_out_sweep_order;
        std::vector< std::vector<int> > 		   	        out_in_sweep_order;
	std::vector<int>						systems_list;
	std::vector<int>						region_1,region_2;
        bool                                                            at_root;
	bool                                                            check_env;
	bool 								get_A_special;
	std::vector< std::vector< std::vector< std::vector<Matrix> > > > A_matrices;
	std::vector< std::vector< std::vector<Matrix> > >                special_A_matrices;
	bool                                                            recordT;
	int 								cluster_root;
	double 								noise;
};

#endif
