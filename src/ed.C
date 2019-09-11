#include"ed.h"
#include"printing_functions.h"
#include<omp.h>

//##
void equatek(vector< complex<double> > &x, vector< complex<double> > &y)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < x.size(); ++i) {y[i]=x[i];}
}

void zscalk(const int64_t size,
	    complex<double> a,
	    vector< complex< double> > &x)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {x[i]*=(a);}

}


void zaxpyk(const int64_t size,
	    complex<double> a,
	    vector< complex< double> > &x,
	    vector< complex< double> > &y)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {y[i]+=(a*x[i]);}

}

complex<double> zdotc(	const int64_t &size,
			const vector< complex< double> > &v1,
			const vector< complex< double> > &v2)
{
	double sumr = 0.0;
	double sumi = 0.0;

	#pragma omp parallel for default(shared) reduction (+ : sumr,sumi)
	for(int64_t i = 0; i < size; ++i)
	{
		complex<double> a=conj(v1[i])*v2[i];
		sumr += real(a);
		sumi += imag(a);
	}

	complex<double> sum=complex<double>(sumr,sumi);

	return sum;

}

void normalize(std::vector< complex<double> > & v) 
{
	int64_t size=v.size();
	complex<double> norminv=1.0/sqrt(real(zdotc(size, v, v)));
	zscalk(size,norminv,v);
}

void orth_wrt_previous_evecs(std::vector< complex<double> > & v, 
		       std::vector< std::vector< complex<double> > > & previous_evecs)
{	
	int64_t size=v.size();
	std::vector< complex<double> > qs;

	for (int i=0; i < previous_evecs.size();i++)
	{
		complex<double> q=conj(zdotc(size, v, previous_evecs[i]));
		qs.push_back(q);
	}
	
	for (int i=0; i < previous_evecs.size();i++)
	{
		zaxpyk(size,-qs[i],previous_evecs[i],v);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void equatek(vector<double> &x, vector<double> &y)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < x.size(); ++i) {y[i]=x[i];}
}

//////////////////////////////////////////////////////////////////////
void dscalk(const int64_t size,
	    double a,
	    vector< double > &x)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {x[i]*=(a);}
}

///////////////////////////////////////////////////////////////////////////////////////////
void daxpyk(const int64_t size,
	    double a,
	    vector< double > &x,
	    vector< double > &y)
{
	#pragma omp parallel for
	for(int64_t i = 0; i < size; ++i) {y[i]+=(a*x[i]);}
}


///////////////////////////////////////////////////////////////////////////////////////////
double       ddotk(	const int64_t &size,
			vector< double > &v1,
			vector< double > &v2)
{
	double sumr = 0.0;
	#pragma omp parallel for default(shared) reduction (+ : sumr)
	for(int64_t i = 0; i < size; ++i)
	{
		double a=v1[i]*v2[i];
		sumr += a;
	}
	return sumr;
}

///////////////////////////////////////////////////////////////////////////////////////////
void orth_wrt_previous_evecs(std::vector< double > &v, 
		       std::vector< std::vector< double > > &previous_evecs)
{	
	int64_t size=v.size();
	std::vector< double > qs;

	for (int i=0; i < previous_evecs.size();i++)
	{
		double q=ddotk(size, v, previous_evecs[i]);
		qs.push_back(q);
	}
	
	for (int i=0; i < previous_evecs.size();i++) daxpyk(size,-qs[i],previous_evecs[i],v);
}

///////////////////////////////////////////////////////////////////////////////////////////
void normalize(std::vector< double > & v) 
{
	int64_t size=v.size();
	double norminv=1.0/sqrt((ddotk(size, v, v)));
	dscalk(size,norminv,v);
}


//////////////////////////////////////////////////////////////////////////////
void ed_with_hints_given(std::vector< std::vector<int> > 		const &map,
		      	 std::vector< std::vector<complex<double> > > 	const &hints,
                         std::vector<double> 			        &eigs,
			 Matrix 					&eigenvecs,
			 bool 					        ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					size=map.size();
   Matrix                               ham(size,size);
   double 				dif;
 
   //if (ipr) {cout<<" Number of states in Matrix diag is "<<size<<endl;}
   for (int i=0;i<size*size;i++) ham[i]=0.0;
   time(&start_0);
   eigs.clear();
   eigenvecs.clear();
   eigs.resize(size);
   eigenvecs.resize(size,size);    
   time (&start_0);
   for (int i=0;i<size;i++) 
   {
   	for (int k=0;k<map[i].size();k++){ham(i,map[i][k])+=real(hints[i][k]);}
   }
   if (size<20) print_real_mat(ham);   
   symmetric_diagonalize(ham,eigs,eigenvecs);
   time(&end_0);
   dif=difftime(end_0,start_0);
}
//////////////////////////////////////////////////////////////////////////////
void ed_with_hints_given(std::vector< std::vector<int> > 		const &map,
		      	 std::vector< std::vector<double> > 	        const &hints,
                         std::vector<double> 			        &eigs,
			 Matrix 					&eigenvecs,
			 bool 					        ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					size=map.size();
   Matrix                               ham(size,size);
   double 				dif;
 
   //if (ipr) {cout<<" Number of states in Matrix diag is "<<size<<endl;}
   for (int i=0;i<size*size;i++) ham[i]=0.0;
   time(&start_0);
   eigs.clear();
   eigenvecs.clear();
   eigs.resize(size);
   eigenvecs.resize(size,size);    
   time (&start_0);
   for (int i=0;i<size;i++) 
   {
   	for (int k=0;k<map[i].size();k++){ham(i,map[i][k])+=(hints[i][k]);}
   }
   if (size<20) print_real_mat(ham);   
   symmetric_diagonalize(ham,eigs,eigenvecs);
   time(&end_0);
   dif=difftime(end_0,start_0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double measure_sz(std::vector<int> const &config)
{ 
  int 		n_up,n_down;
  n_up=(int) count(config.begin(),config.end(),1);
  n_down=config.size()-n_up;
  return double(n_up-n_down)/2.0;
}

/////////////////////////////////////////////////////////////////////////////////////
double compute_spin(double s_s_plus_one)
{return ((sqrt((4.0*s_s_plus_one)+1.0)-1.0)/2.0);}

/////////////////////////////////////////////////////////////////////////////////////
void ed_special(Ham &h, double sz, std::vector<int> &dets, std::vector<double> &gs)
{
   int nsites=h.num_sites;
   int num_ones=int ((nsites/2)+sz+ 1.0e-06);
   constrained_dets(h.num_sites,num_ones,dets);
   int 					size=dets.size();
   std::vector<int> 			config(h.num_sites);
   std::vector< complex<double> > 	hints_list;
   std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list;
   Matrix 				h_mat(size,size);
   Matrix 				eigenvecs(size,size);
   std::vector<double>  eigs(size);
   
   for (int i1=0;i1<size;i1++)
   {
	int i=dets[i1];
        convert_num_to_vec(i,2,h.num_sites,config);
	touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        for (int k=0;k<touched_sites_list.size();k++)
        {
            int j=i;
            for (int l=0;l<touched_sites_list[k].size();l++)
            {j=j-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
            int j1=0;
	    for (int l=0;l<dets.size();l++)
	    {
		if(dets[l]==j) {j1=l;}
	    } 

            if (j>=i){h_mat(i1,j1)=h_mat(i1,j1)+real(hints_list[k]);h_mat(j1,i1)=h_mat(i1,j1);}
        }
   }
   symmetric_diagonalize(h_mat,eigs,eigenvecs);
   cout<<"Lowest eigen-energy = "<<eigs[0]<<endl;
   gs.resize(size);
   for (int i1=0;i1<size;i1++) gs[i1]=eigenvecs(i1,0);
}
//////////////////////////////////////////////////////////////////////////////

void ed_get_eigs(Ham &h, std::vector<double> &eigs)
{
   
   time_t 				start,end;
   int 					size=pow(2,h.num_sites);
   int 					i,j,k,l,m;
   int 					spin;
   double 				dif;
   double 				tmp;
   complex<double> 			tot,c_first_second;
   std::vector<int> 			config(h.num_sites);
   std::vector< complex<double> > 	hints_list;
   std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list;
   std::string 				tmp_string;
   std::string 				eigs_dump_filename;
   std::string 				evecs_dump_filename;
   Matrix 				h_mat(size,size);
   Matrix 				eigenvecs(size,size);
   
   eigs.resize(size);
   
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (k=0;k<touched_sites_list.size();k++)
        {
            j=i;
            for (l=0;l<touched_sites_list[k].size();l++)
            {j=j-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
            
            if (j>=i){h_mat(i,j)=h_mat(i,j)+real(hints_list[k]);h_mat(j,i)=h_mat(i,j);}
        }
   }
   
   time (&start);
   symmetric_diagonalize(h_mat,eigs,eigenvecs);
   time (&end);
   dif=difftime(end,start);

   cout<<"==================================================================="<<endl;
   cout<<"Total time to diagonalize (Exact) was "<<dif<<" seconds"<<endl;
   cout<<"==================================================================="<<endl;
   
   for (j=0;j<size;j++)
   {
      i=0;  	
      while (abs(eigenvecs(i,j))<1.0e-2) {i+=1;}

      //cout<<"Found non zero element i = "<<i<<endl;
        
      convert_num_to_vec(i,2,h.num_sites,config);
      tot=0.0;

      for (int first=0;first<h.num_sites;first++)
      {
       	  for (int second=first+1;second<h.num_sites;second++)
          {  
		touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
		calc_hints_sxsx_sysy(1.0,first,second,config,
			             touched_sites_list,
                                     vals_on_touched_list,
                                     hints_list);
		
		c_first_second=0.0;

		for (k=0;k<touched_sites_list.size();k++)
		{
			m=i;
			for (l=0;l<touched_sites_list[k].size();l++){m=m-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
			c_first_second+=(hints_list[k]*eigenvecs(m,j));
		}
			
		touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
      	  	calc_hints_szsz(1.0,first,second,config,
		      touched_sites_list,vals_on_touched_list,hints_list);
 		
		c_first_second+=(hints_list[0]*eigenvecs(i,j));

		tot+=c_first_second;	  
	}
      }		
     
      tot=tot/eigenvecs(i,j); 
      tmp=(2.0*real(tot))+(0.75*double(h.num_sites));
      //cout<<"S(S+1) = "<<tmp<<endl;
      cout<<"==========================================="<<endl;
      cout<<"Spin   = "<<compute_spin(tmp)<<endl;
      cout<<"S_z    = "<<measure_sz(config)<<endl;
      cout<<"Energy = "<<eigs[j]<<endl;	
      cout<<"==========================================="<<endl;
  }

}
//////////////////////////////////////////////////////////////////////////////
void lanczos_no_sym_any_spin_get_eigs(Ham &h,
                                       int iterations, 
                                       std::vector<double> &eigs)
{
   int base=int(2.0*h.spin+1.0e-6)+1;
   int nsites=h.num_sites;
   int hilbert=pow(base,nsites);
   cout<< " Hilbert = "<<hilbert<<endl;
   iterations=min(iterations,hilbert);
   std::vector<int> 		      config;
   Matrix 			      eigenvecs;
   std::vector< complex<double> >     hints_list;
   std::vector< double >              rehints_list;
   std::vector <std::vector<double> > hints;
   std::vector<int> 		      connected_els;
   std::vector< std::vector<int> >    map;
   std::vector< std::vector<int> >    touched_sites_list,vals_on_touched_list;
   std::vector<int> mapv,inverse_map;

   for (int i=0;i<hilbert;i++)
   {
	mapv.push_back(i);
	inverse_map.push_back(i); 
	touched_sites_list.clear();vals_on_touched_list.clear();
	hints_list.clear();
	rehints_list.clear();
	connected_els.clear();
        convert_num_to_vec(i,base,nsites,config);
	//print_vec(config);
   	h(config,touched_sites_list,vals_on_touched_list,hints_list);
	//cout<<" H done acting on config"<<endl;
	for (int n=0; n< hints_list.size();n++)
	{
		int j=i;
		for (int a=0;a<touched_sites_list[n].size();a++)
		{
			int po=pow(base,nsites-1-touched_sites_list[n][a]);
			j=j+po*(vals_on_touched_list[n][a]-config[touched_sites_list[n][a]]);
		}
		connected_els.push_back(j);
		rehints_list.push_back(real(hints_list[n]));
	}
	map.push_back(connected_els);
	hints.push_back(rehints_list);
   }
   if (hilbert>3000)
   {
   	lanczos_with_hints_given(iterations,1,map,hints,eigs,eigenvecs,true); 
   }
   else
   {
   	ed_with_hints_given(map,hints,eigs,eigenvecs,true); 
   	if (hilbert<100)
	{
		for (int i=0;i<eigenvecs.NRows();i++)
		{
			cout<<boost::format("%3d") %i<<"  "<<boost::format("%+.10f") %eigenvecs(i,0)<<endl;
		}
	}
   }


   std::vector<double> ev;
   for (int i=0;i<eigenvecs.NRows();i++) {ev.push_back(eigenvecs(i,0));}
   for (int i=0;i<hilbert;i++) 
   {
	convert_num_to_vec(i,base,nsites,config);
	for (int j=0;j<nsites;j++) cout<<config[j]-h.spin;
	cout<<"  "<<ev[i]<<endl;
   }	
   //print_triangle_and_hexagon_dms(base, mapv, inverse_map, 
//				   h.num_sites, ev);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//void print_triangle_and_hexagon_dms(int base, std::vector<int> mapv, 
//			            std::vector<int> inverse_map, 
//				    int nsites, std::vector<double> ev)
//{
//   std::vector<int> all_states,which_indices_to_keep;
//   for (int i=0;i<nsites;i++) {all_states.push_back(base);} 
//   
//   // Get the density matrix of a small region (triangle on 12 site) 
//   which_indices_to_keep.push_back(0);
//   which_indices_to_keep.push_back(1);
//   which_indices_to_keep.push_back(2);
//
//   Matrix den_mat;
//   int local_sites=which_indices_to_keep.size();
//   int local_hilbert=pow(base,local_sites);
//   std::vector<double> rdm_eigs(local_hilbert);
//   Matrix              rdm_eigenvecs(local_hilbert,local_hilbert);
//   
//   density_matrix_calc_with_map_general_base(base,
//			 ev, mapv, inverse_map,
//			 all_states, which_indices_to_keep,
//			 den_mat);
//   symmetric_diagonalize(den_mat,rdm_eigs,rdm_eigenvecs);
//   
//   for (int i=0;i<rdm_eigs.size();i++) cout<<rdm_eigs[i]<<endl;
//   
//   // Get the density matrix of a small region (hexagon on 15 site) 
//   which_indices_to_keep.clear();
//   which_indices_to_keep.push_back(0);
//   which_indices_to_keep.push_back(2);
//   which_indices_to_keep.push_back(3);
//   which_indices_to_keep.push_back(4);
//   which_indices_to_keep.push_back(10);
//   which_indices_to_keep.push_back(11);
//
//   local_sites=which_indices_to_keep.size();
//   local_hilbert=pow(base,local_sites);
//   std::vector<double> rdm_eigs_hex(local_hilbert);
//   Matrix              rdm_eigenvecs_hex(local_hilbert,local_hilbert);
//   
//   density_matrix_calc_with_map_general_base(base,
//			 ev, mapv, inverse_map,
//			 all_states, which_indices_to_keep,
//			 den_mat);
//   symmetric_diagonalize(den_mat,rdm_eigs_hex,rdm_eigenvecs_hex);
//   
//   for (int i=0;i<rdm_eigs_hex.size();i++) cout<<rdm_eigs_hex[i]<<endl;
//   
//}

//////////////////////////////////////////////////////////////////////////////
void lanczos_no_sym_get_eigs(Ham &h,
                             int iterations, 
                             std::vector<double> &eigs)
{
   int i,size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   std::vector<int> map,inverse_map;
   std::vector<double> spins;

   for (i=0;i<size;i++){ map.push_back(i);inverse_map.push_back(i);}
   
   lanczos_given_map(h,iterations,map,inverse_map,eigs,spins);
}
//////////////////////////////////////////////////////////////////////////////
void lanczos_requested_sz(Ham &h,
                      int iterations, 
                      std::vector<double> &eigs,
		      double s_z,
		      std::vector<double> &spins,
                      bool measure_s)
{
   int 				i,j;
   int 				n,n_up,n_down;
   int 				size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   std::vector<int> 		map;
   std::vector<int> 		inverse_map(size);
   std::vector<int> 		config(h.num_sites);   
   bool 			change;

   // Initialize
   eigs.clear();spins.clear();
   
   cout<<"==========================================="<<endl;
   cout<<"Starting to make maps"<<endl;
   cout<<"==========================================="<<endl;
   // Make maps
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
  	n_up=(int) count(config.begin(),config.end(),1);
	n_down=h.num_sites-n_up;
	if (abs((double(n_up-n_down)/2.0)-s_z)<1.0e-10) 
	{
   		map.push_back(i);
        	inverse_map[i]=map.size()-1;
  	} 
   }

   cout<<"==========================================="<<endl;
   cout<<"Finished making map in lanczos_requested_sz"<<endl;
   cout<<"==========================================="<<endl;
   
   lanczos_given_map(h,iterations,map,inverse_map,
	             eigs,spins,measure_s);

   change=true;
   n=0;
   // Sort (eigs,spins,szs)
   while (change)
   {	
	   n+=1; 
	   change=false;		
	   for (i=0;i<eigs.size()-n;i++)
	   {
		if (eigs[i]>eigs[i+1])
		{
		   change=true;
		   swap(eigs[i],eigs[i+1]);
		   if (measure_s){swap(spins[i],spins[i+1]);} 	
		}
		
		if (eigs[i]==eigs[i+1])
		{
		  if (measure_s)
		  {
			if (spins[i]>spins[i+1])
			{
				change=true;
				swap(spins[i],spins[i+1]);	
			}
		  }
		}
	   }
   }
}

//////////////////////////////////////////////////////////////////////////////
void lanczos_requested_sz_spin_1(Ham &h,
                      int iterations,std::vector<double> &eigs,
		      double s_z,std::vector<double> &spins,
                      bool measure_s)
{
   int 				n;
   int                          base=int(2.0*h.spin+1.0e-6)+1;
   int 				size=pow(base,h.num_sites);
   iterations=min(iterations,size);
   std::vector<int> 		map;
   std::vector<int> 		inverse_map(size);
   std::vector<int> 		config(h.num_sites);   
   bool 			change;

   // Initialize
   eigs.clear();spins.clear();
   
   cout<<"==========================================="<<endl;
   cout<<"Starting to make maps"<<endl;
   cout<<"==========================================="<<endl;

   // Make maps
   for (int i=0;i<size;i++)
   {
        convert_num_to_vec(i,base,h.num_sites,config);
  	int n_m1=(int) count(config.begin(),config.end(),0);
  	int n_p1=(int) count(config.begin(),config.end(),2);
	double measured_sz=double(n_p1-n_m1);
	if (abs(measured_sz-s_z)<1.0e-10) 
	{
   		map.push_back(i);
        	inverse_map[i]=map.size()-1;
  	} 
   }

   cout<<"==========================================="<<endl;
   cout<<"Finished making map in lanczos_requested_sz"<<endl;
   cout<<"==========================================="<<endl;
   
   //lanczos_given_map(h,iterations,map,inverse_map,
   //	             eigs,spins,false,true,base);
   
   lanczos_given_map_multiple_evecs(h,iterations,12,map,inverse_map,
	             eigs,spins,false,true,base);

   change=true;
   n=0;
   // Sort (eigs,spins,szs)
   while (change)
   {	
	   n+=1; 
	   change=false;		
	   for (int i=0;i<eigs.size()-n;i++)
	   {
		if (eigs[i]>eigs[i+1])
		{
		   change=true;
		   swap(eigs[i],eigs[i+1]);
		   if (measure_s){swap(spins[i],spins[i+1]);} 	
		}
		
		if (eigs[i]==eigs[i+1])
		{
		  if (measure_s)
		  {
			if (spins[i]>spins[i+1])
			{
				change=true;
				swap(spins[i],spins[i+1]);	
			}
		  }
		}
	   }
   }
}


//////////////////////////////////////////////////////////////////////////////
void lanczos_spin_sym(Ham &h,
                      int iterations, 
                      std::vector<double> &eigs,
		      std::vector<double> &szs,
		      std::vector<double> &spins,
                      bool measure_s, bool ipr, 
		      bool only_sz_non_negative)
{
   int 					i,j;
   int 					n,n_up,n_down;
   int 					n_up_max=0;
   int 					size=pow(2,h.num_sites);
   iterations=min(iterations,size);
   double 				sz;
   std::vector< std::vector<int> > 	maps(h.num_sites+1);
   std::vector<int> 			inverse_map(size);
   std::vector<int> 			config(h.num_sites);   
   std::vector<double> 			spins_sz;
   std::vector<double> 			eigs_sz;
   bool 				change;
   bool                                 allowed;
   // Initialize
   eigs.clear();spins.clear();szs.clear();

   // Make maps
   for (i=0;i<size;i++)
   {
        convert_num_to_vec(i,2,h.num_sites,config);
  	n_up=(int) count(config.begin(),config.end(),1);
   	maps[n_up].push_back(i);
	if (n_up>n_up_max) {n_up_max=n_up;}
        inverse_map[i]=maps[n_up].size()-1;
   }

   for (i=0;i<=n_up_max;i++) 
   {
     	     n_down=h.num_sites-i;
     	     sz=double(i-n_down)/2.0;
    
	     allowed=true;
             if (only_sz_non_negative)
	     {
		if (sz>-1.0e-10) {allowed=true;}
		else {allowed=false;}
	     } 
	
             if (allowed)
	     {
		     if (ipr) {cout<<"Lanczos for S_z = "<<sz<<endl;}
		     lanczos_given_map(h,iterations,
				       maps[i],inverse_map,
				       eigs_sz,spins_sz,measure_s,ipr);

		    eigs.insert(eigs.end(),eigs_sz.begin(),eigs_sz.end());
		    if (measure_s) spins.insert(spins.end(),spins_sz.begin(),spins_sz.end());
		    for (j=0;j<eigs_sz.size();j++){szs.push_back(sz);}
		    eigs_sz.clear();spins_sz.clear();	
	     }
   }

   if (ipr)
   {
	   cout<<"==========================================="<<endl;
	   cout<<"Done diagonalizing.... now sorting"<<endl;
	   cout<<"==========================================="<<endl;
   }
   //cout<<"Eigs"<<endl;print_vec(eigs);cout<<"Szs"<<endl;print_vec(szs);
 
   change=true;
   n=0;
   // Sort (eigs,spins,szs)
   /*while (change)
   {	
	   n+=1; 
	   change=false;		
	   for (i=0;i<eigs.size()-n;i++)
	   {
		if (eigs[i]>eigs[i+1])
		{
		   change=true;
		   swap(eigs[i],eigs[i+1]);
		   swap(szs[i],szs[i+1]);
		   if (measure_s){swap(spins[i],spins[i+1]);} 	
		}
		
		if (eigs[i]==eigs[i+1])
		{
		  if (measure_s)
		  {
			if (spins[i]>spins[i+1])
			{
				change=true;
				swap(spins[i],spins[i+1]);	
				swap(szs[i],szs[i+1]);
			}
			if (spins[i]==spins[i+1]) {if (szs[i]>szs[i+1]) {change=true;swap(szs[i],szs[i+1]);} }
		  }
		  else
		  {if (szs[i]>szs[i+1]) {change=true;swap(szs[i],szs[i+1]);}}  	
		}
	   }
   }i*/
   for (j=0;j<eigs.size()-1;j++)
   {	
	   for (i=j+1;i<eigs.size();i++)
	   {
		if (eigs[j]>eigs[i])
		{
		   swap(eigs[j],eigs[i]);
		   swap(szs[j],szs[i]);
		   if (measure_s){swap(spins[j],spins[i]);} 	
		}
		
		if (abs(eigs[j]-eigs[i])<1.0e-8)
		{
		  if (measure_s)
		  {
			if (spins[j]>spins[i])
			{
				swap(spins[j],spins[i]);	
				swap(szs[j],szs[i]);
			}
			if (spins[j]==spins[i]) {if (szs[j]>szs[i]) {swap(szs[j],szs[i]);} }
		  }
		  else
		  {if (szs[j]>szs[i]) {swap(szs[j],szs[i]);}}  	
		}
	   }
   }
}
///////////////////////////////////////////////////////////////////////////////////////////
void lanczos_given_map(Ham &h,
                      int iterations, 
		      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      bool measure_s, bool ipr, int base)
{
   time_t 	      				start,end;
   int 		      				i,it,j,k,l,m;
   int 		      				size=map.size();
   double 	      				spin;
   bool 	      				orth_failed;
   iterations=min(iterations,size);
   eigs.resize(iterations);
   spins.clear();
   double 					dif,tmp;
   double					q,alpha,beta,norm;
   complex<double> 				c_first_second,tot;
   std::vector<int> 				config(h.num_sites),new_configs;
   std::vector<double> 				alphas,betas;
   std::vector<double> 				h_dot_v(size),w(size);
   std::vector<double> 				v_p(size),v_o(size);
   std::vector<double> 				v_p_old(size);
   std::vector<double>				ritz_eigenvec(size);
   std::vector< complex<double> > 		hints_list;
   std::vector< std::vector<int> > 		touched_sites_list,vals_on_touched_list;
   std::vector< std::vector<int> > 		vec_new_configs;
   std::vector< std::vector<double> > 		vs;
   std::vector< std::vector< complex<double> > > vec_hints_list;
   Matrix 					t_mat(iterations,iterations);
   Matrix					t_eigenvecs(iterations,iterations);
  
   //================= 
   // Initializations
   //================= 
   if (ipr) cout<<"Making Hamiltonian" <<endl;
   cout<<"base    ="<<base<<endl;
   cout<<"Hilbert ="<<size<<endl; 
   for (i=0;i<size;i++)
   {
        v_p[i]=uniform_rnd();
        v_o[i]=0.0;w[i]=0.0;
        
        convert_num_to_vec(map[i],base,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
	new_configs.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (k=0;k<touched_sites_list.size();k++)
        {
		j=map[i];
                for (l=0;l<touched_sites_list[k].size();l++)
		{j=j+(double(vals_on_touched_list[k][l]-config[touched_sites_list[k][l]])*pow(base,h.num_sites-1-touched_sites_list[k][l]));}
		new_configs.push_back(inverse_map[j]);
        }
       
        vec_new_configs.push_back(new_configs);     // New configs 
        vec_hints_list.push_back(hints_list);       // Hints
   }
   if (ipr) cout<<"TRACE: Finished making Hamiltonian" <<endl;
       
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;
   betas.push_back(beta);

   //vs.push_back(v_p);

   if (ipr)
   {
   	cout<<"iterations = "<<iterations<<endl;
   	cout<<"vs.size()  = "<<vs.size()<<endl;
   }
 
   for (it=0;it<iterations;it++)
   {
       if (ipr)
       {
       		cout<<"================================================================="<<endl;
       		cout<<"Doing Lanczos iteration = "<<it<<endl;
       }

       if (h.num_sites>24) if (vs.size()>32){vs.erase(vs.begin());}
       
       vs.push_back(v_p);     
       time (&start);
       for (i=0;i<size;i++) // Computing H*v_p - This is the bulk of the operation
       {
            for (k=0;k<vec_new_configs[i].size();k++)
            {w[vec_new_configs[i][k]]+=(real(vec_hints_list[i][k])*v_p[i]);}
       }
       
       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);

       // Now reorthogonalize vectors
       if (vs.size()<iterations)
       {
	       orth_failed=false;
	       for (int i=0;i<vs.size();i++)
	       {
		    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=it;
			if (ipr)
			{
				cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
			}
			orth_failed=true;
		    }
		    else
		    {
			daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
			dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		    }
	       }
	       
	       norm=0.0;
	       if (orth_failed)
	       {
			while (abs(norm)<1.0e-2)
			{
				for (j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				//cout<<"norm="<<norm<<endl;
				v_p=v_p_old;
				for (int i=0;i<vs.size();i++)
				{
					q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
					//cout<<"q (after orth fail )="<<q<<endl;
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
				}
				norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
			}
			dscal(size,1.0/norm,&*v_p.begin(),1);
	       }
       }

       //vs.push_back(v_p);     
       time (&end);
       dif=difftime(end,start);
       
       if (ipr)
       {
       		cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
       		cout<<"================================================================="<<endl;
       }
   }

   if (ipr) {cout<<"vs.size() = "<<vs.size()<<endl;}
   //cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
   
   if (ipr) {cout<<"Time to build  T was "<<dif<<" seconds"<<endl;}
 
   //cout<<"Iterations="<<iterations<<endl; 
   
   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   time (&end);
   dif=difftime(end,start);

   if (ipr) {cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}

   if (measure_s)
   {
	   for (j=0;j<iterations;j++)
	   {
	      if (ipr) {cout<<"Making j= "<<j<<" Ritz eigenvector"<<endl;}

	      for (i=0;i<size;i++)
	      {
		ritz_eigenvec[i]=0.0;
		for (k=0;k<iterations;k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,j);}
	      }
	      
	      //cout<<"Finished computing Ritz eigenvector"<<endl;
	       
	      i=0;  	
	      while (abs(ritz_eigenvec[i])<1.0e-3 and i<map.size()) {i+=1;}

	      //cout<<"Found non zero element i = "<<i<<endl;
	      //cout<<"map.size()="<<map.size()<<endl;

	      convert_num_to_vec(map[i],base,h.num_sites,config);
	      tot=0.0;

	      for (int first=0;first<h.num_sites;first++)
	      {
		  for (int second=first+1;second<h.num_sites;second++)
		  {  
			touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
			calc_hints_sxsx_sysy(1.0,first,second,config,
					     touched_sites_list,
					     vals_on_touched_list,
					     hints_list);
			
			c_first_second=0.0;

			for (k=0;k<touched_sites_list.size();k++)
			{
				m=map[i];
				for (l=0;l<touched_sites_list[k].size();l++){m=m-(pow(-1,vals_on_touched_list[k][l])*pow(2,h.num_sites-1-touched_sites_list[k][l]));}
				c_first_second+=(hints_list[k]*ritz_eigenvec[inverse_map[m]]);
			}
				
			touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
			calc_hints_szsz(1.0,first,second,config,
			      touched_sites_list,vals_on_touched_list,hints_list);
			
			c_first_second+=(hints_list[0]*ritz_eigenvec[i]);

			tot+=c_first_second;	  
		}
	      }		
	     
	      tot=tot/ritz_eigenvec[i]; 
	      tmp=(2.0*real(tot))+(0.75*double(h.num_sites));
	      //cout<<"S(S+1) = "<<tmp<<endl;
	      spin=compute_spin(tmp);
	      spins.push_back(spin);
       	     
	      if (ipr)
	      {
		      cout<<"===================================="<<endl;
		      cout<<"Spin   = "<<spin<<endl;
		      cout<<"S_z    = "<<measure_sz(config)<<endl;
		      cout<<"Energy = "<<eigs[j]<<endl;	
		      cout<<"===================================="<<endl;
	      }
	}

    }

}

////////////////////////////////////////////////////////////////////////////////////////
void lanczos_given_map_return_multiple_evecs(Ham &h,
		      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector< std::vector<double> > &previous_evecs,
                      bool ipr, int base)
{
   previous_evecs.clear();
   ipr=true;
   int how_many_evecs=1;
   time_t 	      				start,end;
   int 		      				size=map.size();
   double 	      				spin;
   bool 	      				orth_failed;
   int iterations=200;                           //# M
   iterations=min(iterations,size);
   how_many_evecs=min(how_many_evecs,size);
   eigs.resize(iterations);
   double 					dif,tmp;
   double					q,alpha,beta,norm;
   complex<double> 				c_first_second,tot;
   std::vector<int> 				config(h.num_sites),new_configs;
   std::vector<double> 				alphas,betas;
   std::vector<double> 				h_dot_v(size),w(size);
   std::vector<double> 				v_p(size),v_o(size);
   std::vector<double> 				v_p_old(size);
   std::vector<double>				ritz_eigenvec(size);
   //std::vector< std::vector<int> > 		vec_new_configs(size);
   std::vector< std::vector<double> > 		vs;
   //std::vector< std::vector< complex<double> > > vec_hints_list(size);
   Matrix 					t_mat(iterations,iterations);
   Matrix					t_eigenvecs(iterations,iterations);
   std::vector<double>                          all_eigs;
   double                                       previous_eig;
   int                                          ncycles=1; 
   //================= 
   // Initializations
   //================= 
   if (ipr) cout<<"Making Hamiltonian (HERE)" <<endl;
   cout<<"base    ="<<base<<endl;
   cout<<"Hilbert ="<<size<<endl; 
   /*#pragma omp parallel for
   for (int i=0;i<size;i++)
   {
   	std::vector<int> 			config(h.num_sites),new_configs;
   	std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list;
   	std::vector< complex<double> > 		hints_list;
        convert_num_to_vec(map[i],base,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
	new_configs.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (int k=0;k<touched_sites_list.size();k++)
        {
		int j=map[i];
                for (int l=0;l<touched_sites_list[k].size();l++)
		{j=j+(double(vals_on_touched_list[k][l]-config[touched_sites_list[k][l]])*pow(base,h.num_sites-1-touched_sites_list[k][l]));}
		new_configs.push_back(inverse_map[j]);
        }
       
        vec_new_configs[i]=new_configs;     // New configs 
        vec_hints_list[i]=hints_list;       // Hints
   }*/
   if (ipr) cout<<"TRACE: Finished making Hamiltonian" <<endl;
   cout<<"how_many_evecs ="<<how_many_evecs<<endl;
   
   /*if (size<10000) 
   {
	ed_with_hints_given(vec_new_configs, vec_hints_list, eigs, t_eigenvecs, true); 
	for (int ne=0;ne<eigs.size();ne++) cout<<boost::format("Energy = %+.15f") %eigs[ne] <<endl;
   	return;
   }*/
 
   for (int num=0;num<how_many_evecs;num++)
   {
   	for (int i=0;i<size;i++) v_p[i]=2.0*uniform_rnd()-1.0;
        previous_eig=1000.0;
	for (int cycle=0;cycle<ncycles;cycle++)
	{
		cout<<"cycle ="<<cycle<<endl;
		cout<<"previous_evecs.size()"<<previous_evecs.size()<<endl;
		for (int i=0;i<previous_evecs.size();i++)
		{
				    //q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    q=ddotk(size,v_p,previous_evecs[i]);
				    cout<<"q ="<<q<<endl;
				    //daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    daxpyk(size,-q,previous_evecs[i],v_p);
				    //double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    double norminv=1.0/sqrt(ddotk(size,v_p,v_p));
				    //dscal(size,norminv,&*v_p.begin(),1);
				    dscalk(size,norminv,v_p);
		}
		//double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
		//dscal(size,norminv,&*v_p.begin(),1);
		double norminv=1.0/sqrt(ddotk(size,v_p,v_p));
		dscalk(size,norminv,v_p);
	 
		# pragma omp prallel for 
   		for (int i=0;i<size;i++){v_o[i]=0.0;w[i]=0.0;}
		vs.clear();
		beta=0.0;
		betas.clear();
		betas.push_back(beta);
		alphas.clear();
		//vs.push_back(v_p);

		if (ipr)
		{
			cout<<"iterations = "<<iterations<<endl;
			cout<<"vs.size()  = "<<vs.size()<<endl;
		}
		
		for (int it=0;it<iterations;it++)
		{
		       if (ipr)
		       {
				cout<<"================================================================="<<endl;
				cout<<"Doing Lanczos iteration = "<<it<<endl;
		       }
		       //if (h.num_sites>24) if (vs.size()>32){vs.erase(vs.begin());}
		       //vs.push_back(v_p);     
		       time (&start);
		       #pragma omp parallel for
		       //good parallelizing
		       for (int i=0;i<size;i++) // Computing H*v_p - This is the bulk of the operation
		       {
			std::vector<int> 			config(h.num_sites),new_configs;
			std::vector< std::vector<int> > 	touched_sites_list,vals_on_touched_list;
			std::vector< complex<double> > 		hints_list;
			convert_num_to_vec(map[i],base,h.num_sites,config);
			touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
			new_configs.clear();

			h(config,touched_sites_list,vals_on_touched_list,hints_list); 
			
			for (int k=0;k<touched_sites_list.size();k++)
			{
				int j=map[i];
				for (int l=0;l<touched_sites_list[k].size();l++)
				{j=j+(double(vals_on_touched_list[k][l]-config[touched_sites_list[k][l]])*pow(base,h.num_sites-1-touched_sites_list[k][l]));}
				new_configs.push_back(inverse_map[j]);
			}
			    //for (int k=0;k<vec_new_configs[i].size();k++)
			    //{w[i]+=(real(vec_hints_list[i][k])*v_p[vec_new_configs[i][k]]);}
			    for (int k=0;k<new_configs.size();k++)
			    {w[i]+=(real(hints_list[k])*v_p[new_configs[k]]);}
		       }
		       
		       //daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
		       daxpyk(size,-beta,v_o,w);
		       //alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
		       alpha=ddotk(size,w,v_p);
		       alphas.push_back(alpha);
		       //daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
		       daxpyk(size,-alpha,v_p,w);
		       v_o=v_p;
		       //beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
		       beta=sqrt(ddotk(size,w,w));
		       v_p=w;
		       //dscal(size,1.0/beta,&*v_p.begin(),1);
		       dscalk(size,1.0/beta,v_p);
		       betas.push_back(beta);
		       //dscal(size,0.0,&*w.begin(),1);
		       dscalk(size,0.0,w);

		       // Orthogonlize v_p with respect to all previous eigenvectors
		       for (int i=0;i<previous_evecs.size();i++)
		       {
				    //q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    q=ddotk(size,v_p,previous_evecs[i]);
				    //daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    daxpyk(size,-q,previous_evecs[i],v_p);
				    //double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    double norminv=1.0/sqrt(ddotk(size,v_p,v_p));
				    //dscal(size,norminv,&*v_p.begin(),1);
				    dscalk(size,norminv,v_p);
		       }	

		  //     *// Now reorthogonalize vectors
		  //     if (vs.size()<iterations)
		  //     {
		  //             orth_failed=false;
		  //             for (int i=0;i<vs.size();i++)
		  //             {
		  //      	    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		  //      	    if (abs(q)>1.0e-8)
		  //      	    {
		  //      		i=it;
		  //      		if (ipr)
		  //      		{
		  //      			cout<<"q (overlap) ="<<q<<endl;
		  //      			cout<<"--------------------------------------------------"<<endl;
		  //      			cout<<"Orthogonalization failed... choosing random vector"<<endl;
		  //      			cout<<"--------------------------------------------------"<<endl;
		  //      		}
		  //      		orth_failed=true;
		  //      	    }
		  //      	    else
		  //      	    {
		  //      		daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
		  //      		dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		  //      	    }
		  //             }
		  //             
		  //             double norm=0.0;
		  //             if (orth_failed)
		  //             {
		  //      		while (abs(norm)<1.0e-2)
		  //      		{
		  //      			for (int j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
		  //      			norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
		  //      			dscal(size,1.0/norm,&*v_p_old.begin(),1);
		  //      			//cout<<"norm="<<norm<<endl;
		  //      			v_p=v_p_old;
		  //      			for (int i=0;i<vs.size();i++)
		  //      			{
		  //      				q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
		  //      				//cout<<"q (after orth fail )="<<q<<endl;
		  //      				daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
		  //      			}
		  //      			norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
		  //      		}
		  //      		dscal(size,1.0/norm,&*v_p.begin(),1);
		  //             }
		  //     }
		       //vs.push_back(v_p);     
		       time (&end);
		       dif=difftime(end,start);
		       
		       if (ipr)
		       {
				cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
				cout<<"================================================================="<<endl;
		       }

			//if (ipr) {cout<<"vs.size() = "<<vs.size()<<endl;}
			//cout<<"Now building T matrix"<<endl;
			   
			time (&start);
			int tmp=alphas.size();
   		        t_mat.resize(tmp,tmp);
   			t_eigenvecs.resize(tmp,tmp);
			for (int j=0;j<tmp*tmp;j++) {t_mat[j]=0.0;t_eigenvecs[j]=0.0;}
			for (int j=0;j<tmp;j++)
			{
				t_mat(j,j)=alphas[j];
				if (j+1<tmp)
				{t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
			}
			time (&end);
			dif=difftime(end,start);
			if (ipr) {cout<<"Time to build  T was "<<dif<<" seconds"<<endl;}
			time (&start); symmetric_diagonalize(t_mat,eigs,t_eigenvecs); time (&end); dif=difftime(end,start);
			if (ipr) {cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}
			if (cycle>0 and it>1) 
			{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
	       			for (int ne=0;ne<eigs.size();ne++) cout<<boost::format("Energy = %+.15f Matrix el = %+.15f") %eigs[ne] %t_eigenvecs(0,ne) <<endl;
				//if (abs(eigs[0]-previous_eig)<1.0e-12) {cout<<"Convergence!!"<<endl;it=iterations;cycle=ncycles;}
				previous_eig=eigs[0];
			}
			else{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
	       			for (int ne=0;ne<eigs.size();ne++) cout<<boost::format("Energy = %+.15f Matrix el = %+.15f") %eigs[ne] %t_eigenvecs(0,ne) <<endl;
				previous_eig=eigs[0];
			    }	

		}
		if (ipr) {cout<<"Making lowest Ritz eigenvector"<<endl;}
		/*#pragma omp parallel for
		for (int i=0;i<size;i++)
		{
			ritz_eigenvec[i]=0.0;
			for (int k=0;k<vs.size();k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,0);}
		}
		//for (int i=0;i<vs.size();i++) {cout<<" Energy = "<<eigs[i]<<"  Matrix element = "<<t_eigenvecs(0,i)<<endl;}
		for (int i=0;i<size;i++) v_p[i]=ritz_eigenvec[i]+(uniform_rnd()*0.00000);	*/
	}
	//previous_evecs.push_back(ritz_eigenvec);
	//cout<<"Eigenvalue number "<<num<<" = "<<boost::format("%+.10f") %eigs[0]<<endl;
	all_eigs.push_back(eigs[0]);
	//if (num==0) evec1=ritz_eigenvec;
	//if (num==1) evec2=ritz_eigenvec;
    }
    eigs=all_eigs;
    cout<<"Done here, now exiting "<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////
void lanczos_given_map_two_evecs(Ham &h,
		      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &evec1,
		      std::vector<double> &evec2,
                      bool ipr, int base)
{
   ipr=true;
   int how_many_evecs=10;
   time_t 	      				start,end;
   int 		      				size=map.size();
   double 	      				spin;
   bool 	      				orth_failed;
   int iterations=30;
   iterations=min(iterations,size);
   eigs.resize(iterations);
   double 					dif,tmp;
   double					q,alpha,beta,norm;
   complex<double> 				c_first_second,tot;
   std::vector<int> 				config(h.num_sites),new_configs;
   std::vector<double> 				alphas,betas;
   std::vector<double> 				h_dot_v(size),w(size);
   std::vector<double> 				v_p(size),v_o(size);
   std::vector<double> 				v_p_old(size);
   std::vector<double>				ritz_eigenvec(size);
   std::vector< complex<double> > 		hints_list;
   std::vector< std::vector<int> > 		touched_sites_list,vals_on_touched_list;
   std::vector< std::vector<int> > 		vec_new_configs;
   std::vector< std::vector<double> > 		previous_evecs;
   std::vector< std::vector<double> > 		vs;
   std::vector< std::vector< complex<double> > > vec_hints_list;
   Matrix 					t_mat(iterations,iterations);
   Matrix					t_eigenvecs(iterations,iterations);
   std::vector<double>                          all_eigs;
   double                                       previous_eig;
   int                                          ncycles=100; 
   //================= 
   // Initializations
   //================= 
   if (ipr) cout<<"Making Hamiltonian" <<endl;
   cout<<"base    ="<<base<<endl;
   cout<<"Hilbert ="<<size<<endl; 
   for (int i=0;i<size;i++)
   {
        convert_num_to_vec(map[i],base,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
	new_configs.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (int k=0;k<touched_sites_list.size();k++)
        {
		int j=map[i];
                for (int l=0;l<touched_sites_list[k].size();l++)
		{j=j+(double(vals_on_touched_list[k][l]-config[touched_sites_list[k][l]])*pow(base,h.num_sites-1-touched_sites_list[k][l]));}
		new_configs.push_back(inverse_map[j]);
        }
       
        vec_new_configs.push_back(new_configs);     // New configs 
        vec_hints_list.push_back(hints_list);       // Hints
   }
   if (ipr) cout<<"TRACE: Finished making Hamiltonian" <<endl;
       
   for (int num=0;num<how_many_evecs;num++)
   {
   	for (int i=0;i<size;i++) v_p[i]=2.0*uniform_rnd()-1.0;
        previous_eig=1000.0;
	for (int cycle=0;cycle<ncycles;cycle++)
	{
		for (int i=0;i<previous_evecs.size();i++)
		{
				    q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    dscal(size,norminv,&*v_p.begin(),1);
		}
		double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
		dscal(size,norminv,&*v_p.begin(),1);
	  
   		for (int i=0;i<size;i++){v_o[i]=0.0;w[i]=0.0;}
		vs.clear();
		beta=0.0;
		betas.clear();
		betas.push_back(beta);
		alphas.clear();
		//vs.push_back(v_p);

		if (ipr)
		{
			cout<<"iterations = "<<iterations<<endl;
			cout<<"vs.size()  = "<<vs.size()<<endl;
		}
		
		for (int it=0;it<iterations;it++)
		{
		       if (ipr)
		       {
				cout<<"================================================================="<<endl;
				cout<<"Doing Lanczos iteration = "<<it<<endl;
		       }
		       //if (h.num_sites>24) if (vs.size()>32){vs.erase(vs.begin());}
		       vs.push_back(v_p);     
		       time (&start);
		       #pragma omp parallel for
		       //good parallelizing
		       for (int i=0;i<size;i++) // Computing H*v_p - This is the bulk of the operation
		       {
			    for (int k=0;k<vec_new_configs[i].size();k++)
			    {w[i]+=(real(vec_hints_list[i][k])*v_p[vec_new_configs[i][k]]);}
		       }
		       
		       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
		       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
		       alphas.push_back(alpha);
		       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
		       v_o=v_p;
		       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
		       v_p=w;
		       dscal(size,1.0/beta,&*v_p.begin(),1);
		       betas.push_back(beta);
		       dscal(size,0.0,&*w.begin(),1);

		       // Orthogonlize v_p with respect to all previous eigenvectors
		       for (int i=0;i<previous_evecs.size();i++)
		       {
				    q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    dscal(size,norminv,&*v_p.begin(),1);
		       }	

		       // Now reorthogonalize vectors
		       if (vs.size()<iterations)
		       {
			       orth_failed=false;
			       for (int i=0;i<vs.size();i++)
			       {
				    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
				    if (abs(q)>1.0e-8)
				    {
					i=it;
					if (ipr)
					{
						cout<<"q (overlap) ="<<q<<endl;
						cout<<"--------------------------------------------------"<<endl;
						cout<<"Orthogonalization failed... choosing random vector"<<endl;
						cout<<"--------------------------------------------------"<<endl;
					}
					orth_failed=true;
				    }
				    else
				    {
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
					dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
				    }
			       }
			       
			       double norm=0.0;
			       if (orth_failed)
			       {
					while (abs(norm)<1.0e-2)
					{
						for (int j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
						norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
						dscal(size,1.0/norm,&*v_p_old.begin(),1);
						//cout<<"norm="<<norm<<endl;
						v_p=v_p_old;
						for (int i=0;i<vs.size();i++)
						{
							q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
							//cout<<"q (after orth fail )="<<q<<endl;
							daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
						}
						norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
					}
					dscal(size,1.0/norm,&*v_p.begin(),1);
			       }
		       }
		       //vs.push_back(v_p);     
		       time (&end);
		       dif=difftime(end,start);
		       
		       if (ipr)
		       {
				cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
				cout<<"================================================================="<<endl;
		       }

			if (ipr) {cout<<"vs.size() = "<<vs.size()<<endl;}
			//cout<<"Now building T matrix"<<endl;
			   
			time (&start);
			int tmp=alphas.size();
   		        t_mat.resize(tmp,tmp);
   			t_eigenvecs.resize(tmp,tmp);
			for (int j=0;j<tmp*tmp;j++) {t_mat[j]=0.0;t_eigenvecs[j]=0.0;}
			for (int j=0;j<tmp;j++)
			{
				t_mat(j,j)=alphas[j];
				if (j+1<tmp)
				{t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
			}
			time (&end);
			dif=difftime(end,start);
			if (ipr) {cout<<"Time to build  T was "<<dif<<" seconds"<<endl;}
			time (&start); symmetric_diagonalize(t_mat,eigs,t_eigenvecs); time (&end); dif=difftime(end,start);
			if (ipr) {cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}
			if (cycle>=0 and it>1) 
			{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
				if (abs(eigs[0]-previous_eig)<1.0e-12) {cout<<"Convergence!!"<<endl;it=iterations;cycle=ncycles;}
				previous_eig=eigs[0];
			}
			else{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
				previous_eig=eigs[0];
			    }	

		}
		if (ipr) {cout<<"Making lowest Ritz eigenvector"<<endl;}
		#pragma omp parallel for
		for (int i=0;i<size;i++)
		{
			ritz_eigenvec[i]=0.0;
			for (int k=0;k<vs.size();k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,0);}
		}
		for (int i=0;i<size;i++) v_p[i]=ritz_eigenvec[i]+(uniform_rnd()*0.000001);	
	}
	previous_evecs.push_back(ritz_eigenvec);
	cout<<"Eigenvalue number "<<num<<" = "<<boost::format("%+.10f") %eigs[0]<<endl;
	all_eigs.push_back(eigs[0]);
	if (num==0) evec1=ritz_eigenvec;
	if (num==1) evec2=ritz_eigenvec;
    }
    eigs=all_eigs;
    cout<<"Done here, now exiting "<<endl;
}



//////////////////////////////////////////////////////////////////////////////
void lanczos_given_map_multiple_evecs(Ham &h,
                      int iterations,
		      int how_many_evecs, 
		      std::vector<int> const &map,
		      std::vector<int> const &inverse_map,	
                      std::vector<double> &eigs,
		      std::vector<double> &spins,
                      bool measure_s, bool ipr, int base)
{
   time_t 	      				start,end;
   int 		      				size=map.size();
   double 	      				spin;
   bool 	      				orth_failed;
   iterations=30;
   iterations=min(iterations,size);
   eigs.resize(iterations);
   spins.clear();
   double 					dif,tmp;
   double					q,alpha,beta,norm;
   complex<double> 				c_first_second,tot;
   std::vector<int> 				config(h.num_sites),new_configs;
   std::vector<double> 				alphas,betas;
   std::vector<double> 				h_dot_v(size),w(size);
   std::vector<double> 				v_p(size),v_o(size);
   std::vector<double> 				v_p_old(size);
   std::vector<double>				ritz_eigenvec(size);
   std::vector< complex<double> > 		hints_list;
   std::vector< std::vector<int> > 		touched_sites_list,vals_on_touched_list;
   std::vector< std::vector<int> > 		vec_new_configs;
   std::vector< std::vector<double> > 		previous_evecs;
   std::vector< std::vector<double> > 		vs;
   std::vector< std::vector< complex<double> > > vec_hints_list;
   Matrix 					t_mat(iterations,iterations);
   Matrix					t_eigenvecs(iterations,iterations);
   std::vector<double>                          all_eigs;
   double                                       previous_eig;
   int                                          ncycles=50; 
   //================= 
   // Initializations
   //================= 
   if (ipr) cout<<"Making Hamiltonian" <<endl;
   cout<<"base    ="<<base<<endl;
   cout<<"Hilbert ="<<size<<endl; 
   for (int i=0;i<size;i++)
   {
        convert_num_to_vec(map[i],base,h.num_sites,config);
        touched_sites_list.clear();vals_on_touched_list.clear();hints_list.clear();
	new_configs.clear();

        h(config,touched_sites_list,vals_on_touched_list,hints_list); 
        
        for (int k=0;k<touched_sites_list.size();k++)
        {
		int j=map[i];
                for (int l=0;l<touched_sites_list[k].size();l++)
		{j=j+(double(vals_on_touched_list[k][l]-config[touched_sites_list[k][l]])*pow(base,h.num_sites-1-touched_sites_list[k][l]));}
		new_configs.push_back(inverse_map[j]);
        }
       
        vec_new_configs.push_back(new_configs);     // New configs 
        vec_hints_list.push_back(hints_list);       // Hints
   }
   if (ipr) cout<<"TRACE: Finished making Hamiltonian" <<endl;
       
   how_many_evecs=min(how_many_evecs,size);
   for (int num=0;num<how_many_evecs;num++)
   {
   	for (int i=0;i<size;i++) v_p[i]=2.0*uniform_rnd()-1.0;
        previous_eig=1000.0;
	for (int cycle=0;cycle<ncycles;cycle++)
	{
		for (int i=0;i<previous_evecs.size();i++)
		{
				    q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    dscal(size,norminv,&*v_p.begin(),1);
		}
		double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
		dscal(size,norminv,&*v_p.begin(),1);
	  
   		for (int i=0;i<size;i++){v_o[i]=0.0;w[i]=0.0;}
		vs.clear();
		beta=0.0;
		betas.clear();
		betas.push_back(beta);
		alphas.clear();
		//vs.push_back(v_p);

		if (ipr)
		{
			cout<<"iterations = "<<iterations<<endl;
			cout<<"vs.size()  = "<<vs.size()<<endl;
		}
		
		for (int it=0;it<iterations;it++)
		{
		       if (ipr)
		       {
				cout<<"================================================================="<<endl;
				cout<<"Doing Lanczos iteration = "<<it<<endl;
		       }
		       //if (h.num_sites>24) if (vs.size()>32){vs.erase(vs.begin());}
		       vs.push_back(v_p);     
		       time (&start);
		       #pragma omp parallel for
		       //good parallelizing
		       for (int i=0;i<size;i++) // Computing H*v_p - This is the bulk of the operation
		       {
			    for (int k=0;k<vec_new_configs[i].size();k++)
			    {w[i]+=(real(vec_hints_list[i][k])*v_p[vec_new_configs[i][k]]);}
		       }
		       
		       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
		       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
		       alphas.push_back(alpha);
		       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
		       v_o=v_p;
		       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
		       v_p=w;
		       dscal(size,1.0/beta,&*v_p.begin(),1);
		       betas.push_back(beta);
		       dscal(size,0.0,&*w.begin(),1);

		       // Orthogonlize v_p with respect to all previous eigenvectors
		       for (int i=0;i<previous_evecs.size();i++)
		       {
				    q=ddot(size,&*v_p.begin(),1,&*previous_evecs[i].begin(),1);
				    daxpy(size,-q,&*previous_evecs[i].begin(),1,&*v_p.begin(),1);
				    double norminv=1.0/sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
				    dscal(size,norminv,&*v_p.begin(),1);
		       }	

		       // Now reorthogonalize vectors
		       if (vs.size()<iterations)
		       {
			       orth_failed=false;
			       for (int i=0;i<vs.size();i++)
			       {
				    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
				    if (abs(q)>1.0e-8)
				    {
					i=it;
					if (ipr)
					{
						cout<<"q (overlap) ="<<q<<endl;
						cout<<"--------------------------------------------------"<<endl;
						cout<<"Orthogonalization failed... choosing random vector"<<endl;
						cout<<"--------------------------------------------------"<<endl;
					}
					orth_failed=true;
				    }
				    else
				    {
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
					dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
				    }
			       }
			       
			       double norm=0.0;
			       if (orth_failed)
			       {
					while (abs(norm)<1.0e-2)
					{
						for (int j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
						norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
						dscal(size,1.0/norm,&*v_p_old.begin(),1);
						//cout<<"norm="<<norm<<endl;
						v_p=v_p_old;
						for (int i=0;i<vs.size();i++)
						{
							q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
							//cout<<"q (after orth fail )="<<q<<endl;
							daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
						}
						norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
					}
					dscal(size,1.0/norm,&*v_p.begin(),1);
			       }
		       }
		       //vs.push_back(v_p);     
		       time (&end);
		       dif=difftime(end,start);
		       
		       if (ipr)
		       {
				cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
				cout<<"================================================================="<<endl;
		       }

			if (ipr) {cout<<"vs.size() = "<<vs.size()<<endl;}
			//cout<<"Now building T matrix"<<endl;
			   
			time (&start);
			int tmp=alphas.size();
   		        t_mat.resize(tmp,tmp);
   			t_eigenvecs.resize(tmp,tmp);
			for (int j=0;j<tmp*tmp;j++) {t_mat[j]=0.0;t_eigenvecs[j]=0.0;}
			for (int j=0;j<tmp;j++)
			{
				t_mat(j,j)=alphas[j];
				if (j+1<tmp)
				{t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
			}
			time (&end);
			dif=difftime(end,start);
			if (ipr) {cout<<"Time to build  T was "<<dif<<" seconds"<<endl;}
			time (&start); symmetric_diagonalize(t_mat,eigs,t_eigenvecs); time (&end); dif=difftime(end,start);
			if (ipr) {cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;}
			if (cycle>=0 and it>1) 
			{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
				if (abs(eigs[0]-previous_eig)<1.0e-12) {cout<<"Convergence!!"<<endl;it=iterations;cycle=ncycles;}
				previous_eig=eigs[0];
			}
			else{
				cout<<"Eigs[0]      = "<<boost::format("%+.10f") %eigs[0]<<endl;
				cout<<"Previous Eig = "<<boost::format("%+.10f") %previous_eig<<endl;
				previous_eig=eigs[0];
			    }	

		}
		if (ipr) {cout<<"Making lowest Ritz eigenvector"<<endl;}
		#pragma omp parallel for
		for (int i=0;i<size;i++)
		{
			ritz_eigenvec[i]=0.0;
			for (int k=0;k<vs.size();k++){ritz_eigenvec[i]+=vs[k][i]*t_eigenvecs(k,0);}
		}
		for (int i=0;i<size;i++) v_p[i]=ritz_eigenvec[i]+(uniform_rnd()*0.000001);	
	}
	previous_evecs.push_back(ritz_eigenvec);
	cout<<"Eigenvalue number "<<num<<" = "<<boost::format("%+.10f") %eigs[0]<<endl;
	all_eigs.push_back(eigs[0]);
    }
    eigs=all_eigs;
}


//////////////////////////////////////////////////////////////////////////////
void lanczos_with_hints_given(int 			      		iterations,
			      int 			      		how_many_eigenvecs, 
		      	      std::vector< std::vector<int> > 		const &map,
		      	      std::vector< std::vector<double> > 	const &hints,
                      	      std::vector<double> 			&eigs,
			      Matrix 					&eigenvecs,
			      bool 					ipr)
{
   time_t 				start,end;
   time_t 				start_0,end_0;
   int 					size=map.size();
   bool 				orth_failed;
   
   time(&start_0);
   iterations=min(iterations,size);
   eigs.resize(iterations);
   
   double 				dif,tmp;
   double				q,alpha,beta,norm;
   std::vector<double> 			alphas,betas;
   std::vector<double> 			h_dot_v(size),w(size);
   std::vector<double> 			v_p(size),v_o(size),v_p_old(size),tmpv;
   std::vector< std::vector<double> >   vs;
   Matrix 				t_mat(iterations,iterations);
   Matrix				t_eigenvecs(iterations,iterations);

   int tid; int chunk=1;

   // Initializations
   //#pragma omp parallel num_threads(8)
   {
       	   //#pragma omp for schedule(dynamic,chunk)
       	   //#pragma omp for
	   for (int ip=0;ip<size;++ip)
	   {
		//tid = omp_get_thread_num();
		//cout<<"Hello from thread = "<< tid;
		v_p[ip]=uniform_rnd();
		//v_p[ip]=1.0;
		v_o[ip]=0.0;w[ip]=0.0;
	   }
	  // #pragma omp barrier
   }   
   
   norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
   dscal(size,1.0/norm,&*v_p.begin(),1);
  
   beta=0.0;betas.push_back(beta);

   for (int it=0;it<iterations;it++)
   {
       if (ipr) cout<<"Doing Lanczos iteration = "<<it<<endl;
       vs.push_back(v_p);     
       
       time (&start);
       //#pragma omp parallel shared(v_p,w,hints,map,size)
       #pragma omp parallel
       {
	       //#pragma omp for schedule(static,chunk)
	       #pragma omp for
	       for (int ip=0;ip<size;ip++)          // Computing H*v_p - This is the bulk of the operation
	       {
		   // if (abs(v_p[i])>1.0e-16)      // Preventing additions of lots of zeros //HJC June 5 2011
		      int mapsize=map[ip].size();
		      for (int kp=0;kp<mapsize;kp++)
		      {
			w[ip]+=(hints[ip][kp]*v_p[map[ip][kp]]);
		      }
	       }
       }

       daxpy(size,-beta,&*v_o.begin(),1,&*w.begin(),1);
       alpha=ddot(size,&*w.begin(),1,&*v_p.begin(),1);
       alphas.push_back(alpha);
       daxpy(size,-alpha,&*v_p.begin(),1,&*w.begin(),1);
       v_o=v_p;
       beta=sqrt(ddot(size,&*w.begin(),1,&*w.begin(),1));
       v_p=w;
       dscal(size,1.0/beta,&*v_p.begin(),1);
       betas.push_back(beta);
       dscal(size,0.0,&*w.begin(),1);

       // Now reorthogonalize vectors
       if (vs.size()<iterations)  // HJC - May 14
       {
	       orth_failed=false;
	       for (int i=0;i<it;i++)
	       {
		    q=ddot(size,&*v_p.begin(),1,&*vs[i].begin(),1);
		    if (abs(q)>1.0e-10)
		    {
			i=it;
			if (ipr)
			{
				cout<<"q (overlap) ="<<q<<endl;
				cout<<"--------------------------------------------------"<<endl;
				cout<<"Orthogonalization failed... choosing random vector"<<endl;
				cout<<"--------------------------------------------------"<<endl;
			}
			orth_failed=true;
		    }
		    else
		    {
			daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
			dscal(size,1.0/(1.0-(q*q)),&*v_p.begin(),1);
		    }
	       }
	      
	       //cout<<vs.size()<<endl;
	       norm=0.0;
	       if (orth_failed)
	       {
			while (abs(norm)<1.0e-2)
			{
				for (int j=0;j<size;j++)  {v_p_old[j]=(1.0+uniform_rnd());}
				norm=sqrt(ddot(size,&*v_p_old.begin(),1,&*v_p_old.begin(),1));
				dscal(size,1.0/norm,&*v_p_old.begin(),1);
				//cout<<"norm="<<norm<<endl;
				v_p=v_p_old;
				for (int i=0;i<it;i++)
				{
					q=ddot(size,&*v_p_old.begin(),1,&*vs[i].begin(),1);
					//cout<<"q (after orth fail )="<<q<<endl;
					daxpy(size,-q,&*vs[i].begin(),1,&*v_p.begin(),1);
				}
				norm=sqrt(ddot(size,&*v_p.begin(),1,&*v_p.begin(),1));
			}
			dscal(size,1.0/norm,&*v_p.begin(),1);
	       }
       } //HJC - May 14

       time (&end);
       dif=difftime(end,start);
       
       if (ipr) cout<<"Time to perform Lanczos iteration "<<it<<" was "<<dif<<" seconds"<<endl;
   }

   if (ipr) cout<<"Now building T matrix"<<endl;
   
   time (&start);
   for (int j=0;j<iterations;j++)
   {
        t_mat(j,j)=alphas[j];
        if (j+1<iterations)
        {t_mat(j,j+1)=betas[j+1];t_mat(j+1,j)=betas[j+1];}
   }
   time (&end);
   dif=difftime(end,start);
  
   if (ipr) cout<<"Time to build T was "<<dif<<" seconds"<<endl;
 
   time (&start);
   symmetric_diagonalize(t_mat,eigs,t_eigenvecs);
   time (&end);
   dif=difftime(end,start);

   if (ipr) cout<<"Time to diagonalize T was "<<dif<<" seconds"<<endl;
	  
   how_many_eigenvecs=min(size,how_many_eigenvecs);
 
   eigenvecs.resize(size,how_many_eigenvecs);

   //if (ipr) cout<<"how_many_eigenvecs="<<how_many_eigenvecs<<endl;//if (ipr) cout<<"size="<<size<<endl;

   for (int j=0;j<how_many_eigenvecs;j++)
   {
	      if (ipr) cout<<"Making "<<j<<"th Ritz eigenvector"<<endl;
	      #pragma omp parallel for
	      for (int i=0;i<size;i++)
	      {
		eigenvecs(i,j)=0.0;
		for (int k=0;k<iterations;k++){eigenvecs(i,j)+=vs[k][i]*t_eigenvecs(k,j);}
	      }
   }      
   //f (ipr) print_real_mat(eigenvecs);
   if (ipr) cout<<"Done computing Ritz eigenvectors"<<endl;
   time(&end_0);
   dif=difftime(end_0,start_0);
   if (ipr) cout<<"Total time for all Lanczos iterations was "<<dif<<" seconds"<<endl;
}

