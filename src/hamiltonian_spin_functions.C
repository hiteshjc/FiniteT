// Hello
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"spin_functions.h"
using namespace std;

void calc_hints_fermion_hop(double t, 
                         int first, int second, std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ S- 
	int a = config[first];
	int b = config[second];
        double sn=1.0;
        int ctr=0;
        for (int i=first+1;i<second;i++)
        {
		if (config[i]==1) ctr=ctr+1;
        }
	if (ctr%2!=0) {sn=-1.0;}

        if (a==1 and b==0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(0);
        	vals_on_touched.push_back(1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
            	hints_list.push_back(complex<double>(-t*sn,0.0));
	}
        if (a==0 and b==1)
	{
        	touched_sites.push_back(second);
        	touched_sites.push_back(first);
        	vals_on_touched.push_back(0);
        	vals_on_touched.push_back(1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
            	hints_list.push_back(complex<double>(-t*sn,0.0));
	}
}

void calc_hints_boson_hop(double t, 
                         int first, int second, std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ S- 
	int a = config[first];
	int b = config[second];
        double sn=1.0;
        int ctr=0;
        for (int i=first+1;i<second;i++)
        {
		if (config[i]==1) ctr=ctr+1;
        }
	//if (ctr%2!=0) {sn=-1.0;}

        if (a==1 and b==0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(0);
        	vals_on_touched.push_back(1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
            	hints_list.push_back(complex<double>(-t*sn,0.0));
	}
        if (a==0 and b==1)
	{
        	touched_sites.push_back(second);
        	touched_sites.push_back(first);
        	vals_on_touched.push_back(0);
        	vals_on_touched.push_back(1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
            	hints_list.push_back(complex<double>(-t*sn,0.0));
	}
}


void calc_hints_toric_fourz(double coupling,
                         std::vector< std::vector<int> > four_site_plaqs, 
                         std::vector<int>                const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> >  &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
	for (int i=0;i<four_site_plaqs.size();i++)
	{ 
		int a = config[four_site_plaqs[i][0]];
	 	int b = config[four_site_plaqs[i][1]];
		int c = config[four_site_plaqs[i][2]];
		int d = config[four_site_plaqs[i][3]];
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
		double hint=0.0;
		hint=coupling*(double(a)-0.5)*(double(b)-0.5)*(double(c)-0.5)*(double(d)-0.5);
            	hints_list.push_back(complex<double>(hint));
	}
}

///////////////////////////////////////////////////////////////////////////////
void calc_hints_toric_fourx(double coupling,
                         std::vector< std::vector<int> >  four_site_vertices, 
                         std::vector<int> const           &config,
                         std::vector< std::vector<int> >  &touched_sites_list,
                         std::vector< std::vector<int> >  &vals_on_touched_list,
                         std::vector< complex<double> >   &hints_list)
{
	for (int i=0;i<four_site_vertices.size();i++)
	{ 
		int a = config[four_site_vertices[i][0]];
	 	int b = config[four_site_vertices[i][1]];
		int c = config[four_site_vertices[i][2]];
		int d = config[four_site_vertices[i][3]];
        	std::vector<int> touched_sites, vals_on_touched;
         	touched_sites.push_back(four_site_vertices[i][0]);   	
         	touched_sites.push_back(four_site_vertices[i][1]);   	
         	touched_sites.push_back(four_site_vertices[i][2]);   	
         	touched_sites.push_back(four_site_vertices[i][3]);   	
		touched_sites_list.push_back(touched_sites);
		vals_on_touched.push_back((a+1)%2);	
	       	vals_on_touched.push_back((b+1)%2);	
	        vals_on_touched.push_back((c+1)%2);	
	        vals_on_touched.push_back((d+1)%2);	
            	vals_on_touched_list.push_back(vals_on_touched);	
		double hint=0.0;
		hint=coupling;
            	hints_list.push_back(complex<double>(hint));
	}
}


///////////////////////////////////////////////////////////////////////////////
void calc_hints_any_spin_all_sip_sjm(double spin, double j_x, double j_bq, 
                         int first, int second, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
        // Apply S+ S- 
	int a = config[first];
	int b = config[second];
        if (a<max_spin_index and b>0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(a+1);
        	vals_on_touched.push_back(b-1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);	
		double hint=0.0;
		hint=0.5*(j_x-j_bq)*splus_fn(spin,a)*sminus_fn(spin,b);
		hint+=j_bq*splus_sz_fn(spin,a)*sminus_sz_fn(spin,b);
		hint+=0.5*j_bq*splus_fn(spin,a)*sminus_sz_fn(spin,b);
		hint-=0.5*j_bq*splus_sz_fn(spin,a)*sminus_fn(spin,b);
            	hints_list.push_back(complex<double>(hint));
		//cout<<"S+ S- hint = "<<hint<<endl;
        	//cout<< "a+1 ="<<a+1<<endl;
        	//cout<< "b-1 ="<<b-1<<endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
// TERM 1 of Ring
void calc_hints_any_spin_all_s1pp_s2m_s3m(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1++ S2- S3-
        if (a<=max_spin_index-2 and b>0 and c>0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(a+2);
        	vals_on_touched.push_back(b-1);
        	vals_on_touched.push_back(c-1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*splus_splus_fn(spin,a)*sminus_fn(spin,b)*sminus_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 2 of Ring
void calc_hints_any_spin_all_s1pm_s2m_s3p(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1+- S2- S3+
        if (a>0 and b>0 and c<=max_spin_index-1)
	{
        	touched_sites.push_back(second);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(b-1);
        	vals_on_touched.push_back(c+1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*splus_sminus_fn(spin,a)*sminus_fn(spin,b)*splus_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 3 of Ring
void calc_hints_any_spin_all_s1pz_s2m_s3z(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1+z S2- S3z
        if (a<=max_spin_index-1 and b>0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(a+1);
        	vals_on_touched.push_back(b-1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*splus_sz_fn(spin,a)*sminus_fn(spin,b)*sz_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 4 of Ring
void calc_hints_any_spin_all_s1mp_s2p_s3m(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1-+ S2p S3m
        if (a<=max_spin_index-1 and b<=max_spin_index-1 and c>0)
	{
        	touched_sites.push_back(second);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(b+1);
        	vals_on_touched.push_back(c-1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*sminus_splus_fn(spin,a)*splus_fn(spin,b)*sminus_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 5 of Ring
void calc_hints_any_spin_all_s1mm_s2p_s3p(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1-- S2p S3p
        if (a>=2 and b<=max_spin_index-1 and c<=max_spin_index-1)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(a-2);
        	vals_on_touched.push_back(b+1);
        	vals_on_touched.push_back(c+1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*sminus_sminus_fn(spin,a)*splus_fn(spin,b)*splus_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 6 of Ring
void calc_hints_any_spin_all_s1mz_s2p_s3z(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1-z S2p S3z
        if (a>0 and b<=max_spin_index-1)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(a-1);
        	vals_on_touched.push_back(b+1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*sminus_sz_fn(spin,a)*splus_fn(spin,b)*sz_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 7 of Ring
void calc_hints_any_spin_all_s1zp_s2z_s3m(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1zp S2z S3m
        if (a<=max_spin_index-1 and c>0)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(a+1);
        	vals_on_touched.push_back(c-1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*sz_splus_fn(spin,a)*sz_fn(spin,b)*sminus_fn(spin,c));
        }
}
///////////////////////////////////////////////////////////////////////////////
// TERM 8 of Ring
void calc_hints_any_spin_all_s1zm_s2z_s3p(double spin, double coupling, 
                         int first, int second, int third, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];
	int c = config[third];

        // Apply S1zm S2z S3p
        if (a>0 and c<=max_spin_index-1)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(third);
        	vals_on_touched.push_back(a-1);
        	vals_on_touched.push_back(c+1);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*sz_sminus_fn(spin,a)*sz_fn(spin,b)*splus_fn(spin,c));
        }
}

///////////////////////////////////////////////////////////////////////////////
void calc_hints_any_spin_all_sipp_sjmm(double spin, double coupling, 
                         int first, int second, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        int max_spin_index=int(2.0*spin+1.0e-06);
	int a = config[first];
	int b = config[second];

        // Apply S++ S-- 
        if (a<max_spin_index-1 and b>1)
	{
        	touched_sites.push_back(first);
        	touched_sites.push_back(second);
        	vals_on_touched.push_back(a+2);
        	vals_on_touched.push_back(b-2);
            	touched_sites_list.push_back(touched_sites);
            	vals_on_touched_list.push_back(vals_on_touched);
            	hints_list.push_back(coupling*splus_splus_fn(spin,a)*sminus_sminus_fn(spin,b));
        }
}
///////////////////////////////////////////////////////////////////////////////
void calc_hints_sxsx_sysy(double coupling, 
                         int first, 
                         int second, 
                         std:: vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ S- 
        if (config[first] != config[second])
        {
            touched_sites.push_back(first);
            touched_sites.push_back(second);
            
            vals_on_touched.push_back(config[second]);
            vals_on_touched.push_back(config[first]);

            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);
            
            hints_list.push_back(coupling*0.5);
        }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_szsz(double coupling, 
                     int first, 
                     int second, 
                     std::vector<int> const &config,
                     std::vector< std::vector<int> > &touched_sites_list,
                     std::vector< std::vector<int> > &vals_on_touched_list,
                     std::vector< complex<double> > &hints_list)

{
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);
            hints_list.push_back(coupling*(double(config[first])-0.5)*(double(config[second])-0.5));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void calc_hints_szsz_all(double coupling, 
                         std::vector< std::vector<int> > const &pairs, 
                         std::vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)

{
            //szsz is diagonal in sz basis
            int first,second;
            double hint=0.0;
            std::vector< std::vector<int> >:: const_iterator p;
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);

            for (p=pairs.begin();p!=pairs.end();p++)
            {
                first=(*p)[0];second=(*p)[1];
                hint+=((double(config[first])-0.5)*(double(config[second])-0.5));
            }
            hints_list.push_back( complex<double> (coupling*hint) );
}
     
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_stag_sz(double hstag, 
                         std::vector<int> const &eta, 
                         std::vector<int> const &config,
                         std::vector< std::vector<int> > &touched_sites_list,
                         std::vector< std::vector<int> > &vals_on_touched_list,
                         std::vector< complex<double> > &hints_list)
{
            //szsz is diagonal in sz basis
            int i;
            double hint=0.0;
            std::vector<int> touched_sites, vals_on_touched;
            
            touched_sites_list.push_back(touched_sites);
            vals_on_touched_list.push_back(vals_on_touched);

            for (i=0;i<eta.size();i++)
            {
                hint+=((double(config[i])-0.5)*(double(eta[i])));
            }
            hints_list.push_back( complex<double> (hstag*hint) );
}
     

///////////////////////////////////////////////////////////////////////////////

void calc_hints_sx(double coupling, 
                   int site, 
                   std:: vector<int> const &config,
                   std::vector< std::vector<int> > &touched_sites_list,
                   std::vector< std::vector<int> > &vals_on_touched_list,
                   std::vector< complex<double> > &hints_list)
{
        std::vector<int> touched_sites, vals_on_touched;
        // Apply S+ + S - / 2 --> Flip a spin! 
        
        touched_sites.push_back(site);
        
        vals_on_touched.push_back((config[site]+1)%2); //Flipped spin half! 

        touched_sites_list.push_back(touched_sites);
        vals_on_touched_list.push_back(vals_on_touched);
        
        hints_list.push_back(coupling*0.5);
}
///////////////////////////////////////////////////////////////////////////////

void compute_c_plus_minus_i(int plus_or_minus,
			    std::vector<double> const &vec_0,
			    std::vector<double> const &vec_1,
   			    std::vector<int> const &maps_0,
   			    std::vector<int> const &inverse_map,
			    std::vector<double> &c_i)
{
      int i,j;
      int bra_config,ket_config;
      int loc;
      int num_sites=c_i.size();
      std::vector<int> config;

      //overlap=0.0;
      for (i=0;i<num_sites;i++) c_i[i]=0.0;
      
      //double norm=0.0;
      // check norms
      //for (int i=0;i<vec_0.size();i++) norm=norm+vec_0[i]*vec_0[i];
      //cout<<"Norm of 0 state"<<norm<<endl;
     
      //norm=0.0; 
      //for (int i=0;i<vec_1.size();i++) norm=norm+vec_1[i]*vec_1[i];
      //cout<<"Norm of 1 state"<<norm<<endl;
      
      if (plus_or_minus==1)
      {
              //cout<<"c_i for S=1 S_z = +1  S=0 S_z=0   FORMULA c_i=<exc|S_i_+|GS>"<<endl;
              for (i=0;i<num_sites;i++)
	      {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			if (config[i]==0)
			{
				bra_config=ket_config+pow(2,num_sites-1-i);
				loc=inverse_map[bra_config];
				c_i[i]+=(vec_1[loc]*vec_0[j]);
			}   
		   }
	      }
      }
      
      if (plus_or_minus==-1)
      {
              //cout<<"c_i for S=1 S_z = -1  S=0 S_z=0   FORMULA c_i=<exc|S_i_-|GS>"<<endl;
              for (i=0;i<num_sites;i++)
	      {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			if (config[i]==1)
			{
				bra_config=ket_config-pow(2,num_sites-1-i);
				loc=inverse_map[bra_config];
				c_i[i]+=(vec_1[loc]*vec_0[j]);
			}   
		   }
	      }
      }
      
      if (plus_or_minus==0)
      {
             //cout<<"c_i for S=1 S_z =  0  S=0 S_z=0   FORMULA c_i=<exc|S_i_z|GS>"<<endl;
             for (i=0;i<num_sites;i++)
	     {
		   for (j=0;j<maps_0.size();j++)
		   {
			ket_config=maps_0[j];
			convert_num_to_vec(ket_config,2,num_sites,config);
			c_i[i]+=(vec_0[j]*vec_1[j]*(double(config[i])-0.5));
			//overlap+=(vec_0[j]*vec_1[j]);
		   }   
	     }
	   
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_sj(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=m;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    
                            if (m!=n)
                            {
				    if (config[m]==0 and config[n]==1)
				    {
					bra_config=ket_config+pow(2,num_sites-1-m);
					bra_config=bra_config-pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					si_sj(m,n)+=(0.5*(vec_1[loc]*vec_0[j]));
				    }

				    if (config[n]==0 and config[m]==1)
				    {
					bra_config=ket_config-pow(2,num_sites-1-m);
					bra_config=bra_config+pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					si_sj(m,n)+=(0.5*(vec_1[loc]*vec_0[j]));
				    }
			    }
                            else
                            {si_sj(m,n)+=(0.5*(vec_1[j]*vec_0[j]));}

			    si_sj(m,n)+=(vec_1[j]*vec_0[j]*(double(config[m])-0.5)*(double(config[n])-0.5));
			    si_sj(n,m)=si_sj(m,n);
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_plus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==0 and config[n]==1) and (m!=n))
			    {
				bra_config=ket_config+pow(2,num_sites-1-m);
				bra_config=bra_config-pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
                            if (m==n and config[m]==1)
		            {
				si_sj(m,n)+=((vec_1[j]*vec_0[j]));
                            }
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_minus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==1 and config[n]==1) and (m!=n))
			    {
				bra_config=ket_config-pow(2,num_sites-1-m);
				bra_config=bra_config-pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_plus(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   Matrix &si_sj)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    if ( (config[m]==1 and config[n]==0) and (m!=n))
			    {
				bra_config=ket_config-pow(2,num_sites-1-m);
				bra_config=bra_config+pow(2,num_sites-1-n);
				loc=inverse_map[bra_config];
				si_sj(m,n)+=((vec_1[loc]*vec_0[j]));
			    }
                            if (m==n and config[m]==0)
		            {
				si_sj(m,n)+=((vec_1[j]*vec_0[j]));
                            }
			}
	      }
      }
}
////////////////////////////////////////////////////////////////////////////////
void compute_si_minus_sj_plus_sk_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
   	 	   std::vector<int> const &inverse_map,
		   std::vector<double> &three_pt)
{
      int 	        bra_config,ket_config;
      int 	   	loc;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      three_pt.clear();
      three_pt.resize(num_sites*num_sites*num_sites);
      three_pt.assign(num_sites*num_sites*num_sites,0.0);
 
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int l=0;l<num_sites;l++)
		       {
			       for (int j=0;j<vec_0.size();j++)
			       {
				    ket_config=maps_0[j];
				    convert_num_to_vec(ket_config,2,num_sites,config);
				    if ((config[m]==1 and config[n]==0) and (m!=n))
				    {
					bra_config=ket_config-pow(2,num_sites-1-m);
					bra_config=bra_config+pow(2,num_sites-1-n);
					loc=inverse_map[bra_config];
					three_pt[num_sites*num_sites*m+num_sites*n+l]+=((vec_1[loc]*vec_0[j]*(double(config[l])-0.5)));
				    }
				    if (m==n and config[m]==0)
				    {
					three_pt[num_sites*num_sites*m+num_sites*n+l]+=((vec_1[j]*vec_0[j])*(double(config[l])-0.5));
				    }
				}
			}
	      }
      }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_z_sj_z(int num_sites,
                   std::vector<double> const &vec_0,
		   std::vector<double> const &vec_1,
   		   std::vector<int> const &maps_0,
		   Matrix &si_sj)
{
      int 	        ket_config;
      int 	   	num_pairs=num_sites*num_sites;
      std::vector<int>  config;

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int m=0;m<num_sites;m++)
      {
	      for (int n=0;n<num_sites;n++)
	      {
		       for (int j=0;j<vec_0.size();j++)
		       {
			    ket_config=maps_0[j];
			    convert_num_to_vec(ket_config,2,num_sites,config);
			    si_sj(m,n)+=((vec_1[j]*vec_0[j]*(double(config[m])-0.5)*(double(config[n])-0.5)));
		       }
	      }
      }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
//// Ron Edit
/////////////////////////////////////////////////////////////////////////////////////////////////////

void compute_sum_siz2(int num_sites,
                   std::vector<double> const &vec_bra,
                   std::vector<double> const &vec_ket,
                   std::vector<int> const &maps_0,
                   double &matrix_el)
{
      int               ket_config;
      int               num_pairs=num_sites*num_sites;
      std::vector<int>  config;
      matrix_el=0.0; 

      for (int j=0;j<vec_bra.size();j++)
      {
              ket_config=maps_0[j];
              convert_num_to_vec(ket_config,3,num_sites,config);
	      double temp=0.0;
	      for (int m=0; m<num_sites; m++) temp+=pow((double(config[m])-1),2.0);
              matrix_el+=((vec_ket[j]*vec_bra[j]*temp));
      }
}
