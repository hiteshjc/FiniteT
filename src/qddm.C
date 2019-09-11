#include"qddm.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//void update_qddm(Matrix &,int,int,Matrix &);
//
////////////////////////////////////////////////////////////////////////////////
//void calc_qddm_pavan_tfim_exact(Ham &h, std::vector<int> kept_indices,
//               			std::vector<double> &eigs, int nlow)
//{
//   cout<<"QDDM for Pavan's spin 1/2 TFIM...."<<endl;
//   int nsites=h.num_sites;
//   int i,j,n,size=pow(2,nsites);
//   std::vector<int> map,inverse_map,all_states;
//   Matrix density_mat,matrix_of_dm,eigenvecs;
//   std::vector<std::vector<double> > evecs;
//   int nlocal=kept_indices.size();
//   int qddm_size=nlow*pow(2,nlocal);
//   int kept_size=pow(2,nlocal);
//   matrix_of_dm.resize(qddm_size,qddm_size);
//   cout <<"QDDM dimension = "<<qddm_size<<endl;
//
//   for (int i=0;i<size;i++)
//   {
//	map.push_back(i);
//	inverse_map.push_back(i);
//   }
//   
//   for (int i=0;i<nlow;i++)
//   {
//	for (int j=0;j<nlow;j++)
//	{
//		cout<<"i,j = "<<i<<", "<<j<<endl;
//		if (i==j)
//		{
//   			density_matrix_calc_with_map(evecs[i],map, inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		else
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs[i], 
//					map,evecs[j],
//					map,inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//   		int shift_r=i*kept_size;int shift_c=j*kept_size;
//   		update_qddm(density_mat,shift_r,shift_c,matrix_of_dm);
//	}
//   }
//
//   cout<<"Matrix of DM = QDDM"<<endl;
//   print_real_mat(matrix_of_dm);
//   eigs.resize(qddm_size);
//   symmetric_diagonalize(matrix_of_dm, eigs, eigenvecs);
//   cout<<"nlow ="<<nlow<<endl;
//   Matrix rho(nlow*nlow,kept_size*kept_size);
//   print_real_mat(rho);
//   
//   // Henley rho[mn,ab]
//   for (int m=0;m<nlow;m++)
//   {
//	for (int n=0;n<nlow;n++)
//	{
//		for (int a=0;a<kept_size;a++)
//		{
//			for (int b=0;b<kept_size;b++)
//		        {
//				rho((m*nlow)+n,(a*kept_size)+b)=matrix_of_dm( (m*kept_size)+a ,(n*kept_size)+b );
//			}		
//		}
//	}
//   }
//   
//   cout<<"RHO"<<endl;
//   print_real_mat(rho);
//   eigs.resize(kept_size*kept_size);
//   // Now take Henley's rho[mn,ab] and make the M Matrix out of it (see Henley, Changlani review paper)
//   Matrix m_matrix(kept_size*kept_size,kept_size*kept_size);
//   for (int i=0;i<kept_size*kept_size*kept_size*kept_size;i++) m_matrix[i]=0.0;
//   
//   for (int alpha=0;alpha<kept_size*kept_size;alpha++)
//   {
//   for (int beta=0;beta<kept_size*kept_size;beta++)
//   {
//   for (int m=0;m<nlow;m++)
//   {
//   for (int n=0;n<nlow;n++)
//   {
//	int mm=(m*nlow)+m;
//	int nn=(n*nlow)+n;
//	int mn=(m*nlow)+n;
//	//m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta)));
//	m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta))-((1.0/double(nlow))*rho(mm,alpha)*rho(nn,beta)));
//   }
//   }
//   }
//   }
//
//   Matrix u(kept_size*kept_size,kept_size*kept_size), vt(kept_size*kept_size,kept_size*kept_size);
//   cout<<"M matrix Rows ="<<m_matrix.NRows()<<endl;
//   cout<<"M matrix Cols ="<<m_matrix.NCols()<<endl;
//   cout<<"kept size     ="<<kept_size<<endl;
//   // Do SVD on M matrix
//   svd(m_matrix, eigs, u, vt);
//   cout<<"Left Eigenvectors of M matrix (u)  "<<endl;
//   print_real_mat(u);
//   cout<<"Right Eigenvectors of M matrix (vt)"<<endl;
//   print_real_mat2(vt);
//   cout<<"Eigenvalues of M matrix (eigs)"<<endl;
//   print_vec_acc(eigs,true);
//   
//   std::vector<double> dominant0,dominant1,dominant2;
//   for (int i=0;i<kept_size*kept_size;i++) dominant0.push_back(u(i,0));
//   for (int i=0;i<kept_size*kept_size;i++) dominant1.push_back(u(i,1));
//   for (int i=0;i<kept_size*kept_size;i++) dominant2.push_back(u(i,2));
//   Matrix dominant0_mat(kept_size,kept_size);
//   Matrix dominant1_mat(kept_size,kept_size);
//   Matrix dominant2_mat(kept_size,kept_size);
//   cout<<"Recasting dominant eigenvector as a matrix....."<<endl;
//   for (int i=0;i<kept_size*kept_size;i++) 
//   {
//        int m=i/kept_size;
//        int n=i-(m*kept_size);
//	cout<<m<<","<<n<<" "<<dominant0[i]<<endl;
//        dominant0_mat(m,n)=dominant0[i];
//        dominant1_mat(m,n)=dominant1[i];
//        dominant2_mat(m,n)=dominant2[i];
//   }
//   cout<<"Printing dominant 0 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant0_mat);
//   cout<<"Printing dominant 1 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant1_mat);
//   cout<<"Printing dominant 2 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant2_mat);
//
//   Matrix tmp,tmp2;
//   cout<<"Matrix mult"<<endl;
//   real_matrix_multiply(dominant0_mat, dominant1_mat,tmp);
//   real_matrix_multiply(dominant1_mat, dominant0_mat,tmp2);
//   for (int i=0;i<tmp.NRows()*tmp.NCols();i++) tmp[i]=tmp[i]-tmp2[i];
//   print_real_mat(tmp);
//}
//
////////////////////////////////////////////////////////////////////////////////
//void calc_qddm_spin1b2(int sz, Ham &h,
//               int iterations,
//	       std::vector<int> kept_indices,
//               std::vector<double> &eigs,
//               int nlow)
//{
//   cout<<"QDDM for spin 1/2...."<<endl;
//   int nsites=h.num_sites;
//   std::vector<double> szs;
//   std::vector<double> spins;
//   int i,j,n,size=pow(2,nsites),n_up_max=0,chosen;
//   iterations=min(iterations,size);
//   std::vector< std::vector<int> > maps(1);
//   std::vector<int> config;
//   std::vector<int> all_states;
//   std::vector<int> inverse_map(size);
//   std::vector<double> spins_0,spins_p_1,spins_m_1;
//   std::vector<double> eigs_0,eigs_p_1,eigs_m_1;
//   std::vector< std::vector<double> > evecs_0;
//   int ket_config,bra_config,loc;
//   int ctr;
//   double max;
//   double overlap=0.0;
//   Matrix density_mat;
//   Matrix matrix_of_dm;
//   Matrix eigenvecs;
//
//   // Initialize
//   eigs.clear();
//   cout<<"Initialize..."<<endl;
//   for (int i=0;i<size;i++)
//   {
//        convert_num_to_vec(i,2,nsites,config);
//  	int n_up=(int) count(config.begin(),config.end(),1);
//  	int n_down=(int) count(config.begin(),config.end(),0);
//	// Sz  = +0 states
//	if (n_up-n_down==2*sz) 
//	{
//   		maps[0].push_back(i);
//        	inverse_map[i]=maps[0].size()-1;
//  	} 
//   }
//   cout<<"Map size ="<<maps[0].size()<<endl;
//   cout<<"Lanczos for spins ..."<<endl;
//   evecs_0.clear();
//   evecs_0.push_back(std::vector<double> ());
//   evecs_0.push_back(std::vector<double> ());
//   lanczos_given_map_return_multiple_evecs(h,maps[0],inverse_map,eigs_0,evecs_0,false,2);
//   eigs.insert(eigs.end(),eigs_0.begin(),eigs_0.begin()+nlow);
//   szs.push_back(0.0);szs.push_back(0.0);
//	
//   cout<<"Done with Lanczos for spin 1/2"<<endl;
//   if (h.pairs1.size()>0)
//   {
//   	cout<<"Mathematica pairs1 for graph"<<endl;
//   	print_mathematica_pairs(h.pairs1);
//   }
//}
//

//////////////////////////////////////////////////////////////////////////////
void calc_qddm_spin1(int sz, Ham &h,
               int iterations,
	       std::vector<int> kept_indices,
               std::vector<double> &eigs,
               int nlow)
{
   cout<<"QDDM for spin 1...."<<endl;
   int nsites=h.num_sites;
   std::vector<double> szs;
   std::vector<double> spins;
   int i,j,n,size=pow(3,nsites),n_up_max=0,chosen;
   iterations=min(iterations,size);
   std::vector< std::vector<int> > maps(1);
   std::vector<int> config;
   std::vector<int> all_states;
   std::vector<int> inverse_map(size);
   std::vector<double> spins_0,spins_p_1,spins_m_1;
   std::vector<double> eigs_0,eigs_p_1,eigs_m_1;
   std::vector< std::vector<double> > evecs_0;
   int ket_config,bra_config,loc;
   int ctr;
   double max;
   double overlap=0.0;
   Matrix density_mat;
   Matrix matrix_of_dm;
   Matrix eigenvecs;

   // Initialize
   eigs.clear();
   cout<<"Initialize..."<<endl;
   for (int i=0;i<size;i++)
   {
        convert_num_to_vec(i,3,nsites,config);
  	int n_up=(int) count(config.begin(),config.end(),0);
  	int n_down=(int) count(config.begin(),config.end(),2);
	// Sz  = +0 states
	if (n_up-n_down==sz) 
	{
   		maps[0].push_back(i);
        	inverse_map[i]=maps[0].size()-1;
  	} 
   }
   cout<<"Map size ="<<maps[0].size()<<endl;
   cout<<"Lanczos for spins ..."<<endl;
   evecs_0.clear();
   evecs_0.push_back(std::vector<double> ());
   evecs_0.push_back(std::vector<double> ());
   lanczos_given_map_return_multiple_evecs(h,maps[0],inverse_map,eigs_0,evecs_0,false,3);
   //eigs.insert(eigs.end(),eigs_0.begin(),eigs_0.begin()+nlow);
   //szs.push_back(0.0);szs.push_back(0.0);
   //     
   //cout<<"Done with Lanczos for spin 1"<<endl;
   //if (h.pairs1.size()>0)
   //{
   //	cout<<"Mathematica pairs1 for graph"<<endl;
   //	print_mathematica_pairs(h.pairs1);
   //}
   //if (h.pairs2.size()>0)
   //{
   //	cout<<"Mathematica pairs2 for graph"<<endl;
   //	print_mathematica_pairs(h.pairs2);
   //	cout<<endl;
   //}

   //cout<<"Energies for Sz = "<<sz<<endl;
   //print_vec_acc(eigs_0,true);
   //for (int i=0;i<nsites;i++) {all_states.push_back(3);}
   //
// // kept_indices.push_back(0);
// // kept_indices.push_back(1);
// // kept_indices.push_back(2);
   //int nlocal=kept_indices.size();

   //int qddm_size=nlow*pow(3,nlocal);
   //int kept_size=pow(3,nlocal);
   //matrix_of_dm.resize(qddm_size,qddm_size);
   //cout <<"QDDM dimension = "<<qddm_size<<endl;

   //for (int i=0;i<nlow;i++)
   //{
   //     for (int j=0;j<nlow;j++)
   //     {
   //     	cout<<"i,j = "<<i<<", "<<j<<endl;
   //     	if (i==j)
   //     	{
   //			density_matrix_calc_with_map(evecs_0[i], maps[0],inverse_map,
   //     				             all_states,kept_indices,density_mat);
   //     	}
   //     	else
   //     	{
   //        		off_diag_density_matrix_calc_with_map(evecs_0[i], 
   //     				maps[0],evecs_0[j],
   //     				maps[0],inverse_map,
   //     				all_states,kept_indices,density_mat);
   //     	}
   //		int shift_r=i*kept_size;int shift_c=j*kept_size;
   //		update_qddm(density_mat,shift_r,shift_c,matrix_of_dm);
   //     }
   //}
   //cout<<"Matrix of DM = QDDM"<<endl;
   //print_real_mat(matrix_of_dm);
   //eigs.resize(qddm_size);
   //symmetric_diagonalize(matrix_of_dm, eigs, eigenvecs);
   //cout<<"nlow ="<<nlow<<endl;
   //Matrix rho(nlow*nlow,kept_size*kept_size);
   //print_real_mat(rho);
   //
   //// Henley rho[mn,ab]
   //for (int m=0;m<nlow;m++)
   //{
   //     for (int n=0;n<nlow;n++)
   //     {
   //     	for (int a=0;a<kept_size;a++)
   //     	{
   //     		for (int b=0;b<kept_size;b++)
   //     	        {
   //     			rho((m*nlow)+n,(a*kept_size)+b)=matrix_of_dm( (m*kept_size)+a ,(n*kept_size)+b );
   //     		}		
   //     	}
   //     }
   //}
   //
   //cout<<"RHO"<<endl;
   //print_real_mat(rho);
   //eigs.resize(kept_size*kept_size);
   //// Now take Henley's rho[mn,ab] and make the M Matrix out of it (see Henley, Changlani review paper)
   //Matrix m_matrix(kept_size*kept_size,kept_size*kept_size);
   //for (int i=0;i<kept_size*kept_size*kept_size*kept_size;i++) m_matrix[i]=0.0;
   //
   //for (int alpha=0;alpha<kept_size*kept_size;alpha++)
   //{
   //for (int beta=0;beta<kept_size*kept_size;beta++)
   //{
   //for (int m=0;m<nlow;m++)
   //{
   //for (int n=0;n<nlow;n++)
   //{
   //     int mm=(m*nlow)+m;
   //     int nn=(n*nlow)+n;
   //     int mn=(m*nlow)+n;
   //     //m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta)));
   //     m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta))-((1.0/double(nlow))*rho(mm,alpha)*rho(nn,beta)));
   //}
   //}
   //}
   //}

   //Matrix u(kept_size*kept_size,kept_size*kept_size), vt(kept_size*kept_size,kept_size*kept_size);
   //cout<<"M matrix Rows ="<<m_matrix.NRows()<<endl;
   //cout<<"M matrix Cols ="<<m_matrix.NCols()<<endl;
   //cout<<"kept size     ="<<kept_size<<endl;
   //// Do SVD on M matrix
   //svd(m_matrix, eigs, u, vt);
   //cout<<"Left Eigenvectors of M matrix (u)  "<<endl;
   //print_real_mat(u);
   //cout<<"Right Eigenvectors of M matrix (vt)"<<endl;
   //print_real_mat2(vt);
   //cout<<"Eigenvalues of M matrix (eigs)"<<endl;
   //print_vec_acc(eigs,true);
   //
   //std::vector<double> dominant0,dominant1,dominant2;
   //for (int i=0;i<kept_size*kept_size;i++) dominant0.push_back(u(i,0));
   //for (int i=0;i<kept_size*kept_size;i++) dominant1.push_back(u(i,1));
   //for (int i=0;i<kept_size*kept_size;i++) dominant2.push_back(u(i,2));
   //Matrix dominant0_mat(kept_size,kept_size);
   //Matrix dominant1_mat(kept_size,kept_size);
   //Matrix dominant2_mat(kept_size,kept_size);
   //cout<<"Recasting dominant eigenvector as a matrix....."<<endl;
   //for (int i=0;i<kept_size*kept_size;i++) 
   //{
   //     int m=i/kept_size;
   //     int n=i-(m*kept_size);
   //     cout<<m<<","<<n<<" "<<dominant0[i]<<endl;
   //     dominant0_mat(m,n)=dominant0[i];
   //     dominant1_mat(m,n)=dominant1[i];
   //     dominant2_mat(m,n)=dominant2[i];
   //}
   //cout<<"Printing dominant 0 eigenvector as a matrix....."<<endl;
   //print_real_mat(dominant0_mat);
   //cout<<"Printing dominant 1 eigenvector as a matrix....."<<endl;
   //print_real_mat(dominant1_mat);
   //cout<<"Printing dominant 2 eigenvector as a matrix....."<<endl;
   //print_real_mat(dominant2_mat);

   //Matrix tmp,tmp2;
   //cout<<"Matrix mult"<<endl;
   //real_matrix_multiply(dominant0_mat, dominant1_mat,tmp);
   //real_matrix_multiply(dominant1_mat, dominant0_mat,tmp2);
   //for (int i=0;i<tmp.NRows()*tmp.NCols();i++) tmp[i]=tmp[i]-tmp2[i];
   //print_real_mat(tmp);
}


////////////////////////////////////////////////////////////////////////////////
//void calc_qddm_spins(Ham &h,
//               int iterations, 
//               std::vector<double> &eigs)
//{
//   cout<<"QDDM spins...."<<endl;
//   int nsites=h.num_sites;
//   std::vector<double> szs;
//   std::vector<double> spins;
//   int i,j,n,size=pow(2,nsites),n_up_max=0,chosen;
//   iterations=min(iterations,size);
//   std::vector< std::vector<int> > maps(3);
//   std::vector<int> config;
//   std::vector<int> all_states,kept_indices;
//   std::vector<int> inverse_map(size);
//   std::vector<double> spins_0,spins_p_1,spins_m_1;
//   std::vector<double> eigs_0,eigs_p_1,eigs_m_1;
//   std::vector< std::vector<double> > evecs_0,evecs_p_1,evecs_m_1;
//   int ket_config,bra_config,loc;
//   int ctr;
//   double max;
//   double overlap=0.0;
//   Matrix density_mat;
//   Matrix matrix_of_dm;
//   Matrix eigenvecs;
//
//   // Initialize
//   eigs.clear();
//   cout<<"Initialize..."<<endl;
//   for (int i=0;i<size;i++)
//   {
//        convert_num_to_vec(i,2,nsites,config);
//  	int n_up=(int) count(config.begin(),config.end(),1);
//	int n_down=h.num_sites-n_up;
//
//	// Sz  = +0 states
//	if (abs((double(n_up-n_down)/2.0)-0.0)<1.0e-10) 
//	{
//   		maps[0].push_back(i);
//        	inverse_map[i]=maps[0].size()-1;
//  	} 
//	
//	// Sz  = +1 states
//	if (abs((double(n_up-n_down)/2.0)-1.0)<1.0e-10) 
//	{
//   		maps[1].push_back(i);
//        	inverse_map[i]=maps[1].size()-1;
//  	} 
//	// Sz  = -1 states
//	if (abs((double(n_up-n_down)/2.0)+1.0)<1.0e-10) 
//	{
//   		maps[2].push_back(i);
//        	inverse_map[i]=maps[2].size()-1;
//  	} 
//   }
//   cout<<"Map size ="<<maps[0].size()<<endl;
//   cout<<"Lanczos for spins ..."<<endl;
//   int nlow=4;
//   lanczos_given_map_evecs(h,iterations,maps[0],inverse_map,eigs_0,spins_0,evecs_0);
//   lanczos_given_map_evecs(h,iterations,maps[1],inverse_map,eigs_p_1,spins_p_1,evecs_p_1);
//   lanczos_given_map_evecs(h,iterations,maps[2],inverse_map,eigs_m_1,spins_m_1,evecs_m_1);
//
//   eigs.insert(eigs.end(),eigs_0.begin(),eigs_0.begin()+2);
//   eigs.insert(eigs.end(),eigs_p_1.begin(),eigs_p_1.begin()+1);
//   eigs.insert(eigs.end(),eigs_m_1.begin(),eigs_m_1.begin()+1);
//   
//   spins.insert(spins.end(),spins_0.begin(),spins_0.begin()+2);
//   spins.insert(spins.end(),spins_p_1.begin(),spins_p_1.begin()+1);
//   spins.insert(spins.end(),spins_m_1.begin(),spins_m_1.begin()+1);
//
//   szs.push_back(0.0);szs.push_back(0.0);szs.push_back(1.0);szs.push_back(-1.0);
//
//   if (h.pairs1.size()>0)
//   {
//   	cout<<"Mathematica pairs1 for graph"<<endl;
//   	print_mathematica_pairs(h.pairs1);
//   }
//   if (h.pairs2.size()>0)
//   {
//   	cout<<"Mathematica pairs2 for graph"<<endl;
//   	print_mathematica_pairs(h.pairs2);
//   	cout<<endl;
//   }
//
//   cout<<"Energies for Sz = +0"<<endl;
//   print_vec_acc(eigs_0,true);
//   cout<<"Energies for Sz = +1"<<endl;
//   print_vec_acc(eigs_p_1,true);
//   cout<<"Energies for Sz = -1"<<endl;
//   print_vec_acc(eigs_m_1,true);
//   //for (int i=0;i<eigs.size();i++) {cout<<eigs[i]<<"  "<<endl;}
//   for (int i=0;i<nsites;i++) {all_states.push_back(2);}
//   
//   kept_indices.push_back(0);
//   kept_indices.push_back(1);
//   kept_indices.push_back(2);
//   //kept_indices.push_back(12);
//   //kept_indices.push_back(6);
//   //kept_indices.push_back(9);
//   //kept_indices.push_back(10);
//   //kept_indices.push_back(11);
//   //kept_indices.push_back(12);
//   //kept_indices.push_back(6);
//   int nlocal=kept_indices.size();
//
//   int qddm_size=nlow*pow(2,nlocal);
//   int kept_size=pow(2,nlocal);
//   matrix_of_dm.resize(qddm_size,qddm_size);
//   cout <<"QDDM dimension = "<<qddm_size<<endl;
//
//   for (int i=0;i<nlow;i++)
//   {
//	for (int j=0;j<nlow;j++)
//	{
//		cout<<"i,j = "<<i<<", "<<j<<endl;
//		if (i==j and i==0)
//		{
//   			density_matrix_calc_with_map(evecs_0[i], maps[0],inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		if (i==j and i==1)
//		{
//   			density_matrix_calc_with_map(evecs_0[i], maps[0],inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		if (i==j and i==2)
//		{
//   			density_matrix_calc_with_map(evecs_p_1[0], maps[1],inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		if (i==j and i==3)
//		{
//   			density_matrix_calc_with_map(evecs_m_1[0], maps[2],inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		if (i==0 and j==1)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[0], 
//					maps[0],evecs_0[1],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==0 and j==2)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[0], 
//					maps[0],evecs_p_1[0],
//					maps[1],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==0 and j==3)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[0], 
//					maps[0],evecs_m_1[0],
//					maps[2],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==1 and j==0)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[1], 
//					maps[0],evecs_0[0],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==1 and j==2)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[1], 
//					maps[0],evecs_p_1[0],
//					maps[1],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==1 and j==3)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[1], 
//					maps[0],evecs_m_1[0],
//					maps[2],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==2 and j==0)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_p_1[0], 
//					maps[1],evecs_0[0],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==2 and j==1)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_p_1[0], 
//					maps[1],evecs_0[1],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==2 and j==3)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_p_1[0], 
//					maps[1],evecs_m_1[0],
//					maps[2],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		
//		if (i==3 and j==0)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_m_1[0], 
//					maps[2],evecs_0[0],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==3 and j==1)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_m_1[0], 
//					maps[2],evecs_0[1],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//		if (i==3 and j==2)
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_m_1[0], 
//					maps[2],evecs_p_1[0],
//					maps[1],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//
//   		int shift_r=i*kept_size;int shift_c=j*kept_size;
//   		update_qddm(density_mat,shift_r,shift_c,matrix_of_dm);
//	}
//   }
//   cout<<"Matrix of DM = QDDM"<<endl;
//   print_real_mat(matrix_of_dm);
//   eigs.resize(qddm_size);
//   symmetric_diagonalize(matrix_of_dm, eigs, eigenvecs);
//   cout<<"nlow ="<<nlow<<endl;
//   //Matrix a(1,3);
//   //print_real_mat(a);
//   Matrix rho(nlow*nlow,kept_size*kept_size);
//   print_real_mat(rho);
//   
//   // Henley rho[mn,ab]
//   for (int m=0;m<nlow;m++)
//   {
//	for (int n=0;n<nlow;n++)
//	{
//		for (int a=0;a<kept_size;a++)
//		{
//			for (int b=0;b<kept_size;b++)
//		        {
//				rho((m*nlow)+n,(a*kept_size)+b)=matrix_of_dm( (m*kept_size)+a ,(n*kept_size)+b );
//			}		
//		}
//	}
//   }
//   
//   cout<<"RHO"<<endl;
//   print_real_mat(rho);
//   eigs.resize(kept_size*kept_size);
//   // Now take Henley's rho[mn,ab] and make the M Matrix out of it (see Henley, Changlani review paper)
//   Matrix m_matrix(kept_size*kept_size,kept_size*kept_size);
//   for (int i=0;i<kept_size*kept_size*kept_size*kept_size;i++) m_matrix[i]=0.0;
//   
//   for (int alpha=0;alpha<kept_size*kept_size;alpha++)
//   {
//   for (int beta=0;beta<kept_size*kept_size;beta++)
//   {
//   for (int m=0;m<nlow;m++)
//   {
//   for (int n=0;n<nlow;n++)
//   {
//	int mm=(m*nlow)+m;
//	int nn=(n*nlow)+n;
//	int mn=(m*nlow)+n;
//	//m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta)));
//	m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta))-((1.0/double(nlow))*rho(mm,alpha)*rho(nn,beta)));
//   }
//   }
//   }
//   }
//
//   Matrix u(kept_size*kept_size,kept_size*kept_size), vt(kept_size*kept_size,kept_size*kept_size);
//   cout<<"M matrix Rows ="<<m_matrix.NRows()<<endl;
//   cout<<"M matrix Cols ="<<m_matrix.NCols()<<endl;
//   cout<<"kept size     ="<<kept_size<<endl;
//   // Do SVD on M matrix
//   svd(m_matrix, eigs, u, vt);
//   cout<<"Left Eigenvectors of M matrix (u)  "<<endl;
//   print_real_mat(u);
//   cout<<"Right Eigenvectors of M matrix (vt)"<<endl;
//   print_real_mat2(vt);
//   cout<<"Eigenvalues of M matrix (eigs)"<<endl;
//   print_vec_acc(eigs,true);
//   
//   std::vector<double> dominant0,dominant1,dominant2;
//   for (int i=0;i<kept_size*kept_size;i++) dominant0.push_back(u(i,0));
//   for (int i=0;i<kept_size*kept_size;i++) dominant1.push_back(u(i,1));
//   for (int i=0;i<kept_size*kept_size;i++) dominant2.push_back(u(i,2));
//   Matrix dominant0_mat(kept_size,kept_size);
//   Matrix dominant1_mat(kept_size,kept_size);
//   Matrix dominant2_mat(kept_size,kept_size);
//   cout<<"Recasting dominant eigenvector as a matrix....."<<endl;
//   for (int i=0;i<kept_size*kept_size;i++) 
//   {
//        int m=i/kept_size;
//        int n=i-(m*kept_size);
//	cout<<m<<","<<n<<" "<<dominant0[i]<<endl;
//        dominant0_mat(m,n)=dominant0[i];
//        dominant1_mat(m,n)=dominant1[i];
//        dominant2_mat(m,n)=dominant2[i];
//   }
//   cout<<"Printing dominant 0 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant0_mat);
//   cout<<"Printing dominant 1 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant1_mat);
//   cout<<"Printing dominant 2 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant2_mat);
//
//   Matrix tmp,tmp2;
//   cout<<"Matrix mult"<<endl;
//   real_matrix_multiply(dominant0_mat, dominant1_mat,tmp);
//   real_matrix_multiply(dominant1_mat, dominant0_mat,tmp2);
//   for (int i=0;i<tmp.NRows()*tmp.NCols();i++) tmp[i]=tmp[i]-tmp2[i];
//   print_real_mat(tmp);
//}
//
////////////////////////////////////////////////////////////////////////////////
//void calc_qddm_fermions(Ham &h,
//               int iterations, 
//               std::vector<double> &eigs)
//{
//   int filling=3;
//   cout<<"Here...."<<endl;
//   int nsites=h.num_sites;
//   int nfermions=nsites/filling;
//   int i,j,n,n_up,n_down,size=pow(2,nsites),n_up_max=0,chosen;
//   iterations=min(iterations,size);
//   std::vector< std::vector<int> > maps(1);
//   std::vector<int> config;
//   std::vector<int> all_states,kept_indices;
//   std::vector<int> inverse_map(size);
//   std::vector<double> spins_0,spins_p_1,spins_m_1;
//   std::vector<double> eigs_0,eigs_p_1,eigs_m_1;
//   std::vector< std::vector<double> > evecs_0,evecs_p_1,evecs_m_1;
//   int ket_config,bra_config,loc;
//   int ctr;
//   double max;
//   double overlap=0.0;
//   Matrix density_mat;
//   Matrix matrix_of_dm;
//   Matrix eigenvecs;
//
//   // Initialize
//   eigs.clear();
//   cout<<"Initialize..."<<endl;
//   for (int i=0;i<size;i++)
//   {
//        convert_num_to_vec(i,2,nsites,config);
//  	int n=(int) count(config.begin(),config.end(),1);
//
//	if (n==nfermions)
//	{
//   		maps[0].push_back(i);
//        	inverse_map[i]=maps[0].size()-1;
//  	} 
//   }
//   cout<<"Map size ="<<maps[0].size()<<endl;
//   cout<<"Lanczos for fermions ..."<<endl;
//   int nlow=3;
//   lanczos_given_map_evecs_for_fermions(h,iterations,maps[0],
//				inverse_map,eigs_0,evecs_0,nlow);
//   eigs.insert(eigs.end(),eigs_0.begin(),eigs_0.end());
//   
//   if (h.pairs1.size()>0)
//   {
//   	cout<<"Mathematica pairs1 for graph"<<endl;
//   	print_mathematica_pairs(h.pairs1);
//   }
//   if (h.pairs2.size()>0)
//   {
//   	cout<<"Mathematica pairs2 for graph"<<endl;
//   	print_mathematica_pairs(h.pairs2);
//   	cout<<endl;
//   }
//
//   cout<<"Energies for nfermions = "<<nfermions<<" particles"<<endl;
//   print_vec_acc(eigs,true);
//   //for (int i=0;i<eigs.size();i++) {cout<<eigs[i]<<"  "<<endl;}
//   for (int i=0;i<nsites;i++) {all_states.push_back(2);}
//   
//   kept_indices.push_back(0);
//   kept_indices.push_back(1);
//   kept_indices.push_back(2);
//   int nlocal=kept_indices.size();
//
//   int qddm_size=nlow*pow(2,nlocal);
//   int kept_size=pow(2,nlocal);
//   matrix_of_dm.resize(qddm_size,qddm_size);
//   cout <<"QDDM dimension = "<<qddm_size<<endl;
//
//   // Calculate <n_0>, <n_1>, <n_2> etc... 
//   double n0_n1=0.0;
//   double n0_n2=0.0;
//   double n1_n2=0.0;
//   double n_0=0.0;
//   double n_1=0.0;
//   double n_2=0.0;
//   for (int i=0;i<maps[0].size();i++)
//   {
//        convert_num_to_vec(maps[0][i],2,nsites,config);
//  	if (config[0]==1) n_0+=(evecs_0[0][i]*evecs_0[0][i]); 
//  	if (config[1]==1) n_1+=(evecs_0[0][i]*evecs_0[0][i]); 
//  	if (config[2]==1) n_2+=(evecs_0[0][i]*evecs_0[0][i]); 
//  	if (config[0]==1 and config[1]==1) n0_n1+=(evecs_0[0][i]*evecs_0[0][i]); 
//  	if (config[0]==1 and config[2]==1) n0_n2+=(evecs_0[0][i]*evecs_0[0][i]); 
//  	if (config[1]==1 and config[2]==1) n1_n2+=(evecs_0[0][i]*evecs_0[0][i]); 
//   }
//
//   cout<<"n0_n1,n0_n2,n1_n2 in GS ="<<n0_n1<<" "<<n0_n2<<" "<<n1_n2<<" "<<endl;
//   cout<<"n_0,n_1,n_2 in GS ="<<n_0<<" "<<n_1<<" "<<n_2<<" "<<endl;
//
//   for (int i=0;i<nlow;i++)
//   {
//	for (int j=0;j<nlow;j++)
//	{
//		cout<<"i,j = "<<i<<", "<<j<<endl;
//		if (i==j)
//		{
//   			density_matrix_calc_with_map(evecs_0[i], maps[0],inverse_map,
//					             all_states,kept_indices,density_mat);
//		}
//		else
//		{
//	   		off_diag_density_matrix_calc_with_map(evecs_0[i], 
//					maps[0],evecs_0[j],
//					maps[0],inverse_map,
//					all_states,kept_indices,density_mat);
//		}
//   		int shift_r=i*kept_size;int shift_c=j*kept_size;
//   		update_qddm(density_mat,shift_r,shift_c,matrix_of_dm);
//	}
//   }
//   cout<<"Matrix of DM = QDDM"<<endl;
//   print_real_mat(matrix_of_dm);
//   eigs.resize(qddm_size);
//   symmetric_diagonalize(matrix_of_dm, eigs, eigenvecs);
//   cout<<"nlow ="<<nlow<<endl;
//   //Matrix a(1,3);
//   //print_real_mat(a);
//   Matrix rho(nlow*nlow,kept_size*kept_size);
//   print_real_mat(rho);
//   
//   // Henley rho[mn,ab]
//   for (int m=0;m<nlow;m++)
//   {
//	for (int n=0;n<nlow;n++)
//	{
//		for (int a=0;a<kept_size;a++)
//		{
//			for (int b=0;b<kept_size;b++)
//		        {
//				rho((m*nlow)+n,(a*kept_size)+b)=matrix_of_dm( (m*kept_size)+a ,(n*kept_size)+b );
//			}		
//		}
//	}
//   }
//   
//   cout<<"RHO"<<endl;
//   print_real_mat(rho);
//   eigs.resize(kept_size*kept_size);
//   // Now take Henley's rho[mn,ab] and make the M Matrix out of it (see Henley, Changlani review paper)
//   Matrix m_matrix(kept_size*kept_size,kept_size*kept_size);
//   for (int i=0;i<kept_size*kept_size*kept_size*kept_size;i++) m_matrix[i]=0.0;
//   
//   for (int alpha=0;alpha<kept_size*kept_size;alpha++)
//   {
//   for (int beta=0;beta<kept_size*kept_size;beta++)
//   {
//   for (int m=0;m<nlow;m++)
//   {
//   for (int n=0;n<nlow;n++)
//   {
//	int mm=(m*nlow)+m;
//	int nn=(n*nlow)+n;
//	int mn=(m*nlow)+n;
//	//m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta)));
//	m_matrix(alpha,beta)+=((rho(mn,alpha)*rho(mn,beta))-((1.0/double(nlow))*rho(mm,alpha)*rho(nn,beta)));
//   }
//   }
//   }
//   }
//
//   Matrix u(kept_size*kept_size,kept_size*kept_size), vt(kept_size*kept_size,kept_size*kept_size);
//   cout<<"M matrix Rows ="<<m_matrix.NRows()<<endl;
//   cout<<"M matrix Cols ="<<m_matrix.NCols()<<endl;
//   cout<<"kept size     ="<<kept_size<<endl;
//   // Do SVD on M matrix
//   svd(m_matrix, eigs, u, vt);
//   cout<<"Left Eigenvectors of M matrix (u)  "<<endl;
//   print_real_mat(u);
//   cout<<"Right Eigenvectors of M matrix (vt)"<<endl;
//   print_real_mat2(vt);
//   cout<<"Eigenvalues of M matrix (eigs)"<<endl;
//   print_vec_acc(eigs,true);
//   
//   std::vector<double> dominant0,dominant1,dominant2;
//   for (int i=0;i<kept_size*kept_size;i++) dominant0.push_back(u(i,0));
//   for (int i=0;i<kept_size*kept_size;i++) dominant1.push_back(u(i,1));
//   for (int i=0;i<kept_size*kept_size;i++) dominant2.push_back(u(i,2));
//   Matrix dominant0_mat(kept_size,kept_size);
//   Matrix dominant1_mat(kept_size,kept_size);
//   Matrix dominant2_mat(kept_size,kept_size);
//   cout<<"Recasting dominant eigenvector as a matrix....."<<endl;
//   for (int i=0;i<kept_size*kept_size;i++) 
//   {
//        int m=i/kept_size;
//        int n=i-(m*kept_size);
//	cout<<m<<","<<n<<" "<<dominant0[i]<<endl;
//        dominant0_mat(m,n)=dominant0[i];
//        dominant1_mat(m,n)=dominant1[i];
//        dominant2_mat(m,n)=dominant2[i];
//   }
//   cout<<"Printing dominant 0 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant0_mat);
//   cout<<"Printing dominant 1 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant1_mat);
//   cout<<"Printing dominant 2 eigenvector as a matrix....."<<endl;
//   print_real_mat(dominant2_mat);
//}
////////////////////////////////////////////////////////////////////////////////
//
//void update_qddm(Matrix &density_mat,int shift_r,int shift_c, Matrix &matrix_of_dm)
//{
//	int i,j;
//	int rows=density_mat.NRows();
//        for (i=0;i<rows;i++)
//	{
//	     for (j=0;j<rows;j++)
//	     {
//		matrix_of_dm(i+shift_r,j+shift_c)=density_mat(i,j);
//		matrix_of_dm(j+shift_c,i+shift_r)=density_mat(i,j);
//	     }
//	}   
//}
//

