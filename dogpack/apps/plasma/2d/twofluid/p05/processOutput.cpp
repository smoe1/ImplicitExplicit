#include "../../lib/defs.h"
#include "../../lib/GetParams.h"
#include "Params.h"
#include <sstream>
#include <string>
#include <cmath>

double trunc(double in)
{
  return fabs(in)<1.0e-99 ? 0.0 : in;
}

// should add support for HDF5.
//
// Also, should have option to use
// H5F_ACC_EXCL rather than H5F_ACC_TRUNC when calling H5Fcreate 
// so that already-written files are not rewritten unless necessary.
//
void ProcessOutput(int argc,char**argv)
{
  string outputdir="output";
  void ParseArguments(int argc,char**argv,string& outputdir);
  ParseArguments(argc,argv,outputdir);
  void GetParams(string inputdir, string outputdir);
  string inputdir = outputdir;
  GetParams(inputdir,outputdir);
  const int      nout     = dogParams.get_nout();
  const double&  tfinal   = dogParams.get_tfinal();
  const int*     method   = dogParams.get_method();
  const int&     meqn     = dogParams.get_meqn();
  const int&     mdim     = dogParams.get_mdim();
  const int&     mx       = dogParams.get_mx();
  const int&     my       = dogParams.get_my();
  const int&     mbc      = dogParams.get_mbc();
  const double&  xlow     = dogParams.get_xlow();
  const double&  ylow     = dogParams.get_ylow();
  const double&  dx       = dogParams.get_dx();
  const double&  dy       = dogParams.get_dy();
  const int&     nstart   = dogParams.get_nstart();
  const int& method1 = method[1];
  const int& maux = method[6]; //iMax(method[6],1);
  int kmax = method[1]*(method[1]+1)/2;
  //
  void getAppParams(string appParams_file);
  getAppParams(inputdir+"/param.data");
  appParams.min_speed=0;
  appParams.max_speed=0;
  const double&  ion_mass = appParams.ion_mass;
  const double&  elc_mass = appParams.elc_mass;
  const double mass_diff = ion_mass-elc_mass;
  const double mass_prod = ion_mass*elc_mass;
  const double cc_clight2 = appParams.cc*appParams.cs_light*appParams.cs_light;
  //const double total_mass = ion_mass+elc_mass;
  //const double reduced_mass = ion_mass*elc_mas/total_mass;
  //
  dTensor3 node(mx+1,my+1,mdim);
  dTensor2 prim_vol(mx,my);
  dTensorBC4 aux(mx,my,maux,kmax,mbc);
  dTensorBC4 q(mx,my,meqn,kmax,mbc);
  dTensorBC4 Lstar(mx,my,meqn,kmax,mbc);
  dTensorBC4 Lstar_pressure(mx,my,meqn,kmax,mbc);
  dTensorBC4 Lstar_gas(mx,my,meqn,kmax,mbc);
  dTensorBC3 smax(mx,my,2,mbc);

  void InitGrid(int method1,double xlow,double ylow,
    double dx,double dy,dTensor3& node, dTensor2& prim_vol);
  InitGrid(method1,xlow,ylow,dx,dy,node,prim_vol);

  const string Xpoint_filename = outputdir+"/Xpoint.dat";
  fstream fs_Xpoint_file (Xpoint_filename.c_str(),ios::in);
  const bool write_Xpoint_file = true; //fs_Xpoint_file.is_open();
  fs_Xpoint_file.close();
  if(write_Xpoint_file) fs_Xpoint_file.open(Xpoint_filename.c_str(),ios::out);
  if(write_Xpoint_file) fs_Xpoint_file << setprecision(16);

  string rec_flx_filename = outputdir+"/rec_flx.dat";
  fstream fs_reconnected_flux_file (rec_flx_filename.c_str(),ios::in);
  const bool write_reconnected_flux_file = true; //!fs_reconnected_flux_file.is_open();
  fs_reconnected_flux_file.close();
  if(write_reconnected_flux_file) fs_reconnected_flux_file.open(rec_flx_filename.c_str(),ios::out);
  if(write_reconnected_flux_file) fs_reconnected_flux_file << setprecision(16);

  const double DT = tfinal/nout;
  //double old_t = 0.;
  //double t = 0.;
  //double dt_half_old = 0.;
  //double dt_half = 0.;
  double E3_old = 0.;
  double E3 = 0.;
  double E3_accum=0.;
  for (int n=nstart; n<=nout; n++)
  {
    cout << "reading in state values for frame " << n << endl;
    double QinitRestart(int nstart, string varname,
        dTensorBC4& q, string outputdir);
    double t = QinitRestart(n, "q", q, outputdir);

    void SampleState(const dTensorBC4& q, dTensor3& qvals,
        int method1, int kmax, int meqn);

    cout << "writing files that only use states" << endl;

    int x_elements=method1*mx;
    int y_elements=method1*my;
    //
    dTensor3 qvals(x_elements,y_elements,meqn);
    SampleState(q, qvals, method1, kmax, meqn);
    
    cout << "saving X-point data" << endl;

    // Get flux through left boundary
    double accum_left_flux = 0.0;
    for (int j=1; j<=my; j++)
    {             
	double Bx_left = q.get(1, j,11,1);
        accum_left_flux += Bx_left;
    }
    double left_flux_2 = accum_left_flux*dogParams.get_dy();

    const double& rho_i    = qvals.get(1,1,1);
    const double& rho_e    = qvals.get(1,1,6);
    const double rho = rho_i + rho_e;
    const double& M3_i     = qvals.get(1,1,4);
    const double& M3_e     = qvals.get(1,1,9);
    const double J3 = M3_i/ion_mass - M3_e/elc_mass;
    E3_old = E3;
    E3                     = qvals.get(1,1,16);
    //dt_half_old = dt_half;
    //dt_half = (t-old_t)/2.;
    cout << "DT=" << DT << endl;
    //E3_accum+=E3_old*DT/2.;
    E3_accum+=E3;//*DT/2.;
     // write values at X-point to a file
    if(write_Xpoint_file) fs_Xpoint_file 
       << setw(24) << scientific << t << " "
       << setw(24) << scientific << E3 << " "
       << setw(24) << scientific << E3_accum << " "
       << setw(24) << scientific << J3 << " "
       << setw(24) << scientific << rho << " "
       << setw(24) << scientific << J3/rho << " "
       << setw(24) << scientific << left_flux_2 << " "
       << endl;

    // compute reconnected flux
    double recon_flux;
    double rght_flux;
    double left_flux;
    double left_lost_flux;
    double total_x_axis_island_flux;
    void calculate_reconnected_flux(
        const dTensorBC4& q,
        double& recon_flux,
        double& rght_flux,
        double& left_flux,
        double& left_lost_flux,
        double& total_x_axis_island_flux);
    calculate_reconnected_flux(
        q,
        recon_flux,
        rght_flux,
        left_flux,
        left_lost_flux,
        total_x_axis_island_flux);
    // reconnected flux info
    //
    if(write_reconnected_flux_file) fs_reconnected_flux_file 
        << setw(24) << scientific << t << " "
        << setw(24) << scientific << recon_flux << " "
        << setw(24) << scientific << rght_flux << " "
        << setw(24) << scientific << left_flux << " "
        << setw(24) << scientific << left_lost_flux << " "
        << setw(24) << scientific << total_x_axis_island_flux << " " << endl;

    // print out state files for ready access
    //
    ostringstream basename;
    basename << setfill('0') << setw(4) << n;

    // write access files
    if(0)
    {
    //string rho_i_filename       = outputdir+"/Q01___."+basename.str()+".dat";
    //string Mi_filename          = outputdir+"/Q02:04."+basename.str()+".dat";
    //string energy_i_filename    = outputdir+"/Q05___."+basename.str()+".dat";
    //string rho_e_filename       = outputdir+"/Q06___."+basename.str()+".dat";
    //string Me_filename          = outputdir+"/Q07:09."+basename.str()+".dat";
    //string energy_e_filename    = outputdir+"/Q10___."+basename.str()+".dat";
    //string B_filename           = outputdir+"/Q11:13."+basename.str()+".dat";
    //string E_filename           = outputdir+"/Q14:16."+basename.str()+".dat";
    //string psi_filename         = outputdir+"/Q17___."+basename.str()+".dat";
    //string phi_filename         = outputdir+"/Q18___."+basename.str()+".dat";

    // check each of these files to see if it has already been created.
    // attempt to open each file for reading.
    //fstream fs_rho_i           (rho_i_filename      .c_str(),ios::in);
    //fstream fs_Mi              (Mi_filename         .c_str(),ios::in);
    //fstream fs_energy_i        (energy_i_filename   .c_str(),ios::in);
    //fstream fs_rho_e           (rho_e_filename      .c_str(),ios::in);
    //fstream fs_Me              (Me_filename         .c_str(),ios::in);
    //fstream fs_energy_e        (energy_e_filename   .c_str(),ios::in);
    //fstream fs_B               (B_filename          .c_str(),ios::in);
    //fstream fs_E               (E_filename          .c_str(),ios::in);
    //fstream fs_psi             (psi_filename        .c_str(),ios::in);
    //fstream fs_phi             (phi_filename        .c_str(),ios::in);
    //
    //const bool write_rho_i           = !fs_rho_i          .is_open();
    //const bool write_Mi              = !fs_Mi             .is_open();
    //const bool write_energy_i        = !fs_energy_i       .is_open();
    //const bool write_rho_e           = !fs_rho_e          .is_open();
    //const bool write_Me              = !fs_Me             .is_open();
    //const bool write_energy_e        = !fs_energy_e       .is_open();
    //const bool write_B               = !fs_B              .is_open();
    //const bool write_E               = !fs_E              .is_open();
    //const bool write_psi             = !fs_psi            .is_open();
    //const bool write_phi             = !fs_phi            .is_open();
    //
    //fs_rho_i          .close();
    //fs_Mi             .close();
    //fs_energy_i       .close();
    //fs_rho_e          .close();
    //fs_Me             .close();
    //fs_energy_e       .close();
    //fs_B              .close();
    //fs_E              .close();
    //fs_psi            .close();
    //fs_phi            .close();

    //if(write_rho_i   )  fs_rho_i   .open(rho_i_filename   .c_str(),ios::out);
    //if(write_Mi      )  fs_Mi      .open(Mi_filename      .c_str(),ios::out);
    //if(write_energy_i)  fs_energy_i.open(energy_i_filename.c_str(),ios::out);
    //if(write_rho_e   )  fs_rho_e   .open(rho_e_filename   .c_str(),ios::out);
    //if(write_Me      )  fs_Me      .open(Me_filename      .c_str(),ios::out);
    //if(write_energy_e)  fs_energy_e.open(energy_e_filename.c_str(),ios::out);
    //if(write_B       )  fs_B       .open(B_filename       .c_str(),ios::out);
    //if(write_E       )  fs_E       .open(E_filename       .c_str(),ios::out);
    //if(write_psi     )  fs_psi     .open(psi_filename     .c_str(),ios::out);
    //if(write_phi     )  fs_phi     .open(phi_filename     .c_str(),ios::out);

    //if(write_rho_i   ) fs_rho_i    << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_Mi      ) fs_Mi       << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_energy_i) fs_energy_i << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_rho_e   ) fs_rho_e    << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_Me      ) fs_Me       << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_energy_e) fs_energy_e << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_B       ) fs_B        << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_E       ) fs_E        << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_psi     ) fs_psi      << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    //if(write_phi     ) fs_phi      << setprecision(16)<<setw(24)<<scientific<<t<<endl;

    // for matlab's sake treat the first coordinate as most rapidly varying
    //
    //for(int j=1; j<=qvals.getsize(2); j++)
    //for(int i=1; i<=qvals.getsize(1); i++)
    //{
    //    const double& rho_i    = qvals.get(i,j,1);
    //    const double& M1_i     = qvals.get(i,j,2);
    //    const double& M2_i     = qvals.get(i,j,3);
    //    const double& M3_i     = qvals.get(i,j,4);
    //    const double& energy_i = qvals.get(i,j,5);
    //    const double& rho_e    = qvals.get(i,j,6);
    //    const double& M1_e     = qvals.get(i,j,7);
    //    const double& M2_e     = qvals.get(i,j,8);
    //    const double& M3_e     = qvals.get(i,j,9);
    //    const double& energy_e = qvals.get(i,j,10);
    //    const double& B1       = qvals.get(i,j,11);
    //    const double& B2       = qvals.get(i,j,12);
    //    const double& B3       = qvals.get(i,j,13);
    //    const double& E1       = qvals.get(i,j,14);
    //    const double& E2       = qvals.get(i,j,15);
    //    const double& E3       = qvals.get(i,j,16);
    //    const double& psi      = qvals.get(i,j,17);
    //    const double& phi      = qvals.get(i,j,18);

    //    // save state variables in individual files for rapid access
    //    //
    //    if(write_rho_i)    fs_rho_i    << setw(24) << scientific << trunc(rho_i   ) << endl;
    //    if(write_Mi)       fs_Mi       << setw(24) << scientific << trunc(M1_i    )
    //                            << " " << setw(24) << scientific << trunc(M2_i    )
    //                            << " " << setw(24) << scientific << trunc(M3_i    ) << endl;
    //    if(write_energy_i) fs_energy_i << setw(24) << scientific << trunc(energy_i) << endl;
    //    if(write_rho_e)    fs_rho_e    << setw(24) << scientific << trunc(rho_e   ) << endl;
    //    if(write_Me)       fs_Me       << setw(24) << scientific << trunc(M1_e    )
    //                            << " " << setw(24) << scientific << trunc(M2_e    )
    //                            << " " << setw(24) << scientific << trunc(M3_e    ) << endl;
    //    if(write_energy_e) fs_energy_e << setw(24) << scientific << trunc(energy_e) << endl;
    //    if(write_B)        fs_B        << setw(24) << scientific << trunc(B1      )
    //                            << " " << setw(24) << scientific << trunc(B2      )
    //                            << " " << setw(24) << scientific << trunc(B3      ) << endl;
    //    if(write_E)        fs_E        << setw(24) << scientific << trunc(E1      )
    //                            << " " << setw(24) << scientific << trunc(E2      )
    //                            << " " << setw(24) << scientific << trunc(E3      ) << endl;
    //    if(write_psi)      fs_psi      << setw(24) << scientific << trunc(psi     ) << endl;
    //    if(write_phi)      fs_phi      << setw(24) << scientific << trunc(phi     ) << endl;
    //}
    //if(write_rho_i)           fs_rho_i          .close();
    //if(write_Mi)              fs_Mi             .close();
    //if(write_energy_i)        fs_energy_i       .close();
    //if(write_rho_e)           fs_rho_e          .close();
    //if(write_Me)              fs_Me             .close();
    //if(write_energy_e)        fs_energy_e       .close();
    //if(write_B)               fs_B              .close();
    //if(write_E)               fs_E              .close();
    //if(write_psi)             fs_psi            .close();
    //if(write_phi)             fs_phi            .close();
    }

    //cout << "NOT calculating flux representation for frame " << n << endl;
    //continue;
    cout << "calculating flux representation for frame " << n << endl;
    // 
    // Extract terms by redefining the state and reconstructing
    // the flux that would obtain. (This is a hack; a more
    // efficient way to do this would be to separate out the
    // boundary integral and internal integral parts of the
    // Riemann solver and use that directly.)
    //
    // To get the divergence of the pressure, construct the
    // fluxes that would obtain if the velocity and EB field
    // were zero.
    //
    dTensorBC4 q_pressure(q); // should overload the assignment operator
    void Reset_q_pressure(dTensorBC4& Q);
    Reset_q_pressure(q_pressure);
    void ConstructL(const int method[],
        const dTensor3& node,
        const dTensor2& prim_vol,
        double t, double dt,
        dTensorBC4& aux, // SetBndValues modifies ghost cells
        dTensorBC4& q,   // SetBndValues modifies ghost cells
        dTensorBC4& Lstar,
        dTensorBC3& smax);
    double dt=0.0; // not actually used by ConstructL
    ConstructL(method,node,prim_vol,t,dt,aux,q_pressure,Lstar_pressure,smax);
    // construct the fluxes that would obtain if the electric
    // and magnetic field were zero
    dTensorBC4 q_gas(q);
    void Reset_q_gas(dTensorBC4& Q);
    Reset_q_gas(q_gas);
    ConstructL(method,node,prim_vol,t,dt,aux,q_gas,Lstar_gas,smax);
    // construct the actual fluxes
    // pretend that we were cleaning the E_field in order to generate
    // divergence constraint error data when we reconstruct the fluxes
    int tmp = appParams.clean_E_field;
    appParams.clean_E_field=1; // hack -- accessed in FluxFunc
    ConstructL(method,node,prim_vol,t,dt,aux,q,Lstar,smax);
    appParams.clean_E_field=tmp;
    
    cout << "evaluating fluxes at sample points" << endl;
    // evaluate the fluxes at the points of interest.
    //
    // Lstar represents the source term minus the divergence of
    // the flux. Evaluate this representation at the points of
    // interest.
    dTensor3 fvals(x_elements,y_elements,meqn);
    SampleState(Lstar, fvals, method1, kmax, meqn);
    //
    dTensor3 fvals_gas(x_elements,y_elements,meqn);
    SampleState(Lstar_gas, fvals_gas, method1, kmax, meqn);
    //
    dTensor3 fvals_pressure(x_elements,y_elements,meqn);
    SampleState(Lstar_pressure, fvals_pressure, method1, kmax, meqn);

    cout << "generating computed data for frame " << n << endl;

    // Compute and print out the terms in Ohm's law
    // in ascii-formatted files
    string Div_J_Flux_filename  = outputdir+"/DFJ"+basename.str()+".dat";
    string E_chk_filename       = outputdir+"/Eck"+basename.str()+".dat";
    string B_cross_u_filename   = outputdir+"/Bxu"+basename.str()+".dat";
    string Hall_filename        = outputdir+"/Hal"+basename.str()+".dat";
    string p_ek_filename        = outputdir+"/pek"+basename.str()+".dat";
    string J_t_filename         = outputdir+"/J_t"+basename.str()+".dat";
    string dtJ_filename         = outputdir+"/dtJ"+basename.str()+".dat";
    string E_t_filename         = outputdir+"/E_t"+basename.str()+".dat";
    string sigma_filename       = outputdir+"/sgm"+basename.str()+".dat";
    string DivB_filename        = outputdir+"/dvB"+basename.str()+".dat";
    string DivE_filename        = outputdir+"/dvE"+basename.str()+".dat";
    // check each of these files to see if it has already been created.
    // attempt to open each file for reading.
    fstream fs_B_cross_u_file  (B_cross_u_filename  .c_str(),ios::in);
    fstream fs_Hall_file       (Hall_filename       .c_str(),ios::in);
    fstream fs_p_ek_file       (p_ek_filename       .c_str(),ios::in);
    fstream fs_J_t_file        (J_t_filename        .c_str(),ios::in);
    fstream fs_dtJ_file        (dtJ_filename        .c_str(),ios::in);
    fstream fs_Div_J_Flux_file (Div_J_Flux_filename .c_str(),ios::in);
    fstream fs_E_chk_file      (E_chk_filename      .c_str(),ios::in);
    fstream fs_E_t_file        (E_t_filename        .c_str(),ios::in);
    fstream fs_sigma_file      (sigma_filename      .c_str(),ios::in);
    fstream fs_DivB_file       (DivB_filename       .c_str(),ios::in);
    fstream fs_DivE_file       (DivE_filename       .c_str(),ios::in);

    const bool write_B_cross_u_file  = !fs_B_cross_u_file .is_open();
    const bool write_Hall_file       = !fs_Hall_file      .is_open();
    const bool write_p_ek_file       = !fs_p_ek_file      .is_open();
    const bool write_J_t_file        = !fs_J_t_file       .is_open();
    const bool write_dtJ_file        = !fs_dtJ_file       .is_open();
    const bool write_Div_J_Flux_file = true; !fs_Div_J_Flux_file.is_open();
    const bool write_E_chk_file      = true; !fs_E_chk_file     .is_open();
    const bool write_E_t_file        = !fs_E_t_file       .is_open();
    const bool write_sigma_file      = !fs_sigma_file     .is_open();
    const bool write_DivB_file       = !fs_DivB_file      .is_open();
    const bool write_DivE_file       = !fs_DivE_file      .is_open();

    fs_B_cross_u_file .close();
    fs_Hall_file      .close();
    fs_p_ek_file      .close();
    fs_J_t_file       .close();
    fs_Div_J_Flux_file.close();
    fs_E_chk_file     .close();
    fs_E_t_file       .close();
    fs_sigma_file     .close();
    fs_DivB_file      .close();
    fs_DivE_file      .close();

    if(write_B_cross_u_file )  fs_B_cross_u_file  .open(B_cross_u_filename  .c_str(),ios::out);
    if(write_Hall_file      )  fs_Hall_file       .open(Hall_filename       .c_str(),ios::out);
    if(write_p_ek_file      )  fs_p_ek_file       .open(p_ek_filename       .c_str(),ios::out);
    if(write_J_t_file       )  fs_J_t_file        .open(J_t_filename        .c_str(),ios::out);
    if(write_dtJ_file       )  fs_dtJ_file        .open(dtJ_filename        .c_str(),ios::out);
    if(write_Div_J_Flux_file)  fs_Div_J_Flux_file .open(Div_J_Flux_filename .c_str(),ios::out);
    if(write_E_chk_file     )  fs_E_chk_file      .open(E_chk_filename      .c_str(),ios::out);
    if(write_E_t_file       )  fs_E_t_file        .open(E_t_filename        .c_str(),ios::out);
    if(write_sigma_file     )  fs_sigma_file      .open(sigma_filename      .c_str(),ios::out);
    if(write_DivB_file      )  fs_DivB_file       .open(DivB_filename       .c_str(),ios::out);
    if(write_DivE_file      )  fs_DivE_file       .open(DivE_filename       .c_str(),ios::out);

    if(write_B_cross_u_file ) fs_B_cross_u_file  << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_Hall_file      ) fs_Hall_file       << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_p_ek_file      ) fs_p_ek_file       << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_J_t_file       ) fs_J_t_file        << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_dtJ_file       ) fs_dtJ_file        << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_Div_J_Flux_file) fs_Div_J_Flux_file << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_E_chk_file     ) fs_E_chk_file      << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_E_t_file       ) fs_E_t_file     << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_sigma_file     ) fs_sigma_file   << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_DivB_file      ) fs_DivB_file    << setprecision(16)<<setw(24)<<scientific<<t<<endl;
    if(write_DivE_file      ) fs_DivE_file    << setprecision(16)<<setw(24)<<scientific<<t<<endl;

    // for matlab's sake treat the first coordinate as most rapidly varying
    //
    for(int j=1; j<=qvals.getsize(2); j++)
    for(int i=1; i<=qvals.getsize(1); i++)
    {
        // Variables
        const double& rho_i    = qvals.get(i,j,1);
        const double& M1_i     = qvals.get(i,j,2);
        const double& M2_i     = qvals.get(i,j,3);
        const double& M3_i     = qvals.get(i,j,4);
        const double& energy_i = qvals.get(i,j,5);
        const double& rho_e    = qvals.get(i,j,6);
        const double& M1_e     = qvals.get(i,j,7);
        const double& M2_e     = qvals.get(i,j,8);
        const double& M3_e     = qvals.get(i,j,9);
        const double& energy_e = qvals.get(i,j,10);
        const double& B1       = qvals.get(i,j,11);
        const double& B2       = qvals.get(i,j,12);
        const double& B3       = qvals.get(i,j,13);
        const double& E1       = qvals.get(i,j,14);
        const double& E2       = qvals.get(i,j,15);
        const double& E3       = qvals.get(i,j,16);
        const double& psi      = qvals.get(i,j,17);
        const double& phi      = qvals.get(i,j,18);

        const double& pi_x = -fvals_pressure.get(i,j,2);
        const double& pi_y = -fvals_pressure.get(i,j,3);
        const double& pe_x = -fvals_pressure.get(i,j,7);
        const double& pe_y = -fvals_pressure.get(i,j,8);

        const double& M1i_t = fvals.get(i,j,2);
        const double& M2i_t = fvals.get(i,j,3);
        const double& M3i_t = fvals.get(i,j,4);

        const double& M1e_t = fvals.get(i,j,7);
        const double& M2e_t = fvals.get(i,j,8);
        const double& M3e_t = fvals.get(i,j,9);

        const double& B1_t = fvals.get(i,j,11);
        const double& B2_t = fvals.get(i,j,12);
        const double& B3_t = fvals.get(i,j,13);

        //const double& E1_t = fvals.get(i,j,14);
        //const double& E2_t = fvals.get(i,j,15);
        //const double& E3_t = fvals.get(i,j,16);

        const double& psi_t = fvals.get(i,j,17);
        // depends on appParams.clean_E_field=1 above
        const double& phi_t = fvals.get(i,j,18);

        const double DivMomFlux_i_1 = -fvals_gas.get(i,j,2);
        const double DivMomFlux_i_2 = -fvals_gas.get(i,j,3);
        const double DivMomFlux_i_3 = -fvals_gas.get(i,j,4);

        const double DivMomFlux_e_1 = -fvals_gas.get(i,j,7);
        const double DivMomFlux_e_2 = -fvals_gas.get(i,j,8);
        const double DivMomFlux_e_3 = -fvals_gas.get(i,j,9);
        
        // compute derived quantities

        const double rho = rho_i + rho_e;
        const double M1 = M1_i + M1_e;
        const double M2 = M2_i + M2_e;
        const double M3 = M3_i + M3_e;
        const double u1 = M1/rho;
        const double u2 = M2/rho;
        const double u3 = M3/rho;
        // J
        const double J1 = M1_i/ion_mass - M1_e/elc_mass;
        const double J2 = M2_i/ion_mass - M2_e/elc_mass;
        const double J3 = M3_i/ion_mass - M3_e/elc_mass;
        // J_t
        const double J1_t = M1i_t/ion_mass - M1e_t/elc_mass;
        const double J2_t = M2i_t/ion_mass - M2e_t/elc_mass;
        const double J3_t = M3i_t/ion_mass - M3e_t/elc_mass;
        // DivCurrFlux
        const double DivCurrFlux_1=DivMomFlux_i_1/ion_mass - DivMomFlux_e_1/elc_mass;
        const double DivCurrFlux_2=DivMomFlux_i_2/ion_mass - DivMomFlux_e_2/elc_mass;
        const double DivCurrFlux_3=DivMomFlux_i_3/ion_mass - DivMomFlux_e_3/elc_mass;
        // J_cross_B
        const double J_cross_B_1 = J2*B3-J3*B2;
        const double J_cross_B_2 = J3*B1-J1*B3;
        const double J_cross_B_3 = J1*B2-J2*B1;

        // gradient of electrokinetic pressure
        const double grad_p_ek_x = pi_x/ion_mass-pe_x/elc_mass;
        const double grad_p_ek_y = pi_y/ion_mass-pe_y/elc_mass;

        const double mass_diff_ovr_rho = mass_diff/rho;
        const double mass_prod_ovr_rho = mass_prod/rho;

        // compute Ohm's law terms

        // B_cross_u
        const double B_cross_u_1 = B2*u3-B3*u2;
        const double B_cross_u_2 = B3*u1-B1*u3;
        const double B_cross_u_3 = B1*u2-B2*u1;

        // Hall
        const double Hall_1 = J_cross_B_1*mass_diff_ovr_rho;
        const double Hall_2 = J_cross_B_2*mass_diff_ovr_rho;
        const double Hall_3 = J_cross_B_3*mass_diff_ovr_rho;

        // electrokinetic pressure term
        const double p_ek_1 = grad_p_ek_x*mass_prod_ovr_rho;
        const double p_ek_2 = grad_p_ek_y*mass_prod_ovr_rho;

        // electrokinetic inertial terms
        //
        const double J_t_term_1 = J1_t*mass_prod_ovr_rho;
        const double J_t_term_2 = J2_t*mass_prod_ovr_rho;
        const double J_t_term_3 = J3_t*mass_prod_ovr_rho;
        //
        const double DivCurrFluxTerm_1 = DivCurrFlux_1*mass_prod_ovr_rho;
        const double DivCurrFluxTerm_2 = DivCurrFlux_2*mass_prod_ovr_rho;
        const double DivCurrFluxTerm_3 = DivCurrFlux_3*mass_prod_ovr_rho;
        //
        const double dtJ_term_1 = J_t_term_1 + DivCurrFluxTerm_1;
        const double dtJ_term_2 = J_t_term_2 + DivCurrFluxTerm_2;
        const double dtJ_term_3 = J_t_term_3 + DivCurrFluxTerm_3;

        // check assumptions of Ohm's law

        // Is displacement current negligible?
        //
        const double& E1_t = fvals.get(i,j,14);
        const double& E2_t = fvals.get(i,j,15);
        const double& E3_t = fvals.get(i,j,16);

        // Does quasineutrality hold?
        const double charge_dens = rho_i/ion_mass - rho_e/elc_mass;

        // Is Ohm's law true?
        const double E1_chk=B_cross_u_1+Hall_1+p_ek_1+J_t_term_1+DivCurrFluxTerm_1;
        const double E2_chk=B_cross_u_2+Hall_2+p_ek_2+J_t_term_2+DivCurrFluxTerm_2;
        const double E3_chk=B_cross_u_3+Hall_3       +J_t_term_3+DivCurrFluxTerm_3;

        // check divergence constraints

        const double DivB = psi_t/cc_clight2;
        const double DivE_mns__qdens_ovr_eps = phi_t/appParams.cc;

        // write results to files in ascii format

        //
        // Ohm's law terms
        //
        if(write_B_cross_u_file) fs_B_cross_u_file
          << setw(24) << scientific << trunc(B_cross_u_1) << " "
          << setw(24) << scientific << trunc(B_cross_u_2) << " "
          << setw(24) << scientific << trunc(B_cross_u_3) << " "
          << endl;
        if(write_Hall_file) fs_Hall_file
          << setw(24) << scientific << trunc(Hall_1) << " "
          << setw(24) << scientific << trunc(Hall_2) << " "
          << setw(24) << scientific << trunc(Hall_3) << " "
          << endl;
        if(write_p_ek_file) fs_p_ek_file
          << setw(24) << scientific << trunc(p_ek_1) << " "
          << setw(24) << scientific << trunc(p_ek_2) << " "
          << setw(24) << scientific << 0e0           << " "
          << endl;
        if(write_dtJ_file) fs_dtJ_file
          << setw(24) << scientific << trunc(dtJ_term_1) << " "
          << setw(24) << scientific << trunc(dtJ_term_2) << " "
          << setw(24) << scientific << trunc(dtJ_term_3) << " "
          << endl;
        if(write_J_t_file) fs_J_t_file
          << setw(24) << scientific << trunc(J_t_term_1) << " "
          << setw(24) << scientific << trunc(J_t_term_2) << " "
          << setw(24) << scientific << trunc(J_t_term_3) << " "
          << endl;
        if(write_Div_J_Flux_file) fs_Div_J_Flux_file 
          << setw(24) << scientific << trunc(DivCurrFluxTerm_1) << " "
          << setw(24) << scientific << trunc(DivCurrFluxTerm_2) << " "
          << setw(24) << scientific << trunc(DivCurrFluxTerm_3) << " "
          << endl;
        if(write_E_chk_file) fs_E_chk_file 
          << setw(24) << scientific << trunc(E1_chk) << " "
          << setw(24) << scientific << trunc(E2_chk) << " "
          << setw(24) << scientific << trunc(E3_chk) << " "
          << endl;
        //
        // Checks (things that ideally are zero)
        //
        // assumptions of Ohm's law
        if(write_E_t_file) fs_E_t_file
          << setw(24) << scientific << trunc(E1_t) << " "
          << setw(24) << scientific << trunc(E2_t) << " "
          << setw(24) << scientific << trunc(E3_t) << " "
          << endl;
        if(write_sigma_file) fs_sigma_file 
          << setw(24) << scientific << charge_dens << endl;
        //
        // divergence constraint error
        //
        if(write_DivB_file) fs_DivB_file
          << setw(24) << scientific << trunc(DivB) << endl;
        if(write_DivE_file) fs_DivE_file
          << setw(24) << scientific << trunc(DivE_mns__qdens_ovr_eps) << endl;
        //
        // write results to files in HDF5 format
    }
    if(write_Div_J_Flux_file) fs_Div_J_Flux_file.close();
    if(write_E_chk_file)      fs_E_chk_file     .close();
    if(write_B_cross_u_file)  fs_B_cross_u_file .close();
    if(write_Hall_file)       fs_Hall_file      .close();
    if(write_p_ek_file)       fs_p_ek_file      .close();
    if(write_J_t_file)        fs_J_t_file       .close();
    if(write_dtJ_file)        fs_dtJ_file       .close();
    if(write_E_t_file)        fs_E_t_file       .close();
    if(write_sigma_file)      fs_sigma_file     .close();
    cout << "finished generating Ohm's law data for frame " << n << endl;
  }
  if(write_reconnected_flux_file) fs_reconnected_flux_file.close();
}

int main(int argc, char* argv[])
{
  ProcessOutput(argc,argv);
  return 0;
}

void Reset_q_pressure(dTensorBC4& Q)
{
    for (int i=1; i<=Q.getsize(1); i++)
    for (int j=1; j<=Q.getsize(2); j++)
    for (int k=1; k<=Q.getsize(4); k++)
    {
        // Variables
        double rho_i    = Q.get(i,j,1, k);
        double M1_i     = Q.get(i,j,2, k);
        double M2_i     = Q.get(i,j,3, k);
        double M3_i     = Q.get(i,j,4, k);
        double energy_i = Q.get(i,j,5, k);
        
        double rho_e    = Q.get(i,j,6, k);
        double M1_e     = Q.get(i,j,7, k);
        double M2_e     = Q.get(i,j,8, k);
        double M3_e     = Q.get(i,j,9, k);
        double energy_e = Q.get(i,j,10,k);
    
        double u1_i     = M1_i/rho_i;
        double u2_i     = M2_i/rho_i;
        double u3_i     = M3_i/rho_i;
    
        double u1_e     = M1_e/rho_e;
        double u2_e     = M2_e/rho_e;
        double u3_e     = M3_e/rho_e;
    
        double thermal_energy_i = energy_i
            - 0.5e0*(u1_i*M1_i + u2_i*M2_i + u3_i*M3_i);
        double thermal_energy_e  = energy_e
            - 0.5e0*(u1_e*M1_e + u2_e*M2_e + u3_e*M3_e);
        
        Q.set(i,j,1, k,  rho_i);
        Q.set(i,j,2, k,  0.0);
        Q.set(i,j,3, k,  0.0);
        Q.set(i,j,4, k,  0.0);
        Q.set(i,j,5, k,  thermal_energy_i);

        Q.set(i,j,6, k, rho_e);
        Q.set(i,j,7, k, 0.0);
        Q.set(i,j,8, k, 0.0);
        Q.set(i,j,9, k, 0.0);
        Q.set(i,j,10,k, thermal_energy_e);

        Q.set(i,j,11,k, 0.0);
        Q.set(i,j,12,k, 0.0);
        Q.set(i,j,13,k, 0.0);
        Q.set(i,j,14,k, 0.0);
        Q.set(i,j,15,k, 0.0);
        Q.set(i,j,16,k, 0.0);
        Q.set(i,j,17,k, 0.0);
        Q.set(i,j,18,k, 0.0);
    }
}

void Reset_q_gas(dTensorBC4& Q)
{
    for (int i=1; i<=Q.getsize(1); i++)
    for (int j=1; j<=Q.getsize(2); j++)
    for (int k=1; k<=Q.getsize(4); k++)
    {
        Q.set(i,j,11,k, 0.0);
        Q.set(i,j,12,k, 0.0);
        Q.set(i,j,13,k, 0.0);
        Q.set(i,j,14,k, 0.0);
        Q.set(i,j,15,k, 0.0);
        Q.set(i,j,16,k, 0.0);
        Q.set(i,j,17,k, 0.0);
        Q.set(i,j,18,k, 0.0);
    }
}
