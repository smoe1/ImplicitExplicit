#include <cmath>
#include <string>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Components.h"
#include "GEM.h"

// within this element compute the basis coefficients of the
// divergence of the magnetic field (the decrease in the number
// of coefficients reflects the decrease in the order of accuracy
// when taking the derivative).
static void compute_coefs_of_divergence_in_element(
  double* divB,
  const double* B1,
  const double* B2,
  double dx,
  double dy,
  int method_order)
{
  // not optimally efficient, but called infrequently
  switch(method_order)
  {
   case 5:
    divB[10] = 2.0*sq3*B1[12]/dx + 6.0*B2[15]*sq7/dy;            
    divB[9]  = 6.0*B1[14]*sq7/dx + 2.0*sq3*B2[11]/dy;
    divB[8]  = 2.0*sq5*B1[13]*sq3/dx + 2.0*sq5*sq7*B2[12]/dy;
    divB[7]  = 2.0*sq5*sq7*B1[11]/dx + 2.0*sq5*B2[13]*sq3/dy;
   case 4:
    divB[6]  = 2.0*sq3*B1[8]/dx + 2.0*sq5*B2[10]*sq7/dy;
    divB[5]  = 2.0*sq5*B1[9]*sq7/dx + 2.0*B2[7]*sq3/dy;
    divB[4]  = 2.0*B1[7]*sq3*sq5/dx + 2.0*sq3*B2[8]*sq5/dy;
   case 3:
    divB[3]  = (2.0*(B1[11]*sq7+sq3*B1[4]))/dx 
              +2.0*sq3*(3*B2[15]+B2[6]*sq5)/dy;
    divB[2]  = 2.0*sq3*(B1[5]*sq5+3*B1[14])/dx 
              +(2.0*(B2[12]*sq7+sq3*B2[4]))/dy;
   case 2:
    divB[1]  = (2.0*(B1[2]*sq3+B1[9]*sq7))/dx 
              +(2.0*(B2[3]*sq3+B2[10]*sq7))/dy;
    break;
   case 1:
      divB[0] = 0.;
      break;
      //eprintf("cannot compute cell divergence for first order");
   default:
      eprintf("unsupported method_order=%d\n",method_order);
  }
}

static void output_local_divergence_error(
    int method_order,
    double dx,
    double dy,
    int n_B1,
    int n_B2,
    const dTensorBC4& q,
    double t)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);

  // Approach #1 (compute L2 norm of global solution using finite difference)
  //
  //#pragma omp parallel for
  assert(method_order<=3);
  double integral_divB_squared=0.0;
  for (int i=1; i<=mx; i++)   
  {
    double inner_integral_divB_squared = 0.;
    for (int j=1; j<=my; j++)
    {
      // compute contribution of element to L2-norm squared of div(B)
      double divB_coef[7];
      for (int k=1; k<=kmax; k++){
        const double dB1 = (q.get(i+1,j,n_B1,k) - q.get(i-1,j,n_B1,k))/(2*dx);
        const double dB2 = (q.get(i,j+1,n_B2,k) - q.get(i,j-1,n_B2,k))/(2*dy);
        divB_coef[k] = dB1 + dB2;
      }
      // corrections (more would be necessary for higher order)
      if(method_order==3)
      {
        divB_coef[1] -= (2*sq5)*divB_coef[5];
      }
      double tmp = 0.;
      for (int k=1; k<=kmax; k++) tmp += pow(divB_coef[k],2);
      inner_integral_divB_squared += tmp;
    }
    integral_divB_squared += inner_integral_divB_squared;
  }
  const double normDivB = sqrt(integral_divB_squared*dx*dy);
    
  // Approach #2: sum L2 norms of solution in interior of each element
  // (exact L2 norm over *interior* of each element)
  double sum_elementDivB_sq=0.0;
  for (int i=1; i<=mx; i++)   
  {
    double B1[16], B2[16]; for(int k=1; k<=15; k++) B1[k]=0,B2[k]=0;
    double divB[11]; for (int k=1; k<=10; k++) divB[k]=0;
    double inner_sum_elementDivB_sq = 0.;
    for (int j=1; j<=my; j++)
    {
      for (int k=1; k<=kmax; k++){
        B1[k]= q.get(i,j,n_B1,k);
        B2[k]= q.get(i,j,n_B2,k);
      }
      compute_coefs_of_divergence_in_element(divB,B1,B2,dx,dy,method_order);
      
      double tmp = 0.;
      for (int k=1; k<=kmax; k++) tmp += pow(divB[k],2);
      inner_sum_elementDivB_sq += tmp;
    }
    sum_elementDivB_sq += inner_sum_elementDivB_sq;
  }
  const double norm_interiorDivB = sqrt(sum_elementDivB_sq*dx*dy);
    
  // write norms to file.
  double matlab_fix_double(double val);
  const char* get_outputdir();
  string fname = string(get_outputdir())+"/divB.dat";
  FILE* file = fopen(fname.c_str(),(t==0 ? "w" : "a"));
  // I see no reason to report double precision.
  fprintf(file, "%e %e %e\n", t,
    matlab_fix_double(normDivB),
    matlab_fix_double(norm_interiorDivB));
  fclose(file);
}

void ConSoln(
    const dTensorBC4& aux,
    const dTensorBC4& q,
    double t)
{
    void OutputMagneticFluxes(
        int n_B1,
        int n_B2,
        const dTensorBC4& q,
        double t);
    OutputMagneticFluxes(
        _B1, _B2, q, t);

    assert(dogParams.get_mcapa()<1); // without capacity function
    void WriteConservation(
        const dTensorBC4& q,
        double t);
    WriteConservation(q, t);

    void Output_center_E3(
      const dTensorBC4& q, double t, int _E3);
    Output_center_E3(q, t, _E3);
    void Output_center_gas_states(
      const dTensorBC4& q, double t, const char* outputfile,
      int species_offset, bool fiveMoment=false);
    const char* get_outputdir();
    string outputfile = string(get_outputdir())+"/xpoint_gas_i.dat";
    Output_center_gas_states(q, t,outputfile.c_str(),0,false);
    void Output_center_cell_state(const dTensorBC4& q, double t);
    Output_center_cell_state(q, t);
    //
    // output the rate of production of entropy at the center
    //
    const int maux = aux.getsize(3);
    if(maux>0)
    {
      assert_eq(maux,1);
      void output_component_at_center(
        const dTensorBC4& q, double t, const char* outputfile,
        int idx);
      outputfile = string(get_outputdir())+"/xpoint_ent_prod_i.dat";
      output_component_at_center(aux, t, outputfile.c_str(),1);
    }

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    output_local_divergence_error(dogParams.get_space_order(), dx, dy,
      _B1, _B2, q, t);
}
