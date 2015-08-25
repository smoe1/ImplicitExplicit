#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "GEMparams.h"
#include "gas05.h"
#include "gas10.h"
#include <string>

using namespace std;

static double cell_value(
  const dTensorBC4& q,
  int i, int j, /* coordinates of cell */
  // expected to be between -1 and 1 (could change type to double)
  int xpos, int ypos, /* position in cell */
  int stateIdx)
{
  double ret=0.;
  switch(dogParams.get_space_order())
  {
    default:
      unsupported_value_error(dogParams.get_space_order());
    case 3: // third-order
      ret += (sq5*0.5)*(3*xpos*xpos-1)*q.get(i,j,stateIdx,5);
      ret += (sq5*0.5)*(3*ypos*ypos-1)*q.get(i,j,stateIdx,6);
      ret += 3*xpos*ypos*q.get(i,j,stateIdx,4);
    case 2: // second-order
      if(xpos) ret += sq3*xpos*q.get(i,j,stateIdx,2);
      if(ypos) ret += sq3*ypos*q.get(i,j,stateIdx,3);
    case 1: // first-order
      ret += q.get(i,j,stateIdx,1);
  }
  return ret;
}

static void get_center_cell_coords(int& xidx, int& yidx)
{
  if(gemParams.x_symmetry_not_enforced())
    xidx = dogParamsCart2.get_mx()/2+1;
  else
    xidx=1;
  if(gemParams.y_symmetry_not_enforced())
    yidx = dogParamsCart2.get_my()/2+1;
  else
    yidx=1;
}

// take the origin value to be the average of the value
// according to all the neighboring cells
static double origin_value(
  int idx,
  const dTensorBC4& q)
{
  double value;
  int i,j; get_center_cell_coords(i,j);
  switch(gemParams.get_enforced_symmetry()&XY)
  {
    case NONE:
      value= 0.25*(cell_value(q,i,j,-1,-1,idx)
                  +cell_value(q,i-1,j-1,1,1,idx)
                  +cell_value(q,i,j-1,-1,1,idx)
                  +cell_value(q,i-1,j,1,-1,idx));
      break;
    case X:
      assert_eq(i,1);
      value= 0.5*(cell_value(q,i,j,-1,-1,idx)
                 +cell_value(q,i,j-1,-1,1,idx));
      break;
    case Y:
      assert_eq(j,1);
      value= 0.5*(cell_value(q,i,j,-1,-1,idx)
                 +cell_value(q,i-1,j,1,-1,idx));
      break;
    case XY:
      assert_eq(i,1); assert_eq(j,1);
      value= cell_value(q,i,j,-1,-1,idx);
      break;
    default:
      invalid_value_error(gemParams.get_enforced_symmetry());
  }
  return value;
}

static double evaluate_at_origin(
  double (*evaluate)(const dTensorBC4& q,
    int i, int j, int xpos, int ypos, int species_offset),
  const dTensorBC4& q, int species_offset)
{
  double value;
  int i,j; get_center_cell_coords(i,j);
  switch(gemParams.get_enforced_symmetry()&XY)
  {
    case NONE:
      value= 0.25*(evaluate(q,i,j,-1,-1,species_offset)
                  +evaluate(q,i-1,j-1,1,1,species_offset)
                  +evaluate(q,i,j-1,-1,1,species_offset)
                  +evaluate(q,i-1,j,1,-1,species_offset));
      break;
    case X:
      assert_eq(i,1);
      value= 0.5*(evaluate(q,i,j,-1,-1,species_offset)
                 +evaluate(q,i,j-1,-1,1,species_offset));
      break;
    case Y:
      assert_eq(j,1);
      value= 0.5*(evaluate(q,i,j,-1,-1,species_offset)
                 +evaluate(q,i-1,j,1,-1,species_offset));
      break;
    case XY:
      assert_eq(i,1); assert_eq(j,1);
      value= evaluate(q,i,j,-1,-1,species_offset);
      break;
    default:
      invalid_value_error(gemParams.get_enforced_symmetry());
  }
  return value;
}

static double evaluate_entropy_per_mass05(
  const dTensorBC4& q,
  int i, int j, /* coordinates of cell */
  // expected to be between -1 and 1 (could change type to double)
  int xpos, int ypos, /* position in cell */
  int species_offset)
{
  using namespace FiveMomentComponentID;
  const double rho = cell_value(q,i,j,xpos,ypos,species_offset+_rho);
  const double M1  = cell_value(q,i,j,xpos,ypos,species_offset+_M1 );
  const double M2  = cell_value(q,i,j,xpos,ypos,species_offset+_M2 );
  const double M3  = cell_value(q,i,j,xpos,ypos,species_offset+_M3 );
  const double nrg = cell_value(q,i,j,xpos,ypos,species_offset+_N  );
  return get_entropy_per_mass05(rho,M1,M2,M3,nrg);
}

static double evaluate_entropy_per_mass10(
  const dTensorBC4& q,
  int i, int j, /* coordinates of cell */
  // expected to be between -1 and 1 (could change type to double)
  int xpos, int ypos, /* position in cell */
  int species_offset)
{
  using namespace TenMomentComponentID;
  const double rho = cell_value(q,i,j,xpos,ypos,species_offset+_rho);
  const double M1  = cell_value(q,i,j,xpos,ypos,species_offset+_M1 );
  const double M2  = cell_value(q,i,j,xpos,ypos,species_offset+_M2 );
  const double M3  = cell_value(q,i,j,xpos,ypos,species_offset+_M3 );
  const double N11 = cell_value(q,i,j,xpos,ypos,species_offset+_N11);
  const double N12 = cell_value(q,i,j,xpos,ypos,species_offset+_N12);
  const double N13 = cell_value(q,i,j,xpos,ypos,species_offset+_N13);
  const double N22 = cell_value(q,i,j,xpos,ypos,species_offset+_N22);
  const double N23 = cell_value(q,i,j,xpos,ypos,species_offset+_N23);
  const double N33 = cell_value(q,i,j,xpos,ypos,species_offset+_N33);
  return get_entropy_per_mass10(rho,M1,M2,M3,N11,N12,N13,N22,N23,N33);
}

static double evaluate_P13(
  const dTensorBC4& q,
  int i, int j, /* coordinates of cell */
  // expected to be between -1 and 1 (could change type to double)
  int xpos, int ypos, /* position in cell */
  int species_offset)
{
  using namespace TenMomentComponentID;
  const double rho = cell_value(q,i,j,xpos,ypos,species_offset+_rho);
  const double M1  = cell_value(q,i,j,xpos,ypos,species_offset+_M1 );
  const double M3  = cell_value(q,i,j,xpos,ypos,species_offset+_M3 );
  const double N13 = cell_value(q,i,j,xpos,ypos,species_offset+_N13);
  const double P13 = N13 - M1*M3/rho;
  
  return P13;
}

static double evaluate_P23(
  const dTensorBC4& q,
  int i, int j, /* coordinates of cell */
  // expected to be between -1 and 1 (could change type to double)
  int xpos, int ypos, /* position in cell */
  int species_offset)
{
  using namespace TenMomentComponentID;
  const int _rho_s = species_offset+_rho;
  const int _M2_s  = species_offset+_M2;
  const int _M3_s  = species_offset+_M3;
  const int _N23_s = species_offset+_N23;

  const double rho = cell_value(q,i,j,xpos,ypos,species_offset+_rho);
  const double M2  = cell_value(q,i,j,xpos,ypos,species_offset+_M2 );
  const double M3  = cell_value(q,i,j,xpos,ypos,species_offset+_M3 );
  const double N23 = cell_value(q,i,j,xpos,ypos,species_offset+_N23);
  const double P23 = N23 - M2*M3/rho;
  
  return P23;
}

// This assumes that the quantity that the "evaluate" function
// is antisymmetric in the direction of the x-axis.
// This is a second-order-accurate estimate, which is the
// best we can expect for a third-order method.
static double x_derivative_at_origin(
  double (*evaluate)(const dTensorBC4& q,
    int i, int j, int xpos, int ypos, int species_offset),
  const dTensorBC4& q, int species_offset)
{
  int i,j; get_center_cell_coords(i,j);
  double dx_val;
  switch(gemParams.get_enforced_symmetry()&XY)
  {
    case NONE:
      dx_val = evaluate(q,i,j,0,0,species_offset)
             - evaluate(q,i-1,j,0,0,species_offset);
             + evaluate(q,i,j-1,0,0,species_offset)
             - evaluate(q,i-1,j-1,0,0,species_offset);
      dx_val *= 0.5;
      break;
    case X:
      assert_eq(i,1);
      dx_val = evaluate(q,i,j,0,0,species_offset)
             + evaluate(q,i,j-1,0,0,species_offset);
      break;
    case Y:
      assert_eq(j,1);
      dx_val = evaluate(q,i,j,0,0,species_offset)
             - evaluate(q,i-1,j,0,0,species_offset);
      break;
    case XY:
      assert_eq(i,1); assert_eq(j,1);
      dx_val = 2.*evaluate(q,i,j,0,0,species_offset);
      break;
    default:
      invalid_value_error(gemParams.get_enforced_symmetry());
  }
  return dx_val/dogParamsCart2.get_dx();
}

// This assumes that the quantity that the "evaluate" function
// is antisymmetric in the direction of the y-axis.
// This is a second-order-accurate estimate, which is the
// best we can expect for a third-order method.
static double y_derivative_at_origin(
  double (*evaluate)(const dTensorBC4& q,
    int i, int j, int xpos, int ypos, int species_offset),
  const dTensorBC4& q, int species_offset)
{
  int i,j; get_center_cell_coords(i,j);
  double dy_val;
  switch(gemParams.get_enforced_symmetry()&XY)
  {
    case NONE:
      dy_val = evaluate(q,i,j,0,0,species_offset)
             - evaluate(q,i,j-1,0,0,species_offset)
             + evaluate(q,i-1,j,0,0,species_offset);
             - evaluate(q,i-1,j-1,0,0,species_offset);
      dy_val *= 0.5;
      break;
    case X:
      assert_eq(i,1);
      dy_val = evaluate(q,i,j,0,0,species_offset)
             - evaluate(q,i,j-1,0,0,species_offset);
      break;
    case Y:
      assert_eq(j,1);
      dy_val = evaluate(q,i,j,0,0,species_offset)
             + evaluate(q,i-1,j,0,0,species_offset);
      break;
    case XY:
      assert_eq(i,1); assert_eq(j,1);
      dy_val = 2.*evaluate(q,i,j,0,0,species_offset);
      break;
    default:
      invalid_value_error(gemParams.get_enforced_symmetry());
  }
  return dy_val/dogParamsCart2.get_dy();
}

void Output_center_gas_states(
  const dTensorBC4& q, double t, const char* outputfile,
  int species_offset, bool fiveMoment)
{
  const int _rho = 1;
  const int _M3  = 4;
  const double rho = origin_value(species_offset+_rho,q);
  const double M3  = origin_value(species_offset+_M3, q);
  
  bool append = (t!=0);
  FILE* file = fopen(outputfile, append?"a":"w");
  fprintf(file, "%24.16e %24.16e %24.16e ", t,rho,M3);
  if(!fiveMoment)
  {
    double P13_x = x_derivative_at_origin(evaluate_P13,q,species_offset);
    double P23_y = y_derivative_at_origin(evaluate_P23,q,species_offset);
    fprintf(file, "%24.16e %24.16e ", P13_x,P23_y);
    double entropy = evaluate_at_origin(evaluate_entropy_per_mass10,q,species_offset);
    fprintf(file, "%24.16e ", entropy);
  }
  else
  {
    fprintf(file, "0 0 ");
    double entropy = evaluate_at_origin(evaluate_entropy_per_mass05,q,species_offset);
    fprintf(file, "%24.16e ", entropy);
  }
  fprintf(file,"\n");
  fclose(file);
}

void output_component_at_center(
  const dTensorBC4& q, double t, const char* outputfile,
  int idx)
{
  const double S = origin_value(idx,q);
  
  bool append = (t!=0);
  FILE* file = fopen(outputfile, append?"a":"w");
  fprintf(file, "%24.16e %24.16e \n", t,S);
  fclose(file);
}

double matlab_fix_double(double val);

static double integrate_along_cell_edge(
  const dTensorBC4& q,
  int xidx, int yidx, int stateIdx,
  int x_or_y, int upper_or_lower)
{
  assert(x_or_y==X||x_or_y==Y);
  assert(abs(upper_or_lower)==1);
  double value=0.;
  // top edge
  switch(dogParams.get_space_order())
    {
    default:
      unsupported_value_error(dogParams.get_space_order());
    case 4:
    case 3: // third-order
      value += sq5*q.get(xidx,yidx,stateIdx,x_or_y==X?5:6);
    case 2: // second-order
      value += upper_or_lower*sq3*q.get(xidx,yidx,stateIdx,(x_or_y==X)?2:3);
    case 1: // first-order
      value += q.get(xidx,yidx,stateIdx,1);
    }
    return value;
}

static double x_edge_integral(
  double& x_axis_accum_flux,
  double& min_x_axis_accum_flux,
  double& total_x_axis_reverse_variation,
  const dTensorBC4& q,
  int stateIdx,
  int xstart, int xend, int yidx,
  int upper_or_lower,
  int direction_of_pos_flux)
{
  assert(abs(direction_of_pos_flux)==1);
  assert_le(xstart,xend);
  for (int i=xstart; i<=xend; i++)
  {             
      const double cell_edge_integral
        = integrate_along_cell_edge(q,i,yidx,stateIdx,Y,upper_or_lower);
        // = q.get(i,yidx,stateIdx,1);
      double flux = direction_of_pos_flux*cell_edge_integral;
      x_axis_accum_flux += flux;
      min_x_axis_accum_flux = Min(min_x_axis_accum_flux,x_axis_accum_flux);
      if(flux<0.0) total_x_axis_reverse_variation -= flux;
  }
}

static double y_edge_integral(
  double& accum_flux,
  const dTensorBC4& q,
  int stateIdx,
  int ystart, int yend,
  int xidx,
  int upper_or_lower,
  int direction_of_pos_flux)
{
  assert(abs(direction_of_pos_flux)==1);
  assert_le(ystart,yend);
  for (int j=ystart; j<=yend; j++)
  {             
      const double cell_edge_integral
        = integrate_along_cell_edge(q,xidx,j,stateIdx,X,upper_or_lower);
        // = q.get(xidx,j,stateIdx,1);
      accum_flux += direction_of_pos_flux*cell_edge_integral;
  }
}

void calculate_boundary_fluxes(
    int n_B1,
    int n_B2,
    const dTensorBC4& q,
    double& bottom_to_rght_flux,
    double& rght_flux,
    double& left_flux,
    double& left_lost_flux,
    double& total_x_axis_island_flux)
{
    int xidx; int yidx; get_center_cell_coords(xidx, yidx);
    int mx   = q.getsize(1);
    int my   = q.getsize(2);

    // 1 = positive quadrant, 0 = negative quadrant
    double x_axis_accum_flux_11=0.0,
           x_axis_accum_flux_01=0.0,
           x_axis_accum_flux_10=0.0,
           x_axis_accum_flux_00=0.0;
    double min_x_acc_flx_11=0.,
           min_x_acc_flx_01=0.,
           min_x_acc_flx_10=0.,
           min_x_acc_flx_00=0.;
    double total_x_axis_island_flux_11=0.,
           total_x_axis_island_flux_01=0.,
           total_x_axis_island_flux_10=0.,
           total_x_axis_island_flux_00=0.;
    // integrate along each edge of each quarter of the domain and average
    // (for perfect symmetry the integrals will all agree)
    x_edge_integral(
      x_axis_accum_flux_11, min_x_acc_flx_11, total_x_axis_island_flux_11,
      q, n_B2, xidx,mx, yidx,-1,1);
    if(gemParams.x_symmetry_not_enforced())
      x_edge_integral(
        x_axis_accum_flux_10, min_x_acc_flx_10, total_x_axis_island_flux_10,
        q, n_B2, 1,xidx-1, yidx,1,-1);
    if(gemParams.y_symmetry_not_enforced())
    {
      x_edge_integral(
        x_axis_accum_flux_01, min_x_acc_flx_01, total_x_axis_island_flux_01,
        q, n_B2, xidx,mx, yidx-1,1,1);
      if(gemParams.x_symmetry_not_enforced())
      {
        x_edge_integral(
          x_axis_accum_flux_00, min_x_acc_flx_00, total_x_axis_island_flux_00,
          q, n_B2, 1,xidx-1, yidx-1,-1,-1);
      }
    }
     double min_x_axis_accum_flux = (
      min_x_acc_flx_11 +
      min_x_acc_flx_01 +
      min_x_acc_flx_10 +
      min_x_acc_flx_00)/gemParams.get_num_quadrants();
    // divide the flux by the number of edges along which we integrated
    const double x_axis_accum_flux = (
      x_axis_accum_flux_11 +
      x_axis_accum_flux_01 +
      x_axis_accum_flux_10 +
      x_axis_accum_flux_00)/gemParams.get_num_quadrants();
    total_x_axis_island_flux = (
      total_x_axis_island_flux_11 +
      total_x_axis_island_flux_01 +
      total_x_axis_island_flux_10 +
      total_x_axis_island_flux_00)/gemParams.get_num_quadrants();
    //printf(
    //  "x_axis_accum_flux_11 = %e\n"
    //  "x_axis_accum_flux_01 = %e\n"
    //  "x_axis_accum_flux_10 = %e\n"
    //  "x_axis_accum_flux_00 = %e\n",
    //   x_axis_accum_flux_11,
    //   x_axis_accum_flux_01,
    //   x_axis_accum_flux_10,
    //   x_axis_accum_flux_00);

    double x_axis_flux = x_axis_accum_flux * dogParamsCart2.get_dx();
    min_x_axis_accum_flux *= dogParamsCart2.get_dx();
    bottom_to_rght_flux = x_axis_flux - min_x_axis_accum_flux;
    total_x_axis_island_flux *= dogParamsCart2.get_dx();

    double accum_left_flux_11 = 0.0,
           accum_left_flux_01 = 0.0,
           accum_left_flux_10 = 0.0,
           accum_left_flux_00 = 0.0;
    y_edge_integral(accum_left_flux_11, q, n_B1, yidx,my, xidx,-1,1);
    if(gemParams.x_symmetry_not_enforced())
      y_edge_integral(accum_left_flux_01, q, n_B1, yidx,my, xidx-1,1,1);
    if(gemParams.y_symmetry_not_enforced())
    {
      y_edge_integral(accum_left_flux_10, q, n_B1, 1,yidx-1, xidx,-1,-1);
      if(gemParams.x_symmetry_not_enforced())
        y_edge_integral(accum_left_flux_00, q, n_B1, 1,yidx-1, xidx-1,1,-1);
    }
    //printf(
    //  "accum_left_flux_11 = %e\n"
    //  "accum_left_flux_01 = %e\n"
    //  "accum_left_flux_10 = %e\n"
    //  "accum_left_flux_00 = %e\n",
    //   accum_left_flux_11,
    //   accum_left_flux_01,
    //   accum_left_flux_10,
    //   accum_left_flux_00);
    double accum_left_flux = (
      accum_left_flux_11 +
      accum_left_flux_01 +
      accum_left_flux_10 +
      accum_left_flux_00)/gemParams.get_num_quadrants();

    double accum_rght_flux_11 = 0.0,
           accum_rght_flux_01 = 0.0,
           accum_rght_flux_10 = 0.0,
           accum_rght_flux_00 = 0.0;
    y_edge_integral(accum_rght_flux_11, q, n_B1, yidx,my, mx,1,1);
    if(gemParams.x_symmetry_not_enforced())
      y_edge_integral(accum_rght_flux_01, q, n_B1, yidx,my, 1,-1,1);
    if(gemParams.y_symmetry_not_enforced())
    {
      y_edge_integral(accum_rght_flux_10, q, n_B1, 1,yidx-1, mx,1,-1);
      if(gemParams.x_symmetry_not_enforced())
        y_edge_integral(accum_rght_flux_00, q, n_B1, 1,yidx-1, 1,-1,-1);
    }
    //printf(
    //  "accum_rght_flux_11 = %e\n"
    //  "accum_rght_flux_01 = %e\n"
    //  "accum_rght_flux_10 = %e\n"
    //  "accum_rght_flux_00 = %e\n",
    //   accum_rght_flux_11,
    //   accum_rght_flux_01,
    //   accum_rght_flux_10,
    //   accum_rght_flux_00);
    double accum_rght_flux = (
      accum_rght_flux_11 +
      accum_rght_flux_01 +
      accum_rght_flux_10 +
      accum_rght_flux_00)/gemParams.get_num_quadrants();

    left_flux = accum_left_flux*dogParamsCart2.get_dy();
    rght_flux = accum_rght_flux*dogParamsCart2.get_dy();

    left_lost_flux = -min_x_axis_accum_flux;
}

double OutputMagneticFluxes(
    int n_B1,
    int n_B2,
    const dTensorBC4& q,
    double t)
{
    // -------------------------
    // RECONNECTED MAGNETIC FLUX
    // -------------------------
    //
    // Taking advantage of the symmetry of the GEM problem,
    // we define reconnected flux in terms of the magnetic flux
    // through the boundaries of the first quadrant, which
    // we here refer to as the left, right, bottom, and top.
    //
    // I now define reconnected flux as the decrease in the
    // integral of the flux that passes rightward through the
    // left boundary.  This ensures that reconnection rate
    // is minus the value of the third component of the
    // electric field at the origin.
    //
    // I previously defined reconnected flux to be the
    // bottom_to_rght_flux, i.e., the integral of the flux that
    // passes upward through the bottom boundary and exits
    // through the right boundary. (This definition was an
    // attempt to ensure that magnetic islands (enclosed by field
    // lines that do not pass through the right boundary) do not
    // contribute to the flux.)
    //
    // Let phi(x) be the accumlation integral of B_y from 0 to x.
    // Then the reconnected flux is phi(L_x/2) - min_{[0,L_x/2]}(phi),
    // assuming that all magnetic islands curl in the same
    // direction as the initial curl of B.
    //
    // We report:
    // (1) the bottom_to_rght_flux flux (the net flux of lines that enters
    //     through the bottom and exits through the right wall),
    // (2) the net flux through the right wall, 
    // (3) the net flux through the left wall, 
    // (4) the flux of lines that enter from the left and exit
    //     through the bottom (i.e. -min(phi)), and
    // (5) the total flux of field lines leaving through the bottom
    //     (i.e. the negative variation of phi); subtracting -min(phi)
    //     measures the strength of the magnetic islands along
    //     the bottom.
    //

    double bottom_to_rght_flux;
    double rght_flux;
    double left_flux;
    double left_lost_flux;
    double total_x_axis_island_flux;
    void calculate_boundary_fluxes(
        int n_B1,
        int n_B2,
        const dTensorBC4& q,
        double& bottom_to_rght_flux,
        double& rght_flux,
        double& left_flux,
        double& left_lost_flux,
        double& total_x_axis_island_flux);
    calculate_boundary_fluxes(
        n_B1,
        n_B2,
        q,
        bottom_to_rght_flux,
        rght_flux,
        left_flux,
        left_lost_flux,
        total_x_axis_island_flux);
    
    const char* get_outputdir();
    string outputfile = string(get_outputdir())+"/recon_flux.dat";
    bool append = (t!=0 || dogParams.get_mrestart()!=0);
    FILE* file = fopen(outputfile.c_str(), append?"a":"w");
    fprintf(file,"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e \n", t,
      bottom_to_rght_flux,
      rght_flux,
      left_flux,
      left_lost_flux,
      total_x_axis_island_flux);
    fclose(file);
    return bottom_to_rght_flux;
}

// generic conservation;
// should move this into library and replace this with something
// appropriate to the plasma context
//
void WriteConservation(
    const dTensorBC4& q,
    double t)
{
    // -----------------
    // CONSERVATION
    // -----------------
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    dTensor1 qsum(meqn);
    const int numelems = mx*my;
    for (int m=1; m<=meqn; m++)
    {
        qsum.set(m,0.0);
      
        for (int i=1; i<=mx; i++)
        for (int j=1; j<=my; j++)
        {            
            const double& qtmp = q.get(i,j,m,1);
        
            qsum.set(m, qsum.get(m) + qtmp/double(numelems) );
        }
    }

  // ideally in the case of a restart
  // we should truncate the recon_flux.dat file
  // at the restart time before we resume appending,
  // but this way we at least don't lose information
  // so the file can be fixed by hand.
  // 
  bool append = (t!=0 || dogParams.get_mrestart()!=0);
  const char* get_outputdir();
  string outputfile = string(get_outputdir())+"/conservation.dat";
  FILE* file = fopen(outputfile.c_str(), append?"a":"w");
  fprintf(file,"%24.16e ",t);
  for (int m=1; m<=meqn; m++)
  {
    fprintf(file,"%24.16e ", matlab_fix_double(qsum.get(m)));
  }
  fprintf(file,"\n");
  fclose(file);
}

static double cell_center_value(int i, int j, /* coordinates of cell */
  int stateIdx,
  const dTensorBC4& q)
{
  double ret=0.;
  switch(dogParams.get_space_order())
  {
    default:
      unsupported_value_error(dogParams.get_space_order());
    case 3: // third-order
      ret -= (sq5*0.5)*(q.get(i,j,stateIdx,5)+q.get(i,j,stateIdx,6));
    case 2: // second-order
    case 1: // first-order
      ret += q.get(i,j,stateIdx,1);
  }
  //assert_almost_eq(ret,cell_value(q,i,j,0,0,stateIdx));
  return ret;
}

void Output_center_E3(
  const dTensorBC4& q, double t,
  int _E3)
{
  const char* get_outputdir();
  string outputfile = string(get_outputdir())+"/xpoint_E3.dat";
  const double E3  = origin_value(_E3, q);
  bool append = (t!=0);
  FILE* file = fopen(outputfile.c_str(), append?"a":"w");
  fprintf(file,"%24.16e %24.16e \n",t,E3);
  fclose(file);
}

// we should really figure out the coefficients that a cell
// at the center would have.
void Output_center_cell_state(const dTensorBC4& q, double t)
{
  int xidx, yidx; get_center_cell_coords(xidx, yidx);
  const char* get_outputdir();
  string outputfile = string(get_outputdir())+"/xpoint_cell_state.dat";
  bool append = (t!=0);
  FILE* file = fopen(outputfile.c_str(), append?"a":"w");
  fprintf(file,"%24.16e ",t);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  for(int m=1;m<=meqn;m++)
  for(int k=1;k<=kmax;k++)
  {
    fprintf(file,"%24.16e ",q.get(xidx,yidx,m,k));
  }
  fprintf(file,"\n");
  fclose(file);
}

static void get_P_cc(double& P13_cc, double& P23_cc,
  const dTensorBC4& q, int xidx, int yidx, int species_offset)
{
  using namespace TenMomentComponentID;
  const double rho_cc = cell_center_value(xidx,yidx,species_offset+_rho,q);
  const double M1_cc  = cell_center_value(xidx,yidx,species_offset+_M1 ,q);
  const double M2_cc  = cell_center_value(xidx,yidx,species_offset+_M2 ,q);
  const double M3_cc  = cell_center_value(xidx,yidx,species_offset+_M3 ,q);
  const double N13_cc = cell_center_value(xidx,yidx,species_offset+_N13,q);
  const double N23_cc = cell_center_value(xidx,yidx,species_offset+_N23,q);
  P13_cc = N13_cc - M1_cc*M3_cc/rho_cc;
  P23_cc = N23_cc - M2_cc*M3_cc/rho_cc;
}

static void get_grad_P_at_origin(double& P13_x, double& P23_y,
  const dTensorBC4& q, int species_offset)
{
  int xidx, yidx; get_center_cell_coords(xidx, yidx);
  double P13_11, P23_11; get_P_cc(P13_11, P23_11, q, xidx, yidx, species_offset);
  // initialize all cells with rotated/mirrored values
  double P13_01 = -P13_11; double P23_10 = -P23_11; // may reset later
  double P13_10 =  P13_11; double P23_01 =  P23_11; // may reset later
  double P13_00 = -P13_11; double P23_00 = -P23_11;
  if(gemParams.x_symmetry_not_enforced())
  {
    get_P_cc(P13_01, P23_01, q, xidx-1, yidx, species_offset);
    // reinitialize mirror cells assuming rotational symmetry
    P13_10 = -P13_01; P23_10 = -P23_01;
  }
  if(gemParams.y_symmetry_not_enforced())
  {
    get_P_cc(P13_10, P23_10, q, xidx, yidx-1, species_offset);
    if(gemParams.x_symmetry_not_enforced())
    {
      get_P_cc(P13_00, P23_00, q, xidx-1, yidx-1, species_offset);
    }
    else
    {
      // reinitialize mirror cells assuming rotational symmetry
      P13_01 = -P13_10; P23_01 = -P23_10;
    }
  }
  P13_x = ((P13_11 - P13_01) + (P13_10 - P13_00))
    *0.5/dogParamsCart2.get_dx();
  P23_y = ((P23_11 - P23_10) + (P13_01 - P13_00))
    *0.5/dogParamsCart2.get_dy();
}

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

void output_local_divergence_error(
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

