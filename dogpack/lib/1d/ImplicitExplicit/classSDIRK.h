#include <vector>

typedef std::vector<std::vector<double> > STDVECARRAY;
typedef std::vector<double> STDVEC;

class SDIRK
{
  public:
    SDIRK(int order,int stages) {
      matrix_ = STDVECARRAY (stages, STDVEC (stages));
      times_  = STDVEC      (stages);
      build_coeffs(order, stages);
    }

    double coeff (int row, int column) {
      return matrix_[row][column];
    }
      
    int n_rows() { return matrix_.size();    }

    int n_cols() { return matrix_[0].size(); }

  private:
    void build_coeffs(int order, int stages);
    STDVECARRAY matrix_;
    STDVEC times_;
};

void SDIRK::build_coeffs(int order, int stages) {
  if (order==3 && stages==4) {
    double gamma=0.4358665215;     // y_{n+1}=Y_4
    // double gamma=0.2928932188;  // y_{n+1}=Y_3
    double gamma_sq=gamma*gamma;
    times_[0]=0.0;
    times_[1]=2*gamma;
    times_[2]=1.0;
    times_[3]=1.0;
    
    // row 0 intentionally zero
    matrix_[1][0]=gamma;
    matrix_[1][1]=gamma;
    matrix_[2][0]=0.25*(-4*gamma_sq+6*gamma-1.0)/gamma;
    matrix_[2][1]=0.25*(-2*gamma+1.0)/gamma;
    matrix_[2][2]=gamma;
    matrix_[3][0]=(6*gamma-1.0)/(12*gamma);
    matrix_[3][1]=-1.0/(24*gamma_sq-12*gamma);
    matrix_[3][2]=(-6*gamma_sq+6*gamma-1.0)/(6*gamma-3.0);
    matrix_[3][3]=gamma;
  }
}
