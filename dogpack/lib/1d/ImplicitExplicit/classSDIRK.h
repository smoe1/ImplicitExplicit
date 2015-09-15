#ifndef SDIRK_H
#define SDIRK_H

#include <vector>

typedef std::vector<std::vector<double> > STDVECARRAY;
typedef std::vector<double> STDVEC;

class SDIRK
{
  public:
    SDIRK(int order,int stages);

    double coeff (int row, int column) {
      return matrix_[row][column];
    }

    double get_time(int row) {
      return times_[row];
    }
      
    int n_rows() { return matrix_.size();    }

    int n_cols() { return matrix_[0].size(); }

  private:
    void build_coeffs(int order, int stages);
    STDVECARRAY matrix_;
    STDVEC times_;
};
 
#endif
