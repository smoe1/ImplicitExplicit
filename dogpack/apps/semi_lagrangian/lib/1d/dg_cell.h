#ifndef DG_CELL_H
#define DG_CELL_H
#include "tensors.h"
class dg_cell
{

    // quadrature weights and points.  NOTE: w1d[0] is never accessed in
    // order to avoid confusion with 1 and 0 based indexing
//  static dTensor1 *w1d[6];
//  static dTensor1 *x1d[6];

    public:

        dg_cell();
        ~dg_cell();

        // function that simplly computes:
        //      phik( k ) = \int_a^b \phi^{(k)} ( xi ) d xi
        void integratePhi(double a, double b, dTensor1& phik );
        double integratePhi(double a, double b, int k );

    private:

        // quadrature weights and points.  NOTE: w1d[0] is never accessed in
        // order to avoid confusion with 1 and 0 based indexing
        void SetW1dX1d(dTensor1 &w1d, dTensor1& x1d);
        void evaluatePolys( const dTensor1& spts, dTensor2& phi );


};
#endif
