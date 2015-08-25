// --------------------------------------------------------------------------
//  IMPLEMENTATION FILE (discontCellConstCoeff.cpp)
//
//  Fill in a description here...
//    This class provides a method for projecting this cell onto 
//     the legendre basis.
//
// --------------------------------------------------------------------------
//		THIS PART CURRENTLY IS UNOPERATIONAL
//		if discont_direction == 0: the integration is broken up into 4 steps.
//			cell1 = lower left cell
//			cell2 = lower right cell
//			cell3 = upper left cell
//			cell4 = upper right cell
// --------------------------------------------------------------------------
//
//
//		discont_direction == 1: the integration is broken up into 2 pieces.
//
//			cell1 = left cell
//			cell2 = right cell
//
//		discont_direction == 2: the integration is broken up into 2 pieces
//
//			cell1 = bottom cell
//			cell2 = upper cell 
//
//     Gaussian quadrature is used to compute the integration.
//
// --------------------------------------------------------------------------

#include "../defs.h"
#include "discontCell.h"

// --------------------------------------------------------------------------
//        Private members of class:
//
//	int mpoints, mpoints1d, kmax
//	double scut, width1, width2
//
//  int discont_direction;		//type of cell constructed.  
//
//	dTensor1 *wgt;		//gauss quad weights
//	dTensor2 *q1, *q2;	//Legendre weights on left and right half of cell
//	
//   q's values on left and right half of cell evaluated at each quadrature
//   point
//	dTensor1 *qvals1, *qvals2;
//	
//	legendre polynomial values evaluated on gauss quad grid points
//	size is (mpoints, kmax) for each of these
//	dTensor2 *phi1, *phi2, *centeredPhi2, *centeredPhi1;
//
//   // gaussian quadrature points necessary for integration
//	dTensor2 *spts1(mpoints,2), *spts2(mpoints,2)     
//			
//   Points are labeled according to the following format:
//		1:mpoints1d - are the points in the 'first row'
//        mpoints1d+1 : 2*mpoints1d - are the hpoints in the 'second row'.
//        2mpoints1d+1 : 3*mpoints1d - are the hpoints in the 'third row'.
//
// --------------------------------------------------------------------------

discontCellConstCoeff::discontCellConstCoeff(int kmax, double scut, int discont_direction)
// Constructor
// POST: Create a discontinuous cell 
{

	//type of discontinuity that will occur
	this->discont_direction = discont_direction;	
	this->kmax = kmax;	//number of polynomials used
	this->scut = scut;	// location of discontinuity

	width1 = scut+1;
	width2 = 1-scut;

	mpoints   = int((1+4*kmax-int(sqrt(1+8*kmax)))/2);
	mpoints1d = int(sqrt(mpoints));
//	scut = new dTensor1(mpoints1d);	//location of cut (located in [-1,1])

	// width of left/right cells inside [-1,1]
//	width1 = new dTensor1(mpoints1d);
//	width2 = new dTensor1(mpoints1d);

	//quick error check:
	if(mpoints1d > 5 || mpoints1d < 1)
	{
		cout << " error constructing cell..." << endl;
		cout << "          kmax = " << kmax << endl;
		cout << "     mpoints1d = " << mpoints1d << endl;
		cout << "     discont_direction = " << discont_direction << endl;
		exit(1);
	}

	///////////////////////////////////////////////////////////////////////////
	// set up gaussian quadrature points on left and right half of cell
	//  find the 1d gaussian quad wgts and points
	w1d = new dTensor1(mpoints1d);
	x1d = new dTensor1(mpoints1d);
	setGaussPoints(); 
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	// the integrated legendre polynomials (integrated in "rows")
	intPhi = new dTensor4(mpoints1d, kmax, kmax, 2);
	///////////////////////////////////////////////////////////////////////////
	
	// integrate the legendre polynomials
	setPhi();

}

discontCellConstCoeff::~discontCellConstCoeff()
// Destructor
// POST: discontCellCoeff no longer exists
{
	delete w1d; delete x1d; 
	delete intPhi;
}

void discontCellConstCoeff::setGaussPoints()
{

    // ---------------------------------
    // Set quadrature weights and points
    // ---------------------------------
    switch ( mpoints1d )
    {
        case 1:
	  w1d->set(1, 2.0 );
	  
	  x1d->set(1, 0.0 );
	  break;

        case 2:
	  w1d->set(1,  1.0 );
	  w1d->set(2,  1.0 );
	  
	  x1d->set(1, -1.0/sq3 );
	  x1d->set(2,  1.0/sq3 );
	  break;
	  
        case 3:
	  w1d->set(1,  5.0/9.0 );
	  w1d->set(2,  8.0/9.0 );
	  w1d->set(3,  5.0/9.0 );
	  
	  x1d->set(1,  -sq3/sq5 );
	  x1d->set(2,  0.0 );
	  x1d->set(3, sq3/sq5 );
	  break;

        case 4:
	  w1d->set(1, (18.0 - sq3*sq10)/36.0 );
	  w1d->set(2, (18.0 + sq3*sq10)/36.0 );
	  w1d->set(3, w1d->get(2) );
	  w1d->set(4, w1d->get(1) );
	  
	  x1d->set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
	  x1d->set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
	  x1d->set(3, -x1d->get(2) );
	  x1d->set(4, -x1d->get(1) );	      
	  break;
	     
        case 5:	     
	  w1d->set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
	  w1d->set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
	  w1d->set(3, 128.0/225.0 );
	  w1d->set(4, w1d->get(2) );
	  w1d->set(5, w1d->get(1) );
	  
	  x1d->set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
	  x1d->set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
	  x1d->set(3,  0.0 );
	  x1d->set(4, -x1d->get(2) );
	  x1d->set(5, -x1d->get(1) );
	  break;
    }

}// end of setGausspts

// function for intPhiprinting information currently stored in the cell
void discontCellConstCoeff::printintPhi() const
{

	cout << " Printing intPhi... " << endl;
	cout << "         kmax = " << kmax << endl;
	for(int n=1; n<=mpoints1d; n++)
	{
	cout << "         scut(n) = " << scut << endl;
	}
	cout << "         mpoints = " << mpoints << endl;
	
	cout << "intPhi(j,k,n,1) has values: " << endl;
	cout.precision(5);
	for(int j=1; j<= mpoints1d; j++ )
	for(int k=1; k<=kmax; k++)
	for(int n=1; n<=kmax; n++)
	{
		cout << "    j = " << j << ":   k = " << k << ":   n = " << n;
		cout << ":    eta(j) = " << x1d->get(j) << 
		":    intPhi = " << intPhi->get(j,k,n,1) << endl;
	}
	
	cout << endl;
	
	cout << "spts2 has values (x,y): " << endl;
	for(int n=1; n<= mpoints; n++ )
	{
//		cout << "("<< spts2->get(n,1) << ", "<<  spts2->get(n,2) << " )" << endl;
	}

}//end of printspts

//function for printing information stored in phivals
void discontCellConstCoeff::printwgts() const
{
	int m;
	cout << "     **** Printing wgts ******* " << endl;
	for(m=1; m<=w1d->getsize(); m++)
	{
		cout << "      w1d->get(m) = " << w1d->get(m) << endl;
	}
	cout << endl << endl;
}


//function for printing information stored in phivals
void discontCellConstCoeff::printphivals() const
{

	return;

}

///////////////////////////////////////////////////////////////////////////////
// 
//	Routine for computing the projection onto the legendre basis.
//
// Parameters:
// 
//   dTensor2(mpoints1d, kmax) q1new = coefficient weights on legendre polys for left
//   half of grid cell
//
//   dTensor2(mpoints1d, kmax) q2new = coefficient weights on legendre polys for right
//   half of grid cell
//
//  Returns:
//
//   dTensor1(kmax) qnew = coefficient weights on legendre polys over entire
//   cell.
//
///////////////////////////////////////////////////////////////////////////////
void discontCellConstCoeff::project(const dTensor2& q1new, 
	const dTensor2& q2new, dTensor1& qnew)
{
	
	double tmp;
	int n, j, k;
	
	for(k = 1; k <= kmax; k++)
	{
		tmp = 0.0;
		for(n = 1; n <= kmax; n++ )
		for(j = 1; j<= mpoints1d; j++ )
		{
			tmp += w1d->get(j)*intPhi->get(j,k,n,1)*q1new.get(j,n) +
				w1d->get(j)*intPhi->get(j,k,n,2)*q2new.get(j,n);
		}
		
		qnew.set( k, tmp / 4.0 );
	
	}//end of setting weights;

}// end of project

///////////////////////////////////////////////////////////////////////////////
// Sets qvals1/qvals2 to the exact value of q on the left/right 
// hand points.  Method is needed to perform the integration
///////////////////////////////////////////////////////////////////////////////
void discontCellConstCoeff::setqValues()
{
	double tmp1, tmp2;
	int m, m1, m2, k;

	m = 0;
	for(m1=1; m1<= mpoints1d; m1++)  // 'row' number
	for(m2=1; m2<= mpoints1d; m2++)  // 'column' number
	{
		m++;
		
		tmp1 = 0.0;
		tmp2 = 0.0;

		for(k=1; k<= kmax; k++)
		{
			// only variation is in the 'vertical' direction, so just need to
			// sample each 'row' for each of the q1, q2.
//			tmp1 += (q1->get(m2, k)) * (phi1->get(m,k)); 
//			tmp2 += (q2->get(m2, k)) * (phi2->get(m,k));
		}
//cout << "  m1 = " << m1 << "   m2 = " << m2 << ":  m = " << m << endl;
		qvals1->set(m, tmp1);
		qvals2->set(m, tmp2);

	}// end of looping over all points
}// end of setting qvalues
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Evalutate the function, phi(k)(xi,eta)
// not pretty, but this will work.
///////////////////////////////////////////////////////////////////////////////
double discontCellConstCoeff::evalPhi(const int& k, const double &xi, 
	const double &eta)
{

	double xi2,xi3,xi4,eta2,eta3,eta4;

	xi2 = xi*xi;
	xi3 = xi*xi2;
	xi4 = xi*xi3;
	eta2 = eta*eta;
	eta3 = eta*eta2;
	eta4 = eta*eta3;

	// Legendre basis functions evaluated at (xi,eta) in the
	// interval [-1,1]x[-1,1].
	switch( k )
	{
		case 15: 
		return ( 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );

		case 14: 
		return(  105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );

		case 13: 
		return(  5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
		
		case 12: 
		return(  sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );

		case 11: 
		return(  sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

		case 10:
		return(  sq7*(2.5*eta3 - 1.5*eta) );
		
		case 9:
		return(  sq7*(2.5*xi3 - 1.5*xi) );
		
		case 8:
		return(   sq3*sq5*xi*(1.5*eta2 - 0.5) );

		case 7: 
		return(   sq3*sq5*eta*(1.5*xi2 - 0.5) );

		case 6:
		return(   sq5*(1.5*eta2 - 0.5) );

		case 5:
		return(   sq5*(1.5*xi2 - 0.5) );
		
		case 4:
		return(   3.0*xi*eta );		    
	
		case 3:
		return(  sq3*eta );

		case 2:
		return(  sq3*xi  );

		case 1: 
		return(  1.0 );

	}

	// still here?
	cout << "bad polynomial value called" << endl;
	exit(1);

}// end of function evalPhi
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
//  Function for evaluating the integral of q(xi,eta) in one direction.
//
//  For discontinuity in the 'x'direction, we have:
//
//  F^(k) ( eta_j ) = \sum_n ( Q1^(n) *intPhi(j, k, 1, n) + 
//                     Q2^(n) * intPhi(j,k,2,n).
//
//  intPhi(j,k,1,n) = \int_{-1}^s \phil^(n) \phi^k
//  intPhi(j,k,2,n) = \int_{s}^1 \phir^(n) \phi^k
//
//  Where j is the index of the 'row'
///////////////////////////////////////////////////////////////////////////////
void discontCellConstCoeff::setPhi()
{

	int k, n, i, j;
	double tmp1, tmp2;

	//x1( i, j, 1 ) = values for phi^(k)  (centered poly)
	//x1( i, j, 2 ) = values for phi^(n)  (off centered poly)
	dTensor3 x1( mpoints1d, mpoints1d, 2 );  // quadrature points for 'left'
	dTensor3 x2( mpoints1d, mpoints1d, 2 );  // quadrature points on 'right' 

	for( j=1; j<= mpoints1d; j++ )
	for( i=1; i<= mpoints1d; i++ )
	{
//		tmp1 =  (scut->get(j) - (-1) ) / 2.0 * x1d->get(i) +
//				(scut->get(j) + (-1) ) / 2.0;
//		tmp2 =  (   1 - scut->get(j) ) / 2.0 * x1d->get(i) +
//				(scut->get(j) + 1 ) / 2.0;

		tmp1 =  (scut - (-1) ) / 2.0 * x1d->get(i) +
				(scut + (-1) ) / 2.0;
		tmp2 =  (   1 - scut ) / 2.0 * x1d->get(i) +
				(scut + 1 ) / 2.0;
		
		x1.set(i, j, 1, tmp1);
		x2.set(i, j, 1, tmp2);
		
		tmp1 -= (scut - 1 );
		tmp2 -= (scut + 1 );

		x1.set(i, j, 2, tmp1);
		x2.set(i, j, 2, tmp2);
	
	}

	for( k=1; k<=kmax; k++)
	for( n=1; n<=kmax; n++)
	for( j=1; j<=mpoints1d; j++)
	{
		tmp1 = 0.0;
		tmp2 = 0.0;
		for( i=1; i<= mpoints1d; i++)
		{

			switch(discont_direction)
			{
			case 1:  //advect in x-direction
			tmp1 += w1d->get(i)*width1 / 2.0 *
					evalPhi(k, x1.get(i, j, 1), x1d->get(j) ) *
					evalPhi(n, x1.get(i, j, 2), x1d->get(j) );
			tmp2 += w1d->get(i)*width2 / 2.0 *
					evalPhi(k, x2.get(i, j, 1), x1d->get(j) ) *
					evalPhi(n, x2.get(i, j, 2), x1d->get(j) );
			break;

			case 2:
			tmp1 += w1d->get(i)*width1 / 2.0 *
					evalPhi(k, x1d->get(j), x1.get(i, j, 1) ) *
					evalPhi(n, x1d->get(j), x1.get(i, j, 2) );
			tmp2 += w1d->get(i)*width2 / 2.0 *
					evalPhi(k, x1d->get(j), x2.get(i, j, 1) ) *
					evalPhi(n, x1d->get(j), x2.get(i, j, 2) );
			break;
			}// end of switch statement
			
		}
//cout << "    j = " << j << ":   k = " << k << ":   n = " << n << endl;
//cout << "    intPhi = " << intPhi->get(j,k,n,1) << endl;

		intPhi->set( j, k, n, 1, tmp1 );
		intPhi->set( j, k, n, 2, tmp2 );

	}

}

