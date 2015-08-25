#include "../../../../../lib/2d/cart/dog_defs.h"
#include "../../../../../lib/tensors.h"
#include "discontCell.h"

//////////////////////////////////////////////////////////////////////////////
// Tester function for discontCell.cpp
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

	int morder = 3;
	int clean_cut = 1;   // == 1 for clean cut, otherwise 0 bad cut.
	int kmax = int(morder*(morder+1)/2);
	dTensor2 q1(morder, kmax); // 'left'  legendre weights
	dTensor2 q2(morder, kmax); // 'right' legendre weights
	dTensor1 qvals1(kmax);
	dTensor1 qvals2(kmax);
	dTensor1 scut(morder);

	if( !clean_cut )
	{   // not a clean cut
		for(int n=1; n<= morder; n++)
		{
			scut.set(n, -1.0);
		}
		scut.set(2, 1.0);
	}
	else
	{  // clean cut.  this occurs at s = 0
		for(int n=1; n<= morder; n++)
		{
			scut.set(n, 0.0);
		}
	}

	// initiate ql and qr to 0
	for(int k=1; k<=kmax; k++)
	for(int l=1; l<=morder; l++)
	{
		q1.set(l, k, 0.0);
		q2.set(l, k, 0.0);
	}
	
	// set some values to be non-zero
	for(int l=1; l<=morder; l++)
	{
		q1.set(l, 1,  1.0);
		q1.set(l, 1,  1.0);
		q2.set(l, 1, 1.0);
		q2.set(l, 1, 1.0);
	}
	
	//discontCell(int kmax, double initscut, double cell_area);
	cout << "constructing cell" << endl;
//	discontCell *cell1 = new discontCell(kmax, scut, 1);
//	discontCell *cell2 = new discontCell(kmax, scut, 2);
	discontCell *cell1 = new discontCell(kmax, 1);
	discontCell *cell2 = new discontCell(kmax, 2);
	cell1->setScut(scut);
	cell1->project(q1, q2, qvals1);
	cell1->setScut(scut);
	cell1->project(q1, q2, qvals1);
//	cell2->project(q1, q2, qvals2);

	//cell->printwgts();
	cout << "calling printphi" << endl;
	cell1->printintPhi();
	cout << "called printphi" << endl;
//	cell2->printintPhi();

//	cell1->printphivals();

	cout << endl << " Printing qvals (in tester file) ... " << endl;
	
	cout << "k=1,2,...,kmax, qvals1(k) has value: " << endl;
	for(int k=1; k<= kmax; k++ )
	{
		cout << "     k = ";
		cout << k << ";     qvals1(k) = " << qvals1.get(k) << endl;
	}
	
	cout << "k=1,2,...,kmax, qvals2(k) has value: " << endl;
	for(int k=1; k<= kmax; k++ )
	{
		cout << "     k = ";
		cout << k << ";     qvals2(k) = " << qvals2.get(k) << endl;
	}
	
	delete cell1;  delete cell2;
	
	return 0;

}
