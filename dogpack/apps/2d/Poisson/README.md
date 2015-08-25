# Unstructured 2D Poisson solvers

The unstructured Poisson examples in this section require three operations:

1. Mesh generation software, MeshGenC++ in order to construct the mesh.
2. MATLAB to compute Sparse Cholesky factorization.
3. DoGPack to run the code.

For example, the following three steps will produce results for the 
annulus\_unst example.

1. Create the mesh 
-------------------

Make sure MeshGenC++ is updated to the latest version). Go to

    cd $DOGPACK/apps/2d/Poisson/annulus_unst/Unstructured\_Mesh

and execute:

    make
    mesh.exe

2. Create the Cholesky factorization
------------------------------------

Each application has an extra folder, called matlab, who purpose is to
create the Cholesky factorization.  We currently only support using
MATLAB (and possible Octave) to create the factorization.  To run this, open 
MATLAB, and change working directory to

    $DOGPACK/apps/2d/Poisson/annulus_unst/matlab 
    
From there, type

    CreateMatrix(morder)

which will call the CreateMatrix.m function.  Valid choices for morder are
1,2,3.

3. Run the code
---------------

After executing steps 1 and 2, you are now ready to run the code.  Go to

    cd $DOGPACK/apps/2d/Poisson/annulus_unst/

and type

    make
    dog.exe

You can plot the results as usual in 
[MATLAB plotting routines](../../../viz/matlab/README.md) 
via "plotdog2.m", or with the 
[Python plotting routines](../../../viz/python/README.md).
