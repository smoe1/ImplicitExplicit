#######################################################################################
##
##   Information you have available to you:
##
##     Basic information:
##                   nframe:  frame number
##               mx, my, mz:  # of grid points
##               dx, dy, dz:  grid spacing
##   mx_old, my_old, mz_old:  # of grid points (original mesh before sub-dividing)
##   dx_old, dy_old, dz_old:  grid spacing (original mesh before sub-dividing)
##                     meqn:  number of equations
##                    meth1:  spatial order of accuracy
##                     time:  time of current slice
##                        m:  solution component that user asked to plot
## [xlow,xhigh,ylow,yhigh,zlow,zhigh]:  min/max values of grid
##
##   Slice information:
##        numslices:  number of slices
##        slicekind:  type of slice (1:z=const, 2:y=const, 3:x=const),
##                    size = (numslices)
##         sliceidx:  index of slice, size = (numslices)
##              ms1:  max value of 1st index in slice, size = (numslices)
##              ms2:  max value of 2nd index in slice, size = (numslices)
##
##   Grid information:
##       (xc_z, yc_z): grid points (cell centers) at constant z, size = (mx,my)
##       (xc_y, zc_y): grid points (cell centers) at constant y, size = (mx,mz)
##       (yc_x, zc_x): grid points (cell centers) at constant x, size = (my,mz)
##       (xl_z, yl_z): grid points (lower left corners) at constant z, size = (mx+1,my+1)
##       (xl_y, zl_y): grid points (lower left corners) at constant y, size = (mx+1,mz+1)
##       (yl_x, zl_x): grid points (lower left corners) at constant x, size = (my+1,mz+1)
##
##   Solution information:
##          qaug:  solution sampled on mesh, with zero padding to
##                 make it compatible with surf and pcolor matlab
##                 plotting tools, size = (mmax+1,mmax+1,meqn,numslices),
##                 where mmax = max([mx,my,mz])
##
#######################################################################################

#----------------------------------------------------------
def plotq3_cart_slices(outputdir,nframe,
                       numslices,slicekind,sliceidx,
                       m,meth1,meqn,ms1,ms2,time,s2d,
                       xlow,xhigh,ylow,yhigh,zlow,zhigh,
                       mx,my,mz,dx,dy,dz,
                       mx_old,my_old,mz_old,dx_old,dy_old,dz_old,
                       points_per_dir,LegVals,
                       xc_z,yc_z,xc_y,zc_y,yc_x,zc_x,
                       xl_z,yl_z,xl_y,zl_y,yl_x,zl_x,qaug):

    import os
    import matplotlib.pyplot as plt
    import numpy as np

    for ns in range(0,numslices):

        if slicekind[ns]==1:            
            zconst = zlow + (sliceidx[ns]-0.5)*dz_old;            
            plt.figure(ns)
            plt.clf()
            plt.gca().set_aspect('equal')
            plt.gca().set_xlim(xl_z[0,0],xl_z[ms1[ns],0])
            plt.gca().set_ylim(yl_z[0,0],yl_z[0,ms2[ns]])
            plt.pcolor(xl_z,yl_z,qaug[:,:,m,ns])
            cbar = plt.colorbar()
            tmp1 = "".join(("Slice of q(",str(m+1),") at z = "))
            tmp2 = "".join((tmp1,str(zconst)))
            tmp3 = "".join((tmp2,(" at time t = ")))
            tmp4 = "".join((tmp3,str(time)))
            title = "".join((tmp4,"     [DoGPack]"))
            plt.title(title)
            plt.draw()            
        elif slicekind[ns]==2:
            yconst = ylow + (sliceidx[ns]-0.5)*dy_old;
            plt.figure(ns)
            plt.clf()
            plt.gca().set_aspect('equal')
            plt.gca().set_xlim(xl_y[0,0],xl_y[ms1[ns],0])
            plt.gca().set_ylim(zl_y[0,0],zl_y[0,ms2[ns]])
            plt.pcolor(xl_y,zl_y,qaug[:,:,m,ns])
            cbar = plt.colorbar()
            tmp1 = "".join(("Slice of q(",str(m+1),") at y = "))
            tmp2 = "".join((tmp1,str(yconst)))
            tmp3 = "".join((tmp2,(" at time t = ")))
            tmp4 = "".join((tmp3,str(time)))
            title = "".join((tmp4,"     [DoGPack]"))
            plt.title(title)
            plt.draw()
        elif slicekind[ns]==3:
            xconst = xlow + (sliceidx[ns]-0.5)*dx_old;
            plt.figure(ns)
            plt.clf()
            plt.gca().set_aspect('equal')
            plt.gca().set_xlim(yl_x[0,0],yl_x[ms1[ns],0])
            plt.gca().set_ylim(zl_x[0,0],zl_x[0,ms2[ns]])
            plt.pcolor(yl_x,zl_x,qaug[:,:,m,ns])
            cbar = plt.colorbar()
            tmp1 = "".join(("Slice of q(",str(m+1),") at x = "))
            tmp2 = "".join((tmp1,str(xconst)))
            tmp3 = "".join((tmp2,(" at time t = ")))
            tmp4 = "".join((tmp3,str(time)))
            title = "".join((tmp4,"     [DoGPack]"))
            plt.title(title)
            plt.draw()


    # Scatter plot -- read-in data from 1D output
    dogpack = os.getenv("DOGPACK");
    outputdir1d = "".join((dogpack,"/apps/1d/euler/spherically_symmetric/output_3d_compare/"))
    TF = os.path.exists(outputdir1d)
    if TF==False:
        print " "
        print " NOTE: If you want to compare 3D scatter plot results to 1D data, "
        print "       you need to run the following code: "
        print " $DOGPACK/apps/1d/euler/spherically_symmetric/dog.exe -o output_3d_compare"
        print " "

    PLOT1D = False
    if TF==True:
        q1d_file = "".join((str(nframe+10000),".dat"))
        q1d_file = "s" + q1d_file[1:]
        q1d_file = "".join(("qpointval",q1d_file))
        q1d_file = "".join((outputdir1d,q1d_file))
        from local_helper import read_q1d_data 
        time1d,r1d,q1d = read_q1d_data(q1d_file)

        if ( (np.abs(time1d-time)<1.0e-8) and (m==0 or m==4)):
            PLOT1D = True
            if m==0:
                m1d = 0
            else:
                m1d = 2
            

    # Scatter plot -- create scatter data from given slices
    nnn = 0;
    for ns in range(0,numslices):
        if slicekind[ns]==1:
            nnn = nnn + (mx)*(my);
        elif slicekind[ns]==2:
            nnn = nnn + (mx)*(mz);
        elif slicekind[ns]==3:
            nnn = nnn + (my)*(mz);
            
    rad  = np.zeros(nnn,float);
    qrad = np.zeros(nnn,float);
    krad = 0;
    
    for ns in range(0,numslices):        

        if slicekind[ns]==1:
            for mi1 in range(0,mx):
                for mi2 in range(0,my):                    
                    rad[krad]  = np.sqrt(xc_z[mi1,mi2]*xc_z[mi1,mi2]+yc_z[mi1,mi2]*yc_z[mi1,mi2]);
                    qrad[krad] = qaug[mi1,mi2,m,ns];
                    krad = krad+1;
        elif slicekind[ns]==2:
            for mi1 in range(0,mx):
                for mi2 in range(0,mz):                    
                    rad[krad]  = np.sqrt(xc_y[mi1,mi2]*xc_y[mi1,mi2]+zc_y[mi1,mi2]*zc_y[mi1,mi2]);
                    qrad[krad] = qaug[mi1,mi2,m,ns];
                    krad = krad+1;
        elif slicekind[ns]==3:
            for mi1 in range(0,my):
                for mi2 in range(0,mz):                    
                    rad[krad]  = np.sqrt(xc_z[mi1,mi2]*xc_z[mi1,mi2]+yc_z[mi1,mi2]*yc_z[mi1,mi2]);
                    qrad[krad] = qaug[mi1,mi2,m,ns];
                    krad = krad+1;

    # Scatter plot -- create plot
    plt.figure(numslices)
    plt.clf()            
    plt.gca().set_xlim(0,1)
    plt.plot(rad,qrad,'ro')            
    tmp1 = "".join(("Scatter of q(",str(m+1),")"))
    tmp3 = "".join((tmp1,(" at time t = ")))
    tmp4 = "".join((tmp3,str(time)))
    title = "".join((tmp4,"     [DoGPack]"))
    plt.title(title)

    print " PLOT1D = ",PLOT1D
    if (PLOT1D == True):
        p1d = plt.plot(r1d,q1d[:,m1d],'k-')
    
    plt.draw()
               
#----------------------------------------------------------
