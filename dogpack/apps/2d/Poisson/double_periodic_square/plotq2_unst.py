
#----------------------------------------------------------
def plotq2_unst(outputdir,nframe,
                m,meqn,NumPhysElems,NumPhysNodes,
                xlow,xhigh,ylow,yhigh,time,
                x,y,tnode,qsoln,
                xmid,ymid,qsoln_elem):
    
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim([xlow,xhigh])
    plt.gca().set_ylim([ylow,yhigh])
    p1=plt.tripcolor(x, y, tnode, qsoln_elem[:,m], shading='flat')
    tmp1 = "".join(("q(",str(m+1),") at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    cbar = plt.colorbar()
    plt.draw()

    # exact solution
    phiex  = np.zeros(NumPhysElems,float);
    const = np.sin(2.4*np.pi)*np.sin(2.4*np.pi);
    E1ex  = np.zeros(NumPhysElems,float);
    E2ex  = np.zeros(NumPhysElems,float);
    for i in range(0,NumPhysElems):
        xm = xmid[i]+0.2;
        ym = ymid[i]+0.2;
        phiex[i] = np.sin(2.0*np.pi*xm)*np.sin(2.0*np.pi*ym)-const;
        E1ex[i] = -2*np.pi*np.cos(2.0*np.pi*xm)*np.sin(2.0*np.pi*ym);
        E2ex[i] = -2*np.pi*np.sin(2.0*np.pi*xm)*np.cos(2.0*np.pi*ym);

    # error
    err1 = np.linalg.norm(phiex-qsoln_elem[:,0])/np.linalg.norm(phiex);
    err2 = np.linalg.norm(E1ex-qsoln_elem[:,1])/np.linalg.norm(E1ex);
    err3 = np.linalg.norm(E2ex-qsoln_elem[:,2])/np.linalg.norm(E2ex);

    print " "
    print " err in phi ="'%16.8e'%err1
    print " err in  E1 ="'%16.8e'%err2
    print " err in  E2 ="'%16.8e'%err3
    print " "

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim([xlow,xhigh])
    plt.gca().set_ylim([ylow,yhigh])
    if m==0:
        p1=plt.tripcolor(x, y, tnode, phiex, shading='flat');
    elif m==1:
        p1=plt.tripcolor(x, y, tnode, E1ex, shading='flat');
    elif m==2:
        p1=plt.tripcolor(x, y, tnode, E2ex, shading='flat');
    tmp1 = "".join(("Exact solution at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    cbar = plt.colorbar()
    plt.draw()

    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('equal')
    plt.gca().set_xlim([xlow,xhigh])
    plt.gca().set_ylim([ylow,yhigh])
    if m==0:
        p1=plt.tripcolor(x, y, tnode, np.abs(phiex-qsoln_elem[:,0]), shading='flat');
    elif m==1:
        p1=plt.tripcolor(x, y, tnode, np.abs(E1ex-qsoln_elem[:,1]), shading='flat');
    elif m==2:
        p1=plt.tripcolor(x, y, tnode, np.abs(E2ex-qsoln_elem[:,2]), shading='flat');
    tmp1 = "".join(("Absolute error at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    cbar = plt.colorbar()
    plt.draw()
    
    
#----------------------------------------------------------
