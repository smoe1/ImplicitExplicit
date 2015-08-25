def qinit( xc, time, direction ):

    import numpy as np

    # Exact solution:
#   width = 2.0*0.25

    if( direction > 0 ):
        xold = ( (xc+1.0 - time)%2. ) - 1.0
    else:
        xold = ( (xc+1.0 + time)%2. ) - 1.0
        
    qex  = 0.*xc
    for (i,x) in enumerate(xold):
        if( abs(x+0.0) < 0.2 ):
            qex[i] = np.cos( (x+0.0)*np.pi/(2.0*0.2) )**6
    return qex

#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln,auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # Exact solution (left and right going characteristics)
    qex = 0.5*( qinit( xc, time, -1 ) + qinit( xc, time, 1 ) )

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,1.1])
    plt.plot(xc,qsoln[:,0],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    # Exact solution (left and right going characteristics)
    qex = 0.5*( qinit( xc, time, 1 ) - qinit( xc, time, -1 ) )
    
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.6,0.6])
    plt.plot(xc,qsoln[:,1],'bo')
    plt.plot(xc,qex, '-r', linewidth=2.0)
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
