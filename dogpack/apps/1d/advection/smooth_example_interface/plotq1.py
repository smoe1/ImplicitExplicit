#----------------------------------------------------------
def plotq1(m, space_order, meqn, mx, time, xc, qsoln, auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # Exact solution:
    width = 2.0*0.25
    xold = (xc - time)%1.
    qex  = 0.*qsoln[:,0]
    for (i,x) in enumerate(xold):
        if( abs(x-0.5) <= 0.5*width ):
            qex[i] = 1.0;#np.cos(np.pi*(x-0.5)/width)**6

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    #plt.gca().set_xlim([xc[0],xc[mx-1]])
    #plt.gca().set_ylim([-1.0,2.1])
    plt.plot(xc,qsoln[:,m],'bo')
    #plt.plot(xc,(abs(qsoln[:,m])+1.0e-16), 'bo', linewidth=2.0)
    #tmp1 = "".join(("q(t,x) at t = ",str(time)))
    tmp1 = "L2 Projection"
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.tight_layout()
    plt.plot([0.5,0.5],[0.0,2.1],'r')
    print mx/5, max(abs(qex-qsoln[:,m])),time
    #plt.savefig('solnalphais3015discon.pdf')
    plt.pause(0.1)
    plt.draw()
    
###----------------------------------------------------------
