
#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln,auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-1.1,1.1])
    plt.plot(xc,qsoln[:,m],'b-')
    tmp1 = "".join(("q(t,x) at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#----------------------------------------------------------
