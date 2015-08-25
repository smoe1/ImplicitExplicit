
#----------------------------------------------------------
def plotq1(m,meth1,meqn,mx,time,xc,qsoln,auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np
    from math import fabs


    rho = qsoln[:,0]

    u = qsoln[:,1]/qsoln[:,0]

    press = qsoln[:,2]-pow(qsoln[:,1],2)/qsoln[:,0]

    qheat = qsoln[:,3] - 3*press*u - rho*pow(u,3)

    r = qsoln[:,4] - 4.0*qheat*u - 6.0*press*pow(u,2) - rho*pow(u,4)     

    qlow  = min(rho)
    qhigh = max(rho)
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1 
      

    plt.figure(1)
    plt.clf()
    plt.plot(xc,rho,'b-')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    tmp1 = "".join(("Density at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    qlow  = min(u)
    qhigh = max(u)
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1 

    plt.figure(2)
    plt.clf()
    plt.plot(xc,u,'b-')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    tmp1 = "".join(("Velocity at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    qlow  = min(press)
    qhigh = max(press)
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1 

    plt.figure(3)
    plt.clf()
    plt.plot(xc,press,'b-')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    tmp1 = "".join(("Pressure at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    qlow  = min(qheat)
    qhigh = max(qheat)
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1 

    plt.figure(4)
    plt.clf()
    plt.plot(xc,qheat,'b-')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    tmp1 = "".join(("Heat at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    qlow  = min(r)
    qhigh = max(r)
    if fabs(qhigh-qlow)>1.0e-7:            
        qeps  = 0.1*(qhigh-qlow)
    else:
        qeps = 0.1 
   
    plt.figure(5)
    plt.clf()
    plt.plot(xc,r,'b-')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([qlow-qeps,qhigh+qeps])
    tmp1 = "".join(("r at t = "))
    tmp2 = "".join((tmp1,str(time)))
    title = "".join((tmp2,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
#----------------------------------------------------------
