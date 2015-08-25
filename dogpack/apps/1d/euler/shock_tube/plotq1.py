
#-----------------------------------------------------------------------------#
def plotq1(m,meth1,meqn,mx,time,xc,qsoln,auxsoln):
    
    import matplotlib.pyplot as plt
    import numpy as np

    # HACK:
    gamma = 1.4
    gp1 = gamma+1.
    gm1 = gamma-1.

    # DENSITY
    m=0
    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.9,3.1])
    plt.plot(xc,qsoln[:,m],'bo')
    tmp1 = "".join(("Density at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()

    # PRESSURE
    m=1

    press = (gamma-1)*(qsoln[:,4]-0.5*(
        qsoln[:,1]**2 + qsoln[:,2]**2 + qsoln[:,3]**2)/qsoln[:,0])
    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([0.8, 3.2])
    plt.plot(xc,press,'bo')
    tmp1 = "".join(("Pressure at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
 
    # Velocity
    m=2
    u = np.zeros(mx,float)
    u = qsoln[:,1] / qsoln[:,0]
        
    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([xc[0],xc[mx-1]])
    plt.gca().set_ylim([-0.1,0.8])
    plt.plot(xc,u,'bo')
    tmp1 = "".join(("Velocity at t = ",str(time)))
    title = "".join((tmp1,"     [DoGPack]"))
    plt.title(title)
    plt.draw()
    
#-----------------------------------------------------------------------------#
# Now, insert the exact solution underneath it all:
#-----------------------------------------------------------------------------#
# location of initial Riemann problem
    xs = 0.5 

    p_star = 1.693387213839243e+00
    u_star = 4.641116216606624e-01

    # two speeds for the left rarefaction:
    sll    = -1.183215956619923e+00
    slr    = -6.262820106271283e-01

    # intermediate values for the density
    rho_sl = 1.993965770327274e+00
    rho_sr = 1.450638447387611e+00

    # shock speed for right hand shock:
    s3 = 1.494009590533840e+00

    plt.figure(1)

    W_l = [3., 0., 3.]
    rho_l = 3.
    rho_r = 1.

    # Region 1:
    plt.plot( [xc[0], xs+time*sll], [rho_l, rho_l],'r-')

    #===============================================================================
    # Rarefaction:
    def LeftFan( x, t, W_l ):
        """ Data in the case of a left rarefaction fan.

        User is responsible for knowing correct bounds of x to supply.
        """

        rho_l     = W_l[0]
        u_l       = W_l[1] 
        press_l   = W_l[2]

        c_l = np.sqrt( abs( gamma*press_l / rho_l ) )

        rho   = rho_l   * ( 2./gp1 + (gm1/(c_l*gp1))*( u_l - (x-xs)/t ) )**(  2./gm1)
        u     = 2/gm1   * ( c_l + gm1/2 *u_l + x/t )
        press = press_l * ( 2./gp1 + (gm1/(c_l*gp1))*( u_l - (x-xs)/t ) )**(2*gamma/gm1)

        return [rho, u, press]

    if( time > 0. ):
        xsample = np.linspace( xs+time*sll, xs+time*slr, 1000 )
        rho_fan, u_fan, press_fan = LeftFan( xsample, time, W_l )
        plt.plot( xsample, rho_fan, 'r-' )

        # Region 2: 
        plt.plot( [xs+time*slr, xs+time*u_star], [rho_sl, rho_sl],'r-')

        # Region 3: 
        plt.plot( [xs+time*u_star, xs+time*u_star], [rho_sl, rho_sr],'r-')
        plt.plot( [xs+time*u_star, xs+time*s3], [rho_sr, rho_sr],'r-')

    else:
        plt.plot( [xs, xs], [rho_l, rho_r], 'r-' )


    # Region 4: 
    plt.plot( [xs+time*s3, xs+time*s3], [rho_sr, rho_r],'r-')
    plt.plot( [xs+time*s3, xc[-1] ], [rho_r, rho_r],'r-')

    plt.draw()


# output from calling riemann_solve_euler.py:
#
#   Pressure and velocity inside the entire the star region:
#       p_star = 1.693387213839243e+00
#       u_star = 4.641116216606624e-01
#    
#   Left Rarefaction: 
#       Left  foot speed  = -1.183215956619923e+00
#       Right foot speed  = -6.262820106271283e-01
#       rho_star_left     = 1.993965770327274e+00 
#    
#   Contact Discontinuity:
#       speed (ustar)     = 4.641116216606624e-01 
#    
#   Right Shock:
#       Right shock speed = 1.494009590533840e+00
#       rho_star_right    = 1.450638447387611e+00 
#    
#   In summary, we have the following values for each region:
#   ('    Region 1: W = ', [3.0, 0, 3.0])
#   ('    Region 2: W = ', [1.993965770327274, 0.4641116216606624, 1.6933872138392425])
#   ('    Region 3: W = ', [1.4506384473876108, 0.4641116216606624, 1.6933872138392425])
#   ('    Region 4: W = ', [1.0, 0, 1.0])


#-----------------------------------------------------------------------------#
