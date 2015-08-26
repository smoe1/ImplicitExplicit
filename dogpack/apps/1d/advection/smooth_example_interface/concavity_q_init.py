import sympy
from sympy import Rational

def QinitFunc( x ):

#   width = Rational(1,2)
#   if ( abs(x) > width/2 ):
#       return 0
#   else:
    width = Rational(1,2)
    return sympy.cos(sympy.pi*x/width)**6

def DQinitFunc( x ):
    return sympy.diff( QinitFunc(x), x )

plt_flag = False
if( plt_flag ):
    from sympy.plotting import plot
    x  = symbols('x')
    p1 = plot( QinitFunc(x) )


