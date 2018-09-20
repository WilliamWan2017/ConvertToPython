#!/usr/bin/env python3

import simpy
import sympy as S
import sys
from src.ode import ODE

step = 0
def NoName(env, cstate=0):

    delta = None               # None to cause failure
    #define constants
    v2=-10
    v1=30
    le=1
    #define continuous variables used in this ha
    theta=0
    x=0
    y=1
    phi=1

    #define location equations
    L1_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(v1*S.sympify('cos(theta(t))')),ttol=0.01,iterations=1000)
    L1_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(v1*S.sympify('sin(theta(t))')),ttol=0.01,iterations=1000)
    L1_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(v1*S.sympify((S.sympify('tan(phi(t))')/le))),ttol=0.01,iterations=1000)
    L1_ode_phi=ODE(env,lvalue= S.sympify('diff(phi(t))'),
                                                  rvalue=S.sympify(v2),ttol=0.01,iterations=1000)
  

    #define location init value
    L1_FT=True
    L2_End_FT=False

    #define location function
    
    
    # The computations in L1
    # Returning state, delta, value, loc1_FT, loc2_FT
    def L1(x,y,theta,phi,L1_FT,L2_End_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y,S.sympify('theta(t)'): theta,S.sympify('phi(t)'): phi}
        # the edge guard take preference
        if ((y >= 1.8 and x <= 2.8) or (y <= 0.8 and x <= 2.8)):
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,x,y,theta,phi))  
            L2_End_FT=True   
            L1_FT=None                                       
            return 1, 0, x,y,theta,phi, L1_FT,L2_End_FT, curr_time
        elif not ((y >= 1.8 and x <= 2.8) or (y <= 0.8 and x <= 2.8)):
            if not L1_FT:
                x = L1_ode_x.compute(vals, curr_time-prev_time)
                y = L1_ode_y.compute(vals, curr_time-prev_time)
                theta = L1_ode_theta.compute(vals, curr_time-prev_time)
                phi = L1_ode_phi.compute(vals, curr_time-prev_time)
                L1_FT = True
            #else:
            L1_FT = False
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,x,y,theta,phi))
            #set a Maximum value for delta 
            dx,dy=9999999,9999999 
            if abs(x- 2.8) > L1_ode_x.vtol:
                dx= min( L1_ode_x.delta(vals, quanta=( 2.8- x),
                                      other_odes=[L1_ode_y,L1_ode_theta,L1_ode_phi]),dx )
            else:
                x=  2.8
                dx= 0
            if abs(y- 0.8) > L1_ode_y.vtol:
                dy= min( L1_ode_y.delta(vals, quanta=( 0.8- y),
                                      other_odes=[L1_ode_x,L1_ode_theta,L1_ode_phi]),dy )
            else:
                y=  0.8
                dy= 0
            if abs(y- 1.8) > L1_ode_y.vtol:
                dy= min( L1_ode_y.delta(vals, quanta=( 1.8- y),
                                      other_odes=[L1_ode_x,L1_ode_theta,L1_ode_phi]),dy )
            else:
                y=  1.8
                dy= 0
            Return_Delta=min(99999,dx,dy)
            return 0, Return_Delta, x,y,theta,phi, L1_FT,L2_End_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in L1')
        # Location End is end state in this example.
    def L2_End(x,y,theta,phi,L1_FT,L2_End_FT,prev_time):
        global step
        print('total steps: ', step)
        # Done
        sys.exit(1)
  # The dictionary for the switch statement.
    switch_case = {
      0:L1,
      1:L2_End
    }
    prev_time = env.now
    while(True):
        (cstate, delta, x,y,theta,phi,L1_FT,L2_End_FT, prev_time) = switch_case[cstate](x,y,theta,phi,L1_FT,L2_End_FT,
                                                            prev_time)
        # This should always be the final statement in this function
        global step
        step += 1
        yield env.timeout(delta)


def main():
    """
    """
    # Need this for this example.
    sys.setrecursionlimit(2000)
    env = simpy.Environment()
    env.process(NoName(env))
    # Run the simulation until all events in the queue are processed.
    # Make it some number to halt simulation after sometime.
    env.run(until=0.07)


if __name__ == '__main__':
    main()




   