#!/usr/bin/env python3

import simpy
import sympy as S
import sys
from src.ode import ODE

step = 0
def HA1(env, cstate=0):

    delta = None               # None to cause failure
    #define constants
    v2=-10
    v3=10
    le=1
    v1=30
    v4=5
    #define continuous variables used in this ha
    y=1
    x=0
    theta=0
    phi=1

    #define location equations
    L2_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(v3*S.sympify('cos(theta(t))')),ttol=0.01,iterations=1000)
    L2_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(v3*S.sympify('sin(theta(t))')),ttol=0.01,iterations=1000)
    L2_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(v3*S.sympify((S.sympify('tan(phi(t))')/le))),ttol=0.01,iterations=1000)
    L2_ode_phi=ODE(env,lvalue= S.sympify('diff(phi(t))'),
                                                  rvalue=S.sympify('v4'),ttol=0.01,iterations=1000)
    Initial0_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(v1*S.sympify('cos(theta(t))')),ttol=0.01,iterations=1000)
    Initial0_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(v1*S.sympify('sin(theta(t))')),ttol=0.01,iterations=1000)
    Initial0_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(v1*S.sympify((S.sympify('tan(phi(t))')/le))),ttol=0.01,iterations=1000)
    Initial0_ode_phi=ODE(env,lvalue= S.sympify('diff(phi(t))'),
                                                  rvalue=S.sympify('v2'),ttol=0.01,iterations=1000)
  

    #define location init value
    L2_FT=False
    LocationEnd_FT=False
    Initial0_FT=True

    #define location function
    
    
    # The computations in Initial0
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Initial0(y,x,theta,phi,L2_FT,LocationEnd_FT,Initial0_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('y(t)'): y,S.sympify('x(t)'): x,S.sympify('theta(t)'): theta,S.sympify('phi(t)'): phi}
        # the edge guard take preference
        if (((y>=0.8) and (x>=3.2)) or ((y<=-0.4) and (x<=2.8))):
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))  
            L2_FT=True   
            Initial0_FT=None                                       
            return 1, 0, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        if x>3:
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))  
            LocationEnd_FT=True   
            Initial0_FT=None                                       
            return 2, 0, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        elif not (((y>=0.8) and (x>=3.2)) or ((y<=-0.4) and (x<=2.8))):
            if not Initial0_FT:
                y = Initial0_ode_y.compute(vals, curr_time-prev_time)
                x = Initial0_ode_x.compute(vals, curr_time-prev_time)
                theta = Initial0_ode_theta.compute(vals, curr_time-prev_time)
                phi = Initial0_ode_phi.compute(vals, curr_time-prev_time)
                Initial0_FT = True
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))
            #set a Maximum value for delta 
            dx,dy=9999999,9999999 
            if abs(x- 2.8) > Initial0_ode_x.vtol:
                dx= min( Initial0_ode_x.delta(vals, quanta=( 2.8- x),
                                      other_odes=[Initial0_ode_y,Initial0_ode_theta,Initial0_ode_phi]),dx )
            else:
                x=  2.8
                dx= 0
            if abs(y- 0.8) > Initial0_ode_y.vtol:
                dy= min( Initial0_ode_y.delta(vals, quanta=( 0.8- y),
                                      other_odes=[Initial0_ode_x,Initial0_ode_theta,Initial0_ode_phi]),dy )
            else:
                y=  0.8
                dy= 0
            if abs(y- 1.8) > Initial0_ode_y.vtol:
                dy= min( Initial0_ode_y.delta(vals, quanta=( 1.8- y),
                                      other_odes=[Initial0_ode_x,Initial0_ode_theta,Initial0_ode_phi]),dy )
            else:
                y=  1.8
                dy= 0
            Return_Delta=min(99999,dx,dy)
            return 0, Return_Delta, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Initial0')
    
    
    # The computations in L2
    # Returning state, delta, value, loc1_FT, loc2_FT
    def L2(y,x,theta,phi,L2_FT,LocationEnd_FT,Initial0_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('y(t)'): y,S.sympify('x(t)'): x,S.sympify('theta(t)'): theta,S.sympify('phi(t)'): phi}
        # the edge guard take preference
        if not (((y>=0.8) and (x>=3.2)) or ((y<=-0.4) and (x<=2.8))):
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))  
            Initial0_FT=True   
            L2_FT=None                                       
            return 0, 0, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        if x>3:
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))  
            LocationEnd_FT=True   
            L2_FT=None                                       
            return 2, 0, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        elif (((y>=0.8) and (x>=3.2)) or ((y<=-0.4) and (x<=2.8))):
            if not L2_FT:
                y = L2_ode_y.compute(vals, curr_time-prev_time)
                x = L2_ode_x.compute(vals, curr_time-prev_time)
                theta = L2_ode_theta.compute(vals, curr_time-prev_time)
                phi = L2_ode_phi.compute(vals, curr_time-prev_time)
                L2_FT = True
            print('%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % (curr_time,y,x,theta,phi))
            #set a Maximum value for delta 
             
            Return_Delta=min(99999,0)
            return 1, Return_Delta, y,x,theta,phi, L2_FT,LocationEnd_FT,Initial0_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in L2')
        # Location End is end state in this example.
    def LocationEnd(y,x,theta,phi,L2_FT,LocationEnd_FT,Initial0_FT,prev_time):
        global step
        print('total steps: ', step)
        # Done
        sys.exit(1)
  # The dictionary for the switch statement.
    switch_case = {
      1:L2,
      2:LocationEnd,
      0:Initial0
    }
    prev_time = env.now
    while(True):
        (cstate, delta, y,x,theta,phi,L2_FT,LocationEnd_FT,Initial0_FT, prev_time) = switch_case[cstate](y,x,theta,phi,L2_FT,LocationEnd_FT,Initial0_FT,
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
    env.process(HA1(env))
    # Run the simulation until all events in the queue are processed.
    # Make it some number to halt simulation after sometime.
    env.run(until=0.07)


if __name__ == '__main__':
    main()




   