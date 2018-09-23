#!/usr/bin/env python3

import simpy
import sympy as S
import sys
from src.ode import ODE

step = 0
def NoName(env, cstate=0):

    delta = None               # None to cause failure
    #define constants
    T1=10
    T2=10
    thM=20
    thm=5
    v1=-1.3
    v2=-2.7
    vr=10.5
    #define continuous variables used in this ha
    theta=11.5
    x=T1
    y=T2

    #define location equations
    Loc0_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc0_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc0_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(vr),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc1_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc1_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc1_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(v1),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc2_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc2_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=10**-10)
    Loc2_ode_theta=ODE(env,lvalue= S.sympify('diff(theta(t))'),
                                                  rvalue=S.sympify(v2),ttol=10**-3,iterations=100,vtol=10**-10)
  

    #define location init value
    Loc0_FT=True
    Loc1_FT=False
    Loc2_FT=False
    Loc3_FT=False

    #define location function
    
    
    # The computations in Loc0
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc0(x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y,S.sympify('theta(t)'): theta}
        # the edge guard take preference
        if theta == thM and x>=T1:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0',curr_time,x,y,theta))  
            Loc1_FT=True   
            Loc0_FT=None                                       
            return 1, 0, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif theta == thM and x<T1 and y<T2:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0',curr_time,x,y,theta))  
            Loc3_FT=True   
            Loc0_FT=None                                       
            return 3, 0, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif theta == thM and y>=T2:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0',curr_time,x,y,theta))  
            Loc2_FT=True   
            Loc0_FT=None                                       
            return 2, 0, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif theta <= thM:
            if not Loc0_FT:
                x = Loc0_ode_x.compute(vals, curr_time-prev_time)
                y = Loc0_ode_y.compute(vals, curr_time-prev_time)
                theta = Loc0_ode_theta.compute(vals, curr_time-prev_time)
                Loc0_FT = True
            #else:
            Loc0_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0',curr_time,x,y,theta))
            #set a Maximum value for delta 
            dtheta=9999999 
            if abs(theta- thM) > Loc0_ode_theta.vtol:
                dtheta= min( Loc0_ode_theta.delta(vals, quanta=( thM- theta),
                                      other_odes=[Loc0_ode_x,Loc0_ode_y]),dtheta )
            else:
                theta=  thM
                dtheta= 0
            Return_Delta=min(99999,dtheta)
            return 0, Return_Delta, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc0')
    
    
    # The computations in Loc1
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc1(x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y,S.sympify('theta(t)'): theta}
        # the edge guard take preference
        if theta==thm:
            x=0
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc1',curr_time,x,y,theta))  
            Loc0_FT=True   
            Loc1_FT=None                                       
            return 0, 0, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif theta>=thm:
            if not Loc1_FT:
                x = Loc1_ode_x.compute(vals, curr_time-prev_time)
                y = Loc1_ode_y.compute(vals, curr_time-prev_time)
                theta = Loc1_ode_theta.compute(vals, curr_time-prev_time)
                Loc1_FT = True
            #else:
            Loc1_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc1',curr_time,x,y,theta))
            #set a Maximum value for delta 
            dtheta=9999999 
            if abs(theta- thm) > Loc1_ode_theta.vtol:
                dtheta= min( Loc1_ode_theta.delta(vals, quanta=( thm- theta),
                                      other_odes=[Loc1_ode_x,Loc1_ode_y]),dtheta )
            else:
                theta=  thm
                dtheta= 0
            Return_Delta=min(99999,dtheta)
            return 1, Return_Delta, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc1')
    
    
    # The computations in Loc2
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc2(x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y,S.sympify('theta(t)'): theta}
        # the edge guard take preference
        if theta == thm:
            y=0
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc2',curr_time,x,y,theta))  
            Loc0_FT=True   
            Loc2_FT=None                                       
            return 0, 0, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif theta>=thm:
            if not Loc2_FT:
                x = Loc2_ode_x.compute(vals, curr_time-prev_time)
                y = Loc2_ode_y.compute(vals, curr_time-prev_time)
                theta = Loc2_ode_theta.compute(vals, curr_time-prev_time)
                Loc2_FT = True
            #else:
            Loc2_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc2',curr_time,x,y,theta))
            #set a Maximum value for delta 
            dtheta=9999999 
            if abs(theta- thm) > Loc2_ode_theta.vtol:
                dtheta= min( Loc2_ode_theta.delta(vals, quanta=( thm- theta),
                                      other_odes=[Loc2_ode_x,Loc2_ode_y]),dtheta )
            else:
                theta=  thm
                dtheta= 0
            Return_Delta=min(99999,dtheta)
            return 2, Return_Delta, x,y,theta, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc2')
        # Location End is end state in this example.
    def Loc3(x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        global step
        print('total steps: ', step)
        # Done
        sys.exit(1)
    
  # The dictionary for the switch statement.
    switch_case = {
      0:Loc0,
      1:Loc1,
      2:Loc2,
      3:Loc3
    }
    prev_time = env.now
    while(True):
        (cstate, delta, x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, prev_time) = switch_case[cstate](x,y,theta,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,
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
    env.run(until=30)
    print('total steps: ', step)

if __name__ == '__main__':
    main()




   