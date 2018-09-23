#!/usr/bin/env python3

import simpy
import sympy as S
import sys
from src.ode import ODE

step = 0
def NoName(env, cstate=0):

    delta = None               # None to cause failure
    #define constants
    #define continuous variables used in this ha
    x=20

    #define location equations
    LOC2_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify('-0.25*x(t) + 125'),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify('(-0.25)')*S.sympify('x(t)'),ttol=10**-3,iterations=100,vtol=0)
  

    #define location init value
    Loc1_FT=False
    LOC2_FT=False
    _FT=False

    #define location function
    
    
    # The computations in Loc1
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc1(x,Loc1_FT,LOC2_FT,_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x}
        # the edge guard take preference
        if x<=18.5:
            
            print('%s %7.4f:%7.4f' % ( 'Loc1',curr_time,x))  
            LOC2_FT=True   
            Loc1_FT=None                                       
            return 1, 0, x, Loc1_FT,LOC2_FT,_FT, curr_time
        elif x>=18.5:
            if not Loc1_FT:
                x = Loc1_ode_x.compute(vals, curr_time-prev_time)
                Loc1_FT = True
            #else:
            Loc1_FT = False
            print('%s %7.4f:%7.4f' % ( 'Loc1',curr_time,x))
            #set a Maximum value for delta 
            dx=9999999 
            if abs(x- 18.5) > Loc1_ode_x.vtol:
                dx= min( Loc1_ode_x.delta(vals, quanta=( 18.5- x),
                                      other_odes=[]),dx )
            else:
                x=  18.5
                dx= 0
            Return_Delta=min(99999,dx)
            return 0, Return_Delta, x, Loc1_FT,LOC2_FT,_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc1')
    
    
    # The computations in LOC2
    # Returning state, delta, value, loc1_FT, loc2_FT
    def LOC2(x,Loc1_FT,LOC2_FT,_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x}
        # the edge guard take preference
        if x>=19.5:
            
            print('%s %7.4f:%7.4f' % ( 'LOC2',curr_time,x))  
            Loc1_FT=True   
            LOC2_FT=None                                       
            return 0, 0, x, Loc1_FT,LOC2_FT,_FT, curr_time
        elif x<=19.5:
            if not LOC2_FT:
                x = LOC2_ode_x.compute(vals, curr_time-prev_time)
                LOC2_FT = True
            #else:
            LOC2_FT = False
            print('%s %7.4f:%7.4f' % ( 'LOC2',curr_time,x))
            #set a Maximum value for delta 
            dx=9999999 
            if abs(x- 19.5) > LOC2_ode_x.vtol:
                dx= min( LOC2_ode_x.delta(vals, quanta=( 19.5- x),
                                      other_odes=[]),dx )
            else:
                x=  19.5
                dx= 0
            Return_Delta=min(99999,dx)
            return 1, Return_Delta, x, Loc1_FT,LOC2_FT,_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in LOC2')
    
  # The dictionary for the switch statement.
    switch_case = {
      0:Loc1,
      1:LOC2
    }
    prev_time = env.now
    while(True):
        (cstate, delta, x,Loc1_FT,LOC2_FT,_FT, prev_time) = switch_case[cstate](x,Loc1_FT,LOC2_FT,_FT,
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
    env.run(until=0.5)
    print('total steps: ', step)

if __name__ == '__main__':
    main()




   