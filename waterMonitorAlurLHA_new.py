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
    x=0
    y=1

    #define location equations
    Loc0_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc0_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc3_ode_x=ODE(env,lvalue= S.sympify('diff(x(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
    Loc3_ode_y=ODE(env,lvalue= S.sympify('diff(y(t))'),
                                                  rvalue=S.sympify(1.00000000000000),ttol=10**-3,iterations=100,vtol=0)
  

    #define location init value
    Loc0_FT=False
    Loc1_FT=False
    Loc2_FT=False
    Loc3_FT=False

    #define location function
    
    
    # The computations in Loc0
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc0(x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y}
        # the edge guard take preference
        if y==10:
            x=0
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc0-1',curr_time,x,y))  
            Loc1_FT=True   
            Loc0_FT=None                                       
            return 1, 0, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif y<=10:
            if not Loc0_FT:
                x = Loc0_ode_x.compute(vals, curr_time-prev_time)
                y = Loc0_ode_y.compute(vals, curr_time-prev_time)
                Loc0_FT = True
            #else:
            Loc0_FT = False
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc0-2',curr_time,x,y))
            #set a Maximum value for delta 
            dy=9999999 
            if abs(y- 10) > Loc0_ode_y.vtol:
                dy= min( Loc0_ode_y.delta(vals, quanta=( 10- y),
                                      other_odes=[Loc0_ode_x]),dy )
            else:
                y=  10
                dy= 0
            Return_Delta=min(99999,dy)
            return 0, Return_Delta, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc0')
    
    
    # The computations in Loc1
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc1(x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y}
        # the edge guard take preference
        if x==2:
            
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc1-1',curr_time,x,y))  
            Loc2_FT=True   
            Loc1_FT=None                                       
            return 2, 0, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif x<=2:
            if not Loc1_FT:
                x = Loc1_ode_x.compute(vals, curr_time-prev_time)
                y = Loc1_ode_y.compute(vals, curr_time-prev_time)
                Loc1_FT = True
            #else:
            Loc1_FT = False
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc1-2',curr_time,x,y))
            #set a Maximum value for delta 
            dx=9999999 
            if abs(x- 2) > Loc1_ode_x.vtol:
                dx= min( Loc1_ode_x.delta(vals, quanta=( 2- x),
                                      other_odes=[Loc1_ode_y]),dx )
            else:
                x=  2
                dx= 0
            Return_Delta=min(99999,dx)
            return 1, Return_Delta, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc1')
    
    
    # The computations in Loc2
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc2(x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y}
        # the edge guard take preference
        if y==5:
            x=0
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc2-1',curr_time,x,y))  
            Loc3_FT=True   
            Loc2_FT=None                                       
            return 3, 0, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif y>=5:
            if not Loc2_FT:
                x = Loc2_ode_x.compute(vals, curr_time-prev_time)
                y = Loc2_ode_y.compute(vals, curr_time-prev_time)
                Loc2_FT = True
            #else:
            Loc2_FT = False
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc2-2',curr_time,x,y))
            #set a Maximum value for delta 
            dy=9999999 
            if abs(y- 5) > Loc2_ode_y.vtol:
                dy= min( Loc2_ode_y.delta(vals, quanta=( 5- y),
                                      other_odes=[Loc2_ode_x]),dy )
            else:
                y=  5
                dy= 0
            Return_Delta=min(99999,dy)
            return 2, Return_Delta, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc2')
    
    
    # The computations in Loc3
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc3(x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('x(t)'): x,S.sympify('y(t)'): y}
        # the edge guard take preference
        if x==2:
            
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc3-1',curr_time,x,y))  
            Loc0_FT=True   
            Loc3_FT=None                                       
            return 0, 0, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        elif x<=2:
            if not Loc3_FT:
                x = Loc3_ode_x.compute(vals, curr_time-prev_time)
                y = Loc3_ode_y.compute(vals, curr_time-prev_time)
                Loc3_FT = True
            #else:
            Loc3_FT = False
            print('%s %7.4f:%7.4f:%7.4f' % ( 'Loc3-2',curr_time,x,y))
            #set a Maximum value for delta 
            dx=9999999 
            if abs(x- 2) > Loc3_ode_x.vtol:
                dx= min( Loc3_ode_x.delta(vals, quanta=( 2- x),
                                      other_odes=[Loc3_ode_y]),dx )
            else:
                x=  2
                dx= 0
            Return_Delta=min(99999,dx)
            return 3, Return_Delta, x,y, Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc3')
    
  # The dictionary for the switch statement.
    switch_case = {
      0:Loc0,
      1:Loc1,
      2:Loc2,
      3:Loc3
    }
    prev_time = env.now
    while(True):
        (cstate, delta, x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT, prev_time) = switch_case[cstate](x,y,Loc0_FT,Loc1_FT,Loc2_FT,Loc3_FT,
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




   