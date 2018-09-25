#!/usr/bin/env python3

import simpy
import sympy as S
import sys
from src.ode import ODE

step = 0
def q(env, cstate=0):

    delta = None               # None to cause failure
    #define constants
    EPIKS=2.0994
    EPIKSO=2.0458
    EPIKWM=65.0
    EPITFI=0.11
    EPITHV=0.3
    EPITO2=6
    EPITS1=2.7342
    EPITS2=16.0
    EPITSI=1.8875
    EPITSO1=30.0181
    EPITSO2=0.9957
    EPITVM2=60.0
    EPITVP=1.4506
    EPITWINF=0.07
    EPITWM1=60
    EPITWP=200.0
    EPIUS=0.9087
    EPIUSO=0.65
    EPIUU=1.55
    EPIUWM=0.03
    EPI_TWM2=15.0
    jfi1=0.0
    jfi3=0.0
    jsi1=0.0
    jsi2=0.0
    stim=1.0
    #define continuous variables used in this ha
    s=0.0
    tau=0.0
    u=0.0
    v=1.0
    w=1.0

    #define location equations
    Loc0_ode_tau=ODE(env,lvalue= S.sympify('diff(tau(t))'),
                                                  rvalue=S.sympify('1.00000000000000'),ttol=10**-3,iterations=100,vtol=0)
    Loc0_ode_u=ODE(env,lvalue= S.sympify('diff(u(t))'),
                                                  rvalue=S.sympify(-jfi + stim - jsi1 - u(t)/EPITO1),ttol=10**-3,iterations=100,vtol=0)
    Loc0_ode_w=ODE(env,lvalue= S.sympify('diff(w(t))'),
                                                  rvalue=S.sympify((-wt + 1.0 - ut/EPITWINF)/(EPITWM1 + (-EPITWM1 + EPITWM2)*(1.0/(exp(((-2)*EPIKWM)*(-EPIUS + ut)) + 1)))),ttol=10**-3,iterations=100,vtol=0)
    Loc0_ode_v=ODE(env,lvalue= S.sympify('diff(v(t))'),
                                                  rvalue=S.sympify((-vt + 1.0)/EPITVM1),ttol=10**-3,iterations=100,vtol=0)
    Loc0_ode_=ODE(env,lvalue= S.sympify(S*dot),
                                                  rvalue=S.sympify((-st + 1/(exp(((-2)*EPIKS)*(-EPIUS + ut)) + 1))/EPITS1),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_tau=ODE(env,lvalue= S.sympify('diff(tau(t))'),
                                                  rvalue=S.sympify('1.00000000000000'),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_u=ODE(env,lvalue= S.sympify('diff(u(t))'),
                                                  rvalue=S.sympify(-jfi3 + stim - 1.0/((-EPITSO1 + EPITSO2)*(exp(u(t)*((-2)*EPIKSO)) + 1))),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_w=ODE(env,lvalue= S.sympify('diff(w(t))'),
                                                  rvalue=S.sympify((-w(t))/EPITWP),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_v=ODE(env,lvalue= S.sympify('diff(v(t))'),
                                                  rvalue=S.sympify((-v(t))/EPITVM2),ttol=10**-3,iterations=100,vtol=0)
    Loc1_ode_s=ODE(env,lvalue= S.sympify('diff(s(t))'),
                                                  rvalue=S.sympify((-s(t) + 1/(exp(((-2)*EPIKS)*(-EPIUS + u(t))) + 1))/EPITS2),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_tau=ODE(env,lvalue= S.sympify('diff(tau(t))'),
                                                  rvalue=S.sympify('1.00000000000000'),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_u=ODE(env,lvalue= S.sympify('diff(u(t))'),
                                                  rvalue=S.sympify(-jfi3 + stim - 1.0/((-EPITSO1 + EPITSO2)*(exp(u(t)*((-2)*EPIKSO)) + 1))),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_w=ODE(env,lvalue= S.sympify('diff(w(t))'),
                                                  rvalue=S.sympify((-w(t))/EPITWP),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_v=ODE(env,lvalue= S.sympify('diff(v(t))'),
                                                  rvalue=S.sympify((-v(t))/EPITVP),ttol=10**-3,iterations=100,vtol=0)
    Loc2_ode_s=ODE(env,lvalue= S.sympify('diff(s(t))'),
                                                  rvalue=S.sympify((-s(t) + 1/(exp(((-2)*EPIKS)*(-EPIUS + u(t))) + 1))/EPITS2),ttol=10**-3,iterations=100,vtol=0)
    loc3_ode_tau=ODE(env,lvalue= S.sympify('diff(tau(t))'),
                                                  rvalue=S.sympify('1.00000000000000'),ttol=10**-3,iterations=100,vtol=0)
    loc3_ode_u=ODE(env,lvalue= S.sympify('diff(u(t))'),
                                                  rvalue=S.sympify(stim + vt*(-EPITHV + ut)*(EPIUU - ut)/EPITFI - (-EPITSO1 + EPI_{TSO2})/(exp(((-2)*EPIKSO)*(-EPIUSO + ut)) + 1) - 1.0/EPITSO1 - st*wt/EPITSI),ttol=10**-3,iterations=100,vtol=0)
    loc3_ode_w=ODE(env,lvalue= S.sympify('diff(w(t))'),
                                                  rvalue=S.sympify((-wt)/EPITWP),ttol=10**-3,iterations=100,vtol=0)
    loc3_ode_v=ODE(env,lvalue= S.sympify('diff(v(t))'),
                                                  rvalue=S.sympify((-vt)/EPITVP),ttol=10**-3,iterations=100,vtol=0)
    loc3_ode_s=ODE(env,lvalue= S.sympify('diff(s(t))'),
                                                  rvalue=S.sympify((-st + 1/(exp(((-2)*EPIKS)*(-EPIUS + ut)) + 1))/EPITS2),ttol=10**-3,iterations=100,vtol=0)
  

    #define location init value
    Loc0_FT=False
    Loc1_FT=False
    Loc2_FT=False
    loc3_FT=False

    #define location function
    
    
    # The computations in Loc0
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc0(tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('tau(t)'): tau,S.sympify('u(t)'): u,S.sympify('w(t)'): w,S.sympify('v(t)'): v,S.sympify('(t)'): ,S.sympify('s(t)'): s}
        # the edge guard take preference
        if u>=0.006:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0-1',curr_time,tau,u,w,v,,s))  
            Loc1_FT=True   
            Loc0_FT=None                                       
            return 1, 0, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        elif u<=0.006:
            if not Loc0_FT:
                tau = Loc0_ode_tau.compute(vals, curr_time-prev_time)
                u = Loc0_ode_u.compute(vals, curr_time-prev_time)
                w = Loc0_ode_w.compute(vals, curr_time-prev_time)
                v = Loc0_ode_v.compute(vals, curr_time-prev_time)
                 = Loc0_ode_.compute(vals, curr_time-prev_time)
                s = Loc0_ode_s.compute(vals, curr_time-prev_time)
                Loc0_FT = True
            #else:
            Loc0_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc0-2',curr_time,tau,u,w,v,,s))
            #set a Maximum value for delta 
            du=9999999 
            if abs(u- 0.006) > Loc0_ode_u.vtol:
                du= min( Loc0_ode_u.delta(vals, quanta=( 0.006- u),
                                      other_odes=[Loc0_ode_tau,Loc0_ode_w,Loc0_ode_v,Loc0_ode_,Loc0_ode_s]),du )
            else:
                u=  0.006
                du= 0
            Return_Delta=min(99999,du)
            return 0, Return_Delta, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc0')
    
    
    # The computations in Loc1
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc1(tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('tau(t)'): tau,S.sympify('u(t)'): u,S.sympify('w(t)'): w,S.sympify('v(t)'): v,S.sympify('(t)'): ,S.sympify('s(t)'): s}
        # the edge guard take preference
        if u>=0.013:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc1-1',curr_time,tau,u,w,v,,s))  
            Loc2_FT=True   
            Loc1_FT=None                                       
            return 2, 0, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        elif u<=0.013:
            if not Loc1_FT:
                tau = Loc1_ode_tau.compute(vals, curr_time-prev_time)
                u = Loc1_ode_u.compute(vals, curr_time-prev_time)
                w = Loc1_ode_w.compute(vals, curr_time-prev_time)
                v = Loc1_ode_v.compute(vals, curr_time-prev_time)
                 = Loc1_ode_.compute(vals, curr_time-prev_time)
                s = Loc1_ode_s.compute(vals, curr_time-prev_time)
                Loc1_FT = True
            #else:
            Loc1_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc1-2',curr_time,tau,u,w,v,,s))
            #set a Maximum value for delta 
            du=9999999 
            if abs(u- 0.013) > Loc1_ode_u.vtol:
                du= min( Loc1_ode_u.delta(vals, quanta=( 0.013- u),
                                      other_odes=[Loc1_ode_tau,Loc1_ode_w,Loc1_ode_v,Loc1_ode_,Loc1_ode_s]),du )
            else:
                u=  0.013
                du= 0
            Return_Delta=min(99999,du)
            return 1, Return_Delta, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc1')
    
    
    # The computations in Loc2
    # Returning state, delta, value, loc1_FT, loc2_FT
    def Loc2(tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('tau(t)'): tau,S.sympify('u(t)'): u,S.sympify('w(t)'): w,S.sympify('v(t)'): v,S.sympify('(t)'): ,S.sympify('s(t)'): s}
        # the edge guard take preference
        if u>=0.3:
            
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc2-1',curr_time,tau,u,w,v,,s))  
            loc3_FT=True   
            Loc2_FT=None                                       
            return 3, 0, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        elif u<=0.3:
            if not Loc2_FT:
                tau = Loc2_ode_tau.compute(vals, curr_time-prev_time)
                u = Loc2_ode_u.compute(vals, curr_time-prev_time)
                w = Loc2_ode_w.compute(vals, curr_time-prev_time)
                v = Loc2_ode_v.compute(vals, curr_time-prev_time)
                 = Loc2_ode_.compute(vals, curr_time-prev_time)
                s = Loc2_ode_s.compute(vals, curr_time-prev_time)
                Loc2_FT = True
            #else:
            Loc2_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'Loc2-2',curr_time,tau,u,w,v,,s))
            #set a Maximum value for delta 
            du=9999999 
            if abs(u- 0.3) > Loc2_ode_u.vtol:
                du= min( Loc2_ode_u.delta(vals, quanta=( 0.3- u),
                                      other_odes=[Loc2_ode_tau,Loc2_ode_w,Loc2_ode_v,Loc2_ode_,Loc2_ode_s]),du )
            else:
                u=  0.3
                du= 0
            Return_Delta=min(99999,du)
            return 2, Return_Delta, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in Loc2')
    
    
    # The computations in loc3
    # Returning state, delta, value, loc1_FT, loc2_FT
    def loc3(tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT,prev_time):
        curr_time=env.now
        vals={S.sympify('tau(t)'): tau,S.sympify('u(t)'): u,S.sympify('w(t)'): w,S.sympify('v(t)'): v,S.sympify('(t)'): ,S.sympify('s(t)'): s}
        # the edge guard take preference
        elif u<=2.0:
            if not loc3_FT:
                tau = loc3_ode_tau.compute(vals, curr_time-prev_time)
                u = loc3_ode_u.compute(vals, curr_time-prev_time)
                w = loc3_ode_w.compute(vals, curr_time-prev_time)
                v = loc3_ode_v.compute(vals, curr_time-prev_time)
                 = loc3_ode_.compute(vals, curr_time-prev_time)
                s = loc3_ode_s.compute(vals, curr_time-prev_time)
                loc3_FT = True
            #else:
            loc3_FT = False
            print('%s %7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f:%7.4f' % ( 'loc3-2',curr_time,tau,u,w,v,,s))
            #set a Maximum value for delta 
            du=9999999 
            if abs(u- 2.0) > loc3_ode_u.vtol:
                du= min( loc3_ode_u.delta(vals, quanta=( 2.0- u),
                                      other_odes=[loc3_ode_tau,loc3_ode_w,loc3_ode_v,loc3_ode_,loc3_ode_s]),du )
            else:
                u=  2.0
                du= 0
            Return_Delta=min(99999,du)
            return 3, Return_Delta, tau,u,w,v,,s, Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, curr_time
        else:
            raise RuntimeError('Reached unreachable branch'
                               ' in loc3')
    
  # The dictionary for the switch statement.
    switch_case = {
      0:Loc0,
      1:Loc1,
      2:Loc2,
      3:loc3
    }
    prev_time = env.now
    while(True):
        (cstate, delta, tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT, prev_time) = switch_case[cstate](tau,u,w,v,,s,Loc0_FT,Loc1_FT,Loc2_FT,loc3_FT,
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
    env.process(q(env))
    # Run the simulation until all events in the queue are processed.
    # Make it some number to halt simulation after sometime.
    env.run(until=30)
    print('total steps: ', step)

if __name__ == '__main__':
    main()




   