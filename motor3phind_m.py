#-------------------------------------------------------------------------------
# Name:        motor3phaseinduc
# Purpose:
#
# Author:      j peter fitzsimmons
#
# Created:     26/11/2016
# Copyright:   (c) j peter fitzsimmons 2016
# Licence:     <MIT>
#
#use equations from Krause
#Add variable load (Tload)

#-------------------------------------------------------------------------------

from scipy import *
import numpy as np
from numpy import linalg as la
import pylab as plt
"""
This implementation is based on the textbook "Analysis of
Electric Machinery by Paul C. Krause, Chapter 4 page 153

"""

def f(x,t,wb,wr,we,Vsd,Vsq,Vrd,Vrq,Rs,Rr,Xls,Xlr,Xss,Xrr,Xm,Dk,Te,D,Tload,J,p):
    # Flux linkages
    Fqs=x[0]
    Fds=x[1]
    Fqr=x[2]
    Fdr=x[3]

    wr=x[4]
    xdot=np.zeros(5)



    xdot[0]=wb*( Vsq-(Fqs*Rs*Xrr/Dk)-we*Fds/wb + Rs*Xm*Fqr/Dk )
    xdot[1]=wb*(Vsd+(we/wb)*Fqs - Rs*Xrr*Fds/Dk  + Rs*Xm*Fdr/Dk )
    xdot[2]=wb*( Vrq-(we-wr)*Fdr/wb -Rr*Xss*Fqr/Dk + Rr*Xm*Fqs/Dk )
    xdot[3]=wb*( Vrd + Rr*Xm*Fds/Dk - Rr*Xss*Fdr/Dk +(we-wr)*Fqr/wb )

    xdot[4]=(p/(2.0*J))*(Te-D*x[4]-Tload)  # rotor speed


    return xdot
#       ########################################
def motorsimulator(fHz,Vp,Lls,Llr,Lm,Rs,Rr,D,J,p,Tload_step,Tload_ramp,Tload,tfinal,delta):

    # Stator angular electrical frequency
    we=4*np.pi*fHz/p;

    # Motor angular electrical base frequency
    wb=2*np.pi*fHz
    ang=2.0*np.pi/3.0 # 120 deg 3-phase angle


    #calculate reactances( impedances)
    Xls=wb*Lls
    Xlr=wb*Llr
    Xm=wb*Lm;
    Xml=1/(1./Xls + 1./Xm + 1./Xlr)
    Xss=Xls+Xm
    Xrr=Xlr+Xm
    Dk=Xss*Xrr-Xm**2


    # transformation matrix constants for stator voltage
    sq3d2=np.sqrt(3.0)/2.0 # sqrt(3)/2
    Kvdq= (2./3.)*np.array(([1.,.5,-.5],[0.,sq3d2,-sq3d2]))

    # transformation matrix for stator 3-phase currrents
    Kiabc=(2./3.)*np.array(([1.,0.0],[-.5,-sq3d2],[-.5,sq3d2]))

    # idq matrix idq=Kdq.Fxx
    Kidq=(1.0/Dk)*np.array( ([Xrr,0,-Xm,0],[0,Xrr,0,-Xm]),dtype=float )

    # Begin the simulation - initialize variables
    k1=np.array([0.,0.,0.,0.,0.])
    k2=np.zeros(5)
    k3=np.zeros(5)
    k4=np.zeros(5)
    x=np.zeros(5,dtype=float)
    xdot=np.zeros(5,dtype=float)
    idq=np.zeros(2)

    Te=0.0
    xa=[]
    i_abc=[]
    Vabc_ar=[]
    Vdq_ar=[]
    idq_ar=[]
    Time=[]
    torque=[]
    t=.00

##    delta=.001 # resolution for ODE solution
##    tfinal=20.0 # length of simulation

    tfinal=int(tfinal/delta)

    # Start simulation

    for N in range(0,tfinal):
        wt=2.0*np.pi*fHz*t
        fi=2.0*np.pi*fHz*t
        wr=x[4]
        wrt=wr*t
        wx=wt-wr


        #  stator voltages
        Va=Vp*np.sin(wt);
        Vb=Vp*np.sin(wt-ang);
        Vc=Vp*np.sin(wt+ang);
        Vabc=np.array(([Va,Vb,Vc]),matrix)
        Vabc_ar.append(Vabc[:])




        # dq stator voltage transformation
        Kwt=np.array(([np.cos(wt), np.sin(wt)],
                     [-np.sin(wt),np.cos(wt)]))
        Vx=np.dot(Kvdq,Vabc)
        Vdq=np.dot(Kwt,Vx)
        Vrd=0.0 # rotor voltage=0,short for squirrel cage
        Vrq=0.0 # rotor voltage
        Vsd=Vdq.item(0)
        Vsq=Vdq.item(1)
        Vdq_ar.append(Vdq[:])


        #kn=delta*xdot runge kutta 4th order
        xdot = f(x,t,wb,wr,we,Vsd,Vsq,Vrd,Vrq,Rs,Rr,Xls,Xlr,Xss,Xrr,Xm,Dk,Te,D,Tload,J,p)
        k1=delta*xdot

        k2 = delta*f(x+0.5*k1,t+0.5*delta,wb,wr,we,Vsd,Vsq,Vrd,Vrq,Rs,Rr,Xls,Xlr,Xss,Xrr,Xm,Dk,Te,D,Tload,J,p)
        k3 = delta*f(x+0.5*k2,t+0.5*delta,wb,wr,we,Vsd,Vsq,Vrd,Vrq,Rs,Rr,Xls,Xlr,Xss,Xrr,Xm,Dk,Te,D,Tload,J,p)
        k4 = delta*f(x+k3,t+delta,wb,wr,we,Vsd,Vsq,Vrd,Vrq,Rs,Rr,Xls,Xlr,Xss,Xrr,Xm,Dk,Te,D,Tload,J,p)
        #x=x+ k1 # Euler solution

        x=x+ (k1+2.0*k2+2.0*k3+k4)/6.0 # runge Kutta 4th order solution

        Fqs=x[0]
        Fds=x[1]
        Fqr=x[2]
        Fdr=x[3]

         # dq currents
        idq=np.dot(Kidq,x[0:4])
        iqs=idq[0]
        ids=idq[1]

        # abc stator currents
        Ktabc=np.array(([np.cos(wt), -np.sin(wt)],
                     [np.sin(wt),np.cos(wt)]))
        isdq=np.array([ids,iqs])
        ix=np.dot(Ktabc,isdq)
        iabc=np.dot(Kiabc,ix)
        # Torque produced by air gap
        Te=(3.0/2.0)*(p/2.0)*(1/wb)*(Fds*iqs-Fqs*ids)
        torque.append(Te)
        Tload=Tload+Tload_ramp
        if t > 10 :
            Tload_ramp=0.0
            Tload=Tload_step # 89.99% efficiency Rs=.345 Tload=13.6 N-M
        xa.append(x[:])  # must use the form x[:] not y otherwise repeats
        i_abc.append(iabc[:])
        idq_ar.append(isdq[:])
        t=t+delta
        Time.append(t)




    print "Tload=",Tload
    xa=np.array(xa,dtype=float)
    i_abc=np.array(i_abc,dtype=float)
    Vabc_ar=np.array(Vabc_ar,dtype=float)
    Vdq_ar=np.array(Vdq_ar,dtype=float)
    slip=(wb-xa[:,4])/wb
    torque=np.array(torque,dtype=float)
    power_out=torque*xa[:,4]  # element multiplication of Torq * speed
    power_in=Vdq_ar*idq_ar
    powerin_comb=power_in[:,0]+power_in[:,1]
    #Eff=power_out/powerin_comb
    print shape(power_in)
    Pin=np.sum(powerin_comb[15000:19999]/500)
    Pout=np.sum(power_out[15000:19999]/500)
    print " Input Power= ", Pin
    print " Output Power =", Pout
    print "Efficiency=", (Pout/Pin)*100,"%"


    plt.figure(1)
    plt.title("Torque vs Speed")
    plt.plot(xa[:,4],torque)
    plt.figure(2)
    plt.title("Speed")
    plt.plot(Time,xa[:,4])
    plt.figure(3)
    plt.title("Torque")
    plt.plot(Time,torque)
    plt.figure(4)
    plt.title("ia(t)")
    plt.plot(Time,i_abc[:,0])
    plt.figure(5)
    plt.title("Slip")
    plt.plot(Time,slip)
    plt.figure(6)
    plt.title("Power Out")
    plt.plot(Time,power_out)
    plt.figure(7)
    plt.title("Power In")
    plt.plot(Time,power_in)
    plt.show()


def main():
##    #Motor Parameters
##    Lls = .0477  # stator self inductance h
##    Llr =.0577   # rotor self inductance h
##    Lm= .012    # mutual inductance between stator and rotor h
##    Rs = 0.335   # stator coil resistance
##    Rr = 0.64    # rotor coil resistance
##    p=6.0       # number of poles
##    fHz=50.0    # AC frequency
##    Vp=380.0  # peak ac voltage on stator
##    D=.001 # Coefficient of drag (friction)
##    J=.28    #Moment of Inertia of rotor


 #Motor Parameters from paper by rahmann
    fHz=60.0    # AC frequency
    p=4.0       # number of poles
    wb=2*np.pi*fHz
    we=4*np.pi*fHz/p;
    Xm=26.13
    Xls=.754
    Xlr=.754
    Lls =Xls/wb  # stator self inductance h
    Llr =Xlr/wb   # rotor self inductance h
    Lm= Xm/we   # mutual inductance between stator and rotor h
    print Lls,Llr,Lm
    Rs = 0.435   # stator coil resistance
    Rr = 0.816   # rotor coil resistance

    Vp=220.0*np.sqrt(2)  # peak ac voltage on stator
    D=.00 # Coefficient of drag (friction)
    J=0.089    #Moment of Inertia of rotor

## #Motor Parameters from paper by rahmann
##    fHz=60.0    # AC frequency
##    wb=2*np.pi*fHz
##    Xm=13.04
##    Xls=.226
##    Xlr=.226
##    Lls =Xls/wb  # stator self inductance h
##    Llr =Xlr/wb   # rotor self inductance h
##    Lm= Xm/wb   # mutual inductance between stator and rotor h
##    print Lls,Llr,Lm
##    Rs = 0.029   # stator coil resistance
##    Rr = 0.816   # rotor coil resistance
##    p=4.0       # number of poles
##    Vp=2300.0*np.sqrt(2)  # peak ac voltage on stator
##    D=.00 # Coefficient of drag (friction)
##    J=63.87    #Moment of Inertia of rotor

###Motor Parameters from paper by rahmann
##    fHz=60.0    # AC frequency
##    wb=2*np.pi*fHz
##    Xm=13.04
##    Xls=.226
##    Xlr=.226
##    Lls =.8e-3  # stator self inductance h
##    Llr =.8e-3   # rotor self inductance h
##    Lm= 34.7e-3  # mutual inductance between stator and rotor h
##    print Lls,Llr,Lm
##    Rs = 0.087   # stator coil resistance
##    Rr = 0.228  # rotor coil resistance
##    p=4.0       # number of poles
##    Vp=460.0  # peak ac voltage on stator
##    D=.00 # Coefficient of drag (friction)
##    J=1.660    #Moment of Inertia of rotor

    # Tload_step is the step change in load torque applied after reaching stable speed
    Tload_step=0.0 #newton meters
    # Tload_ramp sets rate at which torque load is ramped up
    Tload_ramp=.000
    Tload=0.0  # initial load applied
    delta=.001 # resolution for ODE solution
    tfinal=1.0 # length of simulation


    motorsimulator(fHz,Vp,Lls,Llr,Lm,Rs,Rr,D,J,p,Tload_step,Tload_ramp,Tload,tfinal,delta)

if __name__ == '__main__':
    main()
