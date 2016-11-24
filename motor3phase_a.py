from scipy import *
import numpy as np
import pylab as plt
from scipy.integrate import odeint

#  Listing f function Vdqfun originally in scilab then octave now python
#  This is a 3 phase induction motor, see paper by chekz
# Modified to use my own ode solver




def odesolver( tspan, Vm, fHz,Rs,Rr,Ls,Lr,Lm,D,J,p,Tload):

      xdot=np.zeros(5);
      x=[0.0,0.0,0.0,0.0,0.0]
      delta=1e-4
      tstart=int(tspan[0]/delta)
      tfinal=int(tspan[1]/delta)
      Time=0.0
      Time_a=[]
      iabc=[]
      Torq=[]
      xa=[]
      t=0.0

      # Stator angular electrical frequency
      we=4*np.pi*fHz/p;
      # Motor angular electrical base frequency
      wb=2*np.pi*fHz

    # Reverse transform matrix id and iq to ia,ib,ic
      Mdq2abc=np.array([ [1,0],[-.5,-np.sqrt(3)/2.0],[-.5,np.sqrt(3.0)/2.]])
      Mw=np.array([ [np.cos(wb),np.sin(wb)],[np.sin(wb),np.cos(wb)] ])
      Kinv=np.dot(Mdq2abc,Mw)

      Xls=wb*Ls;
      Xlr=wb*Lr;
      Xm=wb*Lm;
      Xml=1/(1/Xls+1/Xm+1/Xlr);

      #   Differential equations:
      #print "Vdq=", Vdq
      for count in range(tstart,tfinal):
         #x=[Fqs;Fds;Fqr;Fdr;wr]
          # Stream Fqs
          Fqs=x[0];
          #   Stream Fds
          Fds=x[1];
          #Strumien Fqr
          Fqr=x[2];
          # Stream Fdr
          Fdr=x[3];
          # Rotor angular electrical speed [rad/sek]
          wr=x[4];
          fi=2*np.pi*fHz*t;
          # q and d-axis stator currents
          iqs=(Fqs-Xml*(Fqs/Xls+Fqr/Xlr))/Xls;
          ids=(Fds-Xml*(Fds/Xls+Fdr/Xlr))/Xls;
          iqrr=(Fqr-Xml*(Fqs/Xls+Fqr/Xlr))/Xlr;
          idrr=(Fdr-Xml*(Fds/Xls+Fdr/Xlr))/Xlr;
          # Transform inverse to obtain the a,b,c stator currents
          Mw=np.array([ [np.cos(fi),np.sin(fi)],[np.sin(fi),np.cos(fi)] ])
          Kinv=np.dot(Mdq2abc,Mw)
          qd=np.transpose(np.array([iqs,ids]) )
          abc=np.dot(Kinv,qd)
          iabc.append(abc[:])
          # Electrical output torque
          Te=(3.0/4.0)*p*(Fds*iqs-Fqs*ids)/wb;
          Torq.append(Te)

          # q and d-axis stator voltages
          Va=Vm*sin(2.0*np.pi*fHz*t);
          Vb=Vm*sin(2.0*np.pi*fHz*t-2.0*np.pi/3.0);
          Vc=Vm*sin(2.0*np.pi*fHz*t+2.0*np.pi/3.0);
        #  *********************************
          #MODEL with a short circuit on stator at t=3.0
##          if  t>=3 :
##              Va=0
##              Vb=0
##              Vc=0
##              #end

          #   The transformation of the voltage to the d-q
          #array([[Va],[Vb],[Vc] ])
          b=array([[Va],[Vb],[Vc] ],matrix)
          #print "shape b=",shape(b)
          c=array([[1,0.5,-0.5],[0,sqrt(3)/2,-sqrt(3)/2]],matrix )
          #print "c=",shape(c)

          Vab=(2.0/3)*dot(c,b)

          #print "shape Vab=",shape(Vab)
          d=array([[cos(fi),sin(fi)],[-sin(fi),cos(fi)] ],matrix )
          Vdq=dot(d,Vab);
          #print "shape Vdq=",shape(Vdq)
          # Vdq(1)=Vq  Vdq(2)=Vd
          Fmq=Xml*(Fqs/Xls+Fqr/Xlr);
          Fmd=Xml*(Fds/Xls+Fdr/Xlr);




          xdot[0]=wb*(Vdq.item(0)-we/wb*Fds+Rs/Xls*(Xml/Xlr*Fqr+(Xml/Xls-1)*Fqs));
          xdot[1]=wb*(Vdq.item(1)+we/wb*Fqs+Rs/Xls*(Xml/Xlr*Fdr+(Xml/Xls-1)*Fds));
          xdot[2]=wb*(-(we-wr)/wb*Fdr+Rr/Xlr*(Xml/Xls*Fqs+(Xml/Xlr-1)*Fqr));
          xdot[3]=wb*((we-wr)/wb*Fqr+Rr/Xlr*(Xml/Xls*Fds+(Xml/Xlr-1)*Fdr));
          xdot[4]=p/(2*J)*(Te-wr*D-Tload);

          x[0]=x[0]+xdot[0]*delta
          x[1]=x[1]+xdot[1]*delta
          x[2]=x[2]+xdot[2]*delta
          x[3]=x[3]+xdot[3]*delta
          x[4]=x[4]+xdot[4]*delta
          xa.append(x[:])  # must use the form x[:] not y otherwise repeats
          t=t+delta
          Time_a.append(t)
      return (xa,Torq,Time_a,iabc)

 # Simulation parameters

tspan=[0.0,5.0] # time span of the solution
fHz=50;
Vm=380;

#*********************************
# MACHINE PARAMETERS INDUCTOR
#*********************************
# Stator resistance [Ohm]
Rs=0.435;
# Rotor resistance [Ohm]
Rr=0.64;
#   Stator winding inductance [H]
Ls=0.0477;
#   Rotor winding inductance [H]
Lr=0.0577;
Lm =0.012;
# Coefficient of drag
D=0.001;
# Moment of inertia
J=0.28;
#   Number of pole pairs
p=6.0;
# Stator angular electrical frequency
we=4*np.pi*fHz/p;
# Motor angular electrical base frequency
wb=2*np.pi*fHz
#   Load torque
Tload=0.0;
# Reactances
Xls=wb*Ls;
Xlr=wb*Lr;
Xm=wb*Lm;
Xml=1/(1/Xls+1/Xm+1/Xlr);

#   Calling the function odesolver
x,Torque,t,iabc=odesolver(tspan,Vm,fHz,Rs,Rr,Ls,Lr,Lm,D,J,p,Tload)
x=np.array(x,dtype=float)
iabc=np.array(iabc,dtype=float)
print shape(iabc)
print shape(x)
print shape(x[:,4])

plt.figure(1)
plt.plot(t,x[:,4])
plt.xlabel("Seconds")
plt.ylabel("radians/sec")
plt.title("Speed (Rotor)")
plt.figure(2)
plt.plot(t,Torque)
plt.title("Torque Newton-Meters")
plt.xlabel("Seconds")
plt.figure(21)
plt.plot(x[:,4],Torque)
plt.title("Torque vs Speed")
plt.ylabel("Newton-Meters")
plt.xlabel("Speed radians/sec")
plt.figure(3)
plt.plot(t,iabc[:,2])
plt.title("ic(t)")
plt.xlabel("Seconds")
plt.ylabel("Amps")
plt.figure(4)
plt.plot(t,iabc[:,1])
plt.title("ib(t)")
plt.xlabel("Seconds")
plt.ylabel("Amps")
plt.figure(5)
plt.plot(t,iabc[:,0])
plt.title("ia(t)")
plt.xlabel("Seconds")
plt.ylabel("Amps")
plt.show()
