from scipy.optimize import minimize,LinearConstraint
import numpy  as math
import matplotlib.pyplot as plt
import matplotlib.axes as axe
import numpy as np
import seaborn as sn
from matplotlib.ticker import FuncFormatter
#Decide decesion variables
# x[0]=x1 Air compressor ratio
# x[1]=x2 The isentropic efficeincy of compressor
# x[2]=x3 The isentropic efficeincy of gas turbine
# x[3]=x4 The oulet temperture of combustion(ie. inlet of turbine) 


#Common parameter
t1=25+273.15
p1=1
p5=1
yc=1.4
yt=1.33
Wel=50000
LHV=50000
Cpg=1.17
Cpa=1.004
c11=39.5
c12=0.9
c21=25.6
c22=0.995
c23=0.018
c24=26.4
c31=266.3
c32=0.92
c33=0.036
c34=54.4
rb=0.95

#Define relationship
def t2(x):
    cal=(yc-1)/yc
    return t1*(1+(((x[0]**cal)-1)/(x[1])))

def p2(x):
    return p1*x[0]

def p4(x):
    return p2(x)*rb

def rT(x):
    return p4(x)/p5

def t5(x):
    cal=(1-yt)/yt
    return x[3]*(1-x[2]*(1-(rT(x)**cal)))

def f(x):
    up=Cpg*(x[3]-t1)-Cpa*(t2(x)-t1)
    down=LHV-Cpg*(x[3]-t1)
    return up/down

def Ga(x):
    return Wel/(Cpg*(1+f(x))*(x[3]-t5(x))-Cpa*(t2(x)-t1))

def Gf(x):
    return f(x)*Ga(x)

def Gg(x):
    return Gf(x)+Ga(x)

def Wc(x):
    return Cpa*Ga(x)*(t2(x)-t1)

def Wt(x):
    return Cpg*Gg(x)*(x[3]-t5(x))

#define cost functions
#for air compressor
def c1(x):
    return x[0]*math.log(x[0])*(c11*Ga(x))/(c12-x[1])

#for combuster
def c2(x):
    left=(c21*Ga(x)/(c22-rb))
    right=1+math.exp(c23*x[3]-c24)
    return left*right

#for gas turbine
def c3(x):
    left=c31*Gg(x)/(c32-x[2])
    right=math.log(rT(x))*(1+math.exp(c33*x[3]-c34))
    return left*right

#Define constrain
def cont1(x):
    return x[3]-t2(x)
def cont2(x):
    return x[3]-t5(x)
def cont3(x):
    return Gf(x)
def cont4(x):
    return Wt(x)-Wc(x)
con1={'type':'ineq','fun':cont1}
con2={'type':'ineq','fun':cont2}
con3={'type':'ineq','fun':cont3}
con4={'type':'ineq','fun':cont4}
cons=([con1,con2,con3,con4])

#initial guess
x0=[12,0.8,0.85,1400]

#boundries
b1=(5,25.0)
b2=(0.5,0.9)
b3=(0.5,0.92)
b4=(1000,2000)
#b5=(0.8,0.994)
bnds=(b1,b2,b3,b4)

#Define a function for thermal-econmoics solver
#cost in $/kJ
def thermaloptimization(CRF,phi,n,cost):
    def z(x):
        return ((CRF*phi)/(n*3600))*(c1(x)+c2(x)+c3(x))
    def cf(x):
        return cost*Gf(x)*LHV
    def objective(x):
        return n*3600*(z(x)+cf(x))
    solution=minimize(objective,x0,method='trust-constr',bounds=bnds,constraints=cons)
    x=solution.x
    ef=Wel/(Gf(x)*LHV)
    print(x)
    print(c1(x),c2(x),c3(x),6500*3600*z(x),6500*3600*cf(x))
    print(objective(x))
    print(ef)
    print('emision',Gf(x)*2.75)
    return ef

#consider only CO2
#Mass flow of pollution is 44/16=2.75*Gf kg/s
#98.42 euro/ton-> 105.15 dollars/ton->0.10515 dollars/kg
def environoptimization(CRF,phi,n,cost,fp,zj):
    def z(x):
        return ((CRF*phi)/(n*3600))*(c1(x)+c2(x)+c3(x))
    def cf(x):
        return cost*Gf(x)*LHV
    def env(x):
        return fp*zj*Gf(x)*2.75
    def objective(x):
        return n*3600*(z(x)+cf(x)+env(x))
    solution=minimize(objective,x0,method='trust-constr',bounds=bnds,constraints=cons)
    x=solution.x
    ef=Wel/(Gf(x)*LHV)
    print(x)
    print(c1(x),c2(x),c3(x),6500*3600*z(x),6500*3600*cf(x),6500*3600*env(x))
    print(objective(x))
    print(ef)
    print(Gf(x)*2.75)
    return ef

def CRFthermalsensitivity(CRF,phi,n,cost,rangel,rangeh):
    result=[]
    CRF_range=np.linspace(CRF*(1+rangel),CRF*(1+rangeh),5)
    range=np.linspace(rangel,rangeh,5)
    for i in CRF_range:
        totoalcost=thermaloptimization(i,phi,n,cost)
        result.append(totoalcost)
    print(result)
    figure=plt.plot()
    plt.plot(range,result,'bo-')
    plt.title('CRF Sensitivity Analysis')
    plt.xlabel('CRF increase in percentage')
    plt.ylabel('Efficiency')
    plt.show()
    return result

def fuelcostthermalsensitivity(CRF,phi,n,cost,rangel,rangeh):
    result=[]
    fuel_range=np.linspace(cost*(1+rangel),cost*(1+rangeh),5)
    range=np.linspace(rangel,rangeh,5)
    for i in fuel_range:
        totoalcost=thermaloptimization(CRF,phi,n,i)
        result.append(totoalcost)
    print(result)
    plt.plot(range,result,'ro-')
    plt.title('Fuel Cost Sensitivity Analysis')
    plt.xlabel('Fuel Cost increase in percentage')
    plt.ylabel('Efficiency')
    plt.show()
    return result



#sensitivity to the cost of CO2
def environsensitivity(CRF,phi,n,cost,fp,zj,rangel,rangeh):
    result=[]
    cost_range=np.linspace(zj*(1+rangel),zj*(1+rangeh),30)
    for i in cost_range:
        totoalcost=environoptimization(CRF,phi,n,cost,fp,i)
        result.append(totoalcost)
    print(result)
    plt.plot(cost_range,result,'go-')
    plt.title('Environmental Sensitivity Analysis')
    plt.xlabel('Cost of Carbon Dioxide [$/kg]')
    plt.ylabel('Efficiency of the plant')
    plt.show()

#thermaloptimization(0.152,1.1,6500,4*(10**(-6)))
#environoptimization(0.152,1.1,6500,4*(10**(-6)),1,0.10515)
#fuel_cost=fuelcostthermalsensitivity(0.152,1.1,6500,4*(10**(-6)),0.5,1)
#CRF_change=CRFthermalsensitivity(0.152,1.1,6500,4*(10**(-6)),0.5,1)

#range=np.linspace(0.5,1,5)
#plt.plot(range,fuel_cost,'ro-',label='Fuel cost sensitivity')
#plt.plot(range,CRF_change,'bo-',label='CRF sensitivity')
#plt.title('Comparison between CRF and Fuel cost change')
#plt.xlabel('Increase amount in percentage')
#plt.ylabel('Effieceincy')
#plt.legend()
#plt.show()
environsensitivity(0.152,1.1,6500,4*(10**(-6)),1,0.10515,-1,20)