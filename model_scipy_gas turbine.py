from scipy.optimize import minimize,LinearConstraint
import numpy as math

#x1=x[0]air compressor pressure ratio 
#x2=x[1]air compressor isoentropic efficiency 
#x3=x[2]Gas turbine isoentropic effiecency
#x4=x[3]Air pre-heating temperature (t3)
#x5=x[4]Turbine Inlet Temperature (t4)

#T0=298.15 y=1.4 for compressor y=1.33 for gas turbine p0=1.013
#Dependent variables
def t2(x):
    x1=x[0]
    x2=x[1]
    return 298.15*(1+((x1**0.286-1)/(x2)))

def p2(x):
    x1=x[0]
    return 1.013*x1

def p3(x):
    return p2(x)*0.95

def p4(x):
    return p3(x)*0.95

#p6=1.066
#p5=1.099

def r(x):
    return p4(x)/1.099

def t5(x):
    x5=x[4]
    x3=x[2]
    return x5*(1-x3*(1-r(x)**(-0.248)))

#Cpg=1.17
#Cpa=1.004
def f(x):
    x5=x[4]
    x4=x[3]
    return (1.17*(x5-298.15)-1.004*(x4-298.15))/(50000-1.17*(x5-298.15))

def t6(x):
    x4=x[3]
    return t5(x)-1.004*(x4-t2(x))/1.17*(1+f(x))

def ga(x):
    x5=x[4]
    return 30000/(1.17*(1+f(x))*(x5-t5(x))-1.004*(t2(x)-298.15))

def gf(x):
    return f(x)*ga(x)

def gg(x):
    return gf(x)+ga(x)

#Gh2o*(h9-h8)=37660
def t7(x):
    return t6(x)-37660/(1.17*gg(x))

#I don't know where did 27384 come from, didn't specifiy it in the slides
def t7p(x):
    return t6(x)-27384/(1.17*gg(x))

def wc(x):
    return 1.004*ga(x)*(t2(x)-298.15)

def wt(x):
    x5=x[4]
    return 1.17*gg(x)*(x5-t5(x))

# Investment cost functions of the component
# First to To express the dependent variables used in 
# the investment cost functions of the components Ci(X), 
# some other equations have to be added:

def deltaTa(x):
    x4=x[3]
    delta=(t6(x)-t2(x))-(t5(x)-x4)
    t=(t6(x)-t2(x))/(t5(x)-x4)
    return delta/math.log(t)

#Ua=0.018
def aa(x):
    return (1.17*gg(x)*(t5(x)-t6(x)))/(0.018*deltaTa(x))

#t9=485.52
#t8=298.15
def deltaTec(x):
    delta=(t7p(x)-485.52)-(t7(x)-298.15)
    t=(t7p(x)-485.52)/(t7(x)-298.15)
    return delta/math.log(t)

def deltaTev(x):
    delta=(t6(x)-485.52)-(t7p(x)-485.52)
    t=(t6(x)-485.52)/(t7p(x)-485.52)
    return delta/math.log(t)

#Cost function for each component

#for air compressor
def c1(x):
    x1=x[0]
    x2=x[1]
    return x1*math.log(x1)*(39.5*ga(x))/(0.9-x2)

#for combustor (rb=0.95)
def c2(x):
    x5=x[4]
    return (25.6*ga(x)/(0.995-0.95))*(1+math.exp(0.018*x5-26.4))

#for gas turbine
def c3(x):
    x3=x[2]
    x5=x[4]
    par1=(266.3*gg(x))/(0.92-x3)
    par2=math.log(r(x))
    par3=1+math.exp(0.036*x5-54.4)
    return par1*par2*par3

#for air preheater
def c4(x):
    return 2290*(aa(x)**0.6)

#for HRSG
#ec=10276
#ev=27384
#steam flow=14
def c5(x):
    par1=3650*(((10276/deltaTec(x))**0.8)+((27384/deltaTev(x))**0.8))
    par2=11820*14
    par3=658*(gg(x)**1.2)
    return par1+par2+par3

#total recovery cost of a component
#CRF=0.182
#maintenance factor=1.06
#N=8000
#tao=3600 
def ztotal(x):
    return ((0.182*1.06)/(8000*3600))*(c1(x)+c2(x)+c3(x)+c4(x)+c5(x))

#operation cost
#unit cost of the fuel=4*(10**-6)
#LHV=50000
def cf(x):
    return 4*(10**-6)*gf(x)*50000

#define the objective x
def objective(x):
    return 8000*3600*(ztotal(x)+cf(x))

#Define constraints
def cont1(x):
    return x[3]-t2(x)
def cont2(x):
    return t5(x)-x[3]
def cont3(x):
    return t6(x)-t2(x)
def cont4(x):
    return t7(x)-373.15
def cont5(x):
    return x[4]-x[3]
def cont6(x):
    return t6(x)-485.52
def cont7(x):
    return t7p(x)-485.52

#initial guess
x0=[8.5,0.8,0.8,900,1450]

#boundries
b1=(5.0,10.0)
b2=(0.5,0.89)
b3=(0.5,0.91)
b4=(500,1200)
b5=(1200,1550)
bnds=(b1,b2,b3,b4,b5)

con1={'type':'ineq','fun':cont1}
con2={'type':'ineq','fun':cont2}
con3={'type':'ineq','fun':cont3}
con4={'type':'ineq','fun':cont4}
con5={'type':'ineq','fun':cont5}
con6={'type':'ineq','fun':cont6}
con7={'type':'ineq','fun':cont7}

cons=([con1,con2,con3,con4,con5,con6,con7])

solution=minimize(objective,x0,method='trust-constr',bounds=bnds,constraints=cons)
x=solution.x
print(x)
print(objective(x))
