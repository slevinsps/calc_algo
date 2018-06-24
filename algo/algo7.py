from math import *

def f(x):
    a0 = 5
    a1 = 3
    a2 = 2
    g = a0*x/(a1+a2*x)
    return g

def f2(x):
    a0 = 5
    a1 = 3
    a2 = 2
    g = (a1*x+a2)/a0
    return g

#print("|{:3d}|{:7.2f}|{:7.2f}|{:9.4f}\t|{:7.3e}\t|{:3d}\t |{:3d}\t|".format(k, nach, kon, nach, f(nach), it, err))

def dif_odnostor(y,kol,h, x):
    odn_pr = []
    for i in range(0, kol-1):   
        odn_pr.append((y[i+1]-y[i])/(x[i+1] - x[i]))
    return odn_pr

def def_centr(y,kol,h):
    centr_pr = []
    for i in range(1, kol-1):   
        centr_pr.append((y[i+1]-y[i-1])/(2*h))
    return centr_pr

def dif_first_last(y,kol,h):
    first = (-3*y[0] + 4*y[1] - y[2])/(2*h)
    last = (-3*y[kol-1] + 4*y[kol-2] - y[kol-3])/(-2*h)
    print("first ",4*y[1])
    return first, last

def dif_runge(y,kol,h):
    if kol < 5:
        return -1
    res = []
 
    for i in range(0, kol-2):
        pr1 = (y[i+1]-y[i])/h
        pr2 = (y[i+2]-y[i])/(2*h)
        res.append(pr1 + (pr1 - pr2)/((2)-1))
    
    return res

def main():
    kol = 7
    h = 0.5
    nach_x = 10
    nach_x_copy = nach_x

    
    x = []
    for i in range(kol):   
        x.append(nach_x)
        nach_x += h
        
    y = []
    for i in x:
        y.append(f(i))
        
    odn_pr = dif_odnostor(y,kol,h, x)
    centr_pr = def_centr(y,kol,h)
    first_pr, last_pr = dif_first_last(y,kol,h)
    runge = dif_runge(y,kol,h)
    '''print("x = ", x)
    print("y = ", y)
    print("odn_pr = ", odn_pr)
    print("centr_pr = ", centr_pr)
    print("first_pr = ", first_pr)
    print("last_pr = ", last_pr)
    print("runge = ", runge)'''

    nach_x = nach_x_copy
    x1 = []
    for i in range(kol):   
        x1.append(1/nach_x)
        nach_x += h
        
    y2 = []
    for i in x1:
        y2.append(f2(i))

    #print(x1)
    #print(y2)
    virav = dif_odnostor(y2,kol,h, x1)
    #print(virav)
    
    print("----------------------------------------------------------------------------")
    print("|   x   |    y   |  одностор. |   центр.  |  y0/yn  |   Рунге |    вырав.  |")
    print("|       |        |   произв.  |  произв.  |         |         | переменные |")
    print("----------------------------------------------------------------------------")
    print("|{:6.3f} |{:6.3f}  |{:9.3f}   |    ---    |{:8.3f} |{:8.3f} |{:8.3f}    |".format(x[0], y[0], odn_pr[0], first_pr, runge[0], virav[0]*(y[0]*y[0]/(x[0]*x[0]))))
    for i in range(1,kol-2):
        print("|{:6.3f} |{:6.3f}  |{:9.3f}   |{:8.3f}   |   ---   |{:8.3f} |{:8.3f}    |".format(x[i], y[i], odn_pr[i], centr_pr[i-1], runge[i], virav[i]*(y[i]*y[i]/(x[i]*x[i]))))
    print("|{:6.3f} |{:6.3f}  |{:9.3f}   |{:8.3f}   |   ---   |   ---   |{:8.3f}    |".format(x[kol-2], y[kol-2], odn_pr[kol-2],  centr_pr[kol-3], virav[kol-2]*(y[kol-2]*y[kol-2]/(x[kol-2]*x[kol-2]))))
    print("|{:6.3f} |{:6.3f}  |    ---     |    ---    |{:8.3f} |   ---   |     ---    |".format(x[kol-1], y[kol-1], last_pr))
    print("----------------------------------------------------------------------------")
main()
    
