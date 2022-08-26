# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 19:35:19 2022

@author: Diana
"""
import numpy as np
import tkinter as tk
import math as m
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def f(t,z):
    return -2*z*(t-2)

def poly2(z1,z2,z3):
    return (1/12)*(23*z1-16*z2+5*z3)

def poly3(z1,z2,z3,z4):
    return (1/24)*(55*z1-59*z2+37*z3-9*z4)

def f2(t):
    return 0.25+5*m.log1p(t)*m.log1p(t)*m.sin(0.5*t)/0.1

def adams(h,b,o):
    k_runge_max=o
    k_max=int(b/h)
    z_prev=0.25
    t_prev=0
    z_massiv=[0.25]
    t_massiv=[0]
    for k in range(k_runge_max):
        kk1=f(t_prev,z_prev)
        kk2=f(t_prev+h/2,z_prev+kk1*(h/2))
        kk3=f(t_prev+h/2,z_prev+kk2*(h/2))
        kk4=f(t_prev+h,z_prev+h*kk3)
        z_delta=(h/6)*(kk1+2*kk2+2*kk3+kk4)
        z_next=z_prev+z_delta
        z_massiv.append(z_next)
        z_prev=z_next
        t_prev=t_prev+h
        t_massiv.append(t_prev)

    if o==2:

        for k in range((k_runge_max+1),k_max):
            z_delta=h*poly2(f(t_massiv[k-1],z_massiv[k-1]),f(t_massiv[k-2],z_massiv[k-2]),f(t_massiv[k-3],z_massiv[k-3]))
            z_next=z_prev+z_delta
            z_massiv.append(z_next)
            z_prev=z_next
            t_prev=t_prev+h
            t_massiv.append(t_prev)

    if o==3:
        for k in range((k_runge_max+1),k_max):
            z_delta=h*poly3(f(t_massiv[k-1],z_massiv[k-1]),f(t_massiv[k-2],z_massiv[k-2]),f(t_massiv[k-3],z_massiv[k-3]),f(t_massiv[k-4],z_massiv[k-4]))
            z_next=z_prev+z_delta
            z_massiv.append(z_next)
            z_prev=z_next
            t_prev=t_prev+h
            t_massiv.append(t_prev)
            
    z2_massiv=[]
    for i in range(len(t_massiv)):
        z2_massiv.append(f2(t_massiv[i]))
        
    intersection_points=[]
    for b in range(len(t_massiv)):
        if(abs(z_massiv[b]-z2_massiv[b])<1.5):
            intersection_points.append(b)
    
    trap_z1=0
    trap_z2=0
    for j in range(intersection_points[0],intersection_points[(len(intersection_points)-1)]+1):
        trap_z1=trap_z1+z_massiv[j]*h
        trap_z2=trap_z2+z2_massiv[j]*h
        
    trap_z1=trap_z1+z_massiv[intersection_points[0]]*h/2+z_massiv[intersection_points[(len(intersection_points)-1)]]*h/2
    trap_z2=trap_z2+z2_massiv[intersection_points[0]]*h/2+z2_massiv[intersection_points[(len(intersection_points)-1)]]*h/2
    ttk.Label(tab2, text = "S of area="+str(round(trap_z2-trap_z1,2))).grid(column = 2, row = 0,padx = 10,pady = 5)
    ttk.Label(tab2, text = "S func2="+str(round(trap_z2,2))).grid(column = 2, row = 1,padx = 10,pady = 5)
    ttk.Label(tab2, text = "S func1="+str(round(trap_z1,2))).grid(column = 2, row = 2,padx = 10,pady = 5)
    ttk.Label(tab2, text = "intersection point 1 = ("+str(t_massiv[intersection_points[0]])+","+str(z_massiv[intersection_points[0]])+")").grid(column = 3, row = 0,padx = 10,pady = 5)
    ttk.Label(tab2, text = "intersection point 2 = ("+str(round(t_massiv[intersection_points[(len(intersection_points)-1)]],2))+","+str(round(z_massiv[intersection_points[(len(intersection_points)-1)]],2))+")").grid(column = 3, row = 1,padx = 10,pady = 5)
    
    tarray = np.array(t_massiv)
    zarray = np.array(z_massiv)
    z2array = np.array(z2_massiv)
    
    data = np.column_stack([tarray, zarray,z2array])
    return data


def runge(h,b):

    k_max=int(b/h)
    z_prev=0.25
    t_prev=0
    z_massiv=[0.25]
    t_massiv=[0]
    for k in range((k_max+1)):
        kk1=f(t_prev,z_prev)
        kk2=f(t_prev+h/2,z_prev+kk1*(h/2))
        kk3=f(t_prev+h/2,z_prev+kk2*(h/2))
        kk4=f(t_prev+h,z_prev+h*kk3)
        z_delta=(h/6)*(kk1+2*kk2+2*kk3+kk4)
        z_next=z_prev+z_delta
        z_massiv.append(z_next)
        z_prev=z_next
        t_prev=t_prev+h
        t_massiv.append(t_prev)
    
    z2_massiv=[]
    for i in range(len(t_massiv)):
        z2_massiv.append(f2(t_massiv[i]))
        
    intersection_points=[]
    for b in range(len(t_massiv)):
        if(abs(z_massiv[b]-z2_massiv[b])<1.5):
            intersection_points.append(b)

    
    trap_z1=0
    trap_z2=0
    for j in range(intersection_points[0],intersection_points[(len(intersection_points)-1)]+1):
        trap_z1=trap_z1+z_massiv[j]*h
        trap_z2=trap_z2+z2_massiv[j]*h
        
    trap_z1=trap_z1+z_massiv[intersection_points[0]]*h/2+z_massiv[intersection_points[(len(intersection_points)-1)]]*h/2
    trap_z2=trap_z2+z2_massiv[intersection_points[0]]*h/2+z2_massiv[intersection_points[(len(intersection_points)-1)]]*h/2
    ttk.Label(tab1, text = "S of area="+str(round(trap_z2-trap_z1,2))).grid(column = 2, row = 0,padx = 10,pady = 5)
    ttk.Label(tab1, text = "S func2="+str(round(trap_z2,2))).grid(column = 2, row = 1,padx = 10,pady = 5)
    ttk.Label(tab1, text = "S func1="+str(round(trap_z1,2))).grid(column = 2, row = 2,padx = 10,pady = 5)
    ttk.Label(tab1, text = "intersection point 1 = ("+str(t_massiv[intersection_points[0]])+","+str(z_massiv[intersection_points[0]])+")").grid(column = 3, row = 0,padx = 10,pady = 5)
    ttk.Label(tab1, text = "intersection point 2 = ("+str(round(t_massiv[intersection_points[(len(intersection_points)-1)]],2))+","+str(round(z_massiv[intersection_points[(len(intersection_points)-1)]],2))+")").grid(column = 3, row = 1,padx = 10,pady = 5)
    tarray = np.array(t_massiv)
    zarray = np.array(z_massiv)
    z2array = np.array(z2_massiv)
 

    data = np.column_stack([tarray, zarray,z2array])
    return data

def replot_runge():
    fig = Figure(figsize = (5, 5),dpi = 100)
    hh=float(ent_hr.get())
    bb=float(ent_br.get())
    zvals=runge(hh,bb)
    plot1 = fig.add_subplot(111)
    ax = fig.gca()
    ax.set_ylim([-300, 200])
    plot1.plot(zvals[:,0],zvals[:,1],label = "z1")
    plot1.plot(zvals[:,0],zvals[:,2],label = "z2")
    canvas = FigureCanvasTkAgg(fig,master = tab1)
    canvas.draw()
    canvas.get_tk_widget().grid(column = 2, row = 6,padx = 10,pady = 5)
    
def replot_adams():
    fig = Figure(figsize = (5, 5),dpi = 100)
    hh=float(ent_ha.get())
    bb=float(ent_ba.get())
    oo=int(ent_oa.get())
    zvals=adams(hh,bb,oo)
    plot1 = fig.add_subplot(111)
    ax = fig.gca()
    ax.set_ylim([-300, 200])
    plot1.plot(zvals[:,0],zvals[:,1],label = "z1")
    plot1.plot(zvals[:,0],zvals[:,2],label = "z2")
    canvas = FigureCanvasTkAgg(fig,master = tab2)
    canvas.draw()
    canvas.get_tk_widget().grid(column = 2, row = 6,padx = 10,pady = 5)



root = tk.Tk()
root.geometry("850x1000")
root.title("Numeric methods")
tabControl = ttk.Notebook(root)
  
tab1 = tk.Frame(tabControl)
tab2 = tk.Frame(tabControl)
  
tabControl.add(tab1, text ='Runge-Kutta method')
tabControl.add(tab2, text ='Adams method')
tabControl.pack(expand = 1, fill ="both")
ttk.Label(tab1, text ="enter h").grid(column = 0, row = 1,padx = 10,pady = 5)
text_hr = tk.StringVar(value="0.1")
ent_hr=ttk.Entry(tab1,width=10,textvariable=text_hr)
ent_hr.grid(column = 1, row = 1,padx = 10,pady = 5)
ttk.Label(tab1, text ="enter b").grid(column = 0, row = 2,padx = 10,pady = 5)
text_br = tk.StringVar(value="8")
ent_br=ttk.Entry(tab1,width=10,textvariable=text_br)
ent_br.grid(column = 1, row = 2,padx = 10,pady = 5)
but_r = ttk.Button(master=tab1,command = replot_runge,text="Calculate").grid(column = 2, row = 3,padx = 10,pady = 5)

ttk.Label(tab2, text ="enter h").grid(column = 0, row = 1,padx = 10,pady = 5)
text_ha = tk.StringVar(value="0.1")
ent_ha=ttk.Entry(tab2,width=10,textvariable=text_ha)
ent_ha.grid(column = 1, row = 1,padx = 10,pady = 5)
ttk.Label(tab2, text ="enter b").grid(column = 0, row = 2,padx = 10,pady = 5)
text_ba = tk.StringVar(value="8")
ent_ba=ttk.Entry(tab2,width=10,textvariable=text_ba)
ent_ba.grid(column = 1, row = 2,padx = 10,pady = 5)
ttk.Label(tab2, text ="enter order").grid(column = 3, row = 1,padx = 10,pady = 5)
text_oa = tk.StringVar(value="3")
ent_oa=ttk.Entry(tab2,width=10,textvariable=text_oa)
ent_oa.grid(column = 3, row = 2,padx = 10,pady = 5)
but_a = ttk.Button(master=tab2,command = replot_adams,text="Calculate").grid(column = 2, row = 3,padx = 10,pady = 5)


root.mainloop()

