#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import copy
import tkinter as tk
import sys


class Point:
    def __init__(self,x1,x2,x3):
        self._x1 = x1
        self._x2 = x2
        self._x3 = x3
    
    #returning coordinates of Point
    def coords(self):
        return [self._x1,self._x2,self._x3]
        
    def __repr__(self):
        return "(x1,x2,x3) = ({},{},{})".format(self._x1,self._x2,self._x3)


class Geometry:
    def __init__(self,points_list):
        coords = [i.coords() for i in points_list]
        self._points = np.array(coords);
    
    
    def __repr__(self):
        return str(self._points)
    
    def naive_algorithm(self,new_points_list):
        coords_list_new = [i.coords() for i in new_points_list]
        new_points = np.array(coords_list_new)
        
        n = len(new_points)
        
        #ABC points
        abc_new = new_points[:n-1]
        abc_old = self._points[:n-1]
        
        #D points
        d_new = new_points[-1]
        d_old = self._points[-1]

        
        #List od deltas
        det_old = []
        det_new = []
        
        #Calculating deltas for old and new points
        for i in range(n-1):
            #old points
            abc_old_cp = copy.deepcopy(abc_old)
            abc_old_cp[i] = d_old
            abc_old_cp = abc_old_cp.transpose()
            det_old.append(np.linalg.det(abc_old_cp))
            
            #new points
            abc_new_cp = copy.deepcopy(abc_new)
            abc_new_cp[i] = d_new
            abc_new_cp = abc_new_cp.transpose()
            det_new.append(np.linalg.det(abc_new_cp))
        
        #check4determinante value
        for i in range(3):
            if det_new[i] == 0:
                print("!!!!!!!!!!!!!!!Determinante is ZERO!!!!!!!!!!!!!!!")
                print()
                sys.exit()

        print("Deltas for old:")
        print(det_old)
        print()
        print("Deltas for new:")
        print(det_new)
        print()
        
        #Calculating matrix G and H
        matrix_old = []
        matrix_new = []
        for i in range(n-1):
            matrix_old.append([j * det_old[i] for j in abc_old[i]])
            matrix_new.append([j*det_new[i] for j in abc_new[i]])
        matrix_old = np.array(matrix_old)
        matrix_new = np.array(matrix_new)
        matrix_old = matrix_old.transpose()
        matrix_new = matrix_new.transpose()
        
        
        G_inv = np.linalg.inv(matrix_old)
        H = matrix_new
        print("G inverse matrix")
        print(G_inv)
        print()
        print("H matrix")
        print(H)
        print()
        
        
        P = H.dot(G_inv)
        
        #Normalizing transformation matrix
        x = 1 / P[0][0]
        m = np.array([[x,0,0],[0,x,0],[0,0,x]])
        P = P.dot(m)

        
        for i in range(n-1):
            for j in range(n-1):
                P[i][j] = round(P[i][j],3)
                    
        print("Transformation Matrix: P = H o G_inv:")
        print(P)

        
    def DLT_algorithm(self,new_points):
        print("DLT")
    
    def DLT_algorithm_modif(self,new_points):
        print("DLT_mod")
    

#Points (mocking)
#A = Point(-3,-1,1)
#B = Point(3,-1,1)
#C = Point(1,1,1)
#D = Point(-1,1,1)

#Ap = Point(-2,-1,1)
#Bp = Point(2,-1,1)
#Cp = Point(2,1,1)
#Dp = Point(-2,1,1)

#Mock arrays
#abcd = [A,B,C,D]
#abcdp = [Ap,Bp,Cp,Dp]


#Testing algorithm
#g = Geometry(abcd)
#g.naive_algorithm(abcdp)


#Initializing window
import tkinter as tk
w = tk.Tk()
w.geometry("600x500")
w.title("Homework")
w.grid()


#Added just for init of grid
lbl_desc = tk.Label(w,text="")
lbl_desc.grid(column = 0,row = 0)

#A - A' Column
lbl_A = tk.Label(w,text="A")
lbl_A.grid(column = 2,row = 1)

txt_A1 = tk.Entry(w,width = 3)
txt_A1.grid(column = 3,row = 1)

txt_A2 = tk.Entry(w,width = 3)
txt_A2.grid(column = 4,row = 1)

txt_A3 = tk.Entry(w,width = 3)
txt_A3.grid(column = 5,row = 1)

lbl_Ap = tk.Label(w,text="A'")
lbl_Ap.grid(column = 6,row = 1)

txt_Ap1 = tk.Entry(w,width = 3)
txt_Ap1.grid(column = 7,row = 1)

txt_Ap2 = tk.Entry(w,width = 3)
txt_Ap2.grid(column = 8,row = 1)

txt_Ap3 = tk.Entry(w,width = 3)
txt_Ap3.grid(column = 9,row = 1)

#B - B' Column
lbl_B = tk.Label(w,text="B")
lbl_B.grid(column = 2,row = 2)

txt_B1 = tk.Entry(w,width = 3)
txt_B1.grid(column = 3,row = 2)

txt_B2 = tk.Entry(w,width = 3)
txt_B2.grid(column = 4,row = 2)

txt_B3 = tk.Entry(w,width = 3)
txt_B3.grid(column = 5,row = 2)

lbl_Bp = tk.Label(w,text="B'")
lbl_Bp.grid(column = 6,row = 2)

txt_Bp1 = tk.Entry(w,width = 3)
txt_Bp1.grid(column = 7,row = 2)

txt_Bp2 = tk.Entry(w,width = 3)
txt_Bp2.grid(column = 8,row = 2)

txt_Bp3 = tk.Entry(w,width = 3)
txt_Bp3.grid(column = 9,row = 2)

#C - C' Column
lbl_C = tk.Label(w,text="C")
lbl_C.grid(column = 2,row = 3)

txt_C1 = tk.Entry(w,width = 3)
txt_C1.grid(column = 3,row = 3)

txt_C2 = tk.Entry(w,width = 3)
txt_C2.grid(column = 4,row = 3)

txt_C3 = tk.Entry(w,width = 3)
txt_C3.grid(column = 5,row = 3)

lbl_Cp = tk.Label(w,text="C'")
lbl_Cp.grid(column = 6,row = 3)

txt_Cp1 = tk.Entry(w,width = 3)
txt_Cp1.grid(column = 7,row = 3)

txt_Cp2 = tk.Entry(w,width = 3)
txt_Cp2.grid(column = 8,row = 3)

txt_Cp3 = tk.Entry(w,width = 3)
txt_Cp3.grid(column = 9,row = 3)

#D - D' Column
lbl_D = tk.Label(w,text="D")
lbl_D.grid(column = 2,row = 4)

txt_D1 = tk.Entry(w,width = 3)
txt_D1.grid(column = 3,row = 4)

txt_D2 = tk.Entry(w,width = 3)
txt_D2.grid(column = 4,row = 4)

txt_D3 = tk.Entry(w,width = 3)
txt_D3.grid(column = 5,row = 4)

lbl_Dp = tk.Label(w,text="D'")
lbl_Dp.grid(column = 6,row = 4)

txt_Dp1 = tk.Entry(w,width = 3)
txt_Dp1.grid(column = 7,row = 4)

txt_Dp2 = tk.Entry(w,width = 3)
txt_Dp2.grid(column = 8,row = 4)

txt_Dp3 = tk.Entry(w,width = 3)
txt_Dp3.grid(column = 9,row = 4)


#Currently selected chooice
selected = 0;

#User choose which algorithm will run
lbl_choose = tk.Label(w,text= "Choose algorithm:")
lbl_choose.grid(column = 10 , row =5 )



radio_naive = tk.Radiobutton(w,text="Naive algorithm",variable = selected,value = 1)
radio_naive.grid(column =10, row =6 )
radio_dlt = tk.Radiobutton(w,text="DLT algorithm",variable = selected,value = 2)
radio_dlt.grid(column =10, row =7 )
radio_dltMod = tk.Radiobutton(w,text="DLT-mod algorithm",variable = selected,value = 3)
radio_dltMod.grid(column =10, row =8 )







def btn_clicked():

    #TODO: ako nista se ne unese obrardi gresku
    A = Point(float(txt_A1.get()),float(txt_A2.get()),float(txt_A3.get()))
    B = Point(float(txt_B1.get()),float(txt_B2.get()),float(txt_B3.get()))
    C = Point(float(txt_C1.get()),float(txt_C2.get()),float(txt_C3.get()))
    D = Point(float(txt_D1.get()),float(txt_D2.get()),float(txt_D3.get()))

    Ap = Point(float(txt_Ap1.get()),float(txt_Ap2.get()),float(txt_Ap3.get()))
    Bp = Point(float(txt_Bp1.get()),float(txt_Bp2.get()),float(txt_Bp3.get()))
    Cp = Point(float(txt_Cp1.get()),float(txt_Cp2.get()),float(txt_Cp3.get()))
    Dp = Point(float(txt_Dp1.get()),float(txt_Dp2.get()),float(txt_Dp3.get()))

    abcd = [A,B,C,D]
    abcdp = [Ap,Bp,Cp,Dp]

    g = Geometry(abcd)
    #TODO: treba if u zavisnosti od selected da odradi odredjeni alg
    g.naive_algorithm(abcdp)





#Button for running choosen algorithm
btn_run = tk.Button(w,text="Run Algorithm",command = btn_clicked)
btn_run.grid(column=12 ,row =9)

#Running main loop of app
w.mainloop()

