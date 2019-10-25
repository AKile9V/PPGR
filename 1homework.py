#!/usr/bin/env python3
# coding: utf-8


import numpy as np
import copy
import tkinter as tk

import sys
import math

def naive_algorithm(old_points,new_points):
    
    new_points = np.array(new_points)
    old_points = np.array(old_points)
    
    n = len(new_points)
    #ABC points
    abc_new = new_points[:n-1]
    abc_old = old_points[:n-1]
    
    #D points
    d_new = new_points[-1]
    d_old = old_points[-1]
    
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
    
    #############Normalizing transformation matrix
    #x = 1 / P[0][0]
    #m = np.array([[x,0,0],[0,x,0],[0,0,x]])
    #P = P.dot(m)

    
    for i in range(n-1):
        for j in range(n-1):
            P[i][j] = round(P[i][j],5)
                
    print("Transformation Matrix: P = H o G_inv:")
    print(P)

##########################

def DLT_algorithm(old_points,new_points):
    
    old_points = np.array(old_points)
    new_points = np.array(new_points)
    n = len(new_points)
    matrixA = []
    for i in range(n):
        cpoint = old_points[i]
        cpointp = new_points[i]
        miniMatrix1 = [0,0,0,-cpointp[2]*cpoint[0],-cpointp[2]*cpoint[1],-cpointp[2]*cpoint[2],cpointp[1]*cpoint[0],cpointp[1]*cpoint[1],cpointp[1]*cpoint[2]]
        miniMatrix2 = [cpointp[2]*cpoint[0],cpointp[2]*cpoint[1],cpointp[2]*cpoint[2],0,0,0,-cpointp[0]*cpoint[0],-cpointp[0]*cpoint[1],-cpointp[0]*cpoint[2]]
        matrixA.append(miniMatrix1)
        matrixA.append(miniMatrix2)

    matrixA = np.array(matrixA)
    print("Matrix A")
    print(matrixA)
    print()

    s,v,d = np.linalg.svd(matrixA)
    lastD = d[-1]

    P = []

    for i in range(3):

        col = [lastD[3*i],lastD[3*i+1],lastD[3*i+2]]
        P.append(col)
    

    d = np.array(d)
    print(d)

    for i in range(3):
        for j in range(3):
            P[i][j] = round(P[i][j],5)


    P = np.array(P)
    print("P:")
    print(P)
    print()

    return P

############################---------DLT-----------##################################


def DLT_algorithm_modif(old_points,new_points):

    n = len(new_points)
    old_points = np.array(old_points)
    new_points = np.array(new_points)
    


    cxp = sum( [ new_points[i][0] for i in range(n)  ]  )/n
    cyp = sum( [ new_points[i][1] for i in range(n)  ]  )/n

    for i in range(n):
        new_points[i][0] -= cxp
        new_points[i][1] -= cyp


    cx = sum( [ old_points[i][0] for i in range(n)  ]  )/float(n)
    cy = sum( [ old_points[i][1] for i in range(n)  ]  )/float(n)

    for i in range(n):
        old_points[i][0] -= cx
        old_points[i][1] -= cy
        

    
    
    lambd = sum([math.sqrt(new_points[i][0]**2 + new_points[i][1]**2) for i in range(n)]) / n
    lambdp = sum([math.sqrt(old_points[i][0]**2 + old_points[i][1]**2) for i in range(n)]) / n

    k = math.sqrt(2)/lambd
    kp = math.sqrt(2)/lambdp

    for i in range(n):
        old_points[i][0] *= k
        old_points[i][1] *= k
        new_points[i][0] *= kp
        new_points[i][1] *= kp
        
    
    matrixPP = DLT_algorithm(old_points,new_points)
    #matrix1 = np.array([[1,0,-cx],[0,1,-cy],[0,0,1]])
    #matrix2 = np.array([[k,0,0],[0,k,0],[0,0,1]])
    #matrixT = matrix1.dot(matrix2)
    
    #matrix1p = np.array([[1,0,-cxp],[0,1,-cyp],[0,0,1]])
    #matrix2p = np.array([[kp,0,0],[0,kp,0],[0,0,1]])
    #matrixTp = matrix1p.dot(matrix2p)
    
    #matrixTp_inv = np.linalg.inv(matrixTp)

    #matrixP = matrixTp_inv.dot(matrixPP)
    #matrixP = matrixP.dot(matrixT)

    print(matrixP)





def test():
    old_points = []
    new_points = []
    
    n = int(input("How many points you wanna input?"))
    a = int(input("Which alg you wanna use: \n 0:Naive\n1:DLT\n2:DLT-modificated"))
    
    if(n == 0):
        return
    
    print("Unos tacaka")
    p = ["A","B","C","D","E","F","G","H"]
    pp = ["A'","B'","C'","D'","E'","F'","G'","H'"]
    for i in range(n):
        print("Unos tacke {}:".format(p[i]))
        x1 = int(input())
        x2 = int(input())
        x3 = int(input())
        old_points.append([x1,x2,x3])
    
    for i in range(n):
        print("Unos tacke {}:".format(pp[i]))
        x1 = int(input())
        x2 = int(input())
        x3 = int(input())
        new_points.append([x1,x2,x3])
    
    
    if(a == 0):
        naive_algorithm(old_points,new_points)
    elif a == 1:
        DLT_algorithm(old_points,new_points)
    elif a == 2:
        DLT_algorithm_modif(old_points,new_points)
    else:
        print("Ne postoji algoritam")
        
        



        



#Points (mocking)
#A = Point(-3,-1,1)
#B = Point(3,-1,1)
#C = Point(1,1,1)
#D = Point(-1,1,1)
#E = Point(1,2,3)
#F = Point(-8,-2,1)

#Ap = Point(-2,-1,1)
#Bp = Point(2,-1,1)
#Cp = Point(2,1,1)
#Dp = Point(-2,1,1)
#Ep = Point(2,1,4)
#Fp = Point(-16,-5,4)


#Mock arrays
#abcd = [A,B,C,D,E,F]
#abcdp = [Ap,Bp,Cp,Dp,Ep,Fp]
#A = Point(1,1,1)
#B = Point(5,2,1)
#C = Point(6,4,1)
#D = Point(-1,7,1)

#Ap = Point(0,0,1)
#Bp = Point(10,0,1)
#Cp = Point(10,5,1)
#Dp = Point(0,5,1)


##Mock arrays
#abcd = [A,B,C,D]
#abcdp = [Ap,Bp,Cp,Dp]





def main():
    #Initializing window
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
    selected = tk.IntVar()

    #User choose which algorithm will run
    lbl_choose = tk.Label(w,text= "Choose algorithm:")
    lbl_choose.grid(column = 10 , row =5 )


    radioFrame = tk.Frame(w)
    radioFrame.grid(column = 10, row = 6)



    radio_naive = tk.Radiobutton(radioFrame,text="Naive algorithm  ",variable = selected,value = 1).pack(side = tk.TOP)
    radio_dlt = tk.Radiobutton(radioFrame,text="DLT algorithm    ",variable = selected,value = 2).pack(side = tk.TOP)
    radio_dltMod = tk.Radiobutton(radioFrame,text="DLT-m algorithm",variable = selected,value = 3).pack(side = tk.TOP)
    
    
    def btn_clicked():


    #try:
        A = [float(txt_A1.get()),float(txt_A2.get()),float(txt_A3.get())]
        B = [float(txt_B1.get()),float(txt_B2.get()),float(txt_B3.get())]
        C = [float(txt_C1.get()),float(txt_C2.get()),float(txt_C3.get())]
        D = [float(txt_D1.get()),float(txt_D2.get()),float(txt_D3.get())]

        Ap = [float(txt_Ap1.get()),float(txt_Ap2.get()),float(txt_Ap3.get())]
        Bp = [float(txt_Bp1.get()),float(txt_Bp2.get()),float(txt_Bp3.get())]
        Cp = [float(txt_Cp1.get()),float(txt_Cp2.get()),float(txt_Cp3.get())]
        Dp = [float(txt_Dp1.get()),float(txt_Dp2.get()),float(txt_Dp3.get())]

        
        abcd = [A,B,C,D]
        abcdp = [Ap,Bp,Cp,Dp]

        if selected.get() == 1:
            naive_algorithm(abcd,abcdp)
        elif selected.get() == 2:
            print("DLT")
        elif selected.get() == 3:
            print("DLT-m")
        else:
            print("Izaberi");

        
    #except:
    #   print("NULL")
    #  print()

    def btn_cclicked():
        txt_A1.delete(0,tk.END)
        txt_A2.delete(0,tk.END)
        txt_A3.delete(0,tk.END)
        txt_Ap1.delete(0,tk.END)
        txt_Ap2.delete(0,tk.END)
        txt_Ap3.delete(0,tk.END)
        txt_B1.delete(0,tk.END)
        txt_B2.delete(0,tk.END)
        txt_B3.delete(0,tk.END)
        txt_Bp1.delete(0,tk.END)
        txt_Bp2.delete(0,tk.END)
        txt_Bp3.delete(0,tk.END)
        txt_C1.delete(0,tk.END)
        txt_C2.delete(0,tk.END)
        txt_C3.delete(0,tk.END)
        txt_Cp1.delete(0,tk.END)
        txt_Cp2.delete(0,tk.END)
        txt_Cp3.delete(0,tk.END)
        txt_D1.delete(0,tk.END)
        txt_D2.delete(0,tk.END)
        txt_D3.delete(0,tk.END)
        txt_Dp1.delete(0,tk.END)
        txt_Dp2.delete(0,tk.END)
        txt_Dp3.delete(0,tk.END)


    #Button for running choosen algorithm
    btn_run = tk.Button(w,text="Run Algorithm",command = btn_clicked)
    btn_run.grid(column=12 ,row =10)

    #Button 4 clearing all cells
    btn_clear = tk.Button(w,text="Clear",command = btn_cclicked)
    btn_clear.grid(column=12 ,row =11)


    test()

    #Running main loop of app
    w.mainloop()
    




if __name__ == "__main__":
    main()

