#!/usr/bin/env python3
# coding: utf-8


import numpy as np
import copy
import sys
import math
import cv2

import time
import tkinter as tk
from PIL import ImageTk, Image
import tkinter.filedialog


# ---------------------------------- NAIVE ----------------------------------

def naive_algorithm(old_points, new_points):
    new_points = np.array(new_points)
    old_points = np.array(old_points)

    n = len(new_points)
    # ABC points
    abc_new = new_points[:n - 1]
    abc_old = old_points[:n - 1]

    # D points
    d_new = new_points[-1]
    d_old = old_points[-1]

    # List od deltas
    det_old = []
    det_new = []

    # Calculating deltas for old and new points
    for i in range(n - 1):
        # originals
        abc_old_cp = copy.deepcopy(abc_old)
        abc_old_cp[i] = d_old
        abc_old_cp = abc_old_cp.transpose()
        det_old.append(np.linalg.det(abc_old_cp))

        # images
        abc_new_cp = copy.deepcopy(abc_new)
        abc_new_cp[i] = d_new
        abc_new_cp = abc_new_cp.transpose()
        det_new.append(np.linalg.det(abc_new_cp))

    # Checking for determinant value
    for i in range(3):
        if det_new[i] == 0:
            print("Error: points are collinear!")
            sys.exit(1)

    # Calculating matrix G, G^-1 and H
    matrix_old = []
    matrix_new = []
    for i in range(n - 1):
        matrix_old.append([j * det_old[i] for j in abc_old[i]])
        matrix_new.append([j * det_new[i] for j in abc_new[i]])

    matrix_old = np.array(matrix_old)
    matrix_new = np.array(matrix_new)
    matrix_old = matrix_old.transpose()
    matrix_new = matrix_new.transpose()
    g_matrix_inv = np.linalg.inv(matrix_old)
    h_matrix = matrix_new
    p_matrix = h_matrix.dot(g_matrix_inv)

    # Normalizing transformation matrix
    # x = 1 / p_matrix[0][0]
    # m = np.array([[x,0,0],[0,x,0],[0,0,x]])
    # p_matrix = p_matrix.dot(m)

    # Round to 9 decimals
    for i in range(3):
        for j in range(3):
            p_matrix[i][j] = round(p_matrix[i][j], 9)

    return p_matrix


# ---------------------------------- DLT ----------------------------------

def dlt_algorithm(old_points, new_points):
    old_points = np.array(old_points)
    new_points = np.array(new_points)
    n = len(new_points)
    matrix_a = []

    # Matrix A[2nx9]
    for i in range(n):
        cpoint = old_points[i]
        cpointp = new_points[i]
        mini_matrix1 = [0, 0, 0, -cpointp[2] * cpoint[0], -cpointp[2] * cpoint[1], -cpointp[2] * cpoint[2],
                        cpointp[1] * cpoint[0], cpointp[1] * cpoint[1], cpointp[1] * cpoint[2]]
        mini_matrix2 = [cpointp[2] * cpoint[0], cpointp[2] * cpoint[1], cpointp[2] * cpoint[2], 0, 0, 0,
                        -cpointp[0] * cpoint[0], -cpointp[0] * cpoint[1], -cpointp[0] * cpoint[2]]
        matrix_a.append(mini_matrix1)
        matrix_a.append(mini_matrix2)

    matrix_a = np.array(matrix_a)

    # SVD(Singular Value Decomposition)
    s, v, d = np.linalg.svd(matrix_a)
    last_d = d[-1]

    # Transforming last column of matrix D to 3x3 P matrix
    p_matrix = []
    for i in range(3):
        col = [last_d[3 * i], last_d[3 * i + 1], last_d[3 * i + 2]]
        p_matrix.append(col)

    # Round to 9 decimals
    for i in range(3):
        for j in range(3):
            p_matrix[i][j] = round(p_matrix[i][j], 9)

    p_matrix = np.array(p_matrix)
    return p_matrix


# ---------------------------------- DLT-M ----------------------------------


def dlt_algorithm_m(old_points, new_points):
    n = len(new_points)
    old_points = np.array(old_points)
    new_points = np.array(new_points)

    # Making 3rd coordinate equal to 1
    for i in range(n):
        # originals
        old_points[i][0] = old_points[i][0] / old_points[i][2]
        old_points[i][1] = old_points[i][1] / old_points[i][2]
        old_points[i][2] = old_points[i][2] / old_points[i][2]

        # images
        new_points[i][0] = new_points[i][0] / new_points[i][2]
        new_points[i][1] = new_points[i][1] / new_points[i][2]
        new_points[i][2] = new_points[i][2] / new_points[i][2]

    # Center of points (G,G')
    cx = sum([old_points[i][0] for i in range(n)]) / float(n)
    cy = sum([old_points[i][1] for i in range(n)]) / float(n)
    cyp = sum([new_points[i][1] for i in range(n)]) / float(n)
    cxp = sum([new_points[i][0] for i in range(n)]) / float(n)

    # Translating points
    for i in range(n):
        new_points[i][0] -= cxp
        new_points[i][1] -= cyp
        old_points[i][0] -= cx
        old_points[i][1] -= cy

    # Homothetic transformation (S,S')
    lambd = 0
    lambdp = 0
    for i in range(n):
        lambd = lambd + math.sqrt(old_points[i][0] ** 2 + old_points[i][1] ** 2)
        lambdp = lambdp + math.sqrt(new_points[i][0] ** 2 + new_points[i][1] ** 2)

    lambd = lambd / float(n)
    lambdp = lambdp / float(n)
    if lambd == 0 or lambdp == 0:
        print("Error: points are collinear!")
        sys.exit(1)
    k = math.sqrt(2) / lambd
    kp = math.sqrt(2) / lambdp

    # Using Homothety on points
    for i in range(n):
        old_points[i][0] *= k
        old_points[i][1] *= k
        new_points[i][0] *= kp
        new_points[i][1] *= kp

    # DLT on new points = P'
    matrix_pp = dlt_algorithm(old_points, new_points)

    # Calculating matrix T = S * G
    matrix1 = np.array([[1, 0, -cx], [0, 1, -cy], [0, 0, 1]])
    matrix2 = np.array([[k, 0, 0], [0, k, 0], [0, 0, 1]])
    matrix_t = matrix2.dot(matrix1)

    # Calculating matrix T' = S' * G'
    matrix1p = np.array([[1, 0, -cxp], [0, 1, -cyp], [0, 0, 1]])
    matrix2p = np.array([[kp, 0, 0], [0, kp, 0], [0, 0, 1]])
    matrix_tp = matrix2p.dot(matrix1p)

    # (T')^-1
    matrix_tp_inv = np.linalg.inv(matrix_tp)
    # P = (T')^-1 * P' * T
    matrix_p = matrix_tp_inv.dot(matrix_pp)
    matrix_p = matrix_p.dot(matrix_t)

    # Round to 9 decimals
    for i in range(3):
        for j in range(3):
            matrix_p[i][j] = round(matrix_p[i][j], 9)

    return matrix_p


# ---------------------------------- CONSOLE ----------------------------------


def console_alghorithms():
    old_points = []
    new_points = []

    # Mocking
    # old_points = [[-3.0, -1.0, 1.0], [3.0, -1.0, 1.0], [1.0, 1.0, 1.0], [-1.0, 1.0, 1.0],
    #               [1.0, 2.0, 3.0], [-8.0, -2.0, 1.0]]
    # new_points = [[-2.0, -1.0, 1.0], [2.0, -1.0, 1.0], [2.0, 1.0, 1.0], [-2.0, 1.0, 1.0],
    #               [2.0, 1.0, 4.0], [-16.0, -5.0, 4.0]]
    # old_points = [[1.0, 1.0, 1.0], [5.0, 2.0, 1.0], [6.0, 4.0, 1.0], [-1.0, 7.0, 1.0]]
    # new_points = [[0.0, 0.0, 1.0], [10.0, 0.0, 1.0], [10.0, 5.0, 1.0], [0.0, 5.0, 1.0]]

    # Points number
    try:
        n = int(input("How many points do you want to input?\n"))
    except ValueError:
        print("Error: Not a number")
        sys.exit(1)

    # Algorithm choice
    while True:
        if n > 4:
            a = input("Which algorithm do you want to use: \n1:DLT\n2:DLT-mod\n").lower()
            if a == "1" or a == "2" or a == "dlt" or a == "dlt-mod":
                break
            else:
                print("Pick algorithm by typing [1|DLT|dlt] or [2|DLT-mod|dlt-mod]")
        elif n == 4:
            a = input("Which algorithm do you want to use: \n0:Naive\n1:DLT\n2:DLT-mod\n").lower()
            if a == "1" or a == "2" or a == "dlt" or a == "dlt-mod" or a == "0" or a == "naive":
                break
            else:
                print("Pick algorithm by typing [0|NAIVE|naive] or [1|DLT|dlt] or [2|DLT-mod|dlt-mod]")
        else:
            print("Error: Value n must be (>=4)")
            sys.exit(1)

    # Points notation, for now 26; TODO: [lowerCase] a,b,c,d...aa,ab,ac,ad...
    p = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q",
         "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
    pp = ["A'", "B'", "C'", "D'", "E'", "F'", "G'", "H'", "I'", "J'", "K'", "L'", "M'", "N'",
          "O'", "P'", "Q'", "R'", "S'", "T'", "U'", "V'", "W'", "Y'", "Z'"]

    # Taking user inputs for originals
    print("Enter [originals]points in format: x:y:z")
    for i in range(n):
        print("Enter point {}:".format(p[i]))
        entered = input().split(":")
        try:
            x1 = float(entered[0])
            x2 = float(entered[1])
            x3 = float(entered[2])
        except (ValueError, IndexError):
            print("Error: Format must be x:y:z")
            sys.exit(1)
        old_points.append([x1, x2, x3])

    # Taking user inputs for images
    print("Enter [images]points in format: x:y:z")
    for i in range(n):
        print("Enter point {}:".format(pp[i]))
        entered = input().split(":")
        try:
            x1 = float(entered[0])
            x2 = float(entered[1])
            x3 = float(entered[2])
        except (ValueError, IndexError):
            print("Error: Format must be x':y':z'")
            sys.exit(1)
        new_points.append([x1, x2, x3])

    # Using chosen algorithm
    if a == "0" or a == "naive":
        naive_p_matrix = naive_algorithm(old_points, new_points)
        print("Using NAIVE algorithm, transformation matrix P is:")
        print(naive_p_matrix)
    elif a == "1" or a == "dlt":
        dlt_p_matrix = dlt_algorithm(old_points, new_points)
        print("Using DLT algorithm, transformation matrix P is:")
        print(dlt_p_matrix)
    elif a == "2" or a == "dlt-m":
        dlt_m_p_matrix = dlt_algorithm_m(old_points, new_points)
        print("Using DLT-M algorithm, transformation matrix P is:")
        print(dlt_m_p_matrix)
    else:
        print("Typo... Try again")

    return


# ---------------------------------- GUI ----------------------------------


def GUI():

    file_types = ["*.jpeg", "*.jpg", "*.bmp", "*.png"]

    try:
        # Open file dialog
        image_file = tkinter.filedialog.askopenfilename(initialdir="~/Desktop", title="Select picture to upload:",
                                                        filetypes=(("Pictures", file_types),))

        # If u want fullscreen:
        # window.attributes("-zoomed", True)

        # Window size = picture dimensions
        ld4win = cv2.imread(image_file)
        width_wc = ld4win.shape[0]
        height_wc = ld4win.shape[1]
        scalcina = 1

        if width_wc>1080 or height_wc>1920:
            if width_wc>height_wc:
                scalcina = 720/width_wc
            else:
                scalcina = 1280/height_wc

        width_wc = int(width_wc*scalcina)
        height_wc = int(height_wc*scalcina)



    except(SystemError, AttributeError):
        print("You must choose file to upload")
        sys.exit(1)


    window = tk.Tk()
    
    window.title("Projective distortion")

    geom_string = "{}x{}".format(height_wc, width_wc)
    window.geometry(geom_string)
    canv = tk.Canvas(window, width=width_wc, height=height_wc, bg='white')
    canv.pack(expand=True, fill="both")

    points = []
    rec_id = []

    def click(eventorigin):
        x0 = float(eventorigin.x)
        y0 = float(eventorigin.y)
        points.append([x0, y0, 1.0])
        
        
        if len(points) != 1:
            rec_id.append(canv.create_rectangle(x0 - 5, y0 - 5, x0 + 5, y0 + 5, outline="#f11", width=2))
            rec_id.append(canv.create_line(points[len(points)-2][0],points[len(points)-2][1],x0,y0,fill="red"));
        else:
            rec_id.append(canv.create_rectangle(x0 - 5, y0 - 5, x0 + 5, y0 + 5, outline="#f11", width=2))
        
        print("#{}Tacka {}:{} je uneta".format(len(points), x0, y0))
        if len(points) == 4:
            canv.unbind("<Button 1>")
            canv.unbind("<Motion>")
            new_points = solve(points, image_file, width_wc, height_wc)
            ext = image_file.split('.')
            new_image = ImageTk.PhotoImage(Image.open("result."+ext[-1]))
            canv.itemconfig(1, image=new_image)
            
            #Drawing points and rectangle on final picture
            for j in new_points:
                canv.create_rectangle(j[0]-5, j[1]-5,j[0]+5,j[1]+5,outline="yellow",width=2)
            for k in range(len(new_points)-1):
                canv.create_line(new_points[k][0],new_points[k][1],new_points[k+1][0],new_points[k+1][1],fill = "yellow")
            canv.create_line(new_points[len(new_points)-1][0],new_points[len(new_points)-1][1],new_points[0][0],new_points[0][1],fill="yellow")
            
            for i in rec_id:
                canv.delete(i)
            canv.delete(text_id)
            canv.update()

    def motion(event):
        x0 = int(event.x)
        y0 = int(event.y)
        canv.itemconfig(text_id,text = "X: {} Y: {}".format(x0,y0))
        canv.update()
        
    


    canv.bind("<Button 1>", click)
    tmp = Image.open(image_file)
    tmp = tmp.resize((height_wc, width_wc), Image.ANTIALIAS)
    # img = ImageTk.PhotoImage(tmp)



    img = ImageTk.PhotoImage(tmp)
    canv.create_image(0, 0, image=img, anchor=tk.NW)

    text_id  = canv.create_text(100,10,fill="red",font="Times 15 italic bold", text="X: {} Y: {}".format(0,0))
    canv.bind("<Motion>",motion)
    

    tk.mainloop()



def sorting_points(points):
    n = len(points)
    distance = []
    for i in range(n):
        distance.append([points[i][0] ** 2 + points[i][1] ** 2, points[i]])

    ret_points = sorted(distance, key=lambda x: x[0])

    old_points = []
    for i in range(n):
        old_points.append(ret_points[i][1])

    return old_points


# added image_file as argument

def solve(points, image_file,width_w, height_w):
    # Sorting points by distance
    old_points = sorting_points(points)

    # Finding images
    new_points = copy.deepcopy(old_points)
    d = (abs(new_points[0][0] - new_points[3][0]) + abs(new_points[1][0] - new_points[2][0])) / 2
    v = (abs(new_points[0][1] - new_points[1][1]) + abs(new_points[3][1] - new_points[2][1])) / 2

    a = (new_points[0][0] + new_points[1][0]) / 2
    b = (new_points[0][1] + new_points[3][1]) / 2

    new_points[0][0] = a / 2
    new_points[0][1] = b / 2
    new_points[1][0] = new_points[0][0]
    new_points[1][1] = new_points[0][1] + v
    new_points[2][0] = new_points[0][0] + d
    new_points[2][1] = new_points[1][1]
    new_points[3][0] = new_points[2][0]
    new_points[3][1] = new_points[0][1]

    #Points that will be returned from function for purpose of drawing
    new_points_return = copy.deepcopy(new_points)

    # Sorting new points TOO
    new_points = sorting_points(new_points)

    M = dlt_algorithm_m(old_points, new_points)
    im = cv2.imread(image_file)
    imf = cv2.resize(im,(height_w,width_w))
    newim = cv2.warpPerspective(imf, M, (imf.shape[1], imf.shape[0]), flags=cv2.INTER_LINEAR)

    # geting extension of sent file
    ext = image_file.split('.')
    # saving file result.ext
    cv2.imwrite("result." + ext[-1], newim)
    
    return new_points_return


# MAIN
def main():
    GUI()
    # console_alghorithms()

    sys.exit()


# GOGOGO
if __name__ == "__main__":
    main()
