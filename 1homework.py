#!/usr/bin/env python3
# coding: utf-8


import numpy as np
import copy
import sys
import math


# import tkinter as tk


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
        # old points
        abc_old_cp = copy.deepcopy(abc_old)
        abc_old_cp[i] = d_old
        abc_old_cp = abc_old_cp.transpose()
        det_old.append(np.linalg.det(abc_old_cp))

        # new points
        abc_new_cp = copy.deepcopy(abc_new)
        abc_new_cp[i] = d_new
        abc_new_cp = abc_new_cp.transpose()
        det_new.append(np.linalg.det(abc_new_cp))

    # check4determinant value
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

    # Calculating matrix G and H
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
    print("G inverse matrix")
    print(g_matrix_inv)
    print()
    print("H matrix")
    print(h_matrix)
    print()

    p_matrix = h_matrix.dot(g_matrix_inv)

    # Normalizing transformation matrix
    # x = 1 / P[0][0]
    # m = np.array([[x,0,0],[0,x,0],[0,0,x]])
    # P = P.dot(m)

    for i in range(n - 1):
        for j in range(n - 1):
            p_matrix[i][j] = round(p_matrix[i][j], 5)

    print("Transformation Matrix: P = H o G_inv:")
    print(p_matrix)


# ---------------------------------- DLT ----------------------------------

def dlt_algorithm(old_points, new_points):
    old_points = np.array(old_points)
    new_points = np.array(new_points)
    n = len(new_points)
    matrix_a = []
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
    print("Matrix A")
    print(matrix_a)
    print()

    s, v, d = np.linalg.svd(matrix_a)
    last_d = d[-1]

    p_matrix = []

    for i in range(3):
        col = [last_d[3 * i], last_d[3 * i + 1], last_d[3 * i + 2]]
        p_matrix.append(col)

    d = np.array(d)
    print(d)

    for i in range(3):
        for j in range(3):
            p_matrix[i][j] = round(p_matrix[i][j], 5)

    p_matrix = np.array(p_matrix)
    print("P:")
    print(p_matrix)
    print()

    return p_matrix


# ----------------------------------DLT-M ----------------------------------


def dlt_algorithm_m(old_points, new_points):
    n = len(new_points)
    old_points = np.array(old_points)
    new_points = np.array(new_points)

    cxp = sum([new_points[i][0] for i in range(n)]) / n
    cyp = sum([new_points[i][1] for i in range(n)]) / n

    for i in range(n):
        new_points[i][0] -= cxp
        new_points[i][1] -= cyp

    cx = sum([old_points[i][0] for i in range(n)]) / float(n)
    cy = sum([old_points[i][1] for i in range(n)]) / float(n)

    for i in range(n):
        old_points[i][0] -= cx
        old_points[i][1] -= cy

    lambd = sum([math.sqrt(old_points[i][0] ** 2 + old_points[i][1] ** 2) for i in range(n)]) / n
    lambdp = sum([math.sqrt(new_points[i][0] ** 2 + new_points[i][1] ** 2) for i in range(n)]) / n

    k = math.sqrt(2) / lambd
    kp = math.sqrt(2) / lambdp

    for i in range(n):
        old_points[i][0] /= k
        old_points[i][1] /= k
        new_points[i][0] /= kp
        new_points[i][1] /= kp

    matrix_pp = dlt_algorithm(old_points, new_points)
    matrix1 = np.array([[1, 0, -cx], [0, 1, -cy], [0, 0, 1]])
    matrix2 = np.array([[k, 0, 0], [0, k, 0], [0, 0, 1]])
    matrix_t = matrix2.dot(matrix1)

    matrix1p = np.array([[1, 0, -cxp], [0, 1, -cyp], [0, 0, 1]])
    matrix2p = np.array([[kp, 0, 0], [0, kp, 0], [0, 0, 1]])
    matrix_tp = matrix2p.dot(matrix1p)

    matrix_tp_inv = np.linalg.inv(matrix_tp)

    matrix_p = matrix_tp_inv.dot(matrix_pp)
    matrix_p = matrix_p.dot(matrix_t)



    print("P:")
    print(matrix_p)

    print("NEKI ISPIS")
    print(matrix_t)


def run_me():
    old_points = [[-3.0, -1.0, 1.0], [3.0, -1.0, 1.0], [1.0, 1.0, 1.0], [-1.0, 1.0, 1.0],
                  [1.0, 2.0, 3.0], [-8.0, -2.0, 1.0]]
    new_points = [[-2.0, -1.0, 1.0], [2.0, -1.0, 1.0], [2.0, 1.0, 1.0], [-2.0, 1.0, 1.0],
                  [2.0, 1.0, 4.0], [-16.0, -5.0, 4.0]]

    n = int(input("How many points you wanna input?"))
    if n == 0:
        return

    a = int(input("Which alg you wanna use: \n0:Naive\n1:DLT\n2:DLT-modificated"))

    #
    # print("Unos tacaka")
    # p = ["A", "B", "C", "D", "E", "F", "G", "H"]
    # pp = ["A'", "B'", "C'", "D'", "E'", "F'", "G'", "H'"]
    # for i in range(n):
    #     print("Unos tacke {}:".format(p[i]))
    #     x1 = int(input())
    #     x2 = int(input())
    #     x3 = int(input())
    #     old_points.append([x1, x2, x3])
    #
    # for i in range(n):
    #     print("Unos tacke {}:".format(pp[i]))
    #     x1 = int(input())
    #     x2 = int(input())
    #     x3 = int(input())
    #     new_points.append([x1, x2, x3])

    if a == 0:
        naive_algorithm(old_points, new_points)
    elif a == 1:
        dlt_algorithm(old_points, new_points)
    elif a == 2:
        dlt_algorithm_m(old_points, new_points)
    else:
        print("Ne postoji algoritam")

def main():
    run_me()


if __name__ == "__main__":
    main()
