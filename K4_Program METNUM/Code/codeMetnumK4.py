#!/usr/bin/env python
# coding: utf-8

# In[49]:


#import library yang digunakan
import sys
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
from numpy import float32, single, double, array, zeros, diag, diagflat, dot
from pprint import pprint
from scipy.linalg import solve
from math import sin
plt.style.use('seaborn-poster')
ipy = get_ipython()
if ipy is not None:
    ipy.run_line_magic('matplotlib','inline')
    
#============================================================================================================ 
def setengah_interval (X1, X2, a, b, c, d):
    print(" =========================================================\n",
            "Akar-Akar Persamaan: Metode Setengah Interval\n",
            "=========================================================\n")
    #===Inisialisasi Persamaan===
    X1 = X1
    X2 = X2
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FXi = (a*(X1*3))+(b*(X1)*2)-(c*(X1))-d
        FXii = (a*(X2*3))+(b*(X2*2))-(c*(X2))-d
        Xt = (X1 + X2)/2
        FXt = (a*(Xt*3))+(b*(Xt*2))-(c*(Xt))-d
        if FXi * FXt > 0:
            X1 = Xt
        elif FXi * FXt < 0:
            X2 = Xt
        else:
            print("Akar Penyelesaian: ", Xt)

        #===Memeriksa Konvergensi===
        if FXt < 0:
            error = FXt * (-1)
        else:
            error = FXt
        if iterasi > 100:
            print("Angka tak hingga")
            break
        print(iterasi, "|", FXi, "|", FXii, "|", Xt, "|", FXt, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", Xt)
    print("Toleransi Error: ", error)


def interpolasi_linier (X1, a, b, c, d):
    print(" =========================================================\n",
          "Akar-Akar Persamaan: Metode Interpolasi Linier\n",
          "=========================================================\n")

    #===Inisiasi Persamaan===
    X1 = X1
    X2 = X1 + 1
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FX1 = (a*(X1*3))+(b*(X1*2))-(c*(X1))-d
        FX2 = (a*(X2*3))+(b*(X2*2))-(c*(X2))-d
        Xt = X2 - ((FX2/(FX2-FX1)))*(X2-X1)
        FXt = (a*(Xt*3))+(b*(Xt*2))-(c*(Xt))-d
        if FXt*FX1 > 0:
            X2 = Xt
            FX2 = FXt
        else:
            X1 = Xt
            FX1 = FXt 
        if FXt < 0:
            error = FXt * (-1)
        else:
            error = FXt

        #===Memeriksa Konvergensi===
        if iterasi > 1000:
            print("Angka tak hingga")
            break         
        print(iterasi, "|", FX1, "|", FX2, "|", Xt, "|", FXt, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", Xt)
    print("Toleransi Error: ", error)


def newton_raphson (X1, a, b, c, d):
    print(" =========================================================\n",
          "Akar-Akar Persamaan: Metode Newton Raphson\n",
          "=========================================================\n")

    #===Inisiasi Persamaan===
    X1 = X1
    iterasi = 0
    akar = 1
    while (akar > 0.0001):
        iterasi += 1
        Fxn = (a*(X1*3))+(b*(X1*2))-(c*(X1))-d
        Fxxn = (a*3*(X1**2))+(b*2*X1)-c
        xnp1 = X1 - (Fxn/Fxxn)
        fxnp1 = (xnp1*3)+(xnp1*2)-(3*xnp1)-3

        #===Memeriksa Konvergensi===
        Ea = ((xnp1-X1)/xnp1)*100
        if Ea < 0.0001:
            X1 = xnp1
            akar = Ea*(-1)
        else:
            akar = xnp1
            print("Nilai akar adalah: ", akar)
            print("Nilai error adalah: ", Ea)
        if iterasi > 100:
            break
        print(iterasi, "|", X1, "|", xnp1, "|", akar, "|", Ea)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", xnp1)
    print("Toleransi Error: ", akar)


def secant (X1, a, b, c, d):
    print(" =========================================================\n",
          "Akar-Akar Persamaan: Metode Secant\n",
          "=========================================================\n")

    #===Inisiasi Persamaan===
    X1 = X1
    X2 = X1 - 1
    error = 1
    iterasi = 0
    while(error > 0.0001):
        iterasi +=1
        FX1 = (a*(X1*3))+(b*(X1*2))-(c*(X1))-d
        FXmin = (a*(X2*3))+(b*(X2*2))-(c*(X2))-d
        X3 = X1 - ((FX1)*(X1-(X2)))/((FX1)-(FXmin))
        FXplus = (a*(X3*3))+(b*(X3*2))-(c*(X3))-d
        if FXplus < 0:
            error = FXplus * (-1)
        else:
            error = FXplus
        if error > 0.0001:
            X2 = X1
            X1 = X3
        else:
            print("Selesai")

        #===Memeriksa Konvergensi===
        if iterasi > 100:
            print("Angka tak hingga")
            break
        print(iterasi, "|", FX1, "|", FXmin, "|", X3, "|", FXplus, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", X3)
    print("Toleransi Error: ", error)


def iterasi (X1, a, b, c, d):
    print(" =========================================================\n",
          "Akar-Akar Persamaan: Metode Iterasi\n",
          "=========================================================\n")

    #===Inisiasi Persamaan===
    X1 = X1
    error = 1
    iterasi = 0
    while (error > 0.0001):
        iterasi +=1
        Fxn = float((a*(X1*3))+(b*(X1*2))-(c*(X1))-d)
        X2 = float(((b*(X1*2))-(c*(X1))-d)/a)**(1/3)
        Ea = ((-(X2-X1)/(X2))*100)
        if Ea < error:
            X1 = X2
            if Ea > 0:
                error = Ea
            else:
                error = Ea*(-1)
        else:
            error = Ea

        #===Memeriksa Konvergensi===#
        if iterasi > 100:
            print("Angka tak hingga")
            break
            print(iterasi, "|", X1, "|", X2, "|", E, "|", error)
    print("Jumlah Iterasi: ", iterasi)
    print("Akar persamaan adalah: ", X2)
    print("Toleransi Error: ", error)
    
def Gauss(A, f):
    A = np.array((A), dtype=float)
    f = np.array((f), dtype=float)
    n = len(f)
    for i in range(0, n - 1):  #Looping untuk kolom matriks
        if np.abs(A[i, i]) == 0:
            for k in range(i + 1, n):
                if np.abs(A[k, i]) > np.abs(A[i, i]):
                    A[[i, k]] = A[[k, i]]  #Tukar antara baris i dan k
                    f[[i, k]] = f[[k, i]]
                    break
        for j in range(i + 1, n):  #Ulangi baris yang ada di bawah diagonal untuk setiap kolom
            m = A[j, i] / A[i, i]
            A[j, :] = A[j, :] - m * A[i, :]
            f[j] = f[j] - m * f[i]
    return A, f
    
    
def GaussJordan(a,n):
    #Step1 ===> Looping untuk pengolahan metode Gauss Jordan
    print('=============Mulai Iterasi===============')
    for i in range(n):
        if a[i][i]==0:
            sys.exit('Dibagi dengan angka nol (proses tidak dapat dilanjutkan)')
        for j in range(n):
            if i !=j:
                ratio=a[j][i]/a[i][i]
                #print('posisi nol di:[',j,i,']', 'nilai rasio:',ratio)
                for k in range(n+1):
                    a[j,k]=a[j][k]-ratio*a[i][k]
                print(a)
                print('============================================')
    # Step 2 ====> Membuat semua variabel(x,y,z,...)==1
    ax=np.zeros((n,n+1))
    for i in range(n):
        for j in range(n+1):
            ax[i,j]=a[i][j]/a[i][i]
    print('===================Akhir Iterasi===================')
    return ax
    
    
def Jacobi (A,B, x_init, epsilon=1e-10, max_interation= 500):
    D = np.diag(np.diag(A))
    LU = A-D
    x = x_init
    for i in range (max_interation):
        D_inv = np.diag(1/np.diag(D))
        x_new = np.dot(D_inv, B - np.dot(LU,x))
        if np.linalg.norm(x_new-x) < epsilon:
            return x_new
            x = x_new
            return x
        
        
def GaussSeidel (A, B, x, n):
    L = np.tril(A)
    U = A - L
    for i in range(n):
        x=np.dot(np.linalg.inv(L), B - np.dot(U,x))
        print (f'Iterasi Ke-{str(i+1).zfill(3)}'),
        print (x)
    return x
    
    
def Trapesium(f, batasbawah, batasatas, pias):
    x = np.linspace(batasatas, batasbawah, pias+1)
    y = f(x)
    y_kanan = y[1:] #Batas y Kanan / atas
    y_kiri = y[:-1] #Batas y Kiri / bawah
    dx = (batasatas - batasbawah)/pias
    T = (dx/2) * np.sum(y_kanan + y_kiri)
    return T
    

def Simpson13(batasbawah, batasatas, interval):
    #Mendefinisikan fungsi
    # Perhitungan jarak
    h = (batasatas - batasbawah) / interval
    # Menjumlahkan hasil integrasi batas atas dengan batas bawah
    integral =  f(batasatas) + f(batasbawah)

    for i in range(1,interval):
        k = batasatas + i*h

        if i%2 == 0:
            integral = integral + a * f(k)
        else:
            integral = integral + b * f(k)
    # Nilai Hasil Akhir 
    integral = integral * 1/3 * h
    return integral
    

def simpson38(batasbawah, batasatas, interval):
    #mendefinisikan fungsi
    h = (batasatas - batasbawah) / interval
    integral = f(batasatas) + f(batasbawah)

    for i in range(1,interval):
        k = batasatas + i*h

        if i%2 == 0:
            integral = integral + a * f(k)
        else:
            integral = integral + b * f(k)
    integral = integral * 3/8 * h
    return integral
    
    
def Euler (a,b,c,d,h,x0,xn,y0):
    h = h
    x0 = x0
    xn = xn
    y0 = y0
    a = a
    b = b
    c = c
    d = d
    x = np.arange(x0, xn + h, h) #Numerical grid
    G = (a*(x*3))+(b*(x**2))+(c*x)+d
    f = lambda x, y: (a*(x*3))+(b*(x**2))+(c*x)+d #Persamaan Differensial Biasa
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(0, len(x) - 1):
        y[i + 1] = y[i] + h*f(x[i], y[i])
    Galat = G-y
    print(Galat)
    Judul = ("\n Grafik Pendekatan PDB dengan Metode Euler \n")
    fig = plt.figure()
    fig.patch.set_facecolor('yellow')
    plt.plot(x, y, 'b--', label='Hasil Pendekatan') 
    plt.plot(x, G, '-g', label='Hasil Analitik')
    plt.title(Judul) # Judul plot
    plt.xlabel('x') # Label sumbu x
    plt.ylabel('F(x)') # Label sumbu y
    plt.grid() # Menampilkan grid
    plt.legend(loc='lower right')
    plt.savefig(input('Masukan Nama Lokasi :')) #format contoh : D:\Cendra Boskanita\Kuliah\Semester 4\METNUM\PRAKTIKUM\Tugas Akhir Praktikum Metnum 2021\K4_Program METNUM\image\HEUN.png
                                                #DAPAT BERUPA .png / .jpeg / .pdf
    
    
def Heun (a,b,c,d,h,x0,xn,y0):
    a = a
    b = b
    c = c
    d = d
    h = h
    x0 = x0
    xn = xn
    y0 = y0
    x = np.arange(x0, xn + h, h)
    G = (a*(x*3))+(b*(x**2))+(c*x)+d
    f = lambda x, y: (a*(x*3))+(b*(x**2))+(c*x)+d #Persamaan Differensial Biasa
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(0, len(x) - 1):
        k1 = f(x[i], y[i])
        k2 = f((x[i]+h), (y[i]+(h*k1)))
        y[i+1] = y[i]+(0.5*h*(k1+k2))
    Galat = G-y
    print(Galat)
    Judul = ("\n Grafik Pendekatan PDB dengan Metode Heun \n")
    fig = plt.figure()
    fig.patch.set_facecolor('yellow')
    plt.plot(x, y, 'b--', label='Hasil Pendekatan') 
    plt.plot(x, G, '-g', label='Hasil Analitik')
    plt.title(Judul) # Judul plot
    plt.xlabel('x') # Label sumbu x
    plt.ylabel('F(x)') # Label sumbu y
    plt.grid() # Menampilkan grid
    plt.legend(loc='lower right')
    plt.savefig(input('Masukan Nama Lokasi :')) #format contoh : D:\Cendra Boskanita\Kuliah\Semester 4\METNUM\PRAKTIKUM\Tugas Akhir Praktikum Metnum 2021\K4_Program METNUM\image\HEUN.png
                                                #DAPAT BERUPA .png / .jpeg / .pdf

        #============================================================================================================    
        #============================================================================================================    
print('HALO SELAMAT DATANG!\nSebelumnya Perkenalkan, Kami dari KELOMPOK 4, Program Studi Oseanografi Angkatan 2019, Universitas Diponegoro menyediakan layanan perhitungan Metode Numerik.\n')
print('Berikut Layanan Perhitungan Metode Numerik yang Kami sediakan.\n',
      "Kode penggunaan metode numerik: \n",
      '0. Anggota\n',
      "1. Metode Setengah Interval \n",
      "2. Metode Interpolasi Linier \n",
      "3. Metode Newton-Raphson \n",
      "4. Metode Secant \n",
      "5. Metode Iterasi \n",
      "6. Metode Gauss \n",
      "7. Metode Gauss Jordan \n",
      "8. Metode Jacobi \n",
      "9. Metode Gauss Seidel \n",
      "10. Metode Trapesium \n",
      "11. Metode Simpson 1/3 \n",
      "12. Metode Simpson 3/8 \n",
      "13. Metode Euler \n",
      "14. Metode Heun \n")
setting = int(input("Masukkan kode penggunaan metode numerik: "))
if (setting == 0):
    print('Prima Riliayunda P		26050119140056	B\n'
          'Cendra Boskanita Petrova	26050119140055	A\n'
          'Amalia Sekar Ayuningtyas	26050119130135	B\n'
          'Fawasatya Bagus Kresna		26050119140038	B\n'
          'Petrik Siano OPL		26050119130125	B\n'
          'Gisela Malya Asoka Anindita	26050119130137	B\n'
          'Wahyu Erfando			26050119130131	B\n'
          'Kenichy Priyoasmoro		26050119130126	B\n'
          'Syifa Arrahmah			26050119140050	B\n'
          'Syifa Agfanita			26050119130132	A\n'
          'Baker Simamora			26050119130142	A\n'
          'Kurniawan Sandres		26050119130133	B\n')
elif (setting == 1):
    X1 = float(input("Masukkan Nilai Pertama Metode Setengah Interval: "))
    X2 = float(input("Masukkan Nilai Kedua Metode Setengah Interval: "))
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    X = setengah_interval(X1, X2, a, b, c, d)
    print(X)
elif (setting == 2):
    X1 = float(input("Masukkan Nilai Pertama Metode Interpolasi Linier: "))
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    X = interpolasi_linier(X1, a, b, c, d)
    print(X)
elif (setting == 3):
    X1 = float(input("Masukkan Nilai Pertama Metode Newton-Raphson: "))
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    X = newton_raphson (X1, a, b, c, d)
    print(X)
elif (setting == 4):
    X1 = float(input("Masukkan Nilai Pertama Metode Secant: "))
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    X = secant (X1, a, b, c, d)
    print(X)
elif (setting == 5):
    X1 = float(input("Masukkan Nilai Pertama Metode Iterasi: "))
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    X = iterasi (X1, a, b, c, d)
    print(X)
elif (setting == 6):
    a1 = float(input('Nilai Variabel A1 :'))
    a2 = float(input('Nilai Variabel A2 :'))
    a3 = float(input('Nilai Variabel A3 :'))
    a4 = float(input('Nilai Variabel A4 :'))
    b1 = float(input('Nilai Variabel B1 :'))
    b2 = float(input('Nilai Variabel B2 :'))
    b3 = float(input('Nilai Variabel B3 :'))
    b4 = float(input('Nilai Variabel B4 :'))
    c1 = float(input('Nilai Variabel C1 :'))
    c2 = float(input('Nilai Variabel C2 :'))
    c3 = float(input('Nilai Variabel C3 :'))
    c4 = float(input('Nilai Variabel C4 :'))
    d1 = float(input('Nilai Variabel D1 :'))
    d2 = float(input('Nilai Variabel D2 :'))
    d3 = float(input('Nilai Variabel D3 :'))
    d4 = float(input('Nilai Variabel D4 :'))
    e1 = float(input('Nilai Konstanta E1 :'))
    e2 = float(input('Nilai Konstanta E2 :'))
    e3 = float(input('Nilai Konstanta E3 :'))
    e4 = float(input('Nilai Konstanta E4 :'))
    f1 = (c1*d1)/(b1-a1)+ e1
    f2 = (c2*d2)/(b2-a2)+ e2
    f3 = (c3*d1)/(b3-a3)+ e3
    f4 = (c4*d1)/(b4-a4)+ e4
    A = np.array([[a1, b1, c1, d1], 
                  [a2, b2, c2, d2], 
                  [a3, b3, c3, d3], 
                  [a4, b4, c4, d4]])
    f = np.array([f1, f2, f3, f4]) #hasil Input ke fungsi
    print('A = \n%s and f = %s \n' % (A, f))
    x = np.linalg.solve(A, f)
    print('Hasil perhitungan dengan metode Gauss adalah x = \n %s \n' % x)
elif (setting == 7):
    a1 = float(input('Nilai Variabel A1 :'))
    a2 = float(input('Nilai Variabel A2 :'))
    a3 = float(input('Nilai Variabel A3 :'))
    a4 = float(input('Nilai Variabel A4 :'))
    b1 = float(input('Nilai Variabel B1 :'))
    b2 = float(input('Nilai Variabel B2 :'))
    b3 = float(input('Nilai Variabel B3 :'))
    b4 = float(input('Nilai Variabel B4 :'))
    c1 = float(input('Nilai Variabel C1 :'))
    c2 = float(input('Nilai Variabel C2 :'))
    c3 = float(input('Nilai Variabel C3 :'))
    c4 = float(input('Nilai Variabel C4 :'))
    d1 = float(input('Nilai Variabel D1 :'))
    d2 = float(input('Nilai Variabel D2 :'))
    d3 = float(input('Nilai Variabel D3 :'))
    d4 = float(input('Nilai Variabel D4 :'))
    e1 = float(input('Nilai Konstanta E1 :'))
    e2 = float(input('Nilai Konstanta E2 :'))
    e3 = float(input('Nilai Konstanta E3 :'))
    e4 = float(input('Nilai Konstanta E4 :'))
    f1 = float((c1*d1)/(b1-a1)+ e1)
    f2 = float((c2*d2)/(b2-a2)+ e2)
    f3 = float((c3*d1)/(b3-a3)+ e3)
    f4 = float((c4*d1)/(b4-a4)+ e4)
    M = np.array([[a1, b1, c1, d1, f1], 
                  [a2, b2, c2, d2, f2], 
                  [a3, b3, c3, d3, f3], 
                  [a4, b4, c4, d4, f4]])
    n = 4
    print('Matriks Persamaan', M) #Menampilkan matrix awal
        
    J = GaussJordan(M,n) #Menampilkan Hasil
    print(f"""Hasil Pengolahan menggunakan metode Gauss Jordan didapatkan hasil sebagai berikut: {J} \n""")
elif (setting == 8):
    a1 = float(input('Nilai Variabel A1 :'))
    a2 = float(input('Nilai Variabel A2 :'))
    a3 = float(input('Nilai Variabel A3 :'))
    a4 = float(input('Nilai Variabel A4 :'))
    b1 = float(input('Nilai Variabel B1 :'))
    b2 = float(input('Nilai Variabel B2 :'))
    b3 = float(input('Nilai Variabel B3 :'))
    b4 = float(input('Nilai Variabel B4 :'))
    c1 = float(input('Nilai Variabel C1 :'))
    c2 = float(input('Nilai Variabel C2 :'))
    c3 = float(input('Nilai Variabel C3 :'))
    c4 = float(input('Nilai Variabel C4 :'))
    d1 = float(input('Nilai Variabel D1 :'))
    d2 = float(input('Nilai Variabel D2 :'))
    d3 = float(input('Nilai Variabel D3 :'))
    d4 = float(input('Nilai Variabel D4 :'))
    e1 = float(input('Nilai Konstanta E1 :'))
    e2 = float(input('Nilai Konstanta E2 :'))
    e3 = float(input('Nilai Konstanta E3 :'))
    e4 = float(input('Nilai Konstanta E4 :'))
    A = np.array([[a1, b1, c1, d1], 
                  [a2, b2, c2, d2], 
                  [a3, b3, c3, d3], 
                  [a4, b4, c4, d4]])
    f1 = float((c1*d1)/(b1-a1)+ e1)
    f2 = float((c2*d2)/(b2-a2)+ e2)
    f3 = float((c3*d1)/(b3-a3)+ e3)
    f4 = float((c4*d1)/(b4-a4)+ e4)
    B = np.array([f1, f2, f3, f4]) #hasil Input ke fungsi
    x_init = np.zeros(len(B)) 
    J = Jacobi(A,B, x_init)
    print(f"""Hasil perhitungan Iterasi jacobi dari matriks persamaan {A} \n dengan hasil {B} \n didapatkan nilai ke-{len(B)} variabel sebagai berikut, x= {J}"""'\n')
elif (setting == 9):
    a1 = float(input('Nilai Variabel A1 :'))
    a2 = float(input('Nilai Variabel A2 :'))
    a3 = float(input('Nilai Variabel A3 :'))
    a4 = float(input('Nilai Variabel A4 :'))
    b1 = float(input('Nilai Variabel B1 :'))
    b2 = float(input('Nilai Variabel B2 :'))
    b3 = float(input('Nilai Variabel B3 :'))
    b4 = float(input('Nilai Variabel B4 :'))
    c1 = float(input('Nilai Variabel C1 :'))
    c2 = float(input('Nilai Variabel C2 :'))
    c3 = float(input('Nilai Variabel C3 :'))
    c4 = float(input('Nilai Variabel C4 :'))
    d1 = float(input('Nilai Variabel D1 :'))
    d2 = float(input('Nilai Variabel D2 :'))
    d3 = float(input('Nilai Variabel D3 :'))
    d4 = float(input('Nilai Variabel D4 :'))
    e1 = float(input('Nilai Konstanta E1 :'))
    e2 = float(input('Nilai Konstanta E2 :'))
    e3 = float(input('Nilai Konstanta E3 :'))
    e4 = float(input('Nilai Konstanta E4 :'))
    A = np.array([[a1, b1, c1, d1], 
                  [a2, b2, c2, d2], 
                  [a3, b3, c3, d3], 
                  [a4, b4, c4, d4]]) #Untuk masuk ke script input, rumus A cara menulisnya : np.array([[2, -6, -7, 4], [6, -7, -7, -1], [-6, 4, 5, -2], [-3, 2, 4, 6]])
    f1 = float((c1*d1)/(b1-a1)+ e1)
    f2 = float((c2*d2)/(b2-a2)+ e2)
    f3 = float((c3*d1)/(b3-a3)+ e3)
    f4 = float((c4*d1)/(b4-a4)+ e4)
    B = np.array([f1, f2, f3, f4]) #hasil Input ke fungsi
    x = np.zeros_like(B)
    n = int(input('Masukan Batas Iterasi:')) #Disarankan batas iterasi dimulai dari 1 sampai 15 agar tidak terlalu banyak
    S = GaussSeidel(A, B, x, n)
    O = solve(A, B)
    print(f'Hasil Gauss Seidel didapatkan nilai tiap variabel {S}')
    print(f'Variabel matriks dengan SciPy {O}\n')
elif (setting == 10):
    #Masukan Koefisien dan Konstanta dari Persamaan yang diketahui
    A = float(input('Koefisien x Pangkat 3 :'))
    B = float(input('Koefisien x Pangkat 2 :'))
    C = float(input('Koefisien x Pangkat 1 :'))
    D = float(input('Konstanta :'))
    #Mendefinisikan fungsi
    def f(x): 
        return A*(x*3) + B(x*2) + C(x**1) - D
    #Inputan batas persamaan
    batasbawah = float(input('Batas Bawah :')) #Input batas bawah
    batasatas = float(input('Batas Atas :')) #Input Batas atas
    pias = int(input('Banyak Pias :')) #pias
    #Mendefinisikan Fungsi dengan Derajat 3 dengan 1 Variabel x
    #Bentuk Fungsi yang dapat dilakukan adalah : f = AX3 + BX2 + CX - D
    x1 : A*(batasbawah*3) + B(batasbawah*2) + C(batasbawah) - D #Nilai luasan dengan batas bawah
    x2 : A*(batasatas*3) + B(batasatas*2) + C(batasatas) - D #Nilai luasan dengan batas atas
    
    #INTERPRETASI TRAPESIUM DENGAN GRAFIK
    x = np.linspace(batasbawah, batasatas, pias+1)
    y = f(x)
    X = np.linspace(batasbawah, batasatas+1, pias)
    Y = f(X)
    plt.plot(X,Y)
    fig = plt.figure()
    fig.patch.set_facecolor('yellow')
    for i in range(pias):
        xs = [x[i],x[i],x[i+1],x[i+1]]
        ys = [0,f(x[i]),f(x[i+1]),0]
        plt.fill(xs,ys,'b',edgecolor='b',alpha=0.2)
    #Format Keterangan Gambar
    plt.title('Trapesium banyak pias {}'.format(pias))
    #format contoh : D:\Cendra Boskanita\Kuliah\Semester 4\METNUM\PRAKTIKUM\Tugas Akhir Praktikum Metnum 2021\K4_Program METNUM\image\HEUN.png
    #DAPAT BERUPA .png / .jpeg / .pdf
    #Unduh
    plt.savefig(input('Masukan Nama Lokasi :')) 
    plt.show()
    #Menampilkan Hasil
    Lt = Trapesium(f, batasbawah, batasatas, pias)
    print("Luas Trapesium Banyak Pias :", (Lt))
elif (setting == 11):
    #Input persamaan
    a = float(input('Masukkan nilai a pada persamaan: '))
    b = float(input('Masukkan nilai b pada persamaan: '))
    c = float(input('Masukkan nilai c pada persamaan: '))
    d = float(input('Masukkan nilai d pada persamaan: '))
    def f(x):
        return (a*(x*3)) + (b*(x**2)) + (c*x) + d
    #Inputan batas persamaan
    batasbawah = float(input('Batas Bawah :')) #Input batas bawah
    batasatas = float(input('Batas Atas :')) #Input Batas atas
    #Input Nilai Interval
    interval = int(input('Interval :')) #Interval
    # Memanggil Metode Simpson13() untuk mendapatkan hasil
    Hasil13 = Simpson13(batasbawah, batasatas, interval)
    print("Hasil Integrasi Persamaan dengan Metode Simpson13 : %0.6f" % (Hasil13) )
elif (setting == 12):
    #Input persamaan
    a = float(input('Masukkan nilai a pada persamaan: '))
    b = float(input('Masukkan nilai b pada persamaan: '))
    c = float(input('Masukkan nilai c pada persamaan: '))
    d = float(input('Masukkan nilai d pada persamaan: '))
    def f(x):
        return (a*(x*3)) + (b*(x**2)) + (c*x) + d
    #Inputan batas persamaan
    batasbawah = float(input('Batas Bawah :')) #Input batas bawah
    batasatas = float(input('Batas Atas :')) #Input Batas atas
    #Input Nilai Interval
    interval = int(input('Interval :')) #Interval
    Hasil38= simpson38(batasbawah, batasatas, interval)
    print("Hasil Integrasi dari Metode Simpson 3/8 adalah : %0.6f" % (Hasil38))
elif (setting == 13):
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    h = float(input("Masukkan nilai h: ")) 
    x0 = float(input("Masukkan nilai x awal: "))
    xn = float(input("Masukkan nilai x akhir: "))
    y0 = float(input("Masukkan nilai y awal: "))
    x = Euler(a, b, c, d, h, x0, xn, y0)
    print(x)
elif (setting == 14):
    a = float(input("Masukkan nilai a: "))
    b = float(input("Masukkan nilai b: "))
    c = float(input("Masukkan nilai c: "))
    d = float(input("Masukkan nilai d: "))
    h = float(input("Masukkan nilai h: ")) 
    x0 = float(input("Masukkan nilai x awal: "))
    xn = float(input("Masukkan nilai x akhir: "))
    y0 = float(input("Masukkan nilai y awal: "))
    x = Heun(a, b, c, d, h, x0, xn, y0)
    print(x)
else:
    print('Tidak Ditemukan. Silahkan coba lagi dengan input menu yang benar ya :)')


# #### 
