# METODE-NUMERIK-2021
Script ini ditujukan untuk memenuhi Tugas Akhir Praktikum Metode Numerik, Prodi Oseanografi Universitas Diponegoro. Di dalam script mencakup perhitungan-perhitungan numerik yang menggunakan berbagai library :

#import library yang digunakan //Apabila belum memiliki library
1. import sys
2. import numpy as np			                                                        //os.system('pip install numpy')
3. from numpy import float32, single, double, array, zeros, diag, diagflat, dot
4. import matplotlib.pyplot as plt		                                            //os.system('pip install matplotlib')
5. from math import sin
6. from IPython import get_ipython	                                              //os.system('pip install ipython')
7. from pprint import pprint
8. from scipy.linalg import solve 		                                            //os.system('pip install scipy')
9. plt.style.use('seaborn-poster')
10. ipy = get_ipython()
11. if ipy is not None:
     ipy.run_line_magic('matplotlib','inline')

#Kegunaan Library:
1. sys : menyediakan akses ke beberapa variabel yang digunakan atau dipelihara oleh interpreter dan ke fungsi yang berinteraksi kuat dengan interpreter.  jenis library ini  selalu tersedia.
2. numpy : menambahkan dukungan untuk array dan matriks multi-dimensi yang besar, bersama dengan koleksi besar fungsi matematika tingkat tinggi untuk beroperasi pada array ini.
3. matplotlib : pustaka lengkap untuk membuat visualisasi statis, animasi, dan interaktif dengan Python.
4. math  : menyediakan akses ke fungsi matematika yang ditentukan oleh standar C.
5. IPython : menyediakan toolkit yang kaya untuk membantu Anda memaksimalkan penggunaan Python secara interaktif.  Komponen utamanya adalah bantuan kernel Jupyter untuk bekerja dengan kode Python di notebook Jupyter dan antarmuka interaktif lainnya.
6. pprint : menyediakan kemampuan untuk "mencetak cukup" struktur data Python arbitrer dalam bentuk yang dapat digunakan sebagai input ke penerjemah.
7. scipy : menyediakan lebih banyak fungsi utilitas untuk pengoptimalan, statistik, dan pemrosesan sinyal.

#Metode Yang Terdapat Dalam Script
Kode penggunaan akar-akar persamaan:
1. Metode Setengah Interval 
2. Metode Interpolasi Linier 
3. Metode Secant 
4. Metode Newton-Raphson 
5. Metode Iterasi
Kode penggunaan Sistem Persamaan Linier dan Matriks:
6. Metode Gauss 
7. Metode Gauss Jordan 
8. Metode Gauss Seidel 
9. Metode Jacobi 
Kode penggunaan Integrasi Numerik:
10. Metode Trapesium 
11. Metode Simpson 1/3 
12. Metode Simpson 3/8
Kode penggunaan Persamaan Diferensial Biasa :
13. Metode Euler 
14. Metode Heun

#Kelompok 4 terdiri dari :
1. Petrik Siano OPL		          26050119130125	Oseanografi B
2. Kenichy Priyoasmoro		      26050119130126	Oseanografi B
3. Wahyu Erfando 		            26050119130131	Oseanografi B
4. Syifa Agfanita			          26050119130132	Oseanografi A
5. Kurniawan Sandres		        26050119130133	Oseanografi B
6. Amalia Sekar Ayuningtyas	    26050119130135	Oseanografi B
7. Gisela Malya Asoka Anindita	26050119130137	Oseanografi B
8. Baker Simamora		            26050119130142	Oseanografi A
9. Fawasatya Bagus Kresna	      26050119140038	Oseanografi B
10. Syifa Arrahmah		          26050119140050	Oseanografi B
11. Cendra Boskanita Petrova	  26050119140055	Oseanografi A
12. Prima Riliayunda P 		      26050119140056	Oseanografi B
 
#Saran:
Semoga praktikum selanjutnya metode numerik dapat memberikan materi dengan bahasa yang mudah dimengerti, pemberian tutor yang efektif agar implementasi yang diharapkan dapat berjalan dengan baik, serta kedepannya hasil pekerjaan ini bisa digunakan sebagai bahan pembelajaran berikutnya.

#Terima kasih kepada:
Bapak Dr. Aris Ismanto, S.Si., M.Si selaku dosen koordinator praktikum mata kuliah metode numerik serta Ibu Dr. Ir. Dwi Haryo Ismunarti, M.Si., Bapak Azis Rifai, S.T., M.Si., dan Ibu Rikha Widiaratih, S.Si., M.Si. selaku dosen mata kuliah metode numerik.
Para asisten praktikum metode numerik yang senantiasa membantu dalam pengerjaan tugas praktikum.
Anggota kelompok 4 yang senantiasa berjuang bersama-sama untuk tugas akhir ini.
Teman-teman dan kakak-kakak baik dari Oseanografi maupun luar departemen yang namanya tidak bisa disebutkan satu per satu.
