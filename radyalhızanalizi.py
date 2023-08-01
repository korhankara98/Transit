import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

data = #csv

time = data['Time']
rv = data['RV']
rv_err = data['err_RV']

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# "Some cancer sinus function things" 
def sinus_func(x, A, omega, phi, C):                                # x: zaman, A: genlik, omega: açısal frekans, phi: faz açısı, C: sinüs eğrisinin orta noktası
    return A * np.sin(omega * x + phi) + C

A_est = (rv.max() - rv.min()) / 2                                   # Genliğin başlangıç tahmini.
omega_est = 1.5                                                     # Açısal frekans tahmini. Burada pi'li ifade kullanmak daha doğru olacaktır. Kolaya kaçtım.
phi_est = 0                                                         # Faz açısı tahmini. 
C_est = rv.mean()                                                   # Radyal hız verilerinin ortalamasına eşit olan orta nokta tahmini

p0 = [A_est, omega_est, phi_est, C_est]                             # Başlangıç tahmin listesi
popt, pcov = curve_fit(sinus_func, time, rv, p0=p0, sigma=rv_err)   # popt: en iyi sinüs fiti tahmini parametreleri, pcov: fit edilen parametrelerin hataları ile ilgili ama tam anlamadım.
x_fit = np.linspace(time.min(), time.max(), 1000)                   # Fit edilen eğri için bağımsız değişkenin -ki burada zaman oluyor- dizi aralığı
y_fit = sinus_func(x_fit, *popt)                                    # Fit edilen eğri için bağımlı değişkenin -ki burada da radyal hız oluyor- dizi aralığı
A_fit, omega_fit, phi_fit, C_fit = popt

hamp = A_fit*(-1)                                                   # "-" çıkıyordu - ile çarptım :))
period = 2 * np.pi / omega_fit                                      # Periyot

print("Half Amp. (K): +-", hamp, "m/s")
print("Period (p):", period, "day")

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

max1_x = (time.min() + phi_fit / omega_fit)/2*(-1)                  # İlk maksimum x değeri. "-" değerinde ve 2 kat fazla çıkıyordu onu da törpüledim :D
max2_x = max1_x + period                                            # İkinci maksimum x değeri
max_y = C_fit + hamp                                                # Maksimum y değeri

#print(max1_x)
#print(max2_x)
#print(max_y)

plt.errorbar(time, rv, yerr=rv_err, fmt='o', capsize=3, markersize=1, label='Data')

plt.plot(x_fit, y_fit, label='Sinüs Fit')

plt.xlabel("Time (Day)")                                           
plt.ylabel("V$_r$ (m/s)")
plt.title("Radial Velocity with Error Bars")

plt.axhline(y=0, ls="--",c="blue")
plt.plot([max1_x, max2_x], [max_y, max_y], 'r:',  label='Period')
plt.plot([max1_x, max1_x], [C_fit, max_y], 'g:', label='Half Amp.')

plt.legend()
plt.show()

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

P = period * 86400                                                   # saniyeye çevirme
G = 6.67430e-11                                                      # m^3 kg^(-1) s^(-2)
M_sun = 1.989e30                                                     # kg
M_star = M_sun*1.2                                                   # 1.25 M_Güneş olarak verilmiş (kg)

a_3 = G * M_star * (P**2) / (4 * np.pi**2)
a_mtr = a_3**(1/3)                                                   # metre cinsinden

au = 149597870700                                                    # 1 AU, metre cinsinden
a_au = a_mtr / au                                                    # AU cinsinden

K = hamp
jp = 1.898e27                                                        # Jüpiter kütlesi
mc = 3.285e23                                                        # Merkür kütlesi
ea = 5.972e24                                                        # Dünya kütlesi
m2 = K*((P * M_star**2)/(2*np.pi*G))**(1/3)                          # sin(90) olduğunda m2 için alt değer belirlenir. 
m2_jp = m2 / jp
m2_mc = m2 / mc
m2_ea = m2 / ea

if m2_jp > 13.5 and m2_mc < 0.6:
    print("Bulduğunuz gezegen limitlerinde büyük olasılıkla hata mevcut. Parametreleri kontrol ediniz.")
else:
    print("Yarı büyük eksen uzunluğu (a):", a_au, "AU")
    print("Gezegenin Kütlesi:", m2_jp, "M jüp")
    print("Gezegenin Kütlesi:", m2_mc, "M merc")
    print("Gezegenin Kütlesi:", m2_ea, "M earth")
    print("Bu gezegen kütle limitleri dahilinde kabul edilebilir.")

#//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////