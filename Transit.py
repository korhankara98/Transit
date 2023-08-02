import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import batman
sqrt = np.sqrt
sin = np.sin

data = #csv

t = np.array(data['BJDTDB'])                                     # BJD-TDB
flux = np.array(data['NF'])                                      # Normalize Akı
flux_err = np.array(data['NF_err'])                              # Flux Error

#Parametreler
A = 0.1
E = 0.8
Teff = 5302                                               #Kelvin
Rjup = 71492                                              #km
Rjupm = 71492*1000                                        #m
Rsun = 695700                                             #km
Mjup = 1.898e27                                           #kg
Msun = 1.989e30                                           #kg
G = 6.67*10**-11                                          #m³/kg*s²
AB = 1.496*10**11                                         #m
ABkm = 1.496*10**8                                        #km
Rstar = 1.04                                              #Rsun
Teff = 5302                                               #K
P = 3.21222*86400                                         #saniye
Py = 3.21222/365.25                                       #yıl
Pg = 3.21222                                              #gün
m2sini = 2.0596                                           #Mjup
A = 0.10
E = 0.8
DF2 = 1-0.987
DF = 1-0.9853
#Soru2
Rplanet = sqrt(DF2*Rstar**2)                            #Rsun cinsinden
Rp_jp = (Rplanet*Rsun)/Rjup                             #Rjüp 

#mainsequence R = M 0.8

M_star = Rstar**(1/0.8)

a_3 = G * (M_star*Msun) * (P**2) / (4 * np.pi**2)
a_mtr = (a_3**(1/3))/AB                                   # AB cinsinden
a=(Py**2)**(1/3)


prmsa = a/((Rstar*Rsun)/ABkm)
print("Params.a değeri:", prmsa)
print("Yarı büyük eksen uzunluğu:", a)
print("Yarı büyük eksen uzunluğu anakol yaklaşımı:", a_mtr)
print("Rplanet:", Rplanet)
print("Mstar:", M_star)
print("Gezegenin Jüpiter cinsinden yarı çapı:", Rp_jp)


params = batman.TransitParams()
params.t0 = 2456706.585              # transitin merkezini belirler
params.per = 3.21222                 # periyot
params.rp = Rplanet/Rstar            # gezegenin yıldıza oranı
params.a = 8.809534397907948/Rstar   # yarı büyük eksen / yıldızın yarıçapı !!!!!Aynı zamanda Soru5'in cevabı!!!!!!
params.inc = 87                     # eğim
params.ecc = 0.                      # eksantriklik
params.w = 90.                       # uzunluk boylamı
params.u = [0.86, -0.6]              # limb darkening coefficients
params.limb_dark = "quadratic"       # limb darkening model

m = batman.TransitModel(params, t)  
model_flux = m.light_curve(params)    


from scipy.optimize import curve_fit
def fit_func(t, t0, rp, a, inc, u1, u2):
    params.t0 = t0
    params.rp = rp
    params.a = a
    params.inc = inc
    params.u = [u1, u2]
    
    m = batman.TransitModel(params, t)
    return m.light_curve(params)

p0 = [params.t0, params.rp, params.a, params.inc, params.u[0], params.u[1]]
popt, pcov = curve_fit(fit_func, t, flux, sigma=flux_err, p0=p0)
final_model_flux = fit_func(t, *popt)

x_coords = []
y_coords = []

#t1, t2, t3, t4 ve transit orta noktasını grafik üzerinde tıklayarak belirleyebilmek için
# def onclick(event):
#     global x_coords, y_coords
#     x = event.xdata
#     y = event.ydata
#     if x is not None and y is not None:
#         x_coords.append(x)
#         y_coords.append(y)
#         if len(x_coords) == 5:  # 5 farklı tıklama için durdurma koşulu, t1,t2,t3,t4 ve transit orta noktası
#             fig.canvas.mpl_disconnect(cid) 
#             plt.close()


fig, ax = plt.subplots()
ax.errorbar(t, flux, yerr=flux_err, fmt='.', label='Observed')
ax.plot(t, final_model_flux, 'r', label='Fit model')
ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Flux')
#cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()

#Soru3
t1=2456706.523
t2=2456706.538
t3=2456706.632
t4=2456706.6478

print("Yaklaşık t1:", t1)
print("Yaklaşık t2:", t2)
print("Yaklaşık t3:", t3)
print("Yaklaşık t4:", t4)

t14 = (t4-t1)
t23 = (t3-t2)
t14dk = t14*24*60
t23dk = t23*24*60
tftt = t23/t14


print("Yaklaşık toplam geçiş süresi(gün kesri):", t14)
print("Yaklaşık tam geçiş süresi(gün kesri):", t23)
print("Yaklaşık toplam geçiş süresi(dakika):", t14dk)
print("Yaklaşık tam geçiş süresi(dakika):", t23dk)



b = sqrt((((1-sqrt(DF))**2)-((tftt)**2)*((1+sqrt(DF))**2))/(1-((tftt)**2)))

print("Etki parametresi (b):", b)



aRy = (sqrt(((1-sqrt(DF))**2)-(b**2)))/(sin((t23*np.pi)/Pg))
aaRy = (aRy * Rstar * Rsun)    #km

print("aRy:", aRy)
print("aRy yarı büyük eksen sağlaması:", aaRy)


irad = np.arccos(b*Rstar*Rsun/aaRy)
i = np.degrees(irad)

print("irad:", irad)
print("i:", i)



m2 = m2sini/sin(irad)
print("m2:", m2)                #Mjüp


Ortgy = (m2*Mjup)/((4/3)*np.pi*(Rp_jp*Rjupm)**3)*0.001      #g/cm³

print("Ortgy:", Ortgy)


Tdng = (((1-A)/E)**(1/4))*Teff*sqrt((Rstar*Rsun)/(2*aaRy))

print("Tdng:", Tdng)


Ortyy = ((3*np.pi*aRy**3)/(G*P**2))*0.001           #g/cm³

print("Ortyy:", Ortyy)
