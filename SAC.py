''' ------------------------------------------------------------------
This module defines the Standard Atmosphere.

The function `get_parameters` takes the input altitude in [km]
and computes temperature, pressure and density at that altitude.

Data was taken from:

https://en.wikipedia.org/wiki/Standard_sea_level
https://en.wikipedia.org/wiki/Standard_gravity
https://en.wikipedia.org/wiki/Gas_constant
https://en.wikipedia.org/wiki/International_Standard_Atmosphere
------------------------------------------------------------------ '''

import numpy as np
from matplotlib import pyplot as plt

# Standard sea level pressure, temperature and air density:
T0 = 288.15     # [K]
p0 = 101325.0   # [Pa]
rho0 = 1.225    # [kg/m3]

# Standard acceleration due to gravity:
g = 9.80665     # [kg*m/s2]

# Specific gas constant for air:
R = 287.058     # [J/(kg*K)]

# Lapse rates and atmospheric zones altitudes:
# TROPOSPHERE .......................................... (0-10.999)km
h_ts = 0        # [m]
a_ts = -0.0065  # [K/m]
# TROPOPAUSE ========================================== (11-19.999)km
h_tp = 11000    # [m]
a_tp = 0        # [K/m] (isothermal)
# STRATOSPHERE ........................................ (20-31.999)km
h_ss1 = 20000   # [m]
a_ss1 = 0.001   # [K/m]
# ..................................................... (32-46.999)km
h_ss2 = 32000   # [m]
a_ss2 = 0.0028  # [K/m]
# STRATOPAUSE ========================================= (47-50.999)km
h_sp = 47000    # [m]
a_sp = 0        # [K/m] (isothermal)
# MESOSPHERE .......................................... (51-70.999)km
h_ms1 = 51000   # [m]
a_ms1 = -0.0028 # [K/m]
# ......................................................... (71-85)km
h_ms2 = 71000   # [m]
a_ms2 = -0.002  # [K/m]
# ===================================================================
h_fin = 85000   # [m]

# Plotting parameters:
# Fonts:
csfont = {'fontname':'Charter', 'fontweight':'regular'}
hfont = {'fontname':'Charter', 'fontweight':'bold'}
ifont = {'fontname':'Charter', 'fontweight':'regular', 'style':'italic'}

# Colours:
temperatureColour = '#0078ff'
pressureColour = '#ff6103'
densityColour = '#18990c'
zonesColor = '#818a8b'
currentAltitudeColour = '#db2727'
font_axes = 10
font_labels = 12
font_title = 18
font_text = 14
step = 1

def get_parameters(altitude):

    # Convert altitude from [km] to [m]:
    altitude = altitude * 1000

    # Temperature, pressure and density at the upper boundaries:
    # Upper boundary of troposphere: ....................................
    T_1 = T0 + a_ts * (h_tp - h_ts)
    p_1 = p0*(T_1/T0)**(-g/(a_ts*R))
    rho_1 = rho0*(T_1/T0)**(-g/(a_ts*R) - 1)
    # Graph plotting data:
    Y1 = np.arange(h_ts, h_tp, step)
    XT1 = T0 + a_ts*(Y1 - h_ts)
    Xp1 = p0*(XT1/T0)**(-g/(a_ts*R))
    Xrho1 = rho0*(XT1/T0)**(-g/(a_ts*R) - 1)

    # Upper boundary of tropopause: .....................................
    T_2 = T_1
    p_2 = p_1 * np.exp(-(g/(R*T_2)) * (h_ss1 - h_tp))
    rho_2 = rho_1 * np.exp(-(g/(R*T_2)) * (h_ss1 - h_tp))
    # Graph plotting data:
    Y2 = np.arange(h_tp, h_ss1, step)
    XT2 = T_1 + a_tp*(Y2 - h_tp)
    Xp2 = p_1 * np.exp(-(g/(R*XT2)) * (Y2 - h_tp))
    Xrho2 = rho_1 * np.exp(-(g/(R*XT2)) * (Y2 - h_tp))

    # Upper boundary of stratosphere (1): ...............................
    T_3 = T_2 + a_ss1 * (h_ss2 - h_ss1)
    p_3 = p_2*(T_3/T_2)**(-g/(a_ss1*R))
    rho_3 = rho_2*(T_3/T_2)**(-g/(a_ss1*R) - 1)
    # Graph plotting data:
    Y3 = np.arange(h_ss1, h_ss2, step)
    XT3 = T_2 + a_ss1*(Y3 - h_ss1)
    Xp3 = p_2*(XT3/T_2)**(-g/(a_ss1*R))
    Xrho3 = rho_2*(XT3/T_2)**(-g/(a_ss1*R) - 1)

    # Upper boundary of stratosphere (2): ...............................
    T_4 = T_3 + a_ss2 * (h_sp - h_ss2)
    p_4 = p_3*(T_4/T_3)**(-g/(a_ss2*R))
    rho_4 = rho_3*(T_4/T_3)**(-g/(a_ss2*R) - 1)
    # Graph plotting data:
    Y4 = np.arange(h_ss2, h_sp, step)
    XT4 = T_3 + a_ss2*(Y4 - h_ss2)
    Xp4 = p_3*(XT4/T_3)**(-g/(a_ss2*R))
    Xrho4 = rho_3*(XT4/T_3)**(-g/(a_ss2*R) - 1)

    # Upper boundary of stratopause: ....................................
    T_5 = T_4
    p_5 = p_4 * np.exp(-(g/(R*T_5)) * (h_ms1 - h_sp))
    rho_5 = rho_4 * np.exp(-(g/(R*T_5)) * (h_ms1 - h_sp))
    # Graph plotting data:
    Y5 = np.arange(h_sp, h_ms1, step)
    XT5 = T_4 + a_sp*(Y5 - h_sp)
    Xp5 = p_4 * np.exp(-(g/(R*XT5)) * (Y5 - h_sp))
    Xrho5 = rho_4 * np.exp(-(g/(R*XT5)) * (Y5 - h_sp))

    # Upper boundary of mezosphere (1): .................................
    T_6 = T_5 + a_ms1 * (h_ms2 - h_ms1)
    p_6 = p_5*(T_6/T_5)**(-g/(a_ms1*R))
    rho_6 = rho_5*(T_6/T_5)**(-g/(a_ms1*R) - 1)
    # Graph plotting data:
    Y6 = np.arange(h_ms1, h_ms2, step)
    XT6 = T_5 + a_ms1*(Y6 - h_ms1)
    Xp6 = p_5*(XT6/T_5)**(-g/(a_ms1*R))
    Xrho6 = rho_5*(XT6/T_5)**(-g/(a_ms1*R) - 1)

    # Upper boundary of mezosphere (2): .................................
    T_7 = T_6 + a_ms2 * (h_fin - h_ms2)
    # Graph plotting data:
    Y7 = np.arange(h_ms2, h_fin, step)
    XT7 = T_6 + a_ms2*(Y7 - h_ms2)
    Xp7 = p_6*(XT7/T_6)**(-g/(a_ms2*R))
    Xrho7 = rho_6*(XT7/T_6)**(-g/(a_ms2*R) - 1)

    # Temperature, pressure and density calculation:
    if altitude >= h_ts and altitude < h_tp:
        print('You are in the troposphere.')
        zone = 'troposphere'
        # In the troposphere:
        T_fin = T0 + a_ts * (altitude - h_ts)
        p_fin = p0*(T_fin/T0)**(-g/(a_ts*R))
        rho_fin = rho0*(T_fin/T0)**(-g/(a_ts*R) - 1)

    elif altitude >= h_tp and altitude < h_ss1:
        print('You are in the tropopause.')
        print('Temperature is constant in this zone.')
        zone = 'tropopause'
        # In the tropopause:
        T_fin = T_1
        p_fin = p_1 * np.exp(-(g/(R*T_fin)) * (altitude - h_tp))
        rho_fin = rho_1 * np.exp(-(g/(R*T_fin)) * (altitude - h_tp))

    elif altitude >= h_ss1 and altitude < h_ss2:
        print('You are in the stratosphere (1).')
        zone = 'stratosphere (1)'
        # In the stratosphere (1):
        T_fin = T_2 + a_ss1 * (altitude - h_ss1)
        p_fin = p_2*(T_fin/T_2)**(-g/(a_ss1*R))
        rho_fin = rho_2*(T_fin/T_2)**(-g/(a_ss1*R) - 1)

    elif altitude >= h_ss2 and altitude < h_sp:
        print('You are in the stratosphere (2).')
        zone = 'stratosphere (2)'
        # In the stratosphere (2):
        T_fin = T_3 + a_ss2 * (altitude - h_ss2)
        p_fin = p_3*(T_fin/T_3)**(-g/(a_ss2*R))
        rho_fin = rho_3*(T_fin/T_3)**(-g/(a_ss2*R) - 1)

    elif altitude >= h_sp and altitude < h_ms1:
        print('You are in the stratopause.')
        print('Temperature is constant in this zone.')
        zone = 'stratopause'
        # In the stratopause:
        T_fin = T_4
        p_fin = p_4 * np.exp(-(g/(R*T_fin)) * (altitude - h_sp))
        rho_fin = rho_4 * np.exp(-(g/(R*T_fin)) * (altitude - h_sp))

    elif altitude >= h_ms1 and altitude < h_ms2:
        print('You are in the mezosphere (1).')
        zone = 'mezosphere (1)'
        # In the mezosphere (1):
        T_fin = T_5 + a_ms1 * (altitude - h_ms1)
        p_fin = p_5*(T_fin/T_5)**(-g/(a_ms1*R))
        rho_fin = rho_5*(T_fin/T_5)**(-g/(a_ms1*R) - 1)

    elif altitude >= h_ms2 and altitude <= h_fin:
        print('You are in the mezosphere (2).')
        zone = 'mezosphere (2)'
        # In the mezosphere (2):
        T_fin = T_6 + a_ms2 * (altitude - h_ms2)
        p_fin = p_6*(T_fin/T_6)**(-g/(a_ms2*R))
        rho_fin = rho_6*(T_fin/T_6)**(-g/(a_ms2*R) - 1)

    print("\nParameters at: " + str(altitude/1000) + " km:\n")
    print("Temperature: " + str(round(T_fin, 2)) + " K")
    print("Pressure: " + str(round(p_fin, 2)) + " Pa")
    print("Density: " + str(round(rho_fin, 5)) + " kg/m^3")

    print("\nPercentage of the sea level values:\n")
    print("Temperature: " + str(round(T_fin/T0*100, 2)) + "%")
    print("Pressure: " + str(round(p_fin/p0*100, 5)) + "%")
    print("Density: " + str(round(rho_fin/rho0*100, 5)) + "%")

    # Plotting:
    figure = plt.figure(figsize=(15, 12))

    figureSubplot = figure.add_subplot(1,3,1)
    plt.plot(XT1, Y1/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT2, Y2/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT3, Y3/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT4, Y4/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT5, Y5/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT6, Y6/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(XT7, Y7/1000, color=temperatureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.scatter(T_fin, altitude/1000, color=currentAltitudeColour, zorder=1)
    plt.hlines([altitude/1000], 150, 300, color=currentAltitudeColour, linestyle='--', linewidth=1.0)
    plt.hlines([h_tp/1000, h_ss1/1000, h_ss2/1000, h_sp/1000, h_ms1/1000, h_ms2/1000], 150, 300, color=zonesColor, linestyle='--', linewidth=1.0)
    plt.text(155, h_ts/1000+1, 'Troposphere')
    plt.text(155, h_tp/1000+1, 'Tropopause')
    plt.text(155, h_ss1/1000+1, 'Stratosphere (1)')
    plt.text(155, h_ss2/1000+1, 'Stratosphere (2)')
    plt.text(155, h_sp/1000+1, 'Stratopause')
    plt.text(155, h_ms1/1000+1, 'Mezosphere (1)')
    plt.text(155, h_ms2/1000+1, 'Mezosphere (2)')
    plt.xlabel('$T(h)$ [$K$]', fontsize=font_labels)
    plt.ylabel('Altitude [km]', fontsize=font_labels)
    plt.xlim([150, 300])
    plt.ylim([0,85])
    plt.grid(True, alpha=0.2)

    for label in (figureSubplot.get_xticklabels()):
        label.set_fontsize(font_axes)

    for label in (figureSubplot.get_yticklabels()):
        label.set_fontsize(font_axes)

    figureSubplot = figure.add_subplot(1,3,2)
    plt.plot(Xp1, Y1/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp2, Y2/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp3, Y3/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp4, Y4/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp5, Y5/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp6, Y6/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xp7, Y7/1000, color=pressureColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.scatter(p_fin, altitude/1000, color=currentAltitudeColour, zorder=1)
    plt.hlines([altitude/1000], 0, p0, color=currentAltitudeColour, linestyle='--', linewidth=1.0)
    plt.hlines([h_tp/1000, h_ss1/1000, h_ss2/1000, h_sp/1000, h_ms1/1000, h_ms2/1000], 0, p0, color=zonesColor, linestyle='--', linewidth=1.0)
    plt.xlabel('$p(h)$ [$Pa$]', fontsize=font_labels)
    plt.xlim([0, p0])
    plt.ylim([0,85])
    plt.yticks([])
    plt.grid(True, alpha=0.2)

    for label in (figureSubplot.get_xticklabels()):
        label.set_fontsize(font_axes)

    for label in (figureSubplot.get_yticklabels()):
        label.set_fontsize(font_axes)

    figureSubplot = figure.add_subplot(1,3,3)
    plt.plot(Xrho1, Y1/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho2, Y2/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho3, Y3/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho4, Y4/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho5, Y5/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho6, Y6/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.plot(Xrho7, Y7/1000, color=densityColour, linestyle='-', linewidth=2.0, zorder=0)
    plt.scatter(rho_fin, altitude/1000, color=currentAltitudeColour, zorder=1)
    plt.hlines([altitude/1000], 0, rho0, color=currentAltitudeColour, linestyle='--', linewidth=1.0)
    plt.hlines([h_tp/1000, h_ss1/1000, h_ss2/1000, h_sp/1000, h_ms1/1000, h_ms2/1000], 0, rho0, color=zonesColor, linestyle='--', linewidth=1.0)
    plt.xlabel(r'$\rho(h)$ [$kg/m^3$]', fontsize=font_labels)
    plt.xlim([0, rho0])
    plt.ylim([0,85])
    plt.yticks([])
    plt.grid(True, alpha=0.2)

    for label in (figureSubplot.get_xticklabels()):
        label.set_fontsize(font_axes)

    for label in (figureSubplot.get_yticklabels()):
        label.set_fontsize(font_axes)

    return(T_fin, p_fin, rho_fin)
