import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#-----------------------------ORIGINAL DATA-----------------------------

#Step 1 - Import antenna pattern information
with open('U7-Outdoor-5GHz.ant', 'r', encoding='utf16') as antenna_file:
    antenna_data = antenna_file.readlines()
antenna_list = [float(line.strip()) for line in antenna_data]

# seperate azimuth and elevation in 1deg increments:
az = antenna_list[0:360]
el = antenna_list[360:720]

theta = np.linspace(0, 2 * np.pi, 360)

#NOTE: This adjustment exists to ensure 0deg is boresight. Remove depending on formatting of file.
temp = az[0:90]
az = az[90:]
az = np.append(az, temp)


#graph azimuth/elevation from original data
az2d = plt.figure();
plt.subplot(projection='polar')
plt.plot(theta, az)
plt.title('Polar Radiation Plot (Azimuth), Normalized Gain [dBi]')
plt.gca().set_yticklabels([])
el2d = plt.figure();
plt.subplot(projection='polar')
plt.plot(theta, el)
plt.title('Polar Radiation Plot (Elevation), Normalized Gain [dBi]')
plt.gca().set_yticklabels([])

#-----------------------------SAMPLED DATA-----------------------------

# Sample every ten degrees to purposely degrade resolution
sampledaz = az[::10]
sampledel = el[::10]
sampledtheta = theta[::10]

#Remove discontinuities from downsampled data 
if sampledaz[0] != sampledaz[-1]:
    sampledaz = np.append(sampledaz, sampledaz[0])
    sampledtheta = np.append(sampledtheta, sampledtheta[0])
    
if sampledel[0] != sampledel[-1]:
    sampledel = np.append(sampledel, sampledel[0])

#graph azimuth/elevation from sampled data
samp_az2d = plt.figure();
plt.subplot(projection='polar')
plt.scatter(sampledtheta, sampledaz)
plt.plot(sampledtheta, sampledaz)
plt.title('Sampled Polar Radiation Plot (Azimuth), Normalized Gain [dBi]')
plt.gca().set_yticklabels([]) 
samp_el2d = plt.figure();
plt.subplot(projection='polar')
plt.scatter(sampledtheta, sampledel)
plt.plot(sampledtheta, sampledel)
plt.title('Sampled Polar Radiation Plot (Elevation), Normalized Gain [dBi]')
plt.gca().set_yticklabels([]) 
#-----------------------------2D INTERPOLATION-----------------------------

# Redefine "sampledtheta" array to be evenly spaced
theta_down = np.linspace(0, 2*np.pi, len(sampledel))


# Post-interpolation resolution is a fraction of 360 (e.g. resolution_factor=0.1 gives 36 point resolution, resolution_factor=1 gives 360 point resolution)
resolution_factor = 1
theta_up = np.linspace(0, 2*np.pi, int(len(el)*resolution_factor))


# Interpolate the sampled data to match the denser range with cubic interpolation
az_denser = interp1d(theta_down, sampledaz, kind='cubic', fill_value="extrapolate")(theta_up)
el_denser = interp1d(theta_down, sampledel, kind='cubic', fill_value="extrapolate")(theta_up)

#graph azimuth/elevation from interpolated data:
intp_az2d = plt.figure();
plt.subplot(projection='polar')
plt.plot(theta_up, az_denser)
plt.title('Interpolated Polar Radiation Plot (Azimuth), Normalized Gain [dBi]')
plt.gca().set_yticklabels([]) 
intp_el2d = plt.figure();
plt.subplot(projection='polar')
plt.plot(theta_up, el_denser)
plt.title('Interpolated Polar Radiation Plot (Elevation), Normalized Gain [dBi]')
plt.gca().set_yticklabels([]) 

#-----------------------------3D INTERPOLATION-----------------------------

def summing_algorithm(g_az,g_el):
    
    azimuth = np.deg2rad(np.arange(0,360,1))
    elevation = np.deg2rad(np.arange(0,180,1))       
    pattern = np.zeros((360,180))
    
    # perform additive scaling to fill 3d gain array:
    for i, az_gain in enumerate(g_az[0:180]):
        for j, el_gain in enumerate(g_el[0:180]):
            pattern[i, j] = az_gain + el_gain

    for i, az_gain in enumerate(g_az[180:360]):
        # flip because we want to go from 0 to 180 but the az-el data comes in 360 deg. format
        for j, el_gain in enumerate(np.flipud(g_el[180:360])):
            pattern[i+180, j] = az_gain + el_gain

    return pattern, azimuth, elevation

def cross_weighting(g_az,g_el):
    # create the angle coordinates necessary for polar plotting:
    azimuth = np.deg2rad(np.arange(0,360,1))
    elevation = np.deg2rad(np.arange(0,180,1))
    
    # create final gain array:
    pattern = np.zeros((360,180))
    
    # perform cross-weighted scaling to make the final gain array:
    #print(10**np.divide(g_el,10))
    
    vert = 10**(np.divide(g_el,10))
    hor = 10**(np.divide(g_az,10))
    k=2
    
    for i, az_gain in enumerate(g_az[0:180]):
        for j, el_gain in enumerate(g_el[0:180]):
            w1 = vert[j] * (1-hor[i])
            w2 = hor[i] * (1-vert[j])
            if (w1 == 0 and w2 == 0):
                pattern[i,j] = 0
            else:
                pattern[i, j] = (az_gain * w1 + el_gain * w2)/(np.sqrt((w1)**k)+(w2)**k)

    for i, az_gain in enumerate(g_az[180:360]):
        # flip because we want to go from 0 to 180 but the az-el data comes in 360 deg. format
        for j, el_gain in enumerate(np.flipud(g_el[180:360])):
            w1 = vert[j] * (1-hor[i+180])
            w2 = hor[i+180] * (1-vert[j])
            if (w1 == 0 and w2 == 0):
                pattern[i+180,j] = 0
            else:
                pattern[i+180, j] = (az_gain * w1 + el_gain * w2)/(np.sqrt((w1)**k)+(w2)**k)

    return pattern, azimuth, elevation

def hybrid(g_az,g_el):
    cross_pat,caz,cel = cross_weighting(g_az, g_el)
    sum_pat,saz,sel = summing_algorithm(g_az,g_el)
    
    sum_pat_dec = 10**(np.divide(sum_pat,10))
    # cross_pat_dec = 10**(np.divide(cross_pat,10))
    
    n = 20
    
    # change  'n' value in **(1/n) for different results
    pattern = sum_pat * (sum_pat_dec**(1/n)) + cross_pat * (1-(sum_pat_dec)**(1/n))
    azimuth = caz
    elevation = cel

    return pattern, azimuth, elevation

    

def spherical_to_cartesian(azimuth,elevation,pattern):
    # account for normalization For plotting, since 0dB and negative scaling doesn't show up well
    if np.any(pattern < 0):
        norm = -np.min(pattern)
    else: 
        norm = 0
    
    x = np.zeros((360,180))
    y = np.zeros((360,180))
    z = np.zeros((360,180))
    gain = np.zeros((360,180))
    
    for i,az in enumerate(azimuth):
        for j,el in enumerate(elevation):
            g = pattern[i,j] + norm
            x[i,j] = (g*np.sin(el)*np.cos(az))
            y[i,j] = (g*np.sin(el)*np.sin(az))
            z[i,j] = g*np.cos(el)
            gain[i,j] = g
            
    return x, y, z, gain
            

def plot_antenna_pattern(az_gain,el_gain,method='cross'):
    if method == 'summing':
        g,az,el = summing_algorithm(az_gain,el_gain)
        name='Antenna Pattern Data, 3D [Summing Algorithim]'
    elif method == 'cross':
         g,az,el = cross_weighting(az_gain,el_gain)
         name='Antenna Pattern Data, 3D [Approximation Algorithim]'
    elif method == 'hybrid':
         g,az,el = hybrid(az_gain,el_gain)
         name='Antenna Pattern Data, 3D [Hybrid Algorithim]'
    else:
         raise Exception("Choose a valid method: 'Summing', 'Cross', 'Hybrid', or leave that field blank")
        
    x,y,z,gain = spherical_to_cartesian(az,el,g)
    
    return g, x, y, z, name


#-----------------------------RUN SCRIPT-----------------------------
g, x, y, z, name = plot_antenna_pattern(az_denser,el_denser,'summing') #NOTE: Pick method here
az_post3d = g[:,0] #Azimuthal slice at elevation = 0deg
el_post3d = np.array(list(g[0,:])+list(np.flipud(g[359,:])))

#Plot 3D Graph after interpolation
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_surface(x,y,z, cmap = 'jet',edgecolor='none')
ax.set_title('Interpolated 3D Antenna pattern')
plt.axis('off')

#Compare original and 2D [Azimuth]
az_2d_cartesian1, ax2 = plt.subplots()
plt.plot(theta, az, label = "Original")
plt.plot(theta_up, az_denser, label = "2D Interpolation")
plt.title('Rectangular Raidiation Plot (Azimuth)')
plt.xlabel('Angle [rad]')
plt.ylabel('Normalized Gain [dBi]')
custom_ticks = [0, np.pi/2, np.pi, (3*np.pi)/2, 2*np.pi]
custom_tick_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
ax2.set_xticks(custom_ticks)
ax2.set_xticklabels(custom_tick_labels)
plt.legend()


#Calculate avg  percent error [2D_AZ]
AvgError2D_Az = (az - az_denser)
AvgError2D_Az = AvgError2D_Az ** 2
AvgError2D_Az = np.mean(AvgError2D_Az)
print(AvgError2D_Az)


#Compare original and 3D [Azimuth]
az_2d_cartesian2, ax3 = plt.subplots()
plt.plot(theta, az, label = "Original")
plt.plot(theta, az_post3d, label = "3D Interpolation")
plt.title('Rectangular Raidiation Plot (Azimuth)')
plt.xlabel('Angle [rad]')
plt.ylabel('Normalized Gain [dBi]')
custom_ticks = [0, np.pi/2, np.pi, (3*np.pi)/2, 2*np.pi]
custom_tick_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
ax3.set_xticks(custom_ticks)
ax3.set_xticklabels(custom_tick_labels)
plt.legend()


#Calculate avg  percent error [3D_AZ]
AvgError3D_Az = (az - az_post3d)
AvgError3D_Az = AvgError3D_Az ** 2
AvgError3D_Az = np.mean(AvgError3D_Az)
print(AvgError3D_Az)


#Compare original and 2D [Elevation]
el_2d_cartesian1, ax4 = plt.subplots()
plt.plot(theta, el, label = "Original Data")
plt.plot(theta_up, el_denser, label = "2D Interpolation")
plt.title('Rectangular Raidiation Plot (Elevation)')
plt.xlabel('Angle [rad]')
plt.ylabel('Normalized Gain [dBi]')
custom_ticks = [0, np.pi/2, np.pi, (3*np.pi)/2, 2*np.pi]
custom_tick_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
ax4.set_xticks(custom_ticks)
ax4.set_xticklabels(custom_tick_labels)
plt.legend()


#Calculate avg  percent error [2D_AZ]
AvgError2D_El = (el - el_denser)
AvgError2D_El = AvgError2D_El ** 2
AvgError2D_El = np.mean(AvgError2D_El)      
print(AvgError2D_El)


#Compare original and 3D [Elevation]
el_2d_cartesian2, ax5 = plt.subplots()
plt.plot(theta, el, label = "Original Data")
plt.plot(theta_up, el_post3d, label = "3D Interpolation")
plt.title('Rectangular Raidiation Plot (Elevation)')
plt.xlabel('Angle [rad]')
plt.ylabel('Normalized Gain [dBi]')
custom_ticks = [0, np.pi/2, np.pi, (3*np.pi)/2, 2*np.pi]
custom_tick_labels = ['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$']
ax5.set_xticks(custom_ticks)
ax5.set_xticklabels(custom_tick_labels)
plt.legend()


#Calculate avg  percent error [3D_AZ]
AvgError3D_El = (el - el_post3d)
AvgError3D_El = AvgError3D_El ** 2
AvgError3D_El = np.mean(AvgError3D_El)
print(AvgError3D_El)


#-----------------------------COALATE GRAPHS-----------------------------
    #Polar Graphs
root = tk.Tk()
root.title('Antenna Pattern Data, Polar')
    
#original resolution, 1 deg increments
canvas_polar1 = FigureCanvasTkAgg(az2d, master=root)
canvas_polar1.get_tk_widget().grid(row=0, column=0)
canvas_polar2 = FigureCanvasTkAgg(el2d, master=root)
canvas_polar2.get_tk_widget().grid(row=1, column=0)

#sampled resolution, 10 deg increments
canvas_polar3 = FigureCanvasTkAgg(samp_az2d, master=root)
canvas_polar3.get_tk_widget().grid(row=0, column=1)
canvas_polar4 = FigureCanvasTkAgg(samp_el2d, master=root)
canvas_polar4.get_tk_widget().grid(row=1, column=1)

#interpolated resolution, 1 deg increments
canvas_polar5 = FigureCanvasTkAgg(intp_az2d, master=root)
canvas_polar5.get_tk_widget().grid(row=0, column=2)
canvas_polar6 = FigureCanvasTkAgg(intp_el2d, master=root)
canvas_polar6.get_tk_widget().grid(row=1, column=2)

    #Cartesian 3D Graph
root = tk.Tk()
root.title(name)
canvas_3d2 = FigureCanvasTkAgg(fig, master=root)
canvas_3d2.get_tk_widget().grid()#row=0, column=0, rowspan=10, columnspan=10)

    #Graph Error
root = tk.Tk()
root.title('Antenna Pattern Data, Error')
canvas_cartesian1 = FigureCanvasTkAgg(az_2d_cartesian1, master=root)
canvas_cartesian1.get_tk_widget().grid(row=0, column=0)
canvas_cartesian2 = FigureCanvasTkAgg(az_2d_cartesian2, master=root)
canvas_cartesian2.get_tk_widget().grid(row=0, column=1)
canvas_cartesian3 = FigureCanvasTkAgg(el_2d_cartesian1, master=root)
canvas_cartesian3.get_tk_widget().grid(row=1, column=0)
canvas_cartesian4 = FigureCanvasTkAgg(el_2d_cartesian2, master=root)
canvas_cartesian4.get_tk_widget().grid(row=1, column=1)

    #display
root.mainloop()
