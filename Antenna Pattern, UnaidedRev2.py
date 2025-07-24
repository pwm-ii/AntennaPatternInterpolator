import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#-----------------------------ORIGINAL DATA-----------------------------
# Import Antenna Pattern

try:
    with open('U7-Outdoor-5GHz.ant', 'r', encoding='utf16') as antenna_file:
        antenna_data = antenna_file.readlines()
        antenna_list = [float(line.strip()) for line in antenna_data]
except FileNotFoundError:
    print("Error: File 'U7-Outdoor-5GHz.ant' not found.")
    exit()
except Exception as e:
    print(f"Error: {e}")
    exit()

# separate azimuth and elevation in 1deg increments:
az = antenna_list[0:360]
el = antenna_list[360:720]
theta = np.linspace(0, 2 * np.pi, 360, endpoint=False)

# ensure 0deg is boresight
#NOTE: Remove depending on formatting of file.
temp = az[0:90]
az = az[90:]
az = np.append(az, temp)

#-----------------------------SAMPLED DATA-----------------------------
#NOTE: This step only serves to create "representative" data of poor resolution. It should be omitted for real-world applications.

sampledaz = np.append(az[::10], az[0])
sampledel = np.append(el[::10], el[0])
sampledtheta = np.linspace(0, 2 * np.pi, len(sampledaz), endpoint=True)

#-----------------------------2D INTERPOLATION-----------------------------
resolution_factor = 1
theta_up = np.linspace(0, 2*np.pi, int(len(el)*resolution_factor))
#this defines the resolution after interpolation as a fraction of 360
# e.g. a resolution_factor of 0.1 gives 36 point resolution, resolution_factor=1 gives 360 point resolution)

az_denser = interp1d(sampledtheta, sampledaz, kind='cubic', fill_value='extrapolate')(theta_up)
el_denser = interp1d(sampledtheta, sampledel, kind='cubic', fill_value='extrapolate')(theta_up)

#-----------------------------ALGORITHMS-----------------------------
def summing_algorithm(g_az, g_el):
    azimuth = np.deg2rad(np.arange(0,360,1))
    elevation = np.deg2rad(np.arange(0,180,1))
    pattern = np.zeros((360,180))
    pattern[:180, :] = g_az[:180, np.newaxis] + g_el[:180]
    pattern[180:, :] = g_az[180:, np.newaxis] + g_el[180:][::-1]
    return pattern, azimuth, elevation

def approx_algorithm(g_az, g_el):
    #create final gain array
    azimuth = np.deg2rad(np.arange(0,360,1))
    elevation = np.deg2rad(np.arange(0,180,1))
    pattern = np.zeros((360,180))
    vert = 10**(g_el/10)
    hor = 10**(g_az/10)
    k = 2
    # first half
    w1 = vert[:180] * (1 - hor[:180, np.newaxis])
    w2 = hor[:180, np.newaxis] * (1 - vert[:180])
    num = g_az[:180, np.newaxis] * w1 + g_el[:180] * w2
    den = np.sqrt(w1**k) + w2**k
    mask = (w1 == 0) & (w2 == 0)
    pattern[:180, :] = np.where(mask, 0, num/den)
    # second half
    w1 = vert[180:][::-1] * (1 - hor[180:, np.newaxis])
    w2 = hor[180:, np.newaxis] * (1 - vert[180:][::-1])
    num = g_az[180:, np.newaxis] * w1 + g_el[180:][::-1] * w2
    den = np.sqrt(w1**k) + w2**k
    mask = (w1 == 0) & (w2 == 0)
    pattern[180:, :] = np.where(mask, 0, num/den)
    return pattern, azimuth, elevation

def hybrid_algorithm(g_az, g_el):
    approx_pat, az1, el1 = approx_algorithm(g_az, g_el)
    sum_pat, _, _ = summing_algorithm(g_az, g_el)
    sum_dec = 10**(sum_pat/10)
    n = 20
    pattern = sum_pat * (sum_dec**(1/n)) + approx_pat * (1 - sum_dec**(1/n))
    return pattern, az1, el1

def spherical_to_cartesian(azimuth, elevation, pattern):
    norm = -np.min(pattern) if np.any(pattern < 0) else 0
    gain = pattern + norm
    el_grid, az_grid = np.meshgrid(elevation, azimuth, indexing='xy')
    x = gain * np.sin(el_grid) * np.cos(az_grid)
    y = gain * np.sin(el_grid) * np.sin(az_grid)
    z = gain * np.cos(el_grid)
    return x, y, z, gain

def plot_antenna_pattern(az_gain, el_gain, method):
    if method == 'Summing':
        g, az, el = summing_algorithm(az_gain, el_gain)
        name = 'Antenna Pattern, 3D [Summing Algorithm]'
    elif method == 'Approximation':
        g, az, el = approx_algorithm(az_gain, el_gain)
        name = 'Antenna Pattern, 3D [Approximation Algorithm]'
    elif method == 'Hybrid':
        g, az, el = hybrid_algorithm(az_gain, el_gain)
        name = 'Antenna Pattern, 3D [Hybrid Algorithm]'
    else:
        raise Exception("Choose a valid method: 'Summing', 'Approximation', or 'Hybrid'")
    x, y, z, gain = spherical_to_cartesian(az, el, g)
    return g, x, y, z, name

#-----------------------------GUI POPUP-----------------------------
root = tk.Tk()
root.withdraw()

def ask_method(parent):
    dlg = tk.Toplevel(parent)
    dlg.title("Method Selection")
    dlg.geometry("500x200")
    tk.Label(dlg, text="Please select an interpolation method.:").pack(pady=5)
    choice = tk.StringVar(value='summing')

    for m in ['Summing', 'Approximation', 'Hybrid']:
        ttk.Radiobutton(dlg, text=m, variable=choice, value=m).pack(anchor='w', padx=20)
    ttk.Button(dlg, text='OK', command=dlg.destroy).pack(pady=10)
    dlg.grab_set(); parent.wait_window(dlg)
    return choice.get()

method = ask_method(root)
print(f"Method Selected: {method}")
root.deiconify()

#-----------------------------BUILD FIGURES-----------------------------
az2d = plt.figure(); 
plt.subplot(projection='polar'); 
plt.plot(theta, az); 
plt.title('Radiation Pattern, Azimuth [dBi]'); 
plt.gca().set_yticklabels([])

el2d = plt.figure(); 
plt.subplot(projection='polar'); 
plt.plot(theta, el); 
plt.title('Radiation Pattern, Elevation [dBi]'); 
plt.gca().set_yticklabels([])

samp_az2d = plt.figure(); plt.subplot(projection='polar'); 
plt.scatter(sampledtheta, sampledaz); 
plt.plot(sampledtheta, sampledaz); 
plt.title('Sampled Radiation Pattern, Azimuth [dBi]'); 
plt.gca().set_yticklabels([])

samp_el2d = plt.figure(); 
plt.subplot(projection='polar'); 
plt.scatter(sampledtheta, sampledel); 
plt.plot(sampledtheta, sampledel); 
plt.title('Sampled Radiation Pattern, Elevation [dBi]'); 
plt.gca().set_yticklabels([])

intp_az2d = plt.figure(); 
plt.subplot(projection='polar'); 
plt.plot(theta_up, az_denser); 
plt.title('Interpolated Radiation Pattern, Azimuth [dBi]'); 
plt.gca().set_yticklabels([])

intp_el2d = plt.figure();
plt.subplot(projection='polar');
plt.plot(theta_up, el_denser);
plt.title('Interpolated Radiation Pattern, Elevation [dBi]');
plt.gca().set_yticklabels([])

g, x, y, z, name = plot_antenna_pattern(az_denser, el_denser, method)

fig = plt.figure(); 
ax = fig.add_subplot(projection='3d'); 
ax.plot_surface(x, y, z, cmap='jet', edgecolor='none'); 
ax.set_title('Interpolated 3D Antenna pattern'); 
plt.axis('off')

az_post3d = g[:,0]
el_post3d = np.concatenate((g[0, :], np.flipud(g[359, :])))
az_2d_cartesian1, ax2 = plt.subplots();

plt.plot(theta, az, label='Original');

plt.plot(theta_up, az_denser, label='2D Interpolation'); 
plt.title('Radiation Pattern, Azimuth (Rectangular)'); 
plt.xlabel('Angle [rad]'); 
plt.ylabel('Normalized Gain [dBi]'); 
custom_ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]; 
custom_tick_labels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']; 
ax2.set_xticks(custom_ticks); 
ax2.set_xticklabels(custom_tick_labels); 
plt.legend()

az_2d_cartesian2, ax3 = plt.subplots(); 
plt.plot(theta, az, label='Original'); 
plt.plot(theta, az_post3d, label='3D Interpolation'); 
plt.title('Radiation Pattern, Azimuth (Rectangular)'); 
plt.xlabel('Angle [rad]'); 
plt.ylabel('Normalized Gain [dBi]'); 
ax3.set_xticks(custom_ticks); 
ax3.set_xticklabels(custom_tick_labels); 
plt.legend()

el_2d_cartesian1, ax4 = plt.subplots(); 
plt.plot(theta, el, label='Original Data'); 
plt.plot(theta_up, el_denser, label='2D Interpolation'); 
plt.title('Radiation Pattern, Elevation (Rectangular)'); 
plt.xlabel('Angle [rad]'); 
plt.ylabel('Normalized Gain [dBi]'); 
ax4.set_xticks(custom_ticks); 
ax4.set_xticklabels(custom_tick_labels); 
plt.legend()

el_2d_cartesian2, ax5 = plt.subplots(); 
plt.plot(theta, el, label='Original Data'); 
plt.plot(theta_up, el_post3d, label='3D Interpolation'); 
plt.title('Radiation Pattern, Elevation (Rectangular)'); 
plt.xlabel('Angle [rad]'); 
plt.ylabel('Normalized Gain [dBi]'); 
ax5.set_xticks(custom_ticks); 
ax5.set_xticklabels(custom_tick_labels); 
plt.legend()

#-----------------------------CALCULATE ERROR-----------------------------
AvgError2D_Az = np.mean((az - az_denser)**2)
AvgError2D_El = np.mean((el - el_denser)**2)
AvgError3D_Az = np.mean((az - az_post3d)**2)
AvgError3D_El = np.mean((el - el_post3d)**2)
print(f"Mean Sq. Error from 2D interpolation (Azimuth): {AvgError2D_Az}")
print(f"Mean Sq. Error from 2D interpolation (Elevation): {AvgError2D_El}")
print(f"Mean Sq. Error from 3D interpolation (Azimuth): {AvgError3D_Az}")
print(f"Mean Sq. Error from 3D interpolation (Elevation): {AvgError3D_El}")

#-----------------------------DISPLAY GRAPHS-----------------------------
root.title("Antenna Pattern Analysis")
root.geometry("1280x720")
root.resizable(False, False)

n = ttk.Notebook(root)
f1, f2, f4, f5, f6, f7 = [ttk.Frame(n) for _ in range(6)]
n.add(f1, text='Original')
n.add(f4, text='Sampled')
n.add(f5, text='Interpolated')
n.add(f2, text='3D Graph')
n.add(f6, text='2D Error')
n.add(f7, text='3D Error')
n.pack(fill="both", expand=True)

for frame in [f1, f2, f4, f5, f6, f7]:
    for r in range(2): frame.grid_rowconfigure(r, weight=1)
    for c in range(2): frame.grid_columnconfigure(c, weight=1)

canvas_polar1 = FigureCanvasTkAgg(az2d, master=f1)
canvas_polar1.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
canvas_polar2 = FigureCanvasTkAgg(el2d, master=f1)
canvas_polar2.get_tk_widget().grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

canvas_polar3 = FigureCanvasTkAgg(samp_az2d, master=f4)
canvas_polar3.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
canvas_polar4 = FigureCanvasTkAgg(samp_el2d, master=f4)
canvas_polar4.get_tk_widget().grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

canvas_polar5 = FigureCanvasTkAgg(intp_az2d, master=f5)
canvas_polar5.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
canvas_polar6 = FigureCanvasTkAgg(intp_el2d, master=f5)
canvas_polar6.get_tk_widget().grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

canvas_3d2 = FigureCanvasTkAgg(fig, master=f2)
canvas_3d2.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

canvas_cartesian1 = FigureCanvasTkAgg(az_2d_cartesian1, master=f6)
canvas_cartesian1.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
canvas_cartesian3 = FigureCanvasTkAgg(el_2d_cartesian1, master=f6)
canvas_cartesian3.get_tk_widget().grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

canvas_cartesian2 = FigureCanvasTkAgg(az_2d_cartesian2, master=f7)
canvas_cartesian2.get_tk_widget().grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
canvas_cartesian4 = FigureCanvasTkAgg(el_2d_cartesian2, master=f7)
canvas_cartesian4.get_tk_widget().grid(row=0, column=1, padx=5, pady=5, sticky="nsew")

root.mainloop()