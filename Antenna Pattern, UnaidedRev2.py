import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sys
import os

#  =================== INTERPOLATION ALGORITHMS =================== 


def summing_algorithm(g_az, g_el):
    azimuth = np.deg2rad(np.arange(0,360,1))
    elevation = np.deg2rad(np.arange(0,180,1))
    pattern = np.zeros((360,180))
    pattern[:180, :] = g_az[:180, np.newaxis] + g_el[:180]
    pattern[180:, :] = g_az[180:, np.newaxis] + g_el[180:][::-1]
    return pattern, azimuth, elevation

def approx_algorithm(g_az, g_el):
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

#  =================== GUI SETUP =================== 
class SetupDialog:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Antenna Analysis Setup")
        self.root.geometry("500x250")
        
        self.filepath = None
        self.method = tk.StringVar(value="Summing")
        self.confirmed = False

        # File Selection Frame
        frame_file = ttk.LabelFrame(self.root, text="Input Data")
        frame_file.pack(fill="x", padx=10, pady=10)
        
        # Label to show selected path
        self.lbl_file = ttk.Label(frame_file, text="No file selected", foreground="gray")
        self.lbl_file.pack(side="left", padx=5, fill='x', expand=True)
        
        # THE BROWSE BUTTON
        btn_browse = ttk.Button(frame_file, text="Browse...", command=self.browse_file)
        btn_browse.pack(side="right", padx=5)

        # Method Selection Frame
        frame_method = ttk.LabelFrame(self.root, text="Interpolation Method")
        frame_method.pack(fill="x", padx=10, pady=5)
        
        modes = ['Summing', 'Approximation', 'Hybrid']
        for m in modes:
            ttk.Radiobutton(frame_method, text=m, variable=self.method, value=m).pack(anchor='w', padx=10)

        # Action Buttons
        btn_frame = ttk.Frame(self.root)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Run Analysis", command=self.on_run).pack(side="left", padx=10)
        ttk.Button(btn_frame, text="Cancel", command=self.root.destroy).pack(side="left", padx=10)

    def browse_file(self):
        # This function opens the OS file explorer
        ftypes = [("Antenna Files", "*.ant"), ("Text Files", "*.txt"), ("All Files", "*.*")]
        filename = filedialog.askopenfilename(initialdir=os.getcwd(), title="Select Antenna Pattern File", filetypes=ftypes)
        if filename:
            self.filepath = filename
            self.lbl_file.config(text=filename, foreground="black")

    def on_run(self):
        if not self.filepath:
            messagebox.showerror("Error", "Please select an antenna file first.")
            return
        self.confirmed = True
        self.root.destroy() # Close setup window to proceed

    def show(self):
        self.root.mainloop()
        return self.filepath, self.method.get(), self.confirmed

#  =================== RUN CODE =================== 
def main():
    # Setup GUI on run
    setup = SetupDialog()
    filepath, method, confirmed = setup.show()

    if not confirmed or not filepath:
        print("Operation cancelled by user.")
        sys.exit()

    # Process Data
    try:
        with open(filepath, 'r', encoding='utf16') as antenna_file:
            antenna_data = antenna_file.readlines()
            antenna_list = [float(line.strip()) for line in antenna_data]
    except UnicodeError:
        # Fallback for UTF-8 if UTF-16 fails
        try:
            with open(filepath, 'r', encoding='utf-8') as antenna_file:
                antenna_data = antenna_file.readlines()
                antenna_list = [float(line.strip()) for line in antenna_data]
        except Exception as e:
            messagebox.showerror("File Error", f"Could not read file:\n{e}")
            sys.exit()
    except Exception as e:
        messagebox.showerror("File Error", f"Could not read file:\n{e}")
        sys.exit()

    # Data Processing
    az = antenna_list[0:360]
    el = antenna_list[360:720]
    theta = np.linspace(0, 2 * np.pi, 360, endpoint=False)

    temp = az[0:90]
    az = az[90:]
    az = np.append(az, temp)

    # Sampled Data
    sampledaz = np.append(az[::10], az[0])
    sampledel = np.append(el[::10], el[0])
    sampledtheta = np.linspace(0, 2 * np.pi, len(sampledaz), endpoint=True)

    # 2D Interpolation
    theta_up = theta
    az_denser = interp1d(sampledtheta, sampledaz, kind='cubic', fill_value='extrapolate')(theta_up)
    el_denser = interp1d(sampledtheta, sampledel, kind='cubic', fill_value='extrapolate')(theta_up)

    # Calculate 3D Pattern
    g, x, y, z, name = plot_antenna_pattern(az_denser, el_denser, method)

    # 3. Build Figures
    plt.close('all') # Clear any existing plots

    # Helper to create polar plots quickly
    def create_polar_plot(theta_data, r_data, title_text):
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(projection='polar')
        ax.plot(theta_data, r_data)
        ax.set_title(title_text, fontsize=10)
        ax.set_yticklabels([])
        return fig

    az2d = create_polar_plot(theta, az, 'Original Azimuth [dBi]')
    el2d = create_polar_plot(theta, el, 'Original Elevation [dBi]')
    
    samp_az2d = plt.figure(figsize=(4,4))
    plt.subplot(projection='polar')
    plt.scatter(sampledtheta, sampledaz, s=10)
    plt.plot(sampledtheta, sampledaz)
    plt.title('Sampled Azimuth [dBi]', fontsize=10)
    plt.gca().set_yticklabels([])

    samp_el2d = plt.figure(figsize=(4,4))
    plt.subplot(projection='polar')
    plt.scatter(sampledtheta, sampledel, s=10)
    plt.plot(sampledtheta, sampledel)
    plt.title('Sampled Elevation [dBi]', fontsize=10)
    plt.gca().set_yticklabels([])

    intp_az2d = create_polar_plot(theta_up, az_denser, 'Interpolated Azimuth [dBi]')
    intp_el2d = create_polar_plot(theta_up, el_denser, 'Interpolated Elevation [dBi]')

    # 3D Plot
    fig_3d = plt.figure(figsize=(5,4))
    ax = fig_3d.add_subplot(projection='3d')
    ax.plot_surface(x, y, z, cmap='jet', edgecolor='none')
    ax.set_title(f'Interpolated 3D ({method})')
    plt.axis('off')

    # Cartesian Error Plots
    az_post3d = g[:,0]
    el_post3d = np.concatenate((g[0, :], np.flipud(g[359, :])))
    
    custom_ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
    custom_tick_labels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

    def create_cart_plot(orig_data, comp_data, comp_label, title_text):
        fig, ax = plt.subplots(figsize=(5,4))
        ax.plot(theta, orig_data, label='Original')
        ax.plot(theta_up, comp_data, label=comp_label)
        ax.set_title(title_text)
        ax.set_xlabel('Angle [rad]')
        ax.set_ylabel('Gain [dBi]')
        ax.set_xticks(custom_ticks)
        ax.set_xticklabels(custom_tick_labels)
        ax.legend()
        return fig

    az_2d_cartesian1 = create_cart_plot(az, az_denser, '2D Interp', 'Azimuth (Rectangular)')
    az_2d_cartesian2 = create_cart_plot(az, az_post3d, '3D Interp', 'Azimuth (Rectangular)')
    el_2d_cartesian1 = create_cart_plot(el, el_denser, '2D Interp', 'Elevation (Rectangular)')
    el_2d_cartesian2 = create_cart_plot(el, el_post3d, '3D Interp', 'Elevation (Rectangular)')

    # Results Displayed
    results_root = tk.Tk()
    results_root.title(f"Results - {method} Method")
    results_root.geometry("1280x720")

    n = ttk.Notebook(results_root)
    f1 = ttk.Frame(n); n.add(f1, text='Original')
    f4 = ttk.Frame(n); n.add(f4, text='Sampled')
    f5 = ttk.Frame(n); n.add(f5, text='Interpolated')
    f2 = ttk.Frame(n); n.add(f2, text='3D Graph')
    f6 = ttk.Frame(n); n.add(f6, text='2D Error')
    f7 = ttk.Frame(n); n.add(f7, text='3D Error')
    n.pack(fill="both", expand=True)

    for frame in [f1, f2, f4, f5, f6, f7]:
        for r in range(2): frame.grid_rowconfigure(r, weight=1)
        for c in range(2): frame.grid_columnconfigure(c, weight=1)

    # Embed Canvases
    FigureCanvasTkAgg(az2d, master=f1).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(el2d, master=f1).get_tk_widget().grid(row=0, column=1, sticky="nsew")
    FigureCanvasTkAgg(samp_az2d, master=f4).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(samp_el2d, master=f4).get_tk_widget().grid(row=0, column=1, sticky="nsew")
    FigureCanvasTkAgg(intp_az2d, master=f5).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(intp_el2d, master=f5).get_tk_widget().grid(row=0, column=1, sticky="nsew")
    FigureCanvasTkAgg(fig_3d, master=f2).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(az_2d_cartesian1, master=f6).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(el_2d_cartesian1, master=f6).get_tk_widget().grid(row=0, column=1, sticky="nsew")
    FigureCanvasTkAgg(az_2d_cartesian2, master=f7).get_tk_widget().grid(row=0, column=0, sticky="nsew")
    FigureCanvasTkAgg(el_2d_cartesian2, master=f7).get_tk_widget().grid(row=0, column=1, sticky="nsew")

    # Console Output for Accuracy
    print(f"Error Analysis for {method}:")
    print(f"2D Az Error: {np.mean((az - az_denser)**2):.4f}")
    print(f"2D El Error: {np.mean((el - el_denser)**2):.4f}")
    print(f"3D Az Error: {np.mean((az - az_post3d)**2):.4f}")
    print(f"3D El Error: {np.mean((el - el_post3d)**2):.4f}")

    results_root.mainloop()

if __name__ == "__main__":
    main()