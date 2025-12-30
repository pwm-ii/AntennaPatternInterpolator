''' 
==============================================================
                                                             
  It's as easy as...                                         
                                                             
   ██████╗ ██╗███████╗                                       
   ██╔══██╗██║██╔════╝                                       
   ██████╔╝██║█████╗                                         
   ██╔═══╝ ██║██╔══╝                                         
   ██║     ██║███████╗                                       
   ╚═╝     ╚═╝╚══════╝                                       
                                                             
   PAUL'S INTERPOLATION ENGINE                               
   (REV 3.5)

 ============================================================== 
 '''

import os
import sys
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm


# ==========================================
# Data Handling / Interpolation
# ==========================================

class AntennaModel:
    def __init__(self):
        self.filepath = None
        
        # Raw Data Containers
        self.raw_az = None
        self.raw_el = None
        self.raw_theta_az = None
        self.raw_theta_el = None
        
        # Normalized Data Containers (normalzied to 1deg steps)
        self.az_norm = None
        self.el_norm = None
        self.theta_norm = np.deg2rad(np.arange(0, 360, 1))
        
        # State
        self.is_loop_closed = False
        
        # Results
        self.pattern_3d = None
        self.x = None
        self.y = None
        self.z = None
        self.method_name = ""

    def load_data(self, filepath, do_autocenter, do_loop_closure):
        self.filepath = filepath
        self.is_loop_closed = do_loop_closure
        try:
            # check encoding
            raw_lines = []
            for encoding in ['utf-16', 'utf-8', 'latin-1']:
                try:
                    with open(filepath, 'r', encoding=encoding) as f:
                        raw_lines = f.readlines()
                    break
                except UnicodeError:
                    continue
            
            if not raw_lines:
                raise ValueError("Could not decode file with standard encodings.")

            data = [float(line.strip()) for line in raw_lines if line.strip()]
            
            if len(data) < 10: 
                raise ValueError("File contains insufficient data.")

            # Assume file is 50% Azimuth, 50% Elevation
            midpoint = len(data) // 2
            self.raw_az = np.array(data[:midpoint])
            self.raw_el = np.array(data[midpoint:])
            
            # --- OPTIONAL: AUTO-CENTERING ---
            if do_autocenter:
                # Shift raw data so max peak (which I am interpreting as main lobe) is at index 0
                # WARNING: !! DO NOT USE IF PEAK OF AZIMUTH AND ELEVATION ARE NOT MEANT TO BE ALLIGNED!!
                if len(self.raw_az) > 0:
                    shift_az = -np.argmax(self.raw_az)
                    self.raw_az = np.roll(self.raw_az, shift_az)
                
                if len(self.raw_el) > 0:
                    shift_el = -np.argmax(self.raw_el)
                    self.raw_el = np.roll(self.raw_el, shift_el)
            
            # --- OPTIONAL: LOOP CLOSURE ---
            # This is to help with exceptionally poor quality antenna resolution
            if do_loop_closure:
                # If the first and last values do not match, append the first to the end.
                if len(self.raw_az) > 0 and self.raw_az[0] != self.raw_az[-1]:
                    self.raw_az = np.append(self.raw_az, self.raw_az[0])

                if len(self.raw_el) > 0 and self.raw_el[0] != self.raw_el[-1]:
                    self.raw_el = np.append(self.raw_el, self.raw_el[0])

            # Generate the angle arrays for the raw data (0 to 2pi)
            # If loop is closed, we use endpoint=True (0 to 360 inclusive match)
            # If loop is open, we use endpoint=False (0 to <360)
            use_endpoint = do_loop_closure
            
            self.raw_theta_az = np.linspace(0, 2*np.pi, len(self.raw_az), endpoint=use_endpoint)
            self.raw_theta_el = np.linspace(0, 2*np.pi, len(self.raw_el), endpoint=use_endpoint)
            
            # Normalize data to 360 points
            self._normalize_inputs()

        except Exception as e:
            raise RuntimeError(f"Failed to load file: {e}")

    def _normalize_inputs(self):
        # Original indices
        # Use the stored state to determine if we are interpolating 
        # from a closed loop (endpoint=True) or open (endpoint=False)
        use_endpoint = self.is_loop_closed
        
        idx_az = np.linspace(0, 360, len(self.raw_az), endpoint=use_endpoint)
        idx_el = np.linspace(0, 360, len(self.raw_el), endpoint=use_endpoint)
        
        # Target grid (0..359)
        target_deg = np.arange(0, 360, 1)

        # Cubic interpolation
        self.az_norm = interp1d(idx_az, self.raw_az, kind='cubic', fill_value="extrapolate")(target_deg)
        self.el_norm = interp1d(idx_el, self.raw_el, kind='cubic', fill_value="extrapolate")(target_deg)

    def run_interpolation(self, method):
        self.method_name = method
        
        g_az = self.az_norm
        
        g_el = self.el_norm 

        if method == 'Summing':
            self.pattern_3d = self._algo_summing(g_az, g_el)
        elif method == 'Approximation':
            self.pattern_3d = self._algo_approx(g_az, g_el)
        elif method == 'Hybrid':
            self.pattern_3d = self._algo_hybrid(g_az, g_el)
        else:
            raise ValueError("Unknown method")

        self._spherical_to_cartesian()

    def _spherical_to_cartesian(self):
        azimuth = np.deg2rad(np.arange(0, 360, 1))
        elevation = np.deg2rad(np.arange(0, 180, 1))
        
        # Normalize gain
        norm_val = -np.min(self.pattern_3d) if np.any(self.pattern_3d < 0) else 0
        r = self.pattern_3d + norm_val
        
        el_grid, az_grid = np.meshgrid(elevation, azimuth, indexing='xy')
        
        self.x = r * np.sin(el_grid) * np.cos(az_grid)
        self.y = r * np.sin(el_grid) * np.sin(az_grid)
        self.z = r * np.cos(el_grid)

    # --- EXPORT 3D PATTERN AS CSV ---
    def export_csv(self, filename):
        # Phi[deg],Theta[deg],dB(DirTotal)
            # Phi: Inner loop (-180 to 180)
            # Theta: Outer loop (Starts at -180)
            # Data is shifted so Peak is at Phi=0

        if self.pattern_3d is None:
            return

        try:
            # Find peak to align Phi=0
            peak_flat_idx = np.argmax(self.pattern_3d)
            peak_az_idx, _ = np.unravel_index(peak_flat_idx, self.pattern_3d.shape)
            
            with open(filename, 'w') as f:
                f.write("Phi[deg],Theta[deg],dB(DirTotal)\n")
                
                rows, cols = self.pattern_3d.shape # rows=360 (Az), cols=180 (El)
                
                # Outer Loop: Theta (Elevation)
                # Maps index 0..179 to -180..-1
                for el_idx in range(cols):
                    theta_deg = el_idx - 180
                    
                    # Inner Loop: Phi (Azimuth)
                    # -180 to 180
                    for phi_deg in range(-180, 181):
                        
                        # Calculate read index (shifting peak to phi=0)
                        actual_idx = (peak_az_idx + phi_deg) % 360
                        
                        val = self.pattern_3d[actual_idx, el_idx]
                        
                        f.write(f"{phi_deg},{theta_deg},{val:.4f}\n")
                        
        except Exception as e:
            raise RuntimeError(f"Export failed: {e}")

    # --- INTERPOLATION ALGORITHMS ---
    def _algo_summing(self, g_az, g_el):
        pattern = np.zeros((360, 180))
    
        # First half (0-180)
        pattern[:180, :] = g_az[:180, np.newaxis] + g_el[:180]
        
        # Second half (180-360) -> Reversed
        pattern[180:, :] = g_az[180:, np.newaxis] + g_el[180:][::-1]
        return pattern

    def _algo_approx(self, g_az, g_el):
        pattern = np.zeros((360, 180))
        
        # Calculate linear scale on FULL arrays first
        vert = 10**(g_el/10)
        hor = 10**(g_az/10)
        k = 2
        
        def calc_segment(db_h, lin_h, db_v, lin_v):
            # w1 = vertical_linear * (1 - horizontal_linear)
            w1 = lin_v * (1 - lin_h[:, np.newaxis])
            # w2 = horizontal_linear * (1 - vertical_linear)
            w2 = lin_h[:, np.newaxis] * (1 - lin_v)
            
            # num = horizontal_dB * w1 + vertical_dB * w2
            num = db_h[:, np.newaxis] * w1 + db_v * w2
            den = np.sqrt(w1**k) + w2**k
            mask = (w1 == 0) & (w2 == 0)
            return np.where(mask, 0, num/den)

        # First half: 0 to 180
        pattern[:180, :] = calc_segment(g_az[:180], hor[:180], g_el[:180], vert[:180])
        
        # Second half: 180 to 360 (Elevation reversed)
        pattern[180:, :] = calc_segment(g_az[180:], hor[180:], g_el[180:][::-1], vert[180:][::-1])
        return pattern

    def _algo_hybrid(self, g_az, g_el):
        # Pass full arrays through to sub-functions
        approx = self._algo_approx(g_az, g_el)
        summing = self._algo_summing(g_az, g_el)
        sum_dec = 10**(summing/10)
        n = 20
        return summing * (sum_dec**(1/n)) + approx * (1 - sum_dec**(1/n))

    def get_reconstruction_error(self):
        # Slice the 3D pattern to get back 2D cuts
        
        # Azimuth: Uses index 0 (was 90)
        az_reconstructed = self.pattern_3d[:, 0] 
        
        # Elevation: Uses index 0 and index 359 reversed (was 0 and 180 reversed)
        el_reconstructed = np.concatenate((self.pattern_3d[0, :], self.pattern_3d[359, :][::-1]))

        n_az = min(len(self.az_norm), len(az_reconstructed))
        n_el = min(len(self.el_norm), len(el_reconstructed))

        mse_az = np.mean((self.az_norm[:n_az] - az_reconstructed[:n_az])**2)
        mse_el = np.mean((self.el_norm[:n_el] - el_reconstructed[:n_el])**2)
        
        return mse_az, mse_el, az_reconstructed, el_reconstructed


# ==========================================
# GRAPHIC USER INTERFACE
# ==========================================

class PlotPanel(tk.Frame):
    def __init__(self, parent, title, projection=None):
        super().__init__(parent)

        self.figure = Figure(figsize=(3, 2), dpi=100)
        self.figure.set_tight_layout(True)

        if projection == '3d':
            self.ax = self.figure.add_subplot(111, projection='3d')
        elif projection == 'polar':
            self.ax = self.figure.add_subplot(111, projection='polar')
        else:
            self.ax = self.figure.add_subplot(111)
            
        self.ax.set_title(title, fontsize=10)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
        # Add right-click context menu for saving plots
        self.canvas.get_tk_widget().bind("<Button-3>", self._show_context_menu)
        
    def _show_context_menu(self, event):
        """Show context menu on right-click"""
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(label="Save Image As...", command=self._save_figure)
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()
    
    def _save_figure(self):
        """Save the figure using file dialog"""
        filetypes = [
            ("PNG Image", "*.png"),
            ("PDF Document", "*.pdf"),
            ("SVG Vector", "*.svg"),
            ("JPEG Image", "*.jpg"),
            ("All Files", "*.*")
        ]
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=filetypes,
            title="Save Plot Image"
        )
        
        if filename:
            try:
                self.figure.savefig(filename, dpi=300, bbox_inches='tight')
                messagebox.showinfo("Success", f"Image saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Save Failed", str(e))

class SetupDialog:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("Antenna Interpolation Setup")
        self.root.geometry("500x350")
        
        self.filepath = None
        self.method = tk.StringVar(value="Summing")
        
        self.var_autocenter = tk.BooleanVar(value=True)
        self.var_loop_closure = tk.BooleanVar(value=True)
        
        self.confirmed = False
        self._build_ui()

    def _build_ui(self):
        frame_file = ttk.LabelFrame(self.root, text="Input Data")
        frame_file.pack(fill="x", padx=10, pady=10)
        self.lbl_file = ttk.Label(frame_file, text="No file selected", foreground="gray")
        self.lbl_file.pack(side="left", padx=5, fill='x', expand=True)
        ttk.Button(frame_file, text="Browse...", command=self._browse).pack(side="right", padx=5)

        frame_method = ttk.LabelFrame(self.root, text="Interpolation Method")
        frame_method.pack(fill="x", padx=10, pady=5)
        for m in ['Summing', 'Approximation', 'Hybrid']:
            ttk.Radiobutton(frame_method, text=m, variable=self.method, value=m).pack(anchor='w', padx=10)

        frame_settings = ttk.LabelFrame(self.root, text="Settings")
        frame_settings.pack(fill="x", padx=10, pady=5)
        
        ttk.Checkbutton(frame_settings, text="Auto-Center Peaks (to 0°)", variable=self.var_autocenter).pack(anchor='w', padx=10)
        ttk.Checkbutton(frame_settings, text="Enforce Loop Closure", variable=self.var_loop_closure).pack(anchor='w', padx=10)

        btn_frame = ttk.Frame(self.root)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Run Interpolation", command=self._on_run).pack(side="left", padx=10)
        ttk.Button(btn_frame, text="Cancel", command=self.root.destroy).pack(side="left", padx=10)

    def _browse(self):
        ftypes = [("Antenna Files", "*.ant"), ("Text Files", "*.txt"), ("All Files", "*.*")]
        fname = filedialog.askopenfilename(initialdir=os.getcwd(), filetypes=ftypes)
        if fname:
            self.filepath = fname
            self.lbl_file.config(text=os.path.basename(fname), foreground="black")

    def _on_run(self):
        if not self.filepath:
            messagebox.showerror("Error", "Please select an antenna file first.")
            return
        self.confirmed = True
        self.root.destroy()

    def show(self):
        self.root.mainloop()
        return (self.filepath, self.method.get(), 
                self.var_autocenter.get(), self.var_loop_closure.get(), 
                self.confirmed)

class ResultsWindow:
    def __init__(self, model):
        self.model = model
        self.root = tk.Tk()
        self.root.title(f"Results - {model.method_name} Method")
        self.root.geometry("800x550")
        
        # Tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True)
        
        self._init_tab_raw()
        self._init_tab_input()
        self._init_tab_3d()
        self._init_tab_error()

    def _init_tab_raw(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Original (Raw)')
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=1)

        p1 = PlotPanel(frame, "Raw Azimuth Samples [dBi]", projection='polar')
        p1.grid(row=0, column=0, sticky="nsew")
        p1.ax.plot(self.model.raw_theta_az, self.model.raw_az, color='g', marker='.', linestyle='-', linewidth=0.5, label='Raw Az')
        
        p2 = PlotPanel(frame, "Raw Elevation Samples [dBi]", projection='polar')
        p2.grid(row=0, column=1, sticky="nsew")
        p2.ax.plot(self.model.raw_theta_el, self.model.raw_el, color='g', marker='.', linestyle='-', linewidth=0.5, label='Raw El')

    def _init_tab_input(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Original (Smoothed)')
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=1)

        p1 = PlotPanel(frame, "Azimuth Cut [dBi]", projection='polar')
        p1.grid(row=0, column=0, sticky="nsew")
        p1.ax.plot(self.model.theta_norm, self.model.az_norm, color='g', label='Input')
        
        p2 = PlotPanel(frame, "Elevation Cut [dBi]", projection='polar')
        p2.grid(row=0, column=1, sticky="nsew")
        p2.ax.plot(self.model.theta_norm, self.model.el_norm, color='g', label='Input')

    def _init_tab_3d(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='3D Pattern')
        frame.columnconfigure(0, weight=1)
        
        frame.rowconfigure(0, weight=0)
        frame.rowconfigure(1, weight=1)

        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=5)
        
        btn_export = ttk.Button(btn_frame, text="Export 3D Pattern as CSV", command=self._export_csv)
        btn_export.pack(anchor="center")

        p3d = PlotPanel(frame, f"Reconstructed Pattern ({self.model.method_name})", projection='3d')
        p3d.grid(row=1, column=0, sticky="nsew")
        
        surf = p3d.ax.plot_surface(
            self.model.x, self.model.y, self.model.z, 
            cmap=cm.jet, edgecolor='none', alpha=0.9
        )
        p3d.figure.colorbar(surf, ax=p3d.ax, shrink=0.5, aspect=5, label='Gain [dBi]')
        p3d.ax.axis('off')

    def _export_csv(self):
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")],
            title="Export 3D Pattern Data"
        )
        if filename:
            try:
                self.model.export_csv(filename)
                messagebox.showinfo("Export Successful", f"Data saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Export Failed", str(e))

    def _init_tab_error(self):
        """Tab 3: Cartesian & Polar Error Cuts (2x2 Grid)"""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Reconstruction Error')
        
        # Layout: 2 Rows, 2 Columns + 1 Row for Stats
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=3) # Rect Plots
        frame.rowconfigure(1, weight=3) # Polar Plots
        frame.rowconfigure(2, weight=1) # Text stats

        # Calculate errors
        mse_az, mse_el, az_rec, el_rec = self.model.get_reconstruction_error()

        ticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]
        labels = [r'$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$']

        # --- ROW 0: RECTANGULAR PLOTS ---
        
        # Rect Azimuth
        p_az_rect = PlotPanel(frame, "Azimuth Error (Rectangular)")
        p_az_rect.grid(row=0, column=0, sticky="nsew")
        p_az_rect.ax.plot(self.model.theta_norm, self.model.az_norm, 'g:', linewidth=2, label='Original')
        p_az_rect.ax.plot(self.model.theta_norm, az_rec, 'k', label='Reconstructed')
        p_az_rect.ax.set_xticks(ticks)
        p_az_rect.ax.set_xticklabels(labels)
        p_az_rect.ax.set_ylabel("Gain [dBi]")
        p_az_rect.ax.legend()

        # Rect Elevation
        p_el_rect = PlotPanel(frame, "Elevation Error (Rectangular)")
        p_el_rect.grid(row=0, column=1, sticky="nsew")
        p_el_rect.ax.plot(self.model.theta_norm, self.model.el_norm, 'g:', linewidth=2, label='Original')
        p_el_rect.ax.plot(self.model.theta_norm, el_rec, 'k', label='Reconstructed')
        p_el_rect.ax.set_xticks(ticks)
        p_el_rect.ax.set_xticklabels(labels)
        p_el_rect.ax.legend()

        # --- ROW 1: POLAR PLOTS ---

        # Polar Azimuth
        p_az_pol = PlotPanel(frame, "Azimuth Error (Polar)", projection='polar')
        p_az_pol.grid(row=1, column=0, sticky="nsew")
        p_az_pol.ax.plot(self.model.theta_norm, self.model.az_norm, 'g:', linewidth=2, label='Original')
        p_az_pol.ax.plot(self.model.theta_norm, az_rec, 'k', label='Reconstructed')

        # Polar Elevation
        p_el_pol = PlotPanel(frame, "Elevation Error (Polar)", projection='polar')
        p_el_pol.grid(row=1, column=1, sticky="nsew")
        p_el_pol.ax.plot(self.model.theta_norm, self.model.el_norm, 'g:', linewidth=2, label='Original')
        p_el_pol.ax.plot(self.model.theta_norm, el_rec, 'k', label='Reconstructed')

        # --- ROW 2: STATS ---
        stats_frame = ttk.Frame(frame)
        stats_frame.grid(row=2, column=0, columnspan=2, sticky="nsew", padx=20)
        
        txt = (f"Method: {self.model.method_name}\n"
               f"----------------------------\n"
               f"Azimuth MSE: {mse_az:.4f}\n"
               f"Elevation MSE: {mse_el:.4f}\n")
        
        lbl = ttk.Label(stats_frame, text=txt, font=("Courier", 12), justify="center")
        lbl.pack(pady=10)

    def show(self):
        self.root.mainloop()

# ==========================================
# RUN CODE
# ==========================================

def main():
    # Run Setup Dialog
    setup = SetupDialog()
    filepath, method, do_center, do_loop, confirmed = setup.show()

    if not confirmed:
        print("--------------------------")
        print("Analysis cancelled.")
        sys.exit()

    # Initialize Model & Process Data
    try:
        model = AntennaModel()
        
        print("--------------------------")
        print(f"Loading: {filepath}...")
        print(f"Settings: Auto-Center={do_center}, Loop-Closure={do_loop}")
        
        # PASS FLAGS TO LOAD_DATA
        model.load_data(filepath, do_center, do_loop)
        
        print("--------------------------")
        print(f"Running {method} Interpolation...")
        model.run_interpolation(method)
        
    except Exception as e:
        tk.messagebox.showerror("Processing Error", str(e))
        sys.exit()

    # Launch Results View
    print("--------------------------")
    print("Launching Visualization...")
    print("--------------------------")
    app = ResultsWindow(model)
    app.show()

if __name__ == "__main__":
    main()


    '''
    ===========================================================================
   Copyright (C) 2025  Paul Mola

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.   
   =========================================================================== 
    '''