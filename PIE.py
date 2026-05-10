''' 
==============================================================
                                                              
  It's as easy as...                                          
                                                              
    $$$$$$$\  $$$$$$\ $$$$$$$$\ 
    $$  __$$\ \_$$  _|$$  _____|
    $$ |  $$ |  $$ |  $$ |      
    $$$$$$$  |  $$ |  $$$$$\    
    $$  ____/   $$ |  $$  __|   
    $$ |        $$ |  $$ |      
    $$ |      $$$$$$\ $$$$$$$$\ 
    \__|      \______|\________|                                       
                                                              
   PAUL'S INTERPOLATION ENGINE

 ============================================================== 
 '''

import os
import sys
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from scipy.interpolate import CubicSpline
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib.colors import Normalize

class AntennaModel:
    
    """
    Antenna pattern model and interpolation engine.

    ----------------
    pattern_3d has shape (360, 180), indexed as [az_index, el_index].
      - Axis 0 (rows)   : Azimuth (Phi),   0–359°, step = ANGULAR_RESOLUTION
      - Axis 1 (columns): Elevation (Theta), 0–179°, step = ANGULAR_RESOLUTION

    Input gain values should be in dB. Absolute gain reference is discarded; and data is peak-normalised.
    """

    ANGULAR_RESOLUTION = 1
    FULL_ROTATION = 360
    HALF_ROTATION = FULL_ROTATION // 2
    ELEVATION_RANGE = HALF_ROTATION
    MIN_DATA_POINTS = 10
    DEFAULT_K = 2.0
    DEFAULT_N = 5.0

    def __init__(self):
        self.filepath = None
        self.raw_az = None
        self.raw_el = None
        self.raw_theta_az = None
        self.raw_theta_el = None
        self.az_norm = None
        self.el_norm = None
        self.theta_norm = np.deg2rad(np.arange(0, self.FULL_ROTATION, self.ANGULAR_RESOLUTION))
        self.is_loop_closed = False
        self.pattern_3d = None
        self.x = None
        self.y = None
        self.z = None
        self.method_name = ""

    def load_data(self, filepath, do_autocenter, do_loop_closure):
            self.filepath = filepath
            self.is_loop_closed = do_loop_closure
            try:
                raw_lines = []
                encodings_to_try = ['utf-8-sig', 'utf-8', 'utf-16', 'latin-1']
                
                for encoding in encodings_to_try:
                    try:
                        with open(filepath, 'r', encoding=encoding) as f:
                            raw_lines = f.readlines()
                        if raw_lines and '\0' in raw_lines[0]:
                            raise UnicodeError("Decoded content contains null bytes, likely wrong encoding.")
                        break
                    except UnicodeError:
                        continue
                
                if not raw_lines:
                    raise ValueError("Could not decode file with standard encodings.")

                data = [float(line.strip()) for line in raw_lines if line.strip()]
                
                if len(data) < self.MIN_DATA_POINTS: 
                    raise ValueError("File contains insufficient data.")

                # Assume file is 50% Azimuth, 50% Elevation
                midpoint = len(data) // 2
                self.raw_az = np.array(data[:midpoint])
                self.raw_el = np.array(data[midpoint:])
                
                if do_autocenter:
                    # Align main lobe (peak datapoint) peaks to 0°. Only use if Az/El peaks should coincide.
                    if len(self.raw_az) > 0:
                        shift_az = -np.argmax(self.raw_az)
                        self.raw_az = np.roll(self.raw_az, shift_az)
                    
                    if len(self.raw_el) > 0:
                        shift_el = -np.argmax(self.raw_el)
                        self.raw_el = np.roll(self.raw_el, shift_el)
                
                if do_loop_closure:
                    # Force first/last continuity for low-quality measurements
                    if len(self.raw_az) > 0 and self.raw_az[0] != self.raw_az[-1]:
                        self.raw_az = np.append(self.raw_az, self.raw_az[0])

                    if len(self.raw_el) > 0 and self.raw_el[0] != self.raw_el[-1]:
                        self.raw_el = np.append(self.raw_el, self.raw_el[0])

                use_endpoint = do_loop_closure
                self.raw_theta_az = np.linspace(0, 2*np.pi, len(self.raw_az), endpoint=use_endpoint)
                self.raw_theta_el = np.linspace(0, 2*np.pi, len(self.raw_el), endpoint=use_endpoint)
                
                self._normalize_inputs()

            except Exception as e:
                raise RuntimeError(f"Failed to load file: {e}")

    def _normalize_inputs(self):
        """
        Resample the raw azimuth and elevation cuts onto a uniform 1°-step grid
        using a Periodic Cubic Spline.
        """
        target_deg = np.arange(0, self.FULL_ROTATION, self.ANGULAR_RESOLUTION)

        def _periodic_spline(raw):
            if self.is_loop_closed:
                # raw[0] == raw[-1]; index spans [0°, 360°] with n points.
                idx = np.linspace(0, self.FULL_ROTATION, len(raw), endpoint=True)
                data = raw
            else:
                # Close the loop by appending raw[0] at 360°.
                idx = np.linspace(0, self.FULL_ROTATION, len(raw) + 1, endpoint=True)
                data = np.append(raw, raw[0])

            cs = CubicSpline(idx, data, bc_type='periodic')
            return cs(target_deg)

        self.az_norm = _periodic_spline(self.raw_az)
        self.el_norm = _periodic_spline(self.raw_el)

    def run_interpolation(self, method, k=DEFAULT_K, n=DEFAULT_N):
        self.method_name = method
        g_az = self.az_norm
        g_el = self.el_norm 

        if method == 'Summing':
            self.pattern_3d = self._algo_summing(g_az, g_el)
        elif method == 'Approximation':
            self.pattern_3d = self._algo_approx(g_az, g_el, k=k)
        elif method == 'Hybrid':
            self.pattern_3d = self._algo_hybrid(g_az, g_el, k=k, n=n)
        else:
            raise ValueError("Unknown method")

        self._spherical_to_cartesian()

    def _spherical_to_cartesian(self):
        azimuth = np.deg2rad(np.arange(0, self.FULL_ROTATION, self.ANGULAR_RESOLUTION))
        elevation = np.deg2rad(np.arange(0, self.ELEVATION_RANGE, self.ANGULAR_RESOLUTION))
        
        norm_val = -np.min(self.pattern_3d) if np.any(self.pattern_3d < 0) else 0
        r = self.pattern_3d + norm_val
        
        el_grid, az_grid = np.meshgrid(elevation, azimuth, indexing='xy')
        
        self.x = r * np.sin(el_grid) * np.cos(az_grid)
        self.y = r * np.sin(el_grid) * np.sin(az_grid)
        self.z = r * np.cos(el_grid)

    def export_csv(self, filename):
        """
        Export 3D pattern with peak aligned to Phi = 0°
        Phi: Azimuth (0 to 359)
        Theta: Elevation (0 to 179) 
        """

        # Check if data exists before proceeding
        if self.pattern_3d is None:
            raise RuntimeError("Export failed: No interpolation data available. Please run interpolation first.")

        try:
            # Find peak to align Phi = 0° to the maximum gain point
            peak_flat_idx = np.argmax(self.pattern_3d)
            peak_az_idx, _ = np.unravel_index(peak_flat_idx, self.pattern_3d.shape)
            
            phi_vals = np.arange(0, self.FULL_ROTATION, self.ANGULAR_RESOLUTION)    
            theta_vals = np.arange(0, self.ELEVATION_RANGE, self.ANGULAR_RESOLUTION)

            # Re-index azimuth so that the peak is at index 0 (Phi=0)
            az_indices = (peak_az_idx + phi_vals) % self.FULL_ROTATION

            az_grid, el_grid = np.meshgrid(az_indices, theta_vals, indexing='xy')
            phi_grid, theta_grid = np.meshgrid(phi_vals, theta_vals, indexing='xy')

            # Extract gains using the re-indexed grid
            gains = self.pattern_3d[az_grid.T, el_grid.T] 
            
            flat_gains = gains.flatten()
            flat_phi = phi_grid.T.flatten()
            flat_theta = theta_grid.T.flatten()

            header = "Phi[deg],Theta[deg],dB10normalize(GainTotal)"
            data_str = [f"{p},{t},{g:.4e}" for p, t, g in zip(flat_phi, flat_theta, flat_gains)]
            
            with open(filename, 'w') as f:
                f.write(f"{header}\n")
                f.write("\n".join(data_str))
                        
        except Exception as e:
            raise RuntimeError(f"Export failed: {e}")

    def _algo_summing(self, g_az, g_el):

        """
        Adds the logarithmic (dB) azimuth and elevation cuts to produce a full
        3D pattern. Assumes pattern separability.

        - Physical elevation (theta) only spans 0–180°, but the input data
          covers 0–360°. The second half (180–359°) is used to fill the back
          hemisphere of the 3D pattern.
        - In the back hemisphere, the elevation profile is geometrically
          mirrored: a ray at azimuth (180+phi) and elevation theta is equivalent
          to the same elevation traversed in the opposite direction.
        - Reversing g_el for the second azimuth half (180–359°) enforces this
          symmetry, ensuring the elevation profile maps correctly across both
          hemispheres.
        """

        pattern = np.zeros((self.FULL_ROTATION, self.ELEVATION_RANGE))
        
        mid = self.HALF_ROTATION
        
        pattern[:mid, :] = g_az[:mid, np.newaxis] + g_el[:mid]
        pattern[mid:, :] = g_az[mid:, np.newaxis] + g_el[mid:][::-1]
        return pattern

    def _algo_approx(self, g_az, g_el, k=DEFAULT_K):

        """
        Summing algorithm with gain-weighted blending (p-norm blend).

        - w1 and w2 are blending weights derived from the linearised gain values
          of the elevation and azimuth cuts respectively. They represent how much
          each cut contributes at a given point in the pattern.
        - The denominator (w1**k + w2**k)**(1/k) is a p-norm (where p = k)
          applied to the weight vector [w1, w2], and is a generalisation of
          familiar norms:
            - k=1  : sum of weights (linear blend of the two cuts)
            - k=2  : Uses Euclidean blending (square root of the sum of squares)
            - k=inf  : The cut with the higher gain at that specific angle will be favored.
        - Dividing the weighted dB sum by this p-norm normalises the result so
          it stays in the same scale as the inputs. Tuning k controls how the
          two cuts compete: low k produces an even blend, high k increasingly
          favours whichever cut has the larger weight.
        - Where both weights are simultaneously zero (both cuts at unity linear
          gain), the output is set to zero. This occurs only at the pattern peak
          and has negligible effect on the overall reconstruction.
        """

        pattern = np.zeros((self.FULL_ROTATION, self.ELEVATION_RANGE))
        
        vert = 10**(g_el/10)
        hor = 10**(g_az/10)
        
        mid = self.HALF_ROTATION
        
        def calc_segment(db_h, lin_h, db_v, lin_v):
            w1 = lin_v * (1 - lin_h[:, np.newaxis])
            w2 = lin_h[:, np.newaxis] * (1 - lin_v)
            num = db_h[:, np.newaxis] * w1 + db_v * w2
            den = (w1**k + w2**k)**(1/k)
            mask = (w1 == 0) & (w2 == 0)
            return np.where(mask, 0, num/den)

        pattern[:mid, :] = calc_segment(g_az[:mid], hor[:mid], g_el[:mid], vert[:mid])
        pattern[mid:, :] = calc_segment(g_az[mid:], hor[mid:], g_el[mid:][::-1], vert[mid:][::-1])
        return pattern

    def _algo_hybrid(self, g_az, g_el, k=DEFAULT_K, n=DEFAULT_N):

        """
        Calculates a 3D pattern by cross-fading between the Summing and 
        Approximation methods based on local gain intensity.

        - The weighting factor 'W' is derived from the Summing method's 
          linearized gain, modified by the power-law exponent (1/n). 
            - High n: Summing' method is dominant (better for
            preserving peak accuracy).
            - Low n: Causes a faster transition to the 'Approximation' method 
              as gain drops (better for smoothing noisy sidelobes or nulls).
        """
        
        approx = self._algo_approx(g_az, g_el, k=k)
        summing = self._algo_summing(g_az, g_el)
        sum_dec = 10**(summing/10)
        return summing * (sum_dec**(1/n)) + approx * (1 - sum_dec**(1/n))


class PlotPanel(tk.Frame):
    FIG_SIZE = (3, 2)
    FIG_DPI = 100
    SAVE_DPI = 300
    TITLE_FONT_SIZE = 10

    def __init__(self, parent, title, projection=None):
        super().__init__(parent)

        self.figure = Figure(figsize=self.FIG_SIZE, dpi=self.FIG_DPI)
        self.figure.set_tight_layout(True)

        if projection == '3d':
            self.ax = self.figure.add_subplot(111, projection='3d')
        elif projection == 'polar':
            self.ax = self.figure.add_subplot(111, projection='polar')
        else:
            self.ax = self.figure.add_subplot(111)
            
        self.ax.set_title(title, fontsize=self.TITLE_FONT_SIZE)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.get_tk_widget().bind("<Button-3>", self._show_context_menu)
        
    def _show_context_menu(self, event):
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(label="Save Image As...", command=self._save_figure)
        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()
    
    def _save_figure(self):
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
                self.figure.savefig(filename, dpi=self.SAVE_DPI, bbox_inches='tight')
                messagebox.showinfo("Success", f"Image saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Save Failed", f"{e}")

class SetupDialog:
    def __init__(self, initial_file=None):
        self.root = tk.Tk()
        self.root.title("Antenna Interpolation Setup")
        self.root.geometry("500x440")
        
        self.filepath = initial_file
        self.method = tk.StringVar(value="Summing")
        self.var_k = tk.StringVar() 
        self.var_n = tk.StringVar() 
        self.var_autocenter = tk.BooleanVar(value=True)
        self.var_loop_closure = tk.BooleanVar(value=True)
        self.confirmed = False
        
        self._build_ui()
        
        if self.filepath:
            self.lbl_file.config(text=os.path.basename(self.filepath), foreground="black")

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

        frame_weights = ttk.LabelFrame(self.root, text="Interpolation Weights")
        frame_weights.pack(fill="x", padx=10, pady=5)
        
        f_k = ttk.Frame(frame_weights)
        f_k.pack(fill="x", padx=5, pady=2)
        ttk.Label(f_k, text="Approximation Weight (k):").pack(side="left")
        ttk.Entry(f_k, textvariable=self.var_k, width=10).pack(side="right")
        
        f_n = ttk.Frame(frame_weights)
        f_n.pack(fill="x", padx=5, pady=2)
        ttk.Label(f_n, text="Hybrid Weight (n):").pack(side="left")
        ttk.Entry(f_n, textvariable=self.var_n, width=10).pack(side="right")

        lbl_warn = ttk.Label(frame_weights, 
                             text="Leave unspecified to use default values.", 
                             foreground="red", font=("Arial", 8, "italic"), justify="center")
        lbl_warn.pack(pady=5)

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
        
        k_val = self.var_k.get().strip()
        n_val = self.var_n.get().strip()
        
        try:
            k_out = float(k_val) if k_val else AntennaModel.DEFAULT_K
            n_out = float(n_val) if n_val else AntennaModel.DEFAULT_N
            if k_out <= 0 or n_out <= 0:
                messagebox.showerror("Input Error", "Weights k and n must be greater than 0.")
                self.confirmed = False
                return None 
            
        except ValueError:
            k_out = AntennaModel.DEFAULT_K
            n_out = AntennaModel.DEFAULT_N

        if not self.confirmed:
            return None

        return (self.filepath, self.method.get(), 
                self.var_autocenter.get(), self.var_loop_closure.get(), 
                k_out, n_out, self.confirmed)

class ResultsWindow:
    def __init__(self, model):
        self.model = model
        self.root = tk.Tk()
        self.root.title(f"Results - {model.method_name} Method")
        self.root.geometry("1100x600") 
        self._init_top_bar()
        
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True)
        
        self._init_tab_raw()
        self._init_tab_input()
        self._init_tab_3d()

    def _init_top_bar(self):
        top_frame = ttk.Frame(self.root)
        top_frame.pack(side="top", fill="x", pady=5, padx=10)
        btn_back = ttk.Button(top_frame, text="< Back to Setup", command=self._on_back)
        btn_back.pack(side="left")

    def _on_back(self):
        prev_file = self.model.filepath
        self.root.destroy()
        main(initial_file=prev_file)

    def _init_tab_raw(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Original (Raw)')
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=1)

        p1 = PlotPanel(frame, "Raw Azimuth Planar Cut [dB]", projection='polar')
        p1.grid(row=0, column=0, sticky="nsew")
        p1.ax.set_theta_zero_location('N')
        p1.ax.set_theta_direction(-1)
        p1.ax.plot(self.model.raw_theta_az, self.model.raw_az, color='g',
                   marker='.', linestyle='none', linewidth=0)

        p2 = PlotPanel(frame, "Raw Elevation Planar Cut [dB]", projection='polar')
        p2.grid(row=0, column=1, sticky="nsew")
        p2.ax.set_theta_zero_location('N')
        p2.ax.set_theta_direction(-1)
        p2.ax.plot(self.model.raw_theta_el, self.model.raw_el, color='g',
                   marker='.', linestyle='none', linewidth=0)

    def _init_tab_input(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Original (Smoothed)')
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=1)

        p1 = PlotPanel(frame, "Azimuth Planar Cut [dB]", projection='polar')
        p1.grid(row=0, column=0, sticky="nsew")
        p1.ax.set_theta_zero_location('N')
        p1.ax.set_theta_direction(-1)
        p1.ax.plot(self.model.theta_norm, self.model.az_norm, color='g', label='Input')
        
        p2 = PlotPanel(frame, "Elevation Planar Cut [dB]", projection='polar')
        p2.grid(row=0, column=1, sticky="nsew")
        p2.ax.set_theta_zero_location('N')
        p2.ax.set_theta_direction(-1)
        p2.ax.plot(self.model.theta_norm, self.model.el_norm, color='g', label='Input')

    def _init_tab_3d(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text='Reconstructed 3D Pattern')
        
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.rowconfigure(0, weight=0)
        frame.rowconfigure(1, weight=1)

        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        btn_export = ttk.Button(btn_frame, text="Export 3D Pattern as CSV", command=self._export_csv)
        btn_export.pack(anchor="center")

        max_gain = np.max(self.model.pattern_3d)
        plot_data = self.model.pattern_3d - max_gain

        # Wrap arrays to close the azimuthal seam
        plot_data_wrapped = np.concatenate([plot_data, plot_data[:1, :]], axis=0)
        x_wrapped = np.concatenate([self.model.x, self.model.x[:1, :]], axis=0)
        y_wrapped = np.concatenate([self.model.y, self.model.y[:1, :]], axis=0)
        z_wrapped = np.concatenate([self.model.z, self.model.z[:1, :]], axis=0)

        vmin = np.min(plot_data_wrapped)
        vmax = np.max(plot_data_wrapped)
        norm = Normalize(vmin=vmin, vmax=vmax)
        surface_colors = cm.nipy_spectral(norm(plot_data_wrapped))

        p3d = PlotPanel(frame, "3D Polar Plot", projection='3d')
        p3d.grid(row=1, column=0, sticky="nsew")
        
        p3d.ax.plot_surface(
            x_wrapped, y_wrapped, z_wrapped, 
            facecolors=surface_colors,
            rstride=2, cstride=2,
            linewidth=0, antialiased=False, shade=False
        )
        
        m = cm.ScalarMappable(cmap=cm.nipy_spectral, norm=norm)
        m.set_array([])
        p3d.figure.colorbar(m, ax=p3d.ax, shrink=0.5, aspect=5, label='Normalized Gain [dB]')
        p3d.ax.axis('off')

        # --- Polar heatmap: Phi (azimuth) as angular axis, Theta (elevation) as radial axis ---
        p_2d = PlotPanel(frame, "2D Polar Heatmap", projection='polar')
        p_2d.grid(row=1, column=1, sticky="nsew")

        # Build coordinate grids
        # phi: 0–360° (azimuth), theta: 0–180° (elevation)
        phi_deg = np.arange(0, self.model.FULL_ROTATION + 1, self.model.ANGULAR_RESOLUTION)
        theta_deg = np.arange(0, self.model.ELEVATION_RANGE, self.model.ANGULAR_RESOLUTION)

        phi_rad = np.deg2rad(phi_deg)
        # plot_data_wrapped has shape (361, 180): axis0=phi, axis1=theta
        phi_grid, theta_grid = np.meshgrid(phi_rad, theta_deg, indexing='ij')

        p_2d.ax.pcolormesh(
            phi_grid, theta_grid, plot_data_wrapped,
            cmap=cm.nipy_spectral,
            vmin=vmin, vmax=vmax,
            shading='auto'
        )

        # Radial axis: theta goes 0 (centre) → 180° (outer edge)
        p_2d.ax.set_theta_zero_location('N')
        p_2d.ax.set_theta_direction(-1)
        p_2d.ax.set_rlabel_position(25)
        p_2d.ax.set_rlim(0, self.model.ELEVATION_RANGE)
        p_2d.ax.set_rticks([0, 30, 60, 90, 120, 150, 180])
        p_2d.ax.set_yticklabels(['0°', '30°', '60°', '90°', '120°', '150°', '180°'], fontsize=6.8)

        # Add a colorbar via the figure
        sm = cm.ScalarMappable(cmap=cm.nipy_spectral, norm=norm)
        sm.set_array([])
        p_2d.figure.colorbar(sm, ax=p_2d.ax, shrink=0.6, aspect=10, pad=0.15, label='Normalized Gain [dB]')

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
                messagebox.showerror("Export Failed", f"{e}")

    def show(self):
        self.root.mainloop()


def main(initial_file=None):
    setup = SetupDialog(initial_file=initial_file)
    result = setup.show()
    
    if result is None:
        sys.exit()
        
    filepath, method, do_center, do_loop, k_val, n_val, confirmed = result

    if not confirmed:
        print("--------------------------")
        print("Analysis cancelled.")
        sys.exit()

    try:
        model = AntennaModel()
        
        print("--------------------------")
        print(f"Loading: {filepath}...")
        print(f"Settings: Auto-Center={do_center}, Loop-Closure={do_loop}")
        
        model.load_data(filepath, do_center, do_loop)
        
        print("--------------------------")
        print(f"Running {method} Interpolation...")
        print(f"Weights: k={k_val}, n={n_val}")
        
        model.run_interpolation(method, k=k_val, n=n_val)
        
    except Exception as e:
        tk.messagebox.showerror("Processing Error", f"{e}")
        sys.exit()

    print("--------------------------")
    print("Launching Visualization...")
    app = ResultsWindow(model)
    app.show()

if __name__ == "__main__":
    main()


# ===========================================================================
# Copyright (C) 2025  Paul Mola
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# ===========================================================================