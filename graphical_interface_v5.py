import tkinter as tk
from tkinter import ttk
import numpy as np
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Π = lambda t: 1.0 if np.abs(t)<.5 else 0.0 # Unitary rectangular window

# def eye_diagram(data, sampling_rate,  time_eye_window_max):
#     data_duration = len(data)/sampling_rate
#     data_time = np.arange(0,data_duration, 1/sampling_rate)
#     window = np.arange(0, time_eye_window_max, 1/sampling_rate)
#     window[-1]=np.nan
#     time_eye_value = np.array([window for i in range(int(np.floor(data_duration/time_eye_window_max)))]).flatten()
#     time_eye_value = np.concatenate((time_eye_value, window[0 : len(data_time) - len(time_eye_value)]))

#     return time_eye_value


class PropOptGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Optical propagation GUI")

        # Frames
        self.control_frame = ttk.Frame(self.master)
        self.control_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ns")

        self.plot_frame = ttk.Frame(self.master)
        self.plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

        # Longueur
        self.fiber_lenght_label = ttk.Label(self.control_frame, text="Fiber lenght [ km ]:")
        self.fiber_lenght_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")

        self.fiber_lenght_slider = ttk.Scale(self.control_frame, from_=0, to=100, orient="horizontal", command=self.update_plot)
        self.fiber_lenght_slider.grid(row=0, column=1, padx=10, pady=10, sticky="ew")

        self.fiber_lenght_value_label = ttk.Label(self.control_frame, text="")
        self.fiber_lenght_value_label.grid(row=0, column=2, padx=10, pady=10, sticky="w")

        # Dispersion coefficient
        self.dispersion_coef_label = ttk.Label(self.control_frame, text="Dispersion [ ps/(km.nm) ]:")
        self.dispersion_coef_label.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        self.dispersion_coef_slider = ttk.Scale(self.control_frame, from_=0, to=100, orient="horizontal", command=self.update_plot)
        self.dispersion_coef_slider.grid(row=1, column=1, padx=10, pady=10, sticky="ew")

        self.dispersion_coef_value_label = ttk.Label(self.control_frame, text="")
        self.dispersion_coef_value_label.grid(row=1, column=2, padx=10, pady=10, sticky="w")

        # Débit
        self.bit_rate_label = ttk.Label(self.control_frame, text="Débit [ kbit/s ]:")
        self.bit_rate_label.grid(row=2, column=0, padx=10, pady=10, sticky="w")

        self.bit_rate_entry = ttk.Entry(self.control_frame)
        self.bit_rate_entry.grid(row=2, column=1, padx=10, pady=10, sticky="ew")
        self.bit_rate_entry.bind("<Return>", self.update_on_enter)

        # self.bit_rate_value_label = ttk.Label(self.control_frame, text="")
        # self.bit_rate_value_label.grid(row=2, column=2, padx=10, pady=10, sticky="w")

        # Facteur de surechantillonage
        self.sursampling_factor_label = ttk.Label(self.control_frame, text="Facteur de suréch :")
        self.sursampling_factor_label.grid(row=3, column=0, padx=10, pady=10, sticky="w")

        self.sursampling_factor_entry = ttk.Entry(self.control_frame)
        self.sursampling_factor_entry.grid(row=3, column=1, padx=10, pady=10, sticky="ew")
        self.sursampling_factor_entry.bind("<Return>", self.update_on_enter)

        self.sursampling_factor_value_label = ttk.Label(self.control_frame, text="")
        self.sursampling_factor_value_label.grid(row=3, column=2, padx=10, pady=10, sticky="w")

        # Wavelenght
        self.lambda_label = ttk.Label(self.control_frame, text="λ [ µm ]:")
        self.lambda_label.grid(row=4, column=0, padx=10, pady=10, sticky="w")

        self.lambda_entry = ttk.Entry(self.control_frame)
        self.lambda_entry.grid(row=4, column=1, padx=10, pady=10, sticky="ew")
        self.lambda_entry.bind("<Return>", self.update_on_enter)

        self.lambda_value_label = ttk.Label(self.control_frame, text="")
        self.lambda_value_label.grid(row=4, column=2, padx=10, pady=10, sticky="w")

        # Max Time
        # self.max_time_label = ttk.Label(self.control_frame, text="Temps Max:")
        # self.max_time_label.grid(row=5, column=0, padx=10, pady=10, sticky="w")

        # self.max_time_entry = ttk.Entry(self.control_frame)
        # self.max_time_entry.grid(row=5, column=1, padx=10, pady=10, sticky="ew")
        # self.max_time_entry.insert(0, "6.28")  # Default value

        # self.max_time_value_label = ttk.Label(self.control_frame, text="")
        # self.max_time_value_label.grid(row=8, column=0, padx=10, pady=10, sticky="w")

        # Max Longueur
        self.max_fiber_len_label = ttk.Label(self.control_frame, text="Longueur max [ km ]:")
        self.max_fiber_len_label.grid(row=5, column=0, padx=10, pady=10, sticky="w")

        self.max_fiber_len_entry = ttk.Entry(self.control_frame)
        self.max_fiber_len_entry.grid(row=5, column=1, padx=10, pady=10, sticky="ew")
        self.max_fiber_len_entry.insert(0, "100.0")  # Default value
        self.max_fiber_len_entry.bind("<Return>", self.update_on_enter)
        self.max_fiber_len_entry.bind("<Return>", self.update_fiber_len_slider)

        # self.max_amplitude_value_label = ttk.Label(self.control_frame, text="")
        # self.max_amplitude_value_label.grid(row=10, column=0, padx=10, pady=10, sticky="w")

        # Create plots
        self.fig, ((self.ax1, self.ax2), (self.ax3, self.ax4)) = plt.subplots(2, 2, figsize=(10, 6))

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initial settings
        self.fiber_lenght_slider.set(10)
        self.dispersion_coef_slider.set(17)
        self.bit_rate_entry.insert(0, "1.0")
        self.lambda_entry.insert(0,"1.55")
        self.sursampling_factor_entry.insert(0, "4.0")

        self.update_plot()

    def update_plot(self, *args):
        L = self.fiber_lenght_slider.get() * u.km 
        D = self.dispersion_coef_slider.get() * u.ps / (u.nm * u.km)
        try:
            B = float(self.bit_rate_entry.get()) * u.kbit / u.s
            F = float(self.sursampling_factor_entry.get())
            #max_time = float(self.max_time_entry.get())
            λ = float(self.lambda_entry.get()) * u.um
            max_fiber_len = float(self.max_fiber_len_entry.get())
        except ValueError:
            B = 1.0 * u.kbit / u.s
            F = 4
            #max_time = 6.28
            λ = 1.55 * u.um
            max_fiber_len = 100.0

                
        # Physical parameters
        β_2 = - D * λ**2 / (2 * np.pi * const.c)    # Group velocity dispersion parameter
        B_nounit = B.to(u.bit/u.s)/(u.bit/u.s)      # Capacity Hertz unit


        I_max = 40 * u.mA # maximum intensity value of diode laser 
        I_nounit = I_max.to(u.A)/(1 * u.A)


        bits_seq = np.array([0,1,0]) # np.random.randint(2, size=70) #     
        bits_sended_lenght = len(bits_seq)

        τ = .375 / B # caracteristic time of the gaussian pulse
        τ_nounit = τ.to(u.s/u.bit)/(1 * u.s/u.bit)


        # Pulse sequence
        duration = 2 * bits_sended_lenght * τ_nounit
        sampling_rate = F * B_nounit 
        time = np.arange(0, duration, 1/sampling_rate)


        Π = lambda t: 1.0 if np.abs(t)<.5 else 0.0 # Unitary rectangular window
        gaussian_pulse = np.sqrt(I_nounit) * np.array([np.sum([bits_seq[k] * Π((ti * B_nounit - k)) for k in range(bits_sended_lenght)]) for ti in time]) 
        #gaussian_pulse = np.sqrt(I_nounit) * np.array([np.sum([bits_seq[k] * np.exp(-((ti-k/B_nounit)/τ_nounit)**2) for k in range(bits_sended_lenght)]) for ti in time])
        #rect_pulse = np.sqrt(I_nounit) * np.array([np.sum([bits_seq[k] * Π((ti * B_nounit - k)) for k in range(bits_sended_lenght)]) for ti in time])



        # Pulse after propagation
        freq = np.fft.fftfreq(len(time), d=1/sampling_rate)

        σ2 = .005 #
        noise = np.random.normal(loc=0, scale=σ2, size=len(gaussian_pulse))
        #ω = 2 * np.pi * const.c/λ
        ω = 2 * np.pi * freq * u.Hz
        fft_gauss_pulse_prop = np.fft.fft(gaussian_pulse + noise) * np.exp(1j * (β_2 * ω**2 * L/2).to(u.dimensionless_unscaled))
        #fft_rect_pulse_prop = np.fft.fft(rect_pulse) * np.exp(1j * (β_2 * ω**2 * L/2).to(u.dimensionless_unscaled))

        gauss_pulse_prop = np.abs(np.fft.ifft(fft_gauss_pulse_prop))

        ## Eye diagram
        #time_eye = eye_diagram(gaussian_pulse, sampling_rate, 4 * 1/B_nounit)

        time_eye_window_max = 4 * 1/B_nounit
        data_duration = len(gaussian_pulse)/sampling_rate
        data_time = np.arange(0,data_duration, 1/sampling_rate)
        window = np.arange(0, time_eye_window_max, 1/sampling_rate)
        window[-1]=np.nan
        time_eye = np.array([window for i in range(int(np.floor(data_duration/time_eye_window_max)))]).flatten()
        time_eye = np.concatenate((time_eye, window[0 : len(data_time) - len(time_eye)]))


        self.ax1.clear()
        self.ax1.plot(time_eye * B_nounit, gaussian_pulse)
        self.ax1.set_xlabel(r"t/Ts")
        self.ax1.set_ylabel(r"Before propagation")
        #self.ax1.set_xlim([0, max_time])
        #self.ax1.set_ylim([0, max_amplitude])
        self.ax1.grid(True)
        #self.ax1.legend()

        # self.ax2.clear()
        # self.ax2.psd(gaussian_pulse, Fs=sampling_rate.value, sides='twosided')
        # self.ax2.set_xlabel(r'$\nu$ (Hz)')
        # self.ax2.set_ylabel(r"P_{xx}")
        # self.ax2.set_xlim([0, max_time])
        # self.ax2.set_ylim([-max_amplitude, max_amplitude])
        #self.ax2.legend()

        # Plot sinus spectrum

        self.ax3.clear()
        self.ax3.plot(time_eye * B_nounit, gauss_pulse_prop)
        self.ax3.set_xlabel(r"t/Ts")
        self.ax3.set_ylabel(r"After propagation")
        self.ax3.grid(True)
        #self.ax3.legend()


        # self.ax4.clear()
        # self.ax4.psd(gauss_pulse_prop, Fs=sampling_rate.value, sides='twosided')
        # self.ax4.set_xlabel(r'$\nu$ (Hz)')
        # self.ax4.set_ylabel(r"P_{xx}")
        #self.ax4.legend()

        self.canvas.draw()

        # Update slider values
        self.fiber_lenght_value_label.config(text=f" {L.value:.2f}")
        self.dispersion_coef_value_label.config(text=f" {D.value:.2f}")
        #self.frequency_value_label.config(text=f" {frequency:.2f}")
        # self.max_time_value_label.config(text=f"Temps Max: {max_time:.2f}")
        # self.max_amplitude_value_label.config(text=f"Amplitude Max: {max_amplitude:.2f}")

        # Update amplitude slider maximum value
        self.fiber_lenght_slider.configure(to=max_fiber_len)

    def update_on_enter(self, event):
        self.update_plot()

    def update_fiber_len_slider(self, event):
        max_fiber_len = float(self.max_fiber_len_entry.get())
        self.fiber_lenght_slider.configure(to=max_fiber_len)
        self.update_plot()

if __name__ == "__main__":
    root = tk.Tk()
    app = PropOptGUI(root)
    root.mainloop()
