import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def update_plot(*args):
    # Récupérer les valeurs des sliders
    amplitude = amplitude_var.get()
    frequency = frequency_var.get()
    phase = phase_var.get()

    # Générer la sinusoïde en fonction des paramètres
    time = np.linspace(0, 1, 1000)
    signal = amplitude * np.sin(2 * np.pi * frequency * time + np.radians(phase))

    # Mettre à jour le graphique
    ax.clear()
    ax.plot(time, signal, color='blue')
    ax.set_title('Sinusoïde Paramétrée')
    ax.set_xlabel('Temps')
    ax.set_ylabel('Amplitude')

    canvas.draw()

# Créer la fenêtre principale
root = tk.Tk()
root.title("Sinusoïde Paramétrée")

# Variables Tkinter pour stocker les valeurs des sliders
amplitude_var = tk.DoubleVar(value=1.0)
frequency_var = tk.DoubleVar(value=1.0)
phase_var = tk.DoubleVar(value=0.0)

# Créer les sliders pour l'amplitude, la fréquence et la phase
amplitude_label = ttk.Label(root, text="Amplitude")
amplitude_slider = ttk.Scale(root, from_=0, to=2, orient='horizontal', length=200, variable=amplitude_var, command=update_plot)

frequency_label = ttk.Label(root, text="Fréquence")
frequency_slider = ttk.Scale(root, from_=1, to=10, orient='horizontal', length=200, variable=frequency_var, command=update_plot)

phase_label = ttk.Label(root, text="Phase")
phase_slider = ttk.Scale(root, from_=0, to=360, orient='horizontal', length=200, variable=phase_var, command=update_plot)

# Positionner les sliders dans la fenêtre
amplitude_label.pack(pady=5)
amplitude_slider.pack(pady=10)
frequency_label.pack(pady=5)
frequency_slider.pack(pady=10)
phase_label.pack(pady=5)
phase_slider.pack(pady=10)

# Initialiser le graphique
fig, ax = plt.subplots(figsize=(5, 3))
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack()

# Mettre à jour le graphique initial
update_plot()

# Démarrer la boucle d'événements Tkinter
root.mainloop()
