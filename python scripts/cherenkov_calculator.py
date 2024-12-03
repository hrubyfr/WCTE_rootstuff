import tkinter as tk
from tkinter import ttk, messagebox
import math

# Constants
c = 299792458  # Speed of light in vacuum, in meters per second (m/s)
u_mass = 105.6
e_mass = 0.511
pi_mass = 139.57

# Define particles and their masses in kg
particles = [
    ("Electron", 0.511),
    ("Muon", 105.6),
    ("Pion", 139.57),
    ("Proton", 937),
    # You can add more particles here
]

def calculate_speed(momentum, mass):
    """Calculate the speed of a particle given its momentum and rest mass."""
    speed = c * math.sqrt(1 - (mass**2) / (momentum**2 + mass**2))
    return speed
def calculate_beta(speed):
    beta = speed / c
    return beta
def calculate_nthr(beta):
    n_thr = 1/beta
    return n_thr
def calculate_speeds(event=None):
    """Calculate and display speeds for all predefined particles."""
    try:
        momentum = float(momentum_entry.get())
    except ValueError:
        messagebox.showerror("Invalid Input", "Please enter a valid number for momentum.")
        return
    
    # Clear previous results from the table
    for row in results_table.get_children():
        results_table.delete(row)
    
    # Calculate and insert each particle's speed into the table
    for name, mass in particles:
        speed = calculate_speed(momentum, mass)
        beta = calculate_beta(speed)
        n_thr = calculate_nthr(beta)
        results_table.insert("", "end", values=(name, f"{mass:.2f}", f"{speed:.2e}", f"{beta:.2f}", f"{n_thr:.4f}"))

def initialize_table():
    """Initialize the table with particle names and masses."""
    for name, mass in particles:
        results_table.insert("", "end", values=(name, f"{mass:.2f}", "", "", ""))


# Create main window
root = tk.Tk()
root.title("Particle Speed Calculator")

# Momentum Input
tk.Label(root, text="Beam Momentum (MeV/c):").grid(row=0, column=0, sticky="e")
momentum_entry = tk.Entry(root)
momentum_entry.grid(row=0, column=1, pady=5)

momentum_entry.bind("<Return>", calculate_speeds)

# Calculate Button
calculate_button = tk.Button(root, text="Calculate Speeds", command=calculate_speeds)
calculate_button.grid(row=1, column=0, columnspan=2, pady=10)

# Results Table
columns = ("name", "mass", "speed", "beta", "n_thresh")
results_table = ttk.Treeview(root, columns=columns, show="headings")
results_table.heading("name", text="Particle Name")
results_table.heading("mass", text="Mass (MeV)")
results_table.heading("speed", text="Speed (m/s)")
results_table.heading("beta", text="Beta")
results_table.heading("n_thresh", text="Threshold index")

# Display the table
results_table.grid(row=2, column=0, columnspan=2, padx=10, pady=10, sticky="nsew")

initialize_table()
# Run the application
root.mainloop()

