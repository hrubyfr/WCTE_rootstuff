import uproot as ur
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk
import numpy as np
import random

# Open the file
file = ur.open("/home/wcte/Desktop/tuples/events_run0532_beamline_newFormat.root")

# File TTree:
# beam_monitor_info
## beamline_event
## beamline_tdc_nhits
## beamline_qdc_nhits
## beamline_tdc_time_id
## beamline_tdc_time
## beamline_qdc_charge_id
## beamline_qdc_charge

# 1D array with event number
beamline_event = file["beam_monitor_info"]["beamline_event"].array(library="np")
# 1D array with TDC number of hits for each event
beamline_tdc_nhits = file["beam_monitor_info"]["beamline_tdc_nhits"].array(library="np")
# 1D array with QDC number of hits for each event
beamline_qdc_nhits = file["beam_monitor_info"]["beamline_qdc_nhits"].array(library="np")

# 1D array of 1D arrays of variable size of time id for each TDC hit
beamline_tdc_time_id = file["beam_monitor_info"]["beamline_tdc_time_id"].array(library="np")
# 1D array of 1D arrays of variable size of time measured for each TDC hit
beamline_tdc_time = file["beam_monitor_info"]["beamline_tdc_time"].array(library="np")
# 1D array of 1D arrays of variable size of charge id for each QDC hit
beamline_qdc_charge_id = file["beam_monitor_info"]["beamline_qdc_charge_id"].array(library="np")
# 1D array of 1D arrays of variable size of time measured for each QDC hit
beamline_qdc_charge = file["beam_monitor_info"]["beamline_qdc_charge"].array(library="np")

beamline_tdc_time = beamline_tdc_time * 25e-3 # ns


# for i in range(4):
#     print(beamline_tdc_nhits[i],len(beamline_tdc_time_id[i]))

print("There was", len(beamline_event), " events in this file.")

# plt.hist(beamline_tdc_nhits, bins=40, alpha=0.7, color='blue', edgecolor='black')
# plt.title('Number of Total TDC hits per event')
# plt.xlabel('Value')
# plt.ylabel('Frequency')


# plt.hist(beamline_qdc_nhits, bins=40, alpha=0.7, color='red', edgecolor='black')
# plt.title('Number of Total QDC hits per event')
# plt.xlabel('Value')
# plt.ylabel('Frequency')


# List of plots
plots = [(beamline_tdc_nhits, 40,'Number of Total TDC hits per event','Total TDC hits', 'Frequency'),
         (beamline_qdc_nhits, 40,'Number of Total QDC hits per event','Total QDC Hits','Frequency')]
current_plot = 0

def get_random_color(): 
    return "#%06x" % random.randint(0, 0xFFFFFF)

def update_plot(index):
    ax.clear()
    color = get_random_color()
    ax.hist(plots[index][0], bins=plots[index][1], alpha=0.7, color=color, edgecolor='black')
    ax.set_title(plots[index][2])
    ax.set_xlabel(plots[index][3])
    ax.set_ylabel(plots[index][4])
    canvas.draw()

def on_next_button_clicked():
    global current_plot
    if current_plot < len(plots) - 1:
        current_plot += 1
        update_plot(current_plot)

def on_back_button_clicked():
    global current_plot
    if current_plot > 0:
        current_plot -= 1
        update_plot(current_plot)

def on_save_button_clicked(): 
    global current_plot 
    filename = f"plot_{current_plot}.png" 
    fig.savefig(filename, dpi=300) 
    print(f"Plot saved as {filename}")

def on_save_all_button_clicked(): 
    for i, plot in enumerate(plots): 
        ax.clear()
        color = get_random_color()
        ax.hist(plot[0], bins=plot[1], alpha=0.7, color=color, edgecolor='black') 
        ax.set_title(plot[2])
        ax.set_xlabel(plots[3])
        ax.set_ylabel(plots[4])
        filename = f"plot_{i}.png" 
        fig.savefig(filename, dpi=300) 
        print(f"Plot {i} saved as {filename}") 
        update_plot(current_plot)

def on_exit_button_clicked(): 
    root.quit()

# Create the main window
root = tk.Tk()
root.title("Interactive Plot Navigation")

# Create a figure and axis
fig, ax = plt.subplots()

# Create a canvas and add it to the Tkinter window
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Create navigation buttons
button_frame = ttk.Frame(root)
button_frame.pack(side=tk.BOTTOM, fill=tk.X)

back_button = ttk.Button(button_frame, text="Back", command=on_back_button_clicked)
back_button.pack(side=tk.LEFT)

next_button = ttk.Button(button_frame, text="Next", command=on_next_button_clicked)
next_button.pack(side=tk.RIGHT)

save_button = ttk.Button(button_frame, text="Save", command=on_save_button_clicked)
save_button.pack(side=tk.RIGHT)

save_all_button = ttk.Button(button_frame, text="Save All", command=on_save_all_button_clicked) 
save_all_button.pack(side=tk.RIGHT)

exit_button = ttk.Button(button_frame, text="Exit", command=on_exit_button_clicked) 
exit_button.pack(side=tk.RIGHT)

# Define button styles 
style = ttk.Style() 
style.configure("Back.TButton", background="lightblue", foreground="black") 
style.configure("Next.TButton", background="lightgreen", foreground="black") 
style.configure("Save.TButton", background="lightcoral", foreground="black") 
style.configure("SaveAll.TButton", background="lightyellow", foreground="black") 
style.configure("Exit.TButton", background="lightgrey", foreground="black")

# Initialize the first plot
update_plot(current_plot)

# Start the Tkinter main loop
root.mainloop()
