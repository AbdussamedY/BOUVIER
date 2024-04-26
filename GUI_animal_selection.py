import tkinter as tk
from tkinter import ttk
import pandas as pd

def selection_change():
    global raw_data_path
    raw_data_path = selected_option.get()
    print(raw_data_path)

def save_figures_change():
    global Saving_boolean
    Saving_boolean = save_figures_var.get()

def show_outputs_change():
    global Showing_boolean
    Showing_boolean = show_outputs_var.get()

root = tk.Tk()
root.title("Select one animal for rotation in pitch dark analysis")


button_frame = ttk.Frame(root)
button_frame.pack(side='right', padx=10, pady=10)  

selected_option = tk.StringVar()


animal_list = pd.read_excel(r'C:\Users\gbouvier\Documents\census.xlsx', sheet_name='Study')['raw_data_path']

for i, animal in enumerate(animal_list):
    radio_button = tk.Radiobutton(root, text=animal, variable=selected_option, value=animal, command=selection_change)
    radio_button.pack(anchor='w')


save_figures_var = tk.BooleanVar()
save_figures_checkbox = ttk.Checkbutton(button_frame, text="Save figures", variable=save_figures_var, command=save_figures_change)
save_figures_checkbox.pack(anchor='w')


show_outputs_var = tk.BooleanVar()
show_outputs_checkbox = ttk.Checkbutton(button_frame, text="Show outputs", variable=show_outputs_var, command=show_outputs_change)
show_outputs_checkbox.pack(anchor='w')


Saving_boolean = False
Showing_boolean = False

root.mainloop()