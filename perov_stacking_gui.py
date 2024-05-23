#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:16:42 2024

@author: mhaefner-chem

Stacking of Perovskites

TODO
seperate fields for element entry
weird doubling error for cif cells (doubled hh cell with double lattice parameter and atoms)
"""

class positions:
    def __init__(self):
        # define the possible layers
        self.layers = {}
        self.layers["A"] = [["A",0/6,0/6],["X",3/6,3/6],["X",3/6,0/6],["X",0/6,3/6]]
        self.layers["B"] = [["A",4/6,2/6],["X",4/6,5/6],["X",1/6,5/6],["X",1/6,2/6]]
        self.layers["C"] = [["A",2/6,4/6],["X",2/6,1/6],["X",5/6,1/6],["X",5/6,4/6]]
    
        # define the interlayer positions
        self.layers["a"] = [["B",0/6,0/6]]
        self.layers["b"] = [["B",4/6,2/6]]
        self.layers["c"] = [["B",2/6,4/6]]

def check_repetition(string):
    # test whether there's repetition in the string
    # check whether there even can be repetition, i.e., it's divisible by 2 or 3
    sg = ["P-6m2 or P-3m","P3m"]
    if len(string)%3 == 0:
        for i in range(len(string)):
            third = int(len(string)/3)
            
            split = [string[0:third-1],string[third:2*third-1],string[2*third:3*third-1]]
            if split[0] == split[1] and split[0] == split[2]:
                sg = ["R-3m","R3m"]
                break
    if len(string)%2 == 0:
        for i in range(len(string)):
            half = int(len(string)/2)
            
            split = [string[0:half-1],string[half:2*half-1]]
            if split[0] == split[1]:
                sg = ["P6_3/mmc","P6_3mc"]
                break
    # CRUDE HACK! confirm that this holds for these states and check for further states
    if string in ["c","ccc"]:
        sg = ["Pm-3m","Pm-3m"]
    return sg

def check_inversion(string)  :      
    triple_string = string*3
    
    lstr = len(string)
    
    k = 0
    invs = []
    for i in range(lstr,lstr*2):
        k += 1
        fw  = triple_string[i+1:i+lstr]
        obw = triple_string[i-lstr+2:i+1]
        bw  = triple_string[i-lstr+1:i]
        bw  = bw[::-1]
        obw = obw[::-1]
        
        # print(triple_string,fw, obw, bw, k)
        # print(len(fw),len(obw),len(bw))
        if fw == bw:
            # print(k,"fwbw")
            if 2*k < lstr*2:
                invs.append(2*k)
            else:
                invs.append(2*k-lstr*2)
        if fw == obw:
            # print(k,"int")
            if 2*k+1 < lstr*2:
                invs.append(2*k+1)
            else:
                invs.append(2*k+1-lstr*2)
    return invs


def get_vectors(parameters):
    # generate vectors from lattice parameters

    vec_c = [0,0,parameters[2]]

    vec_a = [np.sin(parameters[4]*np.pi/180)*parameters[0],0,np.cos(parameters[4]*np.pi/180)*parameters[0]]

    alpha_spread = [np.sin(parameters[3]*np.pi/180)*parameters[1],0,np.cos(parameters[3]*np.pi/180)*parameters[1]]
    vec_b = [np.cos(parameters[5]*np.pi/180)*alpha_spread[0],np.sin(parameters[5]*np.pi/180)*alpha_spread[0],alpha_spread[2]]

    vecs = [vec_a,vec_b,vec_c]
    return vecs

    
def write_cif(lattice_params,seq_code,seq_layer,atom_list,elements):
    n_layers = len(seq_layer)
    
    # get volume
    latt_vec = get_vectors(lattice_params)
    V = np.linalg.norm(np.multiply(latt_vec[0],np.cross(latt_vec[1], latt_vec[2])))

    # generate sum formula
    sum_formula = [elements["A"]+"1",elements["B"]+"1",elements["X"]+"3"]
    
    name = ""
    for item in sum_formula:
        name += item.replace("1","")
    name += "_" + seq_code
    # cif_name = filename #name + ".cif"
    
    def save_structure():
        Files = [('CIF File', '*.cif'),
            ('All Files', '*.*')]
        file = fd.asksaveasfile(initialfile=name,filetypes = Files, defaultextension = Files)
        return file
    
    with save_structure() as cif:
        cif.write("data_"+name+"\n")
        cif.write("\n")
        cif.write("# This structure is a automatically generated perovskite with the layer sequence {}.\n".format(seq_layer))
        cif.write("\n")
        cif.write("# Chemical Data\n")
        cif.write("{:<24} '{} {} {}'\n".format("_chemical_formula_sum",sum_formula[0],sum_formula[1],sum_formula[2]))
        cif.write("\n")
        cif.write("# Unit Cell\n")
        format = "{:<24} {} \n"
        
        cif.write(format.format("_cell_formula_units_Z",round(n_layers/2)))
        cif.write(format.format("_cell_length_a",lattice_params[0]))
        cif.write(format.format("_cell_length_b",lattice_params[1]))
        cif.write(format.format("_cell_length_c",lattice_params[2]))
        cif.write(format.format("_cell_angle_alpha",lattice_params[3]))
        cif.write(format.format("_cell_angle_beta",lattice_params[4]))
        cif.write(format.format("_cell_angle_gamma",lattice_params[5]))
        cif.write(format.format("_cell_volume",V))
        cif.write("\n")
        
        cif.write("# Symmetry Information\n")
        cif.write(format.format("_symmetry_cell_setting","triclinic"))
        cif.write(format.format("_symmetry_space_group_name_H-M","'P 1'"))
        cif.write("loop_\n")
        cif.write(" _symmetry_equiv_pos_as_xyz\n")
        cif.write("x,y,z\n")
        cif.write("\n")
        
        cif.write("# Structural Data\n")
        cif.write("loop_\n")
        cif.write(" _atom_site_label\n")
        cif.write(" _atom_site_type_symbol\n")
        cif.write(" _atom_site_fract_x\n")
        cif.write(" _atom_site_fract_y\n")
        cif.write(" _atom_site_fract_z\n")
        cif.write(" _atom_site_occupancy\n")
        
        atom_counter = {"A":0,"B":0,"X":0}
        for key, element in elements.items():
            for atom in atom_list:
                if atom[0] == element:
                    atom_counter[key] += 1
                    cif.write("  {}{} {} ".format(atom[0],atom_counter[key],atom[0]))
                    for coord in atom[1]:
                        cif.write("{:6.4f} ".format(coord))
                    cif.write("1\n")
                    
        messagebox.showinfo("CIF file saved","CIF file was successfully saved.")
                    
                    



# GUI

# creates the search and request window
class main_window:
    # initializes the base window
    def __init__(self):
        self.invs = []
        self.sg = ["",""]
        self.layer_sequence = ""
        self.jn = ""
        
        self.ABX_orig = ["A","B","X"]
        
        self.root = create_window("450x300+120+120", "PerovGen")
        self.frame_entry_fields()
        self.frame_output()
        self.frame_buttons()
        self.root.mainloop()        
        
     
    # frame holding the buttons for file management and plotting
    def frame_entry_fields(self):
        self._frame_encoded_sequence = tk.Frame(self.root)
        self._frame_encoded_sequence.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        self._frame_ABX = tk.Frame(self.root)
        self._frame_ABX.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        
        
        def entry_fields():
            # label, entry for Jagodszinski notation
            self._label_encoded_sequence = ttk.Label(self._frame_encoded_sequence)
            self._label_encoded_sequence["text"] = "Jagodszinski notation:"
            self._label_encoded_sequence.pack(side=tk.LEFT)
            
            self.encoded_sequence = tk.StringVar()
            self._entry_encoded_sequence = ttk.Entry(
                self._frame_encoded_sequence,
                textvariable=self.encoded_sequence
            )
            self._entry_encoded_sequence.pack(side=tk.LEFT)
            
            # label, entry for elements
            self._label_ABX = [""] * 3
            self._entry_ABX = [""] * 3
            self.ABX = [""] * 3
            for i in range(3):
                self._label_ABX[i] = ttk.Label(self._frame_ABX)
                self._label_ABX[i]["text"] = "Element {}:".format(self.ABX_orig[i])
                self._label_ABX[i].pack(side=tk.LEFT)
                
                self.ABX[i] = tk.StringVar()
                self._entry_ABX[i] = ttk.Entry(
                    self._frame_ABX,
                    textvariable=self.ABX[i],
                    width=4
                )
                
                self._entry_ABX[i].insert(tk.END, self.ABX_orig[i])
                self._entry_ABX[i].pack(side=tk.LEFT)
            
        entry_fields()
        
    def frame_output(self):
        
        message = ""
        if len(self.sg[0]) > 0:
            message += "Full Jagodzinski notation:\n  {}\n\n".format(self.jn)
            message += "Space group (HM notation):\n"
            if len(self.invs) > 0:
                message += "  {}\n\n".format(self.sg[0])
            else:
                message += "  {}\n\n".format(self.sg[1])
        
        # tidies up the layer sequence with Greek letters
        tidy_layer_sequence = ""
        for letter in self.layer_sequence:
            if letter == "a":
                tidy_layer_sequence += "α"
            elif letter == "b":
                tidy_layer_sequence += "β"
            elif letter == "c":
                tidy_layer_sequence += "γ"
            else:
                tidy_layer_sequence += letter
                
        if len(self.layer_sequence) > 0:            
            message += "Sequence of layers and interlayers\n(inversion centers marked with i):\n"
            message += "  {}\n".format(tidy_layer_sequence)
            message += "  "
            for i in range(len(self.layer_sequence)):
                if i in self.invs:
                    message += "i"
                else:
                    message += " "
            message += "\n"
        
        # message += "{} {}".format(self.sg,self.invs)
        
        
        
        self.text_box = tk.Text(self.root, wrap = "word",height=10)
        self.text_box.pack(side=tk.TOP,expand=True,fill=tk.X)
        self.text_box.insert('1.0', message)
        self.text_box.config(state='disabled')
        
    # frame holding the About button
    def frame_buttons(self):
        self._frame_buttons = tk.Frame(self.root)
        self._frame_buttons.pack(side=tk.BOTTOM,fill=tk.X)
        
        sep = ttk.Separator(self._frame_buttons,orient='horizontal')
        sep.pack(side=tk.TOP,fill=tk.X)
        
        # button for structure info
        self._button_info_structure = ttk.Button(self._frame_buttons, text = 'Info Structure', command = lambda : convert_sequence(cif=False))
        self._button_info_structure.pack(side=tk.LEFT,expand=True)
        
        # button for saving the structure
        self._button_save_structure = ttk.Button(self._frame_buttons, text = 'Save Structure', command = lambda : convert_sequence(cif=True))
        self._button_save_structure.pack(side=tk.LEFT,expand=True)
        
        def convert_sequence(cif):
            # initialize positions and elements
            _positions = positions()
            
            
            base_encoded_sequence = self.encoded_sequence.get()
            
            if len(base_encoded_sequence) == 0:
                write = False
                messagebox.showerror("Missing Jagodszinski sequence","No Jagodszinski sequence given! Sequence should only contain c (cubic) and h (hexagonal).")
            
            # shift by one, since notation starts at layer A, but program applies starting at B
            base_encoded_sequence = base_encoded_sequence[1:]+base_encoded_sequence[0]
            
            # check if the encoded sequence is valid
            for letter in base_encoded_sequence:
                if not letter in ["c","h"]:
                    write = False
                    messagebox.showerror("Invalid Jagodszinski sequence","Invalid Jagodszinski sequence! Sequence should only contain c (cubic) and h (hexagonal).")
                    break
                else:
                    write = True
            
            elements = {"A":"A","B":"B","X":"X"}
            
            # if len(elements_input) == 3:
            elements["A"] = self.ABX[0].get()
            elements["B"] = self.ABX[1].get()
            elements["X"] = self.ABX[2].get()
            # elif len(elements_input) != 0:
            #     messagebox.showerror("Invalid ABX","Element input invalid. Defaulting to A, B, and X.")
            #     write = False
            
            
            
            # to be defined by the user
            # encoded_sequence = "hh" #"chcchccchchcchccchchcchccch"
            # elements = {"A":"Ti","B":"Ca","X":"O"}
            
            if write == True:
                                
                def get_layer_sequence(encoded_sequence):
                    # define the initial sequence
                    layer_sequence = "AcB"
                    for char in encoded_sequence:
                    
                        layer_sequence += add_layer(layer_sequence,_positions,char=char)
        
                    layer_sequence += add_layer(layer_sequence,_positions,final=True)
                                        
                    return layer_sequence
                
                for i in range(1,3):
                    encoded_sequence = base_encoded_sequence * i
                    layer_sequence = get_layer_sequence(encoded_sequence)
                    # print(layer_sequence)
                    if not "x" in layer_sequence:
                        break
                
                layer_sequence = layer_sequence[:-4]
                
                self.jn = encoded_sequence[-1] + encoded_sequence[:-1]
                self.layer_sequence = layer_sequence
                
                # if i == 2:
                #     messagebox.showinfo("Invalid Jagodszinski sequence","Invalid initial Jagodszinski sequence! Sequence ends on double A layer and was doubled to {}.".format(encoded_sequence))
                # if i == 3:
                #     messagebox.showinfo("Invalid Jagodszinski sequence","Invalid initial Jagodszinski sequence! Sequence ends on double A layer and was tripled to {}.".format(encoded_sequence))
            
            
            if write == True:
                # check space group
                self.sg = check_repetition(encoded_sequence)
                self.invs = check_inversion(encoded_sequence)
                self.text_box.destroy()
                self.frame_output()
                
                # create a list of atoms with element labels and fractional coordinates
                atom_list = []
                
                i = 0
                n_layers = len(layer_sequence)
                
                for char in layer_sequence:
                    for entry in _positions.layers[char]:
                        atom = [elements[entry[0]],(entry[1],entry[2],i/n_layers)]
                        atom_list.append(atom)
                        
                    i += 1
                
                
                # define the lattice parameters of the cell
                
                a_cubic = 4.0 # empirically set for now
                
                a = a_cubic*np.sqrt(2)
                c = n_layers/2 * a_cubic/np.sqrt(3)
                
                lattice_params = [a,a,c,90,90,120]
                
                if cif == True:
                    write_cif(lattice_params, self.jn, layer_sequence, atom_list, elements)
                
        # function for adding the layers, one by one
        def add_layer(layers, positions, char="c", final=False):
            
            seq_error = False
            
            if final == False:
                if char == "h":
                    new_layer = layers[-3]
                elif char == "c":
                    options = list(positions.layers.keys())
                    for layer in [layers[-3],layers[-1]]:
                        options.remove(layer)
                    new_layer = options[0]
            else:
                new_layer = layers[0]
                
            
                
            if sorted([layers[-1],new_layer]) == sorted(["A","B"]):
                new_interlayer = "c"
            elif sorted([layers[-1],new_layer]) == sorted(["A","C"]):
                new_interlayer = "b"
            elif sorted([layers[-1],new_layer]) == sorted(["B","C"]):
                new_interlayer = "a"
            else:
               # messagebox.showerror("Invalid Jagodszinski sequence","Invalid Jagodszinski sequence! Sequence ends on double A layer.")
                seq_error = True
            
            if seq_error == False:
                if final == False:
                    new_i_l = new_interlayer + new_layer
                else:
                    new_i_l = new_interlayer
            else:
                new_i_l = "x"
            
            return new_i_l
        
        # button for showing the help window
        self._button_help = ttk.Button(self._frame_buttons, text = 'Help', command = lambda : helpbox())
        self._button_help.pack(side=tk.LEFT,expand=True)
        
        def helpbox():
            helpbox = create_window("650x400+120+120", "Help")
            helpbox.config(bg='#AAAAAA')
            
            message ='''This program turns a perovskite stacking sequence and list of elements A, B, and X into a CIF file.
Version 0.8.0

Jagodszinski sequence:
    hhchcchchhhc

Elements:
    abcabc
'''
            text_box = tk.Text(helpbox, wrap = "word")
            text_box.pack(expand=True,fill=tk.X)
            text_box.insert('end', message)
            text_box.config(state='disabled')
        
        # button that displays a window with the program version, license, and brief description
        about_button = ttk.Button(
            self._frame_buttons,
            text='About',
            command = lambda: about()
            )
        about_button.pack(side=tk.LEFT,expand=True)
        
        def about():
            about = create_window("650x400+120+120", "About PerovGen")
            about.config(bg='#AAAAAA')
            message ='''PerovGen turns a perovskite stacking sequence and a list of elements A, B, and X into a CIF file.
Version 0.8.0

MIT License
Copyright (c) 2024 mhaefner-chem
Contact: michael.haefner@uni-bayreuth.de

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
            text_box = tk.Text(about, wrap = "word")
            text_box.pack(expand=True,fill=tk.X)
            text_box.insert('end', message)
            text_box.config(state='disabled')

# function that ensures that the created windows do not become bigger than the screen
def window_size_limiter(avail_wxh,req_wxh,req_offset_xy):

    actual_wxh = [0,0]
    actual_offsets = [0,0]
    
    # check whether window fits on the current screen with and without offsets
    for i in range(len(avail_wxh)):
        if req_wxh[i] > avail_wxh[i]:
            actual_wxh[i] = avail_wxh[i]
            print("Caution, requested window doesn't fit the screen!")
        elif req_wxh[i] + req_offset_xy[i] > avail_wxh[i]:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = avail_wxh[i] - req_wxh[i]
            print("Caution, requested offset would move window off the screen!")
        else:
            actual_wxh[i] = req_wxh[i]
            actual_offsets[i] = req_offset_xy[i]
    
    return actual_wxh,actual_offsets

# function that creates a new window
def create_window(dimensions="500x350+100+100", title = "Tkinter Hello World", icon = ""):
   
    w = int(dimensions.split("x")[0])
    h = dimensions.split("x")[1]
    h = int(h.split("+")[0])
    
    offset_x = int(dimensions.split("+")[1])
    offset_y = int(dimensions.split("+")[2])
    
    # initializes the Tk root window
    window = tk.Tk()
    
    # gets screen properties and centers in upper third
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()
    
    offset_x = int(screen_width/3 - w / 3)
    offset_y = int(screen_height/3 - h / 3)
    
    # makes sure the window stays within bounds
    actual_wxh, actual_offsets = window_size_limiter([screen_width,screen_height],[w,h], [offset_x,offset_y])
    
    # set a title
    window.title(title)
    
    # specify geometry and max and min measurements
    window.geometry(f"{actual_wxh[0]}x{actual_wxh[1]}+{actual_offsets[0]}+{actual_offsets[1]}")
    window.minsize(10,10)
    window.maxsize(screen_width,screen_height)
    if icon != "":
        window.iconbitmap(icon)
    
    return window

# beginning of main program
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter import messagebox
if __name__ == "__main__":
    
    # make window
    main = main_window()
    