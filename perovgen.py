#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 13:16:42 2024

@author: mhaefner-chem

Stacking of Perovskites, working

"""


def symmops_data(n_sg):
    
     symmops = {"156":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z")],
"160":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z"),
("x+2/3,y+1/3,z+1/3"),
("-y+2/3,x-y+1/3,z+1/3"),
("-x+y+2/3,-x+1/3,z+1/3"),
("-y+2/3,-x+1/3,z+1/3"),
("-x+y+2/3,y+1/3,z+1/3"),
("x+2/3,x-y+1/3,z+1/3"),
("x+1/3,y+2/3,z+2/3"),
("-y+1/3,x-y+2/3,z+2/3"),
("-x+y+1/3,-x+2/3,z+2/3"),
("-y+1/3,-x+2/3,z+2/3"),
("-x+y+1/3,y+2/3,z+2/3"),
("x+1/3,x-y+2/3,z+2/3")],
"164":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("y,x,-z"),
("x-y,-y,-z"),
("-x,-x+y,-z"),
("-x,-y,-z"),
("y,-x+y,-z"),
("x-y,x,-z"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z")],
"166":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("y,x,-z"),
("x-y,-y,-z"),
("-x,-x+y,-z"),
("-x,-y,-z"),
("y,-x+y,-z"),
("x-y,x,-z"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z"),
("x+2/3,y+1/3,z+1/3"),
("-y+2/3,x-y+1/3,z+1/3"),
("-x+y+2/3,-x+1/3,z+1/3"),
("y+2/3,x+1/3,-z+1/3"),
("x-y+2/3,-y+1/3,-z+1/3"),
("-x+2/3,-x+y+1/3,-z+1/3"),
("-x+2/3,-y+1/3,-z+1/3"),
("y+2/3,-x+y+1/3,-z+1/3"),
("x-y+2/3,x+1/3,-z+1/3"),
("-y+2/3,-x+1/3,z+1/3"),
("-x+y+2/3,y+1/3,z+1/3"),
("x+2/3,x-y+1/3,z+1/3"),
("x+1/3,y+2/3,z+2/3"),
("-y+1/3,x-y+2/3,z+2/3"),
("-x+y+1/3,-x+2/3,z+2/3"),
("y+1/3,x+2/3,-z+2/3"),
("x-y+1/3,-y+2/3,-z+2/3"),
("-x+1/3,-x+y+2/3,-z+2/3"),
("-x+1/3,-y+2/3,-z+2/3"),
("y+1/3,-x+y+2/3,-z+2/3"),
("x-y+1/3,x+2/3,-z+2/3"),
("-y+1/3,-x+2/3,z+2/3"),
("-x+y+1/3,y+2/3,z+2/3"),
("x+1/3,x-y+2/3,z+2/3")],
"186":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("-x,-y,z+1/2"),
("y,-x+y,z+1/2"),
("x-y,x,z+1/2"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z"),
("y,x,z+1/2"),
("x-y,-y,z+1/2"),
("-x,-x+y,z+1/2")],
"187":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("x,y,-z"),
("-y,x-y,-z"),
("-x+y,-x,-z"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z"),
("-y,-x,-z"),
("-x+y,y,-z"),
("x,x-y,-z")],
"194":[
("x,y,z"),
("-y,x-y,z"),
("-x+y,-x,z"),
("-x,-y,z+1/2"),
("y,-x+y,z+1/2"),
("x-y,x,z+1/2"),
("y,x,-z"),
("x-y,-y,-z"),
("-x,-x+y,-z"),
("-y,-x,-z+1/2"),
("-x+y,y,-z+1/2"),
("x,x-y,-z+1/2"),
("-x,-y,-z"),
("y,-x+y,-z"),
("x-y,x,-z"),
("x,y,-z+1/2"),
("-y,x-y,-z+1/2"),
("-x+y,-x,-z+1/2"),
("-y,-x,z"),
("-x+y,y,z"),
("x,x-y,z"),
("y,x,z+1/2"),
("x-y,-y,z+1/2"),
("-x,-x+y,z+1/2")],
"221":[
("x,y,z"),
("-x,-y,z"),
("-x,y,-z"),
("x,-y,-z"),
("z,x,y"),
("z,-x,-y"),
("-z,-x,y"),
("-z,x,-y"),
("y,z,x"),
("-y,z,-x"),
("y,-z,-x"),
("-y,-z,x"),
("y,x,-z"),
("-y,-x,-z"),
("y,-x,z"),
("-y,x,z"),
("x,z,-y"),
("-x,z,y"),
("-x,-z,-y"),
("x,-z,y"),
("z,y,-x"),
("z,-y,x"),
("-z,y,x"),
("-z,-y,-x"),
("-x,-y,-z"),
("x,y,-z"),
("x,-y,z"),
("-x,y,z"),
("-z,-x,-y"),
("-z,x,y"),
("z,x,-y"),
("z,-x,y"),
("-y,-z,-x"),
("y,-z,x"),
("-y,z,x"),
("y,z,-x"),
("-y,-x,z"),
("y,x,z"),
("-y,x,-z"),
("y,-x,-z"),
("-x,-z,y"),
("x,-z,-y"),
("x,z,y"),
("-x,z,-y"),
("-z,-y,x"),
("-z,y,-x"),
("z,-y,-x"),
("z,y,x")]}
     
     return symmops[n_sg]







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
    sg = ["P-3m1 (164)","P3m1 (156)","P-6m2 (187)"]
    triplestring = string*3
    if len(string)%3 == 0:
        for i in range(len(string)):
            third = int(len(string)/3)
            
            split = [triplestring[i:third+i],triplestring[third+i:2*third+i],triplestring[2*third+i:3*third+i]]
            if split[0] == split[1] and split[0] == split[2]:
                sg = ["R-3m (166)","R3m (160)","ERR"]
                break
    if len(string)%2 == 0:
        for i in range(len(string)):
            half = int(len(string)/2)
            
            split = [triplestring[i:half+i],triplestring[half+i:2*half+i]]
            if split[0] == split[1]:
                sg = ["P6_3/mmc (194)","P6_3mc (186)","ERR"]
                break
    # special case: triple cubic
    if not "h" in string and len(string)%3 == 0:
        sg = ["Pm-3m (221)","Pm-3m (221)","ERR"]
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
        

        if fw == bw:
            if 2*k < lstr*2:
                invs.append(2*k)
            else:
                invs.append(2*k-lstr*2)
        if fw == obw:
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
    
    
def write_cif(lattice_params,seq_code,seq_layer,atom_list,elements,sg):
    n_layers = len(seq_layer)
    
    # get volume
    latt_vec = get_vectors(lattice_params)
    V = np.linalg.norm(np.multiply(latt_vec[0],np.cross(latt_vec[1], latt_vec[2])))

    # generate sum formula
    sum_formula = [elements["A"]+"1",elements["B"]+"1",elements["X"]+"3"]
    
    # determine setting
    # print(sg.split(")")[0].split[1])
    n_sg = sg[0].split(")")[0]
    hm_sg = n_sg.split("(")[0]
    hm_sg = hm_sg.strip()
    n_sg = n_sg.split("(")[1]
    if int(n_sg) < 168:
        setting = "trigonal"
    elif int(n_sg) < 200:
        setting = "hexagonal"
    else:
        setting = "cubic"
    
    symmops = symmops_data(n_sg)
        
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
    
    file = save_structure()
    if file != None:
        with file as cif:
            cif.write("data_"+name+"\n")
            cif.write("\n")
            cif.write("# This structure is an automatically generated perovskite with the layer sequence {}.\n".format(seq_layer))
            cif.write("\n")
            cif.write("# Chemical Data\n")
            cif.write("{:<24} '{} {} {}'\n".format("_chemical_formula_sum",sum_formula[0],sum_formula[1],sum_formula[2]))
            cif.write("\n")
            cif.write("# Unit Cell\n")
            format = "{:<24} {} \n"
            
            cif.write(format.format("_cell_formula_units_Z",round(n_layers/2)))
            cif.write(format.format("_cell_length_a","{:10.6f}".format(lattice_params[0])))
            cif.write(format.format("_cell_length_b","{:10.6f}".format(lattice_params[1])))
            cif.write(format.format("_cell_length_c","{:10.6f}".format(lattice_params[2])))
            cif.write(format.format("_cell_angle_alpha","{:10.6f}".format(lattice_params[3])))
            cif.write(format.format("_cell_angle_beta","{:10.6f}".format(lattice_params[4])))
            cif.write(format.format("_cell_angle_gamma","{:10.6f}".format(lattice_params[5])))
            cif.write(format.format("_cell_volume","{:10.3f}".format(V)))
            cif.write("\n")
            
            cif.write("# Symmetry Information\n")
            cif.write(format.format("_symmetry_cell_setting",setting))
            cif.write(format.format("_symmetry_space_group_name_H-M","'{}'".format(hm_sg)))
            cif.write(format.format("_space_group_IT_number","{}".format(n_sg)))
            cif.write("loop_\n")
            cif.write(" _symmetry_equiv_pos_as_xyz\n")
            for symmop in symmops:
                cif.write("{}\n".format(symmop))
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
                    

def symmetry_analysis(cell):
    import spglib
    sg = spglib.get_spacegroup(cell)
    symmetry = spglib.get_symmetry(cell)
    stdcell = spglib.standardize_cell(cell)
    # print(symmetry)
    # print(spglib.get_symmetry_dataset(cell))
    # print(spglib.get_symmetry(stdcell))
    # print(stdcell)

    eq_atoms = list(set(symmetry["equivalent_atoms"]))
    return sg, stdcell, eq_atoms



# GUI

# creates the search and request window
import os
class main_window:
    
    # initializes the base window
    def __init__(self):
        
        self.version = "0.9.5"
        
        # attempt to import spglib for in-depth symmetry analysis
        try:
            import spglib
            self.has_spglib = True
        except:
            self.has_spglib = False
        
        self.invs = []
        self.sg = ["",""]
        self.layer_sequence = ""
        self.jn = ""
        
        self.ABX_label = ["A","B","X"]
        self.ABX_orig = ["Ca","Ti","O"]
        
        if self.has_spglib == True:
            window_name = "PerovGen"
        else:
            window_name = "PerovGen"
        
        if os.path.isfile("PerovGen.ico"):
            self.icon = "PerovGen.ico"
        else:
            self.icon = ""
        self.root = create_window("450x350+120+120", window_name,self.icon)
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
                self._label_ABX[i]["text"] = "Element {}:".format(self.ABX_label[i])
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
        
        # check for -6m2 symmety in layer to distinguish between P-6m2 and P-3m1
        has_6m2 = False
        for inv in self.invs:
            if inv%2 == 0 and self.layer_sequence[inv] in ["A","B","C"]:
                has_6m2 = True
            else:
                has_6m2 = False
                break
            
        
        message = ""
        if len(self.sg[0]) > 0:
            message += "Full Jagodzinski notation:\n  {}\n\n".format(self.jn)
            message += "Space group (HM notation):\n"
            if len(self.invs) > 0:
                if "P-3m1" in self.sg[0] and has_6m2 == True:
                    message += "  {}".format(self.sg[2])
                    self.sg[0] = self.sg[2]
                else:
                    message += "  {}".format(self.sg[0])
            else:
                message += "  {}".format(self.sg[1])
                self.sg[0] = self.sg[1]
                
            if self.has_spglib == True:
                message += ", SG from spglib: {}".format(self.sg[3])
            message += "\n\n"
            
        
        
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
            message += "Sequence of layers and interlayers\n(layers with symmetry centers marked with s):\n"
            message += "  {}\n".format(tidy_layer_sequence)
            message += "  "
            for i in range(len(self.layer_sequence)):
                if i in self.invs:
                    message += "s"
                else:
                    message += " "
            message += "\n"
        else:
            if self.has_spglib == True:
                message = "Program ready, spglib loaded for symmetry analysis."
            else:
                message = "Program ready, no spglib found. Using internal symmetry analysis."
        
        
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
                return
                
            # shift by one, since notation starts at layer A, but program applies starting at B
            base_encoded_sequence = base_encoded_sequence[1:]+base_encoded_sequence[0]
            
            # check if the encoded sequence is valid
            for letter in base_encoded_sequence:
                if not letter in ["c","h"]:
                    write = False
                    messagebox.showerror("Invalid Jagodszinski sequence","Invalid Jagodszinski sequence! Sequence should only contain c (cubic) and h (hexagonal).")
                    break
                    return
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
            
            
            
            # reduce encoded sequence to base unit if it contains repeated sequences
            if len(base_encoded_sequence) > 1:
                for i in range(int(len(base_encoded_sequence)),0,-1):
                    if len(base_encoded_sequence)%i == 0:
                        for k in range(i):
                            reduceable = True
                            for j in range(int(len(base_encoded_sequence)/i)):
                                start = j*i + k
                                end = (j+1)*i + k
                                if end > len(base_encoded_sequence):
                                    end = end - len(base_encoded_sequence)
                                # print("compare:",base_encoded_sequence[0:i],base_encoded_sequence[start:end])
                                if not base_encoded_sequence[0:i] == base_encoded_sequence[start:end]:
                                    reduceable = False
                            if reduceable == True:
                                base_encoded_sequence = base_encoded_sequence[0:i]
            # print("final:",base_encoded_sequence)
            
            if write == True:
                
                
                                
                def get_layer_sequence(encoded_sequence):
                    # define the initial sequence
                    layer_sequence = "AcB"
                    for char in encoded_sequence:
                    
                        layer_sequence += add_layer(layer_sequence,_positions,char=char)
        
                    layer_sequence += add_layer(layer_sequence,_positions,final=True)
                                        
                    return layer_sequence
                
                
                for i in range(1,4):
                    encoded_sequence = base_encoded_sequence * i
                    layer_sequence = get_layer_sequence(encoded_sequence)
                    # print(layer_sequence, layer_sequence[-4:-1])
                    if layer_sequence[-4:-1] == "AcB":
                        break
                    
                # check symmetry
                self.sg = check_repetition(encoded_sequence)
                self.invs = check_inversion(encoded_sequence)
                
                # adjust for fact that AcB is preset
                layer_sequence = layer_sequence[:-4]
                
                self.jn = encoded_sequence[-1] + encoded_sequence[:-1]
                self.layer_sequence = layer_sequence
            
            if write == True:
                
                
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
                latt_vec = get_vectors(lattice_params)
                
                pos = []
                numbers = []
                
                for atom in atom_list:
                    
                    pos.append(atom[1])
                    if atom[0] == elements["A"]:
                        numbers.append(1)
                    elif atom[0] == elements["B"]:
                        numbers.append(2)
                    elif atom[0] == elements["X"]:
                        numbers.append(3)
                
                cell = (latt_vec,pos,numbers)
                
                
                
                # check symmetry with spglib
                if self.has_spglib == True:
                    spglib_sym, spglib_cell, spglib_ineq = symmetry_analysis(cell)
                    self.sg.append(spglib_sym)
                    # print(spglib_cell)
                    # print(spglib_ineq)
                    atom_list = []
                    for i in range(len(spglib_cell[1])):
                        if i in spglib_ineq:
                            atom = [self.ABX[spglib_cell[2][i]-1].get(),spglib_cell[1][i]]
                            atom_list.append(atom)
                            
                    def vec_length(vec):
                        length = 0
                        for i in vec:
                            length += i**2
                        length = np.sqrt(length)
                        return length
                    
                    def vec_angle(vec_a,vec_b):
                        l_a = vec_length(vec_a)
                        l_b = vec_length(vec_b)
                        
                        dot = 0
                        for i in range(3):
                            dot += vec_a[i] * vec_b[i]
                            
                        angle = np.arccos(dot/(l_a * l_b))*180.0/np.pi
                        return angle
                
                    latt_vec = spglib_cell[0]
                    for i in range(3):
                        lattice_params[i] = vec_length(latt_vec[i])
                    
                    lattice_params[3] = abs(vec_angle(latt_vec[1],latt_vec[2]))
                    lattice_params[4] = abs(vec_angle(latt_vec[0],latt_vec[2]))
                    lattice_params[5] = abs(vec_angle(latt_vec[0],latt_vec[1]))
                    
                    # print(lattice_params)
                    
                    
                    
                self.text_box.destroy()
                self.frame_output()
                
                
                if cif == True:
                    write_cif(lattice_params, self.jn, layer_sequence, atom_list, elements, self.sg)
                
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
                new_layer = "B" #layers[0]
    
                
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
            helpbox = create_window("650x400+120+120", "Help",self.icon)
            helpbox.config(bg='#AAAAAA')
            
            message ='''PerovGen does a basic symmetry analysis on a given Jagodzinski sequence and a list of elements A, B, and X and turns them into a symmetrized CIF file.
Version {}

Jagodszinski notation:
Specifies the ordering of the layers of the perovskite.
"c" stands for cubic ordering, yielding the stacking pattern ABC.
"h" stands for hexagonal ordering, yielding the stacking pattern ABA.

Elements:
Specifies the elements put into sites A, B, and X of a general perovskite ABX3,with A being the in-layer cation, B the interstitial cation, and X the anion.
    
Notes:
A Jagodzinski sequence has to be provided, but no elements need to be given.
The program automatically reduces or multiplies the Jagodzinski sequence to yield a valid perovskite structure, so the sequence shown by the program may differ from the input sequence.
Furthermore, the final structure written into the CIF file will deviate from the sequence due to the standardization of the unit cell.
    
'''.format(self.version)
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
            about = create_window("650x400+120+120", "About PerovGen",self.icon)
            about.config(bg='#AAAAAA')
            message ='''PerovGen does a basic symmetry analysis on a given Jagodzinski sequence and a list of elements A, B, and X and turns them into a symmetrized CIF file.
Version {}

LICENSE:
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

CREDITS:
The icon for this software was created with VESTA (version 3.5.8)
K. Momma and F. Izumi, 
"VESTA 3 for three-dimensional visualization of crystal, volumetric and morphology data," 
J. Appl. Crystallogr., 44, 1272-1276 (2011).

The binaries were compiled with spglib for symmetry analysis.
Copyright (c) 2024, Spglib team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the <organization> nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT
HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''.format(self.version)
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
    