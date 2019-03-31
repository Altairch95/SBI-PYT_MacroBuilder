#!/usr/bin/env python
# coding=utf-8
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from MB.MacroB import *
import sys
import os
import threading


class StdRedirector():
    """Class that redirects the stdout and stderr to the GUI console"""
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        """Updates the console widget with the stdout and stderr output"""
        self.text_space.config(state=NORMAL)
        self.text_space.insert("end", string)
        self.text_space.see("end")
        self.text_space.config(state=DISABLED)

class MB(Frame):
    def show_about(self):
        messagebox.showinfo(title="About", message="This program has been done in the PYT subject by Natàlia Segura, "
                                                   "Altaïr CHinchilla and Pau Badia. For more information visit:\n"
                                                   "https://github.com/NataliaSegura/Macrocomplex-builder")
    def quit(self):
        """EXists the application"""
        if messagebox.askyesno("Quit", "Are you sure you want to exit?"):
            Frame.quit(self)

    def get_dir(self):
        """Gets the input directory path"""
        directory.set(filedialog.askdirectory())

    def get_template(self):
        """Gets the template path"""
        template_path.set(filedialog.askopenfilename(title="Select a template structure file", filetypes=[
            ("PDB files", "*.pdb"), ("mmCIF files", "*.cif")]))

    def create_menu(self):
        """Creates the menu of the app"""
        self.menubar = Menu(self)

        # CREATE THE FILEMENU
        filemenu = Menu(self.menubar)
        filemenu.add_command(label="Select Directory", command=self.get_dir)
        filemenu.add_separator()
        filemenu.add_command(label="QUIT", command=self.quit)

        # CREATE THE HELP MENU
        helpmenu = Menu(self.menubar)
        helpmenu.add_command(label="About", command=self.show_about)

        self.menubar.add_cascade(label="File", menu=filemenu)
        self.menubar.add_cascade(label="Help", menu=helpmenu)
        self.master.config(menu=self.menubar)

    def create_options(self):
        """Creates the label and frame for the option's widget"""
        self.options = LabelFrame(self, text="Options")
        self.create_options_frame()
        self.options.grid(row=0, column=0, stick="we", columnspan=2)

    def clear_template(self):
        """Clears the template's path"""
        template_path.set("")

    def create_options_frame(self):
        """Creates the contents of the options widget"""
        frame = Frame(self.options)
        label_output = Label(frame, text="Output name")
        label_num_chains = Label(frame, text="Max number of chains")
        label_num_models = Label(frame, text="Number of models")
        label_current_dir = Label(frame, text="Selected directory:")
        label_current_dir_path = Label(frame, textvariable=directory)
        self.label_dirty = Checkbutton(frame, text="dirty mode (generates tmp files)", onvalue=True, offvalue=False, variable=dirty)
        self.label_verbose = Checkbutton(frame, text="Verbose (Show all steps)", onvalue=True, offvalue=False, variable=verbose)
        label_template = Label(frame, text="Stechiometry by template (Optional)")
        entry_template = Button(frame, text="Browse", command=self.get_template)
        clear_template = Button(frame, text="Clear", command=self.clear_template)
        clear_template.grid(row=4, column=1, sticky="e")
        label_template_path = Label(frame, textvariable=template_path)
        entry_output = Entry(frame, textvariable=output_name)
        self.entry_max_chains = Spinbox(frame, from_=100, to=1000)
        self.entry_num_models = Spinbox(frame, from_=1, to=100)
        self.entry_run = Button(frame, text="RUN", command=self.thread_MB)
        label_stech = Label(frame, text="Stechiometry by string (Optional)")
        entry_label_stech = Entry(frame, textvariable=stech_string)
        label_stech.grid(row=5, column=0, sticky="w")
        entry_label_stech.grid(row=5, column=1, columnspan=2, sticky="w")
        label_current_dir.grid(row=0, column=0, sticky="w")
        label_current_dir_path.grid(row=0, column=1, sticky="w", columnspan=3)
        label_output.grid(row=1, column=0, sticky="w")
        entry_output.grid(row=1, column=1, sticky="w")
        label_num_chains.grid(row=2, column=0, sticky="w")
        self.entry_max_chains.grid(row=2, column=1, sticky="w")
        label_num_models.grid(row=3, column=0, sticky="w")
        self.entry_num_models.grid(row=3, column=1, sticky="w")
        self.label_dirty.grid(row=1, column=2, sticky="w", columnspan=3)
        self.label_verbose.grid(row=2, column=2, sticky="w", columnspan=3)
        label_template.grid(row=4, column=0, sticky="w")
        entry_template.grid(row=4, column=1, sticky="w")
        label_template_path.grid(row=4, column=2)
        self.entry_run.grid(row=2, column=6, sticky="w")
        frame.grid(row=0)

    def create_sequence_dict(self):
        """Creates the label and frame of sequence dictionary widget"""
        self.seq_dict_frame = LabelFrame(self, text="Sequences")
        self.create_sequence_dict_frame()
        self.seq_dict_frame.grid(row=1, column=0, stick="wn")

    def create_sequence_dict_frame(self):
        """Creates contents of the seqquence dictionary widget """
        seq_frame = Frame(self.seq_dict_frame)
        scrollbar_v = Scrollbar(seq_frame, orient=VERTICAL)
        scrollbar_h = Scrollbar(seq_frame, orient=HORIZONTAL)
        self.seq_listbox = Text(seq_frame, yscrollcommand=scrollbar_v.set, xscrollcommand=scrollbar_h.set, width=90,
                                height=6, state=DISABLED, wrap="none", background="black", foreground="white")
        scrollbar_v.config(command=self.seq_listbox.yview)
        scrollbar_h.config(command=self.seq_listbox.xview)
        scrollbar_v.pack(side=RIGHT, fill=Y)
        scrollbar_h.pack(side=BOTTOM, fill=X)
        self.seq_listbox.pack(side=LEFT, expand=True, fill=BOTH)
        seq_frame.grid(row=0)

    def create_estequiometry(self):
        """Creates the label and the frame of stechometry widget"""
        self.estequiometry_frame = LabelFrame(self, text="Model composition")
        self.create_estequiometry_frame()
        self.estequiometry_frame.grid(row=2, column=0, stick="wn")

    def create_estequiometry_frame(self):
        """Creates the content of the stechometry widget"""
        frame = Frame(self.estequiometry_frame)
        scrollbar = Scrollbar(frame, orient=HORIZONTAL)
        self.estequiometry_listbox = Text(frame, xscrollcommand=scrollbar.set, width=90, height=2,
                                          state=DISABLED, background="black", foreground="white")
        # scrollbar.config(command=self.seq_listbox.yview)
        scrollbar.pack(side=BOTTOM, fill=X)
        self.estequiometry_listbox.pack(side=LEFT, expand=True, fill=BOTH)
        frame.grid(row=0)

    def create_console(self):
        """Generates the label and console widget"""
        self.console_frame = LabelFrame(self, text="Console", padx=5, pady=5)
        self.create_console_frame()
        self.console_frame.grid(row=1, column=1, rowspan=3, stick="ns")

    def create_console_frame(self):
        """Generates the console widget"""
        frame = Frame(self.console_frame)
        scrollbar = Scrollbar(frame, orient=VERTICAL)
        self.console = Text(frame, yscrollcommand=scrollbar.set, width=60, height=30, state=DISABLED, wrap='word',
                            background="black", foreground="white")
        scrollbar.grid(row=0,column=1, stick="ns")
        scrollbar.config(command=self.console.yview)
        self.console.grid(row=0,column=0, stick="ns")
        sys.stdout = StdRedirector(self.console)
        sys.stderr = StdRedirector(self.console)
        frame.grid(row=0)

    def update_image(self):
        """Updates the image widget with the first macrocomplex model"""
        try:
            cwd = os.getcwd()
            best_model_name = cwd +"/"+ output_name.get()+"_1.cif"
            image_name = "%s/%s.png" % (cwd,output_name.get()+"_1")
            os.system("pymol %s -c -d 'hide all;show ribbon;util.cbc;png %s, width=300, height=300'" % (best_model_name, image_name))
            image = PhotoImage(file=image_name)
            self.model_image.create_image(150, 150, anchor=CENTER, image=image)
            self.model_image.image = image
        except:
            sys.stderr.write("Pymol couldn't create the image. Please, check if pymol is installed")

    def create_image(self):
        """Generates the image label"""
        self.image_frame = LabelFrame(self, text="Structure")
        self.create_image_frame()
        self.image_frame.grid(row=3, column=0, stick="n")

    def create_image_frame(self):
        """Generates the image frame"""
        frame = Frame(self.image_frame)
        self.model_image = Canvas(frame, width=300, height=300, background="black")
        self.model_image.pack()
        frame.grid(row=0, column=1)


    def update_estequiometry(self, best_model):
        """Updates the stecheometry widget"""
        stq_dict = generate_model_profile(best_model)
        self.estequiometry_listbox.config(state=NORMAL)
        self.estequiometry_listbox.delete(1.0, END)
        for key in sorted(stq_dict.keys()):
            self.estequiometry_listbox.insert("end", "%s:%s " % (key, stq_dict[key]))
        self.estequiometry_listbox.config(state=DISABLED)


    def run_MB(self):
        """Using MacroB functions, builds the macrocomplex"""
        self.in_pdbmodels = read_pdbs(directory.get()+"/", verbose.get())
        if self.in_pdbmodels:
            self.seq_dict = unify_ids(self.in_pdbmodels)
            self.interaction_dict = get_interaction_dict(self.in_pdbmodels, verbose.get())
            self.update_seq_dict()
            update_interactions_dict(self.interaction_dict, verbose)
            self.stech_dict = None
            if template_path.get() != "":
                self.stech_dict = get_template_stech_dict(template_path.get(), self.seq_dict, verbose=verbose)
            elif stech_string.get():
                self.stech_dict = get_string_stech_dict(stech_string.get())
            out_pdbmodels = main_loop(int(self.entry_num_models.get()), output_name.get(), self.interaction_dict, verbose.get(),
                                      int(self.entry_max_chains.get()), dirty.get(), self.stech_dict)
            self.update_estequiometry(out_pdbmodels[0])
            save_results(out_pdbmodels, output_name.get())
            self.update_image()
        else:
            sys.stderr.write("No valid input PDB given in %s/ directory. Please, inputs should be "
                             "pdbs that contain only 2 chains interacting" % directory.get())
        self.entry_run.config(state='normal')

    def thread_MB(self):
        """Starts a new thread to run the modeling apart from the GUI"""
        if os.path.isdir(directory.get() + "/"):
            self.console.config(state=NORMAL)
            self.console.delete(1.0, END)
            self.console.config(state=DISABLED)
            self.estequiometry_listbox.config(state=NORMAL)
            self.estequiometry_listbox.delete(1.0, END)
            self.estequiometry_listbox.config(state=DISABLED)
            self.seq_listbox.config(state=NORMAL)
            self.seq_listbox.delete(1.0, END)
            self.seq_listbox.config(state=DISABLED)
            self.model_image.delete("all")
            self.entry_run.config(state='disabled')
            if template_path.get() and stech_string.get():
                sys.stderr.write("For stechometry, select only one of the options, not both.")
            else:
                t = threading.Thread(target=self.run_MB)
                t.daemon = True # close pipe if GUI process exits
                t.start()
        else:
            sys.stderr.write("Directory %s/ doesn't exists, please select a valid directory" % directory.get())

    def createWidgets(self):
        """Creates all widgets"""
        self.create_menu()
        self.create_console()
        self.create_options()
        self.create_sequence_dict()
        self.create_estequiometry()
        self.create_image()
        self.grid(row=0)

    def update_seq_dict(self):
        """Updates the sequence dictionary widget"""
        self.seq_listbox.config(state=NORMAL)
        self.seq_listbox.delete(1.0, END)
        for key in self.seq_dict:
            self.seq_listbox.insert("end", "%s: %s\n" % (self.seq_dict[key], key))
        self.seq_listbox.config(state=DISABLED)



    def __init__(self, master=None, **kwargs):
        """Initalizates the app"""
        Frame.__init__(self, master, **kwargs)
        self.master.wm_title("Macrocomplex Builder")
        self.master.resizable(width=True, height=True)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        self.createWidgets()


if not os.path.exists('tmp'):
    os.mkdir('tmp')
if not os.path.exists('models'):
    os.mkdir('models')
root = Tk()
root.geometry("1300x800")
directory = StringVar()
template_path = StringVar()
template_path.set("")
output_name = StringVar()
max_number_chains = IntVar()
verbose = BooleanVar()
stech_string = StringVar()
dirty = BooleanVar()
directory.set('None')
app = MB(master=root, padx=10, pady=10)
root.mainloop()
