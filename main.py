# Tkinter for python 2, tkinter for python 3
import tkinter as tk
import tkinter.scrolledtext as tk_scrolled_text
from tkinter import N, S, E, W
import pwa
import tkinter.ttk as ttk

class GUI(tk.Tk):
    def __init__(self, prog=None):
        super().__init__()
        self.style = ttk.Style()

        self.program = prog
        # labels and configurations
        self.title("Pairwise Alignment Tool")
        label_font = ('Courier New', 10, '')

        # input boxes
        self.input = tk_scrolled_text.ScrolledText(self, width=100, height=10, undo=True, wrap='word')
        self.input.grid(row=1, column=1, columnspan=10, padx=(5,5), pady=(5,0), sticky=E)

        # execution output box
        self.output = tk_scrolled_text.ScrolledText(self, width=100, height=20, undo=True, font=label_font, wrap='word')
        self.output.grid(row=2, column=1, columnspan=10, padx=(5,5), pady=(5,5), sticky=E)

        self.typ = tk.StringVar(value="PROTEIN")
        type_options_p = ttk.Radiobutton(self, text='Align Peptides', variable=self.typ, value=pwa.AlignmentType.PROTEIN.name)
        type_options_p.grid(row=3, column=6, sticky=W)
        type_options_n = ttk.Radiobutton(self, text='Align Nucleotides', variable=self.typ, value=pwa.AlignmentType.NUCLEOTIDE.name)
        type_options_n.grid(row=4, column=6, sticky=W)

        # options
        self.met = tk.StringVar(value="LOCAL")
        method_options_l = ttk.Radiobutton(self, text='Local Alignment', variable=self.met, value=pwa.AlignmentMethod.LOCAL.name)
        method_options_l.grid(row=3, column=5, sticky=W)
        method_options_g = ttk.Radiobutton(self, text='Global Alignment', variable=self.met, value=pwa.AlignmentMethod.GLOBAL.name)
        method_options_g.grid(row=4, column=5, sticky=W)

        # execution button
        self.button_align = ttk.Button(self, text="Align")
        self.button_align.grid(row=5, column=3, columnspan=5, padx=(5,5), pady=(5,5), sticky="ew")

        self.status = ttk.Label(text="waiting")
        self.status.grid(row=6, column=0, columnspan=10, padx=(10,0), pady=(0,5), sticky=W)

        # labels
        input_label = ttk.Label(self)
        input_label.configure(text='Input:', font=label_font)
        input_label.grid(row=1, column=0, padx=(5,0), sticky=W)

        output_label = ttk.Label(self)
        output_label.configure(text='Output:', font=label_font)
        output_label.grid(row=2, column=0, padx=(5,0), sticky=W)

        # event bindings
        self.button_align.bind("<Button-1>", self.submission)
        self.bind("<Return>", self.submission)

    def submission(self, event=None):
        self.button_align.config(state="disabled")
        self.status.config(text="Working on alignment. It may take a while depending on inputs...\n")
        self.status.update()

        alignment_type = None
        alignment_method = None
        if self.typ.get() == "PROTEIN":
            alignment_type = pwa.AlignmentType.PROTEIN
        else:
            alignment_type = pwa.AlignmentType.NUCLEOTIDE
        if self.met.get() == "LOCAL":
            alignment_method = pwa.AlignmentMethod.LOCAL
        else:
            alignment_method = pwa.AlignmentMethod.GLOBAL

        output = None
        try:
            prog = self.program(alignment_type, alignment_method, self.input.get('1.0', 'end-1c'))
            prog.do_alignment()
            output = prog.get_aligned_seq_formatted()
        except:
            self.status.config(text="Error")
            self.status.update()

        if output:
            self.status.config(text="Done")
            for line in reversed(output):
                self.output.insert('1.0', '\n')
                self.output.insert('1.0', line)
        else:
            self.status.config(text="Inputs not recognized")

        self.button_align.config(state="normal")

if __name__ == '__main__':
    code = pwa.PWA

    gui = GUI(code)
    gui.mainloop()
