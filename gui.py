import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import gnt



class Sliders:
    def __init__(self, root, row, col, start, end, iter):
        return tk.Scale(
            root,
            from_=start,
            to=end,
            resolution=iter,
            orient="horizontal",
        ).grid(row=row, column=col, command=print("1"))


class MainWindow:
    def __init__(self, root, color):
        self.color = color
        self.root = root
        self.root.resizable(0, 0)
        self.root.geometry("400x600")
        root.columnconfigure(0, weight=1)
        self.root.configure(bg=self.color)
        '''
        LABELS
        '''
        population_label = tk.Label(
            self.root,
            text="Πληθυσμός: ",
            fg="#000000",
            font="Courier 14 bold",
            bg=self.color,
        ).grid(row=0, column=0, sticky=tk.W)
        
        generations_label = tk.Label(
            self.root, text="Γενιές:", fg="black", font="Courier 14 bold", bg=self.color
        ).grid(row=1, column=0, sticky=tk.W)
        
        pm_label = tk.Label(
            self.root,
            text="Πιθανότητα Μετάλλαξης: ",
            fg="black",
            font="Courier 14 bold",
            bg=self.color,
        ).grid(row=2, column=0, sticky=tk.W)
        
        pc_label = tk.Label(
            self.root,
            text="Πιθανότητα Διασταύρωσης: ",
            fg="black",
            font="Courier 14 bold",
            bg=self.color,
        ).grid(row=3, column=0, sticky=tk.W)
        
        cp_label = tk.Label(
            self.root,
            text="Σημεία διασταύρωσης: ",
            fg="black",
            font="Courier 14 bold",
            bg=self.color,
        ).grid(row=4, column=0, sticky=tk.W)
        
        bits_label = tk.Label(
            self.root, text="Bits: ", fg="black", font="Courier 14 bold", bg=self.color
        ).grid(row=5, column=0, sticky=tk.W)
        '''
        SLIDERS
        '''
        self.pop_slider = tk.Scale(
            root,
            from_=2,
            to=500,
            resolution=2,
            orient="horizontal",
        ).grid(
            row=0,
            column=1,
        )
        
        self.generation_slider = tk.Scale(
            root,
            from_=2,
            to=1000,
            resolution=2,
            orient="horizontal",
        ).grid(
            row=1,
            column=1,
        )

        self.pm_slider = tk.Scale(
            root,
            from_=0,
            to=1,
            resolution=0.01,
            orient="horizontal",
        ).grid(
            row=2,
            column=1,
        )

        self.pc_slider = tk.Scale(
            root,
            from_=0,
            to=1,
            resolution=0.01,
            orient="horizontal",
        ).grid(
            row=3,
            column=1,
        )

        self.cp_slider = tk.Scale(
            root,
            from_=1,
            to=38,
            resolution=1,
            orient="horizontal",
        ).grid(
            row=4,
            column=1,
        )

        self.bits_slider = tk.Scale(
            root,
            from_=2,
            to=40,
            resolution=1,
            orient="horizontal",
        ).grid(
            row=5,
            column=1,
        )

        '''BUTTONS'''
        self.run_button = tk.Button(
            self.root,
            text="Εκτέλεση",
            width=10,
            font="none 14"
        ).grid(row=7, column=0)

        exit_button = tk.Button(
            self.root,
            text="Έξοδος",
            width=10,
            font="none 14",
            command=self.root.destroy,
        ).grid(row=7, column=1)


        fig = plt.Figure(figsize=(4.5, 3), facecolor="#efebe9")

        plt.rcParams.update({"font.size": 8})
        ax = fig.add_subplot()
        ax.set_facecolor("#efebe9")
        ax.set_xlabel("generations")
        ax.set_ylabel("fitness")
        canvas = FigureCanvasTkAgg(
            fig,
            master=self.root,
        )
        
        canvas.get_tk_widget().grid(row=6, column=0, columnspan=2)
        
        self.root.mainloop()
    
    def run(self):
        ga = gnt.GeneticAlgorithm(self.pop_slider.get())
        
        


def main():
    root = tk.Tk()
    window = MainWindow(root, "#efebe9")


main()

