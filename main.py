# Tkinter for python 2, tkinter for python 3
import tkinter as tk


class MainWindow:
    def __init__(self, parent):
        self.parent = parent
        self.frame = tk.Frame(self.parent)
        self.button1 = tk.Button(self.frame, text ="Test")
        self.button1.pack()
        self.frame.pack()


# Entry point
def main():
    root = tk.Tk()
    app = MainWindow(root)
    root.mainloop()


if __name__ == '__main__':
    main()