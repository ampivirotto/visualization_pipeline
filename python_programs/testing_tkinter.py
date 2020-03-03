#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      tuk32868
#
# Created:     27/02/2020
# Copyright:   (c) tuk32868 2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from tkinter import *
from tkinter import filedialog

def browse_button():
    global folder_path
    filename = filedialog.askdirectory()
    folder_path.set(filename)
    print(filename)


def main():
    root = Tk()
    root.withdraw()
    folder_selected = filedialog.askdirectory()

if __name__ == '__main__':
    main()
