# import libraries
import streamlit as st
import tkinter as tk
from tkinter import filedialog

def folder_picker():
    # Set up tkinter
    root = tk.Tk()
    root.withdraw()

    # Make folder picker dialog appear on top of other windows
    root.wm_attributes('-topmost', 1)

    # Folder picker button
    st.title('Folder Picker')
    st.write('Please select a folder:')
    clicked = st.button('Folder Picker')
    if clicked:
        dirname = st.text_input('Selected folder:', filedialog.askdirectory(master=root))
