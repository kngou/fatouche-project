#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
from tkinter import *
from tkinter.messagebox import *



def firstPage (message):
	window = Tk()
	label = Label(window, text = message)
	label.pack()
	window.mainloop()

# def getInput():
#     showinfo("Alerte", entree.get())
# 	value = StringVar() 
# 	value.set("Valeur")
# 	entree = Entry(fenetre, textvariable=value, width=30)
# 	entree.pack()



def EndOfProgram():
	if askyesno('New amorce list', 'Would you like to design, your amorce from output exon list ?'):
		showinfo('Titre 2', 'Tant pis...')
	else:
		showinfo('Quit', 'End of the program, have a nice day :D')
	
	Button(text='Action', command=EndOfProgram).pack()
# EndOfProgram()