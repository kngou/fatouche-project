from cx_Freeze import setup, Executable

 

exe = Executable(

    script="main.py",

    base="Win32GUI",

    )

 

setup(

    name = "programFatiha",

    version = "0.1",

    description = "the first version of her program",

    executables = [exe]

    )