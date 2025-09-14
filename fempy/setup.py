from setuptools import setup, find_packages

setup(
    name="fempy",  # Replace with your library's name
    version="0.2.0",  # Version of your library
    packages=find_packages(),  # Automatically find all packages
    install_requires=[
        "dill==0.3.8",
        "gmsh==4.13.1",
        "matplotlib==3.8.4",
        "numpy==2.1.3",
        "pandas==2.2.3",
        "scipy==1.14.1",
        "setuptools==65.5.0",
        "sympy==1.13.1",
        "tqdm==4.66.5",
        "vtk==9.3.1"
    ],
    entry_points={
        "console_scripts": [
            "tovtk=fempy.solution.tovtk:main",  # Entry point for the tovtk script
            "mtxviewer=fempy.mtxviewer.app:main",  # Entry point for the mtxviewer script
            "generatebeaminp=fempy.generate.generate_beam_inp:main",  # Entry point for the generatebeaminp script
        ],
    },
)
