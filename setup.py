from setuptools import setup


with open('README.md') as f:
    long_description = f.read()

setup(
    name='BioModTool',
    version='1.0.0',    
    description='Package to generate biomass objective functions to update genome-scale metabolic models.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="",
    authors='Clemence Dupont Thibert, Sylvaine Roy, Gilles Curien, Maxime Durot',
    corresponding_author_email='maxime.durot@totalenergies.com',
    license='LGPL',
    packages=['BioModTool'],
    install_requires=['cobra', 
                      'numpy',
                      'pandas',
                      'openpyxl',          
                     ],
    )   
