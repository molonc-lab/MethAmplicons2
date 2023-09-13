from setuptools import setup, find_packages
from setuptools.command.install import install
import os

# Directory containing setup.py
setup_dir = os.path.abspath(os.path.dirname(__file__))

# Path to requirements.txt
requirements_path = os.path.join(setup_dir, "requirements.txt")

REQUIREMENTS = [i.strip() for i in open(requirements_path).readlines()]

setup(
    name='methamplicons',
    version='0.1',
    description='CLI tool for plotting targeted bisulfite sequencing',
    author='Brett Liddell, Olga Kondrashova', 
    author_email='brett.liddell@qimrberghofer.edu.au, olga.kondrashova@qimrberghofer.edu.au',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=REQUIREMENTS,
    entry_points={
        'console_scripts': 'methamplicons = methamplicons.main:main'
    },
)
