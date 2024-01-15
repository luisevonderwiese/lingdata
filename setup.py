from setuptools import find_packages
from setuptools import setup

setup(
    name='lingdata',
    version='0.0.1',
    install_requires=['requests', 'pandas', 'numpy', 'ete3', 'biopython', 'PyGithub', 'getpass4', 'termcolor', 'argparse'],
    packages=find_packages('.'),
    package_dir={'': '.'},
    url='https://github.com/luisevonderwiese/lingdata',
    license='GNU',
    author='Luise Häuser',
    author_email='luise.haeuser@h-its.org',
    description='Database interface for language phylogenies',
    entry_points={'console_scripts': ['lingdata = lingdata.lingdata:main']}
)
