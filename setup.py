from setuptools import setup, find_packages

setup(
    name='raresim',
    version='0.2.1',
    description='Raresim library',
    author='Ryan Barnard',
    author_email='rbarnard1107@gmail.com',
    packages=find_packages(exclude=['tests']),
    install_requires=[
        'numpy<=2.1.0',
        'numba',
        'pandas',
        'scipy',
    ],
    entry_points={
        'console_scripts': ['raresim=raresim.cli:main'],
    },
    #scripts=['bin/raresim']
)
