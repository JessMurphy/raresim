from setuptools import setup, find_packages

setup(
    name='raresim',
    version='0.2',
    description='Raresim library',
    author='Ryan Barnard',
    author_email='rbarnard1107@gmail.com',
    packages=find_packages(),
    install_requires=[
        'numpy',
    ],
    entry_points={
        'console_scripts': ['raresim=raresim.cli.cli:main'],
    },
)
