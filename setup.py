"""setup.py -- setup script

Setuptools config
"""

from setuptools import setup


setup(
    name='feps',
    version='0.0.0',
    packages=['feps'],
    install_requires=[
        'pytest-cov',
        'bigfloat',
        'pytest',
        'coverage',
        'click',
        'requests'
    ],
    entry_points={
       'console_scripts': [
           'feps = feps.main:cli',
           'seps = feps.epitopeenergies:epitope_energies'
       ]
    }
)
