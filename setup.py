#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', 'gfapy>=1.0.0']

setup_requirements = [ ]

test_requirements = ['gfapy', 'IPython' ]

setup(
    author="Josiah Seaman, Toshiyuki Yokoyama, Simon Heumos",
    author_email='graph-genome-browser@googlegroups.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Python Boilerplate contains all the boilerplate you need to create a Python package.",

    # we have to change that at some point!
    ####
    entry_points={
        'console_scripts': [
            'test_python_template=test_python_template.cli:main',
        ],
    },
    ####
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='vgbrowser',
    name='vgbrowser',
    packages=find_packages(include=['vgbrowser']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/graph-genome/vgbrowser',
    version='0.1.0',
    zip_safe=False,
)