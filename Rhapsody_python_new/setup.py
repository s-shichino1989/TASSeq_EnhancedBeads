#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools import find_packages


def main():
    description = 'RhapsodyPython'

    setup(
        name='RhapsodyPython',
        version='0.0.1',
        author='SS',
        packages=find_packages(),
        description=description,
        long_description=description,
        zip_safe=False,
        include_package_data=True,
        install_requires=[],
        tests_require=[],
        setup_requires=[],
    )


if __name__ == '__main__':
    main()
