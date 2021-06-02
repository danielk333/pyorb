import setuptools
__version__ = '0.4.3'

with open('README.rst', 'r') as fh:
    long_description = fh.read()


with open('requirements', 'r') as fh:
    pip_req = fh.read().split('\n')
    pip_req = [x.strip() for x in pip_req if len(x.strip()) > 0]


setuptools.setup(
    name='pyorb',
    version=__version__,
    long_description=long_description,
    url='https://github.com/danielk333/pyorb',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU-GPLv3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.0',
    install_requires=pip_req,
    packages=setuptools.find_packages(),
    extras_require={
        "dev":  [
            "pytest",
            "sphinx",
            'sphinx-gallery',
            'coverage',
        ],
    },
    # metadata to display on PyPI
    author='Daniel Kastinen',
    author_email='daniel.kastinen@irf.se',
    description='Kepler orbit functions in Python',
    license='GNU-GPLv3',
)
