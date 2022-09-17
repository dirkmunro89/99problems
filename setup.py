from setuptools import setup, find_packages

_py3_packages = find_packages(include=['./'])

# additional requirements besides `install_requires` for development
#_dev = ['yapf', 'flake8', 'tox', 'pytest', 'pytest-cov']

# additional requirements to generate package documentation
#_docs = ['sphinx', 'sphinx_rtd_theme', 'sphinxcontrib-napoleon']

# optional requirements to link with other libraries
#_opt = ['cvxopt']

# combined dependencies
#_all = list(set(_dev) | set(_docs) | set(_opt))

setup(
    name="Python 3",
    version="0.0.1",
    description='General environment for Python 3',
    author_email='dirkmunro8@gmail.com',
    packages=_py3_packages,
    keywords=['Python 3'],
    python_requires='>=3.8, <4',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'cvxopt',
        'pytest',
        'osqp'
    ]
#    extras_require={
#        'all': _all,
#        'dev': _dev,
#        'doc': _docs,
#    },
)
