import setuptools


setuptools.setup(
    name='Gaussium',
    version='1.0.0',
    description='A basic quantum chemistry program written in Python3',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'}
)
