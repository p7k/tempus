from setuptools import setup

setup(
    name='tempus',
    version='0.1',
    description='tempus biofx tech challenge',
    url='http://github.com/p7k/tempus',
    author='Pavel Katsev',
    author_email='pkatsev@gmail.com',
    packages=['tempus'],
    entry_points={'console_scripts': ['tempus = tempus.__main__:cli']},
    zip_safe=True,
    install_requires=['docopt-ng', 'hgvs', 'pyvcf'],
    test_requires=['pytest'],
    python_requires='>=3.7.0'
)
