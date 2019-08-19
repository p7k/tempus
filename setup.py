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
    install_requires=[
        'docopt-ng==0.7.2',
        'flatten-dict==0.1.0',
        'hgvs==1.3.0.post0',
        'more-itertools==7.2.0',
        'pyvcf==0.6.8',
    ],
    tests_require=['pytest'],
    python_requires='>=3.7.0'
)
