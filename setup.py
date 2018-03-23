from distutils.core import setup

setup(
    name='phievo',
    version='1.1',
    packages=['phievo','phievo.Networks','phievo.AnalysisTools','phievo.Populations_Types','phievo.Networks','phievo.Networks.PlotGraph','phievo.ConfigurationTools'],
    #package_dir={'phievo': 'phievo'},
    package_data={'phievo': ['CCodes/*',"ConfigurationTools/default_cfiles/*"]},
    license='GNU LESSER GENERAL PUBLIC LICENSE',
    description='In silico gene network evolution algorithm.',
    author='P. Francois, A. Henry, M. Hemery',
    author_email='paul.francois2@mcgill.ca',
    url='https://phievo.github.io/',
    install_requires=[
        "matplotlib>=1.5.3",
        "networkx==1.11",
        "numpy>=1.11.3",
        "scipy>=0.18.1",
        "Sphinx>=1.5.3",
    ],
    long_description=open('README.md').read(),
)
