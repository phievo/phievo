from distutils.core import setup

setup(
    name='phievo',
    version='1.0',
    packages=['phievo','phievo.Networks','phievo.AnalysisTools','phievo.Populations_Types','phievo.Networks',''],
    license='GNU LESSER GENERAL PUBLIC LICENSE',
    description='In silico gene network evolution algorithm.',
    author='P. Francois, A. Henry, M. Hemery',
    author_email='paulf@physics.mcgill.ca',
    url='https://bitbucket.org/onnetworkevolution/phievo',
    
    long_description=open('README.md').read(),
)
