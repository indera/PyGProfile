import os
from setuptools import setup, find_packages

setup(name='pygprofile',
      version='0.0.1',
      packages=find_packages(),
      description='PyGProfile is a Python module to profile gene lists using Gene Ontology, Pathways and phenotypes annotation',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      package_data={"pygprofile":['data/*']},
      scripts=['scripts/pygp_go.py', 'scripts/pygp_path.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

PyGProfile is a Python module for GO enrichment and pathway analysis

 Contact
=============

If you have any questions or comments about PyGProfile, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
