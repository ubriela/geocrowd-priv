__author__="ubriela"
__date__ ="$Sep 14, 2013 11:03:36 AM$"

from setuptools import setup,find_packages

setup (
  name = 'PrivateGeoCrowd',
  version = '0.1',
  packages = find_packages(),

  # Declare your packages' dependencies here, for eg:
  install_requires=['foo>=3'],

  # Fill in these to make your Egg ready for upload to
  # PyPI
  author = 'ubriela',
  author_email = '',

  summary = 'Just another Python package for the cheese shop',
  url = '',
  license = '',
  long_description= 'Long description of the package',

  # could also include long_description, download_url, classifiers, etc.

  
)