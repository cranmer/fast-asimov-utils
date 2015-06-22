from setuptools import setup, find_packages

setup(
  name = 'fasimov',
  version = '0.0.1',
  description = 'fast Asimov statistical tools',
  url = '',
  author = 'Kyle Cranmer',
  author_email = 'kyle.cranmer@nyu.edu',
  packages = find_packages(),
  include_package_data = True,
  install_requires = [
    'scipy',
    'numpy',
  ],
  dependency_links = [ ]
)