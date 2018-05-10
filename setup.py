from setuptools import setup
from setuptools.extension import Extension

setup(name='grizli_aws',
      version='0.1',
      description='Grizli processing on AWS',
      author='Gabriel Brammer',
      author_email='brammer@stsci.edu',
      url='https://github.com/grizli-project/grizli-aws',
      packages=['grizli_aws'],
      scripts=['scripts/grizli_extract_and_fit.sh',
               'scripts/grizli_extract_and_fit.py',
               'scripts/sync_extractions_TO_s3',
               'scripts/sync_extractions_FROM_s3'],
      install_requires=['boto3>=1.7.16']
     )
