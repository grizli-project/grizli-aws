from distutils.core import setup
import glob

setup(name='grizli_aws',
      version='0.1',
      description='Grizli processing on AWS',
      author='Gabriel Brammer',
      author_email='brammer@stsci.edu',
      url='https://github.com/grizli-project/grizli-aws',
      packages=['grizli_aws'],
      scripts=glob.glob('scripts/*'),
      requires=['boto3']
     )
