#from distutils.core import setup
from setuptools import setup
import glob
import subprocess

args = 'git describe --tags'
p = subprocess.Popen(args.split(), stdout=subprocess.PIPE)
version = p.communicate()[0].decode("utf-8").strip()

#version = '0.1.0'

version_str = """# git describe --tags
__version__ = "{0}"\n""".format(version)

fp = open('grizli_aws/version.py','w')
fp.write(version_str)
fp.close()
print('Git version: {0}'.format(version))

setup(name='grizli_aws',
      version=version,
      description='Grizli processing on AWS',
      author='Gabriel Brammer',
      author_email='brammer@stsci.edu',
      url='https://github.com/grizli-project/grizli-aws',
      packages=['grizli_aws'],
      scripts=glob.glob('scripts/*'),
      install_requires=['boto3']
     )
