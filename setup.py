from distutils.core import setup

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
               'scripts/fit_redshift_lambda.py',
               'scripts/sync_extractions_FROM_s3',
               'scripts/grizli_run_all_synced.sh',
               'scripts/grizli_check_bad_full.sh'],
      requires=['boto3']
     )
