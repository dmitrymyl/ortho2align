from setuptools import setup

setup(name='ortho2align',
      version='0.9',
      description='Sequence alignment tool based on syntenic protein neighbourhood.',
      url='http://github.com/dmitrymyl/ortho2align',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='GPL-3.0',
      packages=['ortho2align'],
      install_requires=['numpy',
                        'pebble',
                        'scipy',
                        'sortedcontainers',
                        'tqdm'],
      python_requires='~=3.6',
      entry_points={'console_scripts':
                    ['ortho2align=ortho2align.cli_scripts:ortho2align']},
      zip_safe=False)
