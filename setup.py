from setuptools import setup

setup(name='ortho2align',
      version='0.7',
      description='Sequence alignment tool based on syntenic protein neighbourhood.',
      url='http://github.com/dmitrymyl/ortho2align',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='MIT',
      packages=['ortho2align'],
      install_requires=['pandas',
                        'tqdm',
                        'sortedcontainers',
                        'matplotlib',
                        'scipy',
                        'numpy'],
      python_requires='~=3.6',
      entry_points={'console_scripts':
                    ['ortho2align=ortho2align.cli_scripts:ortho2align']},
      zip_safe=False)
