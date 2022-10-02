from setuptools import setup


setup(name='ortho2align',
      version='1.0.5',
      description='ortho2align is a lncRNA ortholog discovery CLI tool based on syntenic regions and statistical assessment of alignment nonrandomness.',
      url='http://github.com/dmitrymyl/ortho2align',
      author='Dmitry Mylarshchikov',
      author_email='dmitrymyl@gmail.com',
      license='GPL-3.0',
      packages=['ortho2align'],
      install_requires=['numpy',
                        'pebble',
                        'scipy',
                        'sortedcontainers~=2.1.0',
                        'tqdm'],
      python_requires='>=3.7',
      entry_points={'console_scripts':
                    ['ortho2align=ortho2align.cli_scripts:ortho2align']},
      zip_safe=False)
