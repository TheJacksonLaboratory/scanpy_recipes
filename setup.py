from setuptools import setup, find_packages
from pathlib import Path
#import versioneer

package_name = 'scanpy_recipes'

req_path = Path('requirements.txt')
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

with open('README.md', encoding='utf-8') as fin:
    readme = fin.read()

author = 'Bill Flynn'

setup(
    name=package_name,
    version="0.1",#versioneer.get_version(),
    #cmdclass=versioneer.get_cmdclass(),
    description='Generalized and aggregated analysis functions for working with Anndata and Scanpy',
    long_description=readme,
    url='http://github.com/TheJacksonLaboratory/scanpy_recipes',
    author=author,
    author_email='bill.flynn@jax.org',
    license='',
    install_requires=requires,
    packages=find_packages(),
    # `package_data` does NOT work for source distributions!!!
    # you also need MANIFTEST.in
    # https://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    package_data={"scanpy_recipes": ["data/*.csv"]},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
