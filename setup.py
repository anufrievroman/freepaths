from pathlib import Path
import re
import setuptools

setup_dir = Path(__file__).resolve().parent

version = re.search( r'__version__ = "(.*)"', Path(setup_dir, 'freepaths/__main__.py').open().read())
if version is None:
    raise SystemExit("Could not determine version to use")
version = version.group(1)

setuptools.setup(
    name='freepaths',
    author='Roman Anufriev',
    author_email='anufriev.roman@protonmail.com',
    url='https://github.com/anufrievroman/freepaths',
    description='Phonon Monte Carlo simulator',
    long_description=Path(setup_dir, 'README.md').open().read(),
    long_description_content_type='text/markdown',
    license='GPL',
    entry_points={
        "console_scripts": [
            "freepaths = freepaths.__main__:run"
        ]
    },
    install_requires=['numpy', 'matplotlib', 'scipy', 'imageio', 'colorama'],
    version=version,
    python_requires='~=3.11',
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Environment :: Console",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Utilities",
    ],
    packages=["freepaths"],
    include_package_data=True,
)
