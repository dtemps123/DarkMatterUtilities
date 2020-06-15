from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='dmutils',
      version='1.0',
      description='A plugin for calculating dark matter physics',
	  long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/dtemps123/DarkMatterUtilities',
      author='Dylan Temples',
      author_email='dtemps123@gmail.com',
      license='MIT',
      packages=setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      python_requires='>=3.6',
      zip_safe=False)