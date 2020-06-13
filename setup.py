from setuptools import setup, find_packages
import glob

setup(
    name='mmic',
    version='0.1',
    url='https://github.com/langfzac/mmic',
    packages=find_packages(),
    package_data={'mmic': ['utils/*']},
    license='LICENSE',
    install_requires=["numpy","pytipsy","astropy"],
    include_package_data=True
)
