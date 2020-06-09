from setuptools import setup, find_packages
import glob

setup(
    name='mmic',
    url='https://github.com/langfzac/mmic',
    packages=find_packages(),
    package_data={'mmic': ['utils/*']},
    install_requires=["numpy","pytipsy","astropy"],
    include_package_data=True
)