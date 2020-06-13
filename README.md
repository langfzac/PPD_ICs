# mmic

Initial conditions for the moving-mesh version of ChaNGa, specifically for protoplanetary disks.

## Getting Started

After installing, you can find a detailed tutorial notebook in the examples directory.

### Dependencies
Assuming a recent version of python is installed (i.e. 3.7+), we will need a few packages to get started. Follow the links below to find installation instructions.
- [Astropy](https://www.astropy.org)
- [Numpy](https://numpy.org)

We'll also need [PyTipsy](https://github.com/bwkeller/pytipsy), to write out our ICs as tipsy files. The main repository is currently only compatible with python 2.7. For now, we'll use my fork which is on python 3.7+.
    
    git clone git@github.com:langfzac/pytipsy
    cd pytipsy/
    pip install .
    
#### Optional (for visualizations)
Currently, the [yt](https://yt-project.org) package is best equiped for MANGA visualizations. The most recent, stable version will do. 

### Install mmic
Now, simply run the following to install mmic!

    git clone git@github.com:langfzac/mmic
    cd mmic/
    pip install .
