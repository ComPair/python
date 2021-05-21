python scripts to manipulate the output of MEGAlib

For a description of the tools ad how to make performance plots, see https://confluence.slac.stanford.edu/display/COMPair/AMEGO+Python+Tools

New:

Use `Prep_cosima_ring.py` and `Prep_revan_ring.py` to simulate FarFieldAreaSources, with the directions of the incoming photons randomly distributed on a ring centered around the boresight, with the ring radius given by the theta values and the ring width set to one. Edit the scripts to adjust theta values and energies if necessary. 

usage: 

```
Prep_cosima_ring.py [-h] [-g G] [-o O]

optional arguments:
  -h, --help  show this help message and exit
  -g G        Geometry file
  -o O        Output file dir (Default: current directory)
```

```
Prep_revan_ring.py [-h] [-g G] [-f F] [-o O] [-b B]

optional arguments:
  -h, --help  show this help message and exit
  -g G        Geometry file
  -f F        Source file dir (Default: current directory)
  -o O        Output file dir (Default: current directory); will create
              symbolic links to input files if output directory is different
              from input
  -b B        Base revan config file, default: revan_AMEGO_X.cfg
```


Use `doAnalysis.py` for all-in-one analysis (using EventAnalysis and FigureOfMeritPlotter).

```
usage: doAnalysis.py [-h] [-d DIR] [-t TAG]

optional arguments:
  -h, --help         show this help message and exit
  -d DIR, --dir DIR  Directory with .sim and .tra files. (Default: current
                     directory)
  -t TAG, --tag TAG  Tag for output txt files. (Default: xxx)
```
