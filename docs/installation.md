# Installing PyVIF

The installation is pretty simple but some dependencies need to be installed.
We propose to use Anaconda with Bioconda channel.

You neet to add `conda-forge` and `bioconda` channels:
```
conda config --add channels conda-forge
conda config --add channels bioconda
```

And the just download the repository by copying and pasting these commandes:
```
git clone https://gitlab.curie.fr/ddesvill/pyvif.git
cd pyvif
conda create --name pyvif --file requirements.txt
source activate pyvif
python setup.py install
```

Run the test to check if installation works well:
```
$ pytest
============================= test session starts ==============================
platform linux -- Python 3.6.3, pytest-3.6.1, py-1.5.3, pluggy-0.6.0
rootdir: /bioinfo/users/ddesvill/git/pyvif, inifile:
plugins: cov-2.5.1
collected 10 items                                                             

tests/test_bamtools.py .                                                 [ 10%]
tests/test_bp_finder.py ...                                              [ 40%]
tests/test_cli.py ..                                                     [ 60%]
tests/test_paftools.py ....                                              [100%]

========================== 10 passed in 16.15 seconds ==========================
```

Now, you can [run PyVIF](usage.md) !
