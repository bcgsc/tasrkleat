# Usage

There is only one script that matters for usage,
[`search_hexamer.py`](https://github.com/bcgsc/tasrkleat/blob/master/hexamer_search/search_hexamer.py),
and it is designed to be used as a Python module.

It requires [pysam](https://github.com/pysam-developers/pysam) and
[Biopython](http://biopython.org/). If haven't had a virtual environment with
them installed already, please use [Miniconda](https://conda.io/miniconda.html)
or
[pip](https://pypi.python.org/pypi/pip)+[virtualenv](https://pypi.python.org/pypi/virtualenv)
to create one and install the required packages. Then, copy over the script, and

```
import search_hexamer

search_hexamer.search(*args, **kwargs)
```

See the
[source code](https://github.com/bcgsc/tasrkleat/blob/master/hexamer_search/search_hexamer.py)
for the API.


# What does it do?

Given an cleavage site, search for the existence of one of the following 16
polyadenylation signal (PAS) hexamers (ordered by strength, 1 is the strongest)
in a window (default to 50 bp) upstream of it. If multiple PAS hexamers exist,
the strongest one will be selected.

1. AATAAA
2. ATTAAA
3. AGTAAA
4. TATAAA
5. CATAAA
6. GATAAA
7. AATATA
8. AATACA
9. AATAGA
10. AAAAAG
11. ACTAAA
12. AAGAAA
13. AATGAA
14. TTTAAA
15. AAAACA
16. GGGGCT

The coordinate system used in this patch is 0-based because
[`pysam`](https://github.com/pysam-developers/pysam) uses 0-based coordinates.

Besides, FYI, personal experience tells that UCSC, IGV, GTF are 1-based. For the
coordinate system used by more file formats, please see
https://www.biostars.org/p/84686/.


# Development

To run tests

```
python test_search_hexamer.py
```
