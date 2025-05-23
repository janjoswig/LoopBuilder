LoopBuilder: model missing loops in protein structures
======================================================

This package provides lightweight infrastructure to organise the
loop building process by leveraging arbitrary external tools. The basic
design principle follows a simple protocol:

  * Use a PDBx/mmCIF file containing a partial protein structure as input
  * Choose a `Builder` (e.g. `PDBFixerBuilder`)
  * Optionally, define the missing `Segment`s to be modelled or let the builder find them
  * Run `build` specifying a target number of `SegmentModels` to generate
    * The builder will try to generate models for each segment one after each other
    * Optionally, `Scorer`s and `Filter`s can be used to accept or reject models
    * The successful models will be saved to PDBx/mmCIF files, one for each segment
    * A `segments.csv` file describes the modelled segments
    * A `models.csv` file describes the generated models including their scores


Example
-------

For a working example see the Jupyter notebook under [sandbox/PDBFixer/explore.ipynb](https://github.com/janjoswig/LoopBuilder/blob/main/sandbox/explore.ipynb).

Installation
------------

Install the package via `pip` from the developement repository and optionally
satisfy dependencies by creating a `conda` environment:

```bash
$ git clone https://github.com/janjoswig/LoopBuilder.git
$ cd LoopBuilder
$ conda env create -f env.yml -n LoopBuilder
$ pip install .
```

Additional tools may need to
be installed, depending on what should be used by `LoopBuilder` under the hood.
Here is a list of possible sources:

  * MolProbity: [Docker image](https://hub.docker.com/r/francecosta/molprobity)

