GABI - GESTALT Analysis using Bayesian Inference in BEAST2:
-------

GABI (GESTALT Analysis using Bayesian Inference) is a [BEAST 2](http://www.beast2.org/) package to estimate time-scaled lineage trees and cell population dynamics from alignments of barcodes that accumulate Cas9-mediated indels, such as GESTALT,  as described [here](https://doi.org/10.1126/science.aaf7907).
GABI implements a phylogenetic likelihood adapted from [GAPML](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9387344/), adapted for time-scaled lineage tree inference by incoporating a molecular clock rate.
The preprint describing GABI's implementation and example applications can be found [here](https://doi.org/10.1098/rstb.2023.0318).

Installation
-------
GABI requires an installation of BEAST 2.7 which can be obtained from https://www.beast2.org/. GABI can then be installed using BEAUti with the following instructions:

- Open BEAUti.
- In the `File` menu, select `Manage Packages`.
- Click the `Package repositories` button at the bottom of the dialog box.
- Click on`Add URL` and enter the following URL:
  https://raw.githubusercontent.com/azwaans/GABI/refs/heads/master/package.xml. Close this window.
- `GABI` should now appear as a package in the list of available packages. Select it and click the `Install/Upgrade` button.
- GABI should now be available for use in BEAST 2 and Beauti. Close and restart Beauti to start generating an xml.

Running GABI on your data
-------

Start from a csv file, where the rows are cells, the columns are the GSM formatted barcode sequences. Then use the `convert_cell_id_barcode_to_beast_input.R` script under `scripts/` to generate a `[your_filename].gestalt` file that can be read by BEAUTi.

Open BEAUTi, click on File -> Template -> GABI to open the GABI template. Then, click on the "plus" button in the bottom left corner and read your `[your_filename].gestalt` into BEAUTi. Now, you can specify your preferred Site Model, Clock Model etc. For general information on setting up a Bayesian phylogenetic analysis, refer to the [BEAST tutorials](https://taming-the-beast.org/).

You can find an example xml file in `examples/gestalt_model.xml` file.

To run the analysis from the command line, run the following command:
`beast -seed 1 gestalt_model.xml`

You can also run the analysis using the GUI. To do that, open BEAST. When asked for the input file, choose your xml and hit the "Run" button.

License
-------

GABI is free software.  It is distributed under the terms of version 3
of the GNU General Public License.  A copy of this license should
be found in the file COPYING located in the root directory of this repository.
If this file is absent for some reason, it can also be retrieved from
https://www.gnu.org/licenses.