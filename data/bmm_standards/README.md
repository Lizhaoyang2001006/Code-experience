# bmm-standards

NIST Public Data Resource:
A collection of X-ray Absorption Spectroscopy data of metals, stable oxides, and
other common compounds measured at NIST's Beamline for Materials Measurement.
Version 1.0.0
DOI: https://doi.org/10.18434/mds2-4032


Contact:
  Bruce Ravel
    bruce.ravel@nist.gov
	bravel@bnl.gov

Description:

X-ray Absorption Spectroscopy (XAS) data and associated metadata.  Data includes
XAS and XRF data serialized as text using the XAS Data Interchange specification
(https://github.com/XraySpectroscopy/XAS-Data-Interchange). These data were
measured at NIST's Beamline for Materials Measurement
(https://www.bnl.gov/nsls2/beamlines/beamline.php?r=6-bm) at the National
Synchrotron Light Source II (https://www.bnl.gov/nsls2).

The aim of this collection is to provide a chemically and structurally diverse
set of example XAS spectra measured on elements across the periodic table and
accessible at BMM. As such, these data are collected using a single photon
delivery system and a consistent set of detectors and detector readout systems.

The physical samples from which these data were collected are maintained in a
library at BMM and made available for use by visitors to the beamline.


--------------
Data Use Notes
--------------

This data is publicly available according to the NIST statements of
copyright, fair use and licensing; see
https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications

You may cite the use of this data as follows:

> Bruce Ravel (2025), A collection of X-ray Absorption Spectroscopy data of
> metals, stable oxides, and other common compounds measured at NIST's Beamline
> for Materials Measurement., Version 1.0.0, National Institute of Standards and
> Technology, https://doi.org/10.18434/mds2-4032 (Accessed: [give download date])

-------------
Data Overview
-------------

The "Data/" folder contains subfolders for each element represented in the collection.
Inside those folders are column ASCII data files using the XDI formatting specification
for serializing metadata.  The filnames all follow the format of "Element-Edge-Material.xdi"

The "tmpl/" folder contains HTML templates for formatting the web presentation of these data.

The "tiled/" folder contains configuration and example code for converting to a Tiled
representation of this data collection.  Tiled is the data management service used at NBSLS-II:
https://blueskyproject.io/tiled/

The remaining files are used to generate the web page representation of the collection.

Many of the physical samples in this collection were provided by Dr. Martin Stennett of
the University of Sheffield.  Many other users of BMM have made contributions to this library
over the years of BMM's operation.

---------------------------
Building the periodic table
---------------------------

The simplest way to build the web page is to use a python virtual
environment.

IN the main folder, do

    python3 -m venv ./my-venv
	
Then install the dependencies locally by doing:

    ./my-venv/bin/pip install mendeleev
    ./my-venv/bin/pip install xraylarch
	
The build by doing:

    ./my-venv/bin/python ./json2html.py 



---------------
Version History
---------------

1.0.0 (this version)
  initial release

