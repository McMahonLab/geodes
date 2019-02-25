---
title: "GEODES"
author: "Alex Linz"
output: html_document
---

__Diel transcriptomics of freshwater lakes__

GEODES stands for "Gene Expression in Oligotrophic, Dystrophic, and Eutrophic Systems."

This dataset is the subject of the manuscript "Time-series metatranscriptomes reveal conserved patterns and ecological interactions between phototrophic and heterotrophic microbes in diverse freshwater systems." Alexandra M. Linz, Frank O. Aylward, Stefan Bertilsson, and Katherine McMahon. Submitted 2019.


Want more background about GEODES? Check out these blog posts! https://uwmcmahonlab.wordpress.com/2016/05/04/introducing-geodes/ & https://uwmcmahonlab.wordpress.com/2016/11/14/geodes-data-preview/


__Where do I get the data?__

Are you in the McMahon Lab or one of our collaborators?

- Go to /lakes_data/Metatranscriptomes/geodes/ on our storage drive.
- If you can't access the storage drive, see this page: https://github.com/McMahonLab/McMahon_Fileshare

Not in the lab but curious about our data?

- Download our raw sequencing data from http://genome.jgi.doe.gov/Diecycetabolisms/Diecycetabolisms.info.html . Most metatranscriptomes became public on Aug 1, 2017, but some samples are earlier and others later.
- Go to our GEODES folder on the OSF https://osf.io/9gr62/) - our large processed files are here.
- Shoot us an email! trina.mcmahon@wisc.edu or amlinz@wisc.edu

__BE WARNED__ We're not making the raw data or intermediate files available on Github because it would far exceed our space limit. The quality filtered, compressed metatranscriptome files total 185GB! Download at your own risk.

__Repo Structure__

- protocols
	- Field and lab protocols
- environmental_data
	- for_database
		- Files formatted for upload to the McMahon Lab SQL database
	- for_humans
		- Excel files that are easily interpreted by people
- bioinformatics_workflow
  - executable_files
    - Bash scripts that actually process the data
  - lab_notebooks
    - My notes from writing the CHTC workflow. May contain sass.
  - programs
    - Pre-built tarballs of the programs used for use on CHTC
  - R_processing
    - Table formatting once the data was small enough to load into R
  - submit_files
    - Tells CHTC how to run your jobs
  - submit_node_scripts
    - Extremely minor tasks to run on the submit node itself that need to run between parallel tasks - like making a list of files to queue or moving files around
    - Full_GEODES_Workflow.Rmd/.html
      - 1st attempt at workflow. Now out of date but has some useful info on why I chose to do things the way I did. Will be updated on 2nd interation of workflow
    - GEODES_bioinformatics_instructions.Rmd
      - How to run the current workflow on CHTC
- sample_data
    - Metadata about our samples, such as which were sequenced and how many metatranscriptomic reads were in each
- Manuscript
    - Files associated with our publication, as well as previous drafts
- README.md
    - This file
