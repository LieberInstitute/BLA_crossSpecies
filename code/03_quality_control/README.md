# BLA cross-species project: 03_quality_control

This directory contains the scripts used to run empty droplets and quality control procedures for human, macaque, and baboon datasets. <empty_droplets_{species}> files should be ran first to detect emplet droplets. <remove_droplets_add_PerCellQC_{species}> should next be run to drop empty droplets as well as run standard quality control procedures with the <scuttle> package. 

The species order does not matter.