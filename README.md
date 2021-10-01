# epitope_spotter
  
Toying around with epitope selection from a given protein based on MHC binding, posttranslational modifications and sequence homology

## Usage

Development is in process!

Before building the dockr image, please obtain a standalone version of NetMHC from cbs.dtu.dk and copy the archive to the `third_party` folder.

To run rstudio inside the container, issue 
```
docker run --rm -p 127.0.0.1:8787:8787 -v $PWD:/home/rstudio/local_files -e USERID=$UID -e PASSWORD=[SOME PASWORD] szabogtamas/epitope_spotter
```