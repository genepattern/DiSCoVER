from edjuaro/pdx-hts:3.0
MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C



USER root

# Adding things to /build/
#RUN mkdir /build #it already exists on edjuaro/pdx-hts:3.0

RUN mkdir /module
ADD companion_script.py /module/companion_script.py
ADD call_discover.py /module/call_discover.py
ADD null_module.sh /module/null_module.sh

# using the same user as the parent container as described in:
# https://hub.docker.com/r/jupyter/datascience-notebook/dockerfile
#USER $NB_UID

# build using this:
# docker build -t genepattern/discover:1.1 .
