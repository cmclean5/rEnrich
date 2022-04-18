# rEnrich
Test enrichment of network clusters given node annotation data


## Compile

1) Install GSL:
1.1) Using mac can run: > brew install gsl

2) Set environment variables (replacing your GSL path in GSL_HOME):
2.1) in .bashrc:
      ## GSL environment variables setup
      export GSL_HOME="/opt/homebrew/Cellar/gsl/2.7.1"
      export GSL_CFLAGS="${GSL_HOME}/include"
      export GSL_LIBS="${GSL_HOME}/lib"
      export GSL_CONFIG="${GSL_HOME}/bin/gsl-config"
      export PATH="${GSL_HOME}/bin:${PATH}"
      export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_HOME}/lib"

3) Run autoconfig.ac:
3.1) > autoconf

4) Build and install:
4.1) > R CMD build rEnrich
4.2) > R CMD INSTALL rEnrich_1.0.tar.gz

5) Run example:
5.1) > cd rEnrich/example/
5.2) > R
5.3) > source('example.R')
