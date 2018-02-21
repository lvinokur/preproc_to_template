FROM ubuntu:14.04
MAINTAINER Lea Vinokur <lea.vinokur@gmail.com>
ARG DEBIAN_FRONTEND=noninteractive
# Core system capabilities required
RUN apt-get update && apt-get install -y git python libeigen3-dev zlib1g-dev wget bsdtar software-properties-common curl tar perl-modules
RUN curl -sL https://deb.nodesource.com/setup_4.x | bash -
RUN apt-get install -y nodejs

# Now that we have software-properties-common, can use add-apt-repository to get to g++ version 5, which is required by JSON for Modern C++
RUN add-apt-repository ppa:ubuntu-toolchain-r/test
RUN apt-get update && apt-get install -y g++-5

# Neurodebian
RUN wget -O- http://neuro.debian.net/lists/trusty.au.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9 && apt-get update


# fsl
RUN apt-get install -y ants
RUN apt-get install -y fsl-5.0-core
RUN apt-get install -y fsl-first-data
RUN apt-get install -y fsl-mni152-templates


#MRtrix3 setup
ENV CXX=/usr/bin/g++-5
RUN git clone https://github.com/MRtrix3/mrtrix3.git mrtrix3 && cd mrtrix3 && git checkout dev && python configure -nogui && NUMBER_OF_PROCESSORS=1 python build && git describe --tags > /mrtrix3_version

#RUN echo $'FailOnWarn: 1\n' > /etc/mrtrix.conf
# Environment variables setup
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSL_DIR="${FSLDIR}"
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLMULTIFILEQUIT=TRUE
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV PATH=/opt/freesurfer/bin:/opt/freesurfer/mni/bin:/usr/lib/fsl/5.0:/usr/lib/ants:/mrtrix3/bin:/opt/eddy:$PATH
ENV PYTHONPATH=/mrtrix3/lib

ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLMULTIFILEQUIT=TRUE
# Note: Would prefer NIFTI, but need to stick to compressed for now due to FSL Ubuntu not honoring this variable. May be able to revert once fsl.checkFirst() is merged in.
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV PATH=/usr/lib/fsl/5.0:/usr/lib/ants:/mrtrix3/bin:$PATH
ENV PYTHONPATH=/mrtrix3/lib







RUN mkdir /in
RUN mkdir /out
RUN mkdir /hsm
RUN mkdir /scratch
RUN mkdir /vlsci

COPY preproc_to_template.py /code/run.py
RUN chmod 775 /code/run.py

ENTRYPOINT ["/code/run.py"]
