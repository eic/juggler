FROM  eicweb.phy.anl.gov:4567/containers/eic_container/eic:2.3.0

LABEL maintainer="Whitney Armstrong <warmstrong@anl.gov>" \
      name="juggler" \
      group="EIC/juggler" \
      march="native" \
      basedist="debian" \
      base="eic_container/eic" \
      version="2.2.2"

RUN cd /tmp \
      && git clone https://gitlab.cern.ch/hepmc/HepMC3.git  \
      && mkdir HepMC3/build && cd HepMC3/build \
      && cmake ../. -DHEPMC3_ENABLE_ROOTIO=ON -DCMAKE_INSTALL_PREFIX=/usr/local -DHEPMC3_BUILD_EXAMPLES=OFF \
      && make -j20 \
      && make install  \
      && cd /tmp && rm -rf /tmp/HepMC3

RUN cd /tmp \
      && wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8303.tgz \
      && tar -zxf pythia8303.tgz && cd pythia8303 \
      && ./configure --prefix=/usr/local --enable-shared --enable-64bit --with-root --with-hepmc3 \
      && make -j30 \
      && make install  \
      && cd /tmp && rm -rf /tmp/pythia8303


RUN cd /tmp             \
      && git clone  https://eicweb.phy.anl.gov/EIC/NPDet.git \
      && mkdir -p NPDet/build && cd NPDet/build \
      && cmake ../. -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=/usr/local \
      && make -j30  && make -j4 install \
      && cd /tmp && rm -rf /tmp/NPDet 

RUN cd /tmp             \
      && git clone  https://eicweb.phy.anl.gov/EIC/eicd.git \
      && mkdir -p eicd/build && cd eicd/build \
      && cmake ../. -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=/usr/local \
      && make -j30  && make -j4 install \
      && cd /tmp && rm -rf /tmp/eicd 

RUN cd /tmp             \
      && git clone  https://eicweb.phy.anl.gov/EIC/juggler.git \
      && mkdir -p juggler/build && cd juggler/build \
      && cmake ../. -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=/usr/local \
      && make -j30  && make -j4 install \
      && cd /tmp && rm -rf /tmp/juggler 

