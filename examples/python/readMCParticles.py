#!/usr/bin/env python
#
# This is just a simple test script to check whether the python bindings are
# actually working as intended

import pyLCIO as lcio
import sys

reader = lcio.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(sys.argv[1])

for event in reader:
    mcs = event.getCollection('MCParticle')
    for mc in mcs:
        print(mc.getEnergy())
