import numpy as np
import os
import argparse
import struct
from astropy import units as u
from astropy import constants as const

path = os.path.abspath(os.path.dirname(__file__))

class ICs:
    '''
    Parent class for handling units and the computational volume.
    '''
    def __init__(self,units='CGS'):
        self._set_units(units)
        
    def _set_volume(self):
        # Need to make adjustable
        # Volume parameters (with default values)
        self.X = 25.*self.AU.value 
        self.Y = 25.*self.AU.value
        self.Z = 25.*self.AU.value
        self.rhoAtm = (1e-26 * self.vDens).value  # Volume Density
        self.tempAtm = 600                        # Temperature in Kelvin
    
    def _set_units(self,units='cgs'):
        _units = ['cgs','CGS']
        assert units in _units, 'Invalid units. Try cgs.'
        if units == 'cgs' or units == 'CGS':
            self.G = const.G.cgs
            self.kb = const.k_B.cgs
            self.Msol = const.M_sun.cgs
            self.AU = const.au.cgs
            self.mp = const.m_p.cgs
            self.sDens = u.g / (u.cm * u.cm)
            self.vDens = u.g / (u.cm * u.cm * u.cm)
        else:
            # Add more units later
            print('Invalid units.')
            
    def _generate_mesh_points(self,method='glass'):
        """Create a volume of mesh generating points"""    
    
        if method == 'grid':
            pos = self._makeGridTile()
        elif method == 'glass':
            pos = self._makeGlassTile()
        else:
            raise ValueError('Invalid tiling method. Try glass.')
    
        N = 12
        spacing = 1.0
        numParticles = int(pos.size/3)
        posArray = np.zeros([int(N*N*N*numParticles),3])
    
        totNumParticles = posArray.size/3
        offset = np.zeros(3)

        for i in range(N) :
            offset[0] = 1.*i - 0.5*N;

            for j in range(N) :
                offset[1] = 1.*j - 0.5*N;

                for k in range(N) :
                    offset[2] = 1.*k - 0.5*N;
                    pStart = int(i*N*N+j*N + k)*numParticles
                    pEnd = pStart+numParticles
                    posArray[pStart:pEnd,:] = spacing*(pos[0:numParticles,:] + offset[:])
   
        self.points = posArray/N*2
    
    def _makeGlassTile(self,filename=path+'/utils/d2300.glass'):
        """Make glass tile of particles to fill volume"""

        f = open(filename,'r')
    
        nf77,v0,v1,v2,nf77e = np.fromfile(f,dtype=np.int32,count=5)
        array = np.fromfile(f,dtype=np.int32,count=402)
        nobj=swap_endian(array[302]) # get the number of gas particles

        # read type
        nf77 = np.fromfile(f,dtype=np.int32,count=1)
        types = np.fromfile(f,dtype=np.float32,count=nobj)
        nf77e = np.fromfile(f,dtype=np.int32,count=1)

        # read mass
        nf77 = np.fromfile(f,dtype=np.int32,count=1)
        mass = np.fromfile(f,dtype=np.float32,count=nobj)
        nf77e = np.fromfile(f,dtype=np.int32,count=1)

        # read position
        nf77 = np.fromfile(f,dtype=np.int32,count=1)
        pos = np.fromfile(f,dtype=np.float32,count=nobj*3)
        nf77e = np.fromfile(f,dtype=np.int32,count=1)

        # read velocity
        nf77 = np.fromfile(f,dtype=np.int32,count=1)
        vel = np.fromfile(f,dtype=np.float32,count=nobj*3)
        nf77e = np.fromfile(f,dtype=np.int32,count=1)

        for i in range(pos.size) :
            pos[i] = swap_endianf(pos[i])
    
        f.close()

        return pos.reshape([nobj,3])

    def _makeGridTile(self):
        """Make grid tile of particles to fill volume."""

        pos = np.zeros( [4096,3])
        N = 16
        ipart = 0
        for i in range(N) :
            for j in range(N) :
                for k in range(N) :
                    noise = (np.random.rand(3) - 0.5) * 0.025/N
                    x = (1.*i + 0.5)/N + noise[0]
                    y = (1.*j + 0.5)/N + noise[1]
                    z = (1.*k + 0.5)/N + noise[2]
                    pos[ipart,0] = x
                    pos[ipart,1] = y
                    pos[ipart,2] = z
                    ipart = ipart + 1

        return pos

# Some utilites for makeGlassTile().
def swap_endian(i):
    return struct.unpack("<I", struct.pack(">I",i))[0]
def swap_endianf(i):
    return struct.unpack("<f", struct.pack(">f",i))[0]
