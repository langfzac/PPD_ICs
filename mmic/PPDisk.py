from .mmic import ICs
import numpy as np
from pytipsy import wtipsy

class PPDisk(ICs):
    '''
    Class for holding protoplanetary disk parameters and generating initial conditions.
    
    '''
    def __init__(self,units='CGS'):
        
        self._set_units(units)
        self._set_volume()
        
        # Disk Parameters (with default values)
        self.rhoPower = -1.0                       # Radial density profile power
        self.tempPower = -0.5                      # Temperature profile power
        self.rout = 1.0*self.AU.value              # Outer radius in AU
        self.rin = 0.3*self.AU.value               # Inner radius of the disk in AU
        self.Mstar = 1./3.*self.Msol.value         # Mass of the star
        self.T0 = 150.                             # Central temperature
        self.sigma0 = (48000 * self.sDens).value   # Central surface density
    
    def makeIC(self,filename='snapshot.std',output=True):
        '''Generate tipsy snapshot'''
        try:
            self.points
        except:
            print('Generating initial mesh points...')
            self._generate_mesh_points()
        try:
            self.pos
        except: 
            print('Building disk...')
            self.makeDisk()
        
        narr = np.array([])
        zarr = np.zeros(1)
        zsarr = np.zeros(self.ngas)
        n = self.ngas + 1
        self.header={'time':0.,'n':n,'ndim':3,'ngas':self.ngas,'ndark':1,'nstar':0}
        
        # Gas quantities
        rho = np.ones(self.ngas)
        h = self.rin * 0.1 * np.ones(self.ngas)
        x = self.pos[:,0]; y = self.pos[:,1]; z = self.pos[:,2]
        vx = self.vel[:,0]; vy = self.vel[:,1]; vz = self.vel[:,2]
        temp = self.T;
        
        # Dark matter quantities (should probably let this be set by user)
        dmeps = np.array([2.0*self.rin])
        
        self.catg = {'mass':self.mass,'x':x,'y':y,'z':z,'vx':vx,'vy':vy,'vz':vz,'dens':rho,'tempg':temp,'h':h,'zmetal':zsarr,'phi':zsarr}
        self.catd = {'mass':[self.Mstar], 'x':zarr,'y':zarr,'z':zarr,'vx':zarr,'vy':zarr,'vz':zarr,'eps':dmeps,'phi':zarr}
        self.cats = {'mass':narr,'x':narr,'y':narr,'z':narr,'vx':narr,'vy':narr,'vz':narr,'metals':narr,'tform':narr,'eps':narr,'phi':narr}
        
        print('Writing file...')
        wtipsy(filename,self.header,self.catg,self.catd,self.cats)
        print('Done.\nWritten to: {0}'.format(filename))
    
    def makeDisk(self,doBC = False): 
        """Populate disk particles with rho, temp, vel."""
        
        self.doBC = doBC
        posArray = self.points
        xmax,ymax,zmax = self.X,self.Y,self.Z
        #normalize posArray to span from -xmax to xmax (in the box)
        rout = self.rout
        posArray *= self.AU.value
        posArray *= 2.
        ro = np.sqrt(posArray[:,0]**2 + posArray[:,1]**2)
        z = np.abs(posArray[:,2])
        scale = self._scaling(ro, z, rout)
     
        cellVolume = posArray[:,0].max()*posArray[:,1].max()*posArray[:,2].max()*8*np.ones(posArray.shape[0])/posArray.shape[0]
        posArray[ro > rout] *= scale[ro > rout]
        posArray[ro <= rout] *= scale[ro <= rout]
        cellVolume *= scale[:,0]*scale[:,1]*scale[:,2]
    
        # cut out stuff not in the box
        mask = np.abs(posArray[:,0]) < xmax
        mask = np.logical_and(mask, np.abs(posArray[:,1]) < ymax)
        mask = np.logical_and(mask, np.abs(posArray[:,2]) < zmax)
        posArray = posArray[mask]
        cellVolume = cellVolume[mask]

        # Populate disk particles
        r = np.sqrt(posArray[:,0]**2 + posArray[:,1]**2) 
        z = posArray[:,2]
        rho, T, vorb, h = self._disk_profile(z, r, mu = 2.)

        # Rotate counter-clockwise
        R = np.sqrt(r*r + z*z)
        vArray = np.zeros(posArray.shape)
        vArray[:,1] = vorb[:]*posArray[:,0]/R[:]
        vArray[:,0] = -vorb[:]*posArray[:,1]/R[:]

        # Populate physical cell quantities
        rhoMax = rho.max()
    
        mArray = rho*cellVolume 
        #mmin = mArray[mArray != 0.0].min()
        #mArray[mArray == 0.0] = mmin
        eArray = T

        assert (mArray.size == eArray.size),"Mass and Temp arrays not the same size"
    
        print('Disk Created')
        disk_p = mArray[T < self.tempAtm]
        disk_mass = (np.sum(disk_p) / self.Msol.value)
        print("Disk Mass: {0} Msol".format(disk_mass))
        print("Disk Particles: {0}".format(len(disk_p)))
        print('Total Particles: {0}'.format(len(eArray)))

        self.pos = posArray 
        self.vel = vArray
        self.T = eArray
        self.mass = mArray
        self.rhoMax = rhoMax
        self.ngas = len(mArray)
        
    def _disk_profile(self, z, r, mu = 2.) :
        """Generate disk profiles and particle velocities."""

        G,kB,Msol,AU,mp = self.G.value,self.kb.value,self.Msol.value,self.AU.value,self.mp.value
        n,p,rin,rout,Mstar,T0,sigma0,rhoAtm,tempAtm = self.rhoPower,self.tempPower,self.rin,self.rout,self.Mstar,self.T0,self.sigma0,self.rhoAtm,self.tempAtm
    
        Sigma = sigma0*(r/AU)**(n)
        T = T0*(r/AU)**(p)
        cs = np.sqrt(kB*T/(mu*mp))
        vKep = np.sqrt(G*Mstar/r**3)
        h = cs/vKep
        rhoMid = (1.0/np.sqrt(2.0*np.pi))*Sigma/h
        rho = np.maximum(rhoMid*np.exp(-z*z/(2.0*h*h)), rhoAtm)

        # two smooth cutoffs
        L = 0.3*rout
        rho[r<rin] = rhoAtm#np.maximum(rhoAtm, rho*(r/rin)**27)[r<rin]
        rho[r>=rout] = np.maximum(rhoAtm, (rho*np.exp((-(r-rout)**2)/L**2)))[r>=rout]
        #rho[r<=0.08*rout] = rhoAtm
    
        tempMask = rho == rhoAtm
        T[tempMask] = tempAtm
    
        # For inner BC
        if self.doBC:
            rsph = np.sqrt(r*r + z*z)
            BC = rsph < rin
            T[BC] = 50.
            rho[BC] = 1e-30
    
        # G = 1 in code units
        # Calculate vorb separately, in code units
        vorb = np.sqrt(Mstar/np.sqrt(r*r + z*z))
        if self.doBC: vorb[BC] = 0. 
        #vorb[T == tempAtm] = 0.0
        #rho[rho == rhoAtm] = 0.0 
    
        #print("Number of disk particles: {0}".format(T[T!=tempAtm].size))

        return rho, T, vorb, h    
        
    def _scaling(self, r, z, rout, hr = 0.1) :
        r_sph = np.sqrt(z*z + r*r)
        z0 = hr*r
        shape = np.ones( [z.size, 3]) 
        shape *= np.minimum( np.maximum((r_sph/rout)**10,1.), 1e10)[:,np.newaxis]
        rltrout = r_sph < rout
        shape[rltrout,2] = np.minimum(np.maximum((z/z0), 1.), 1e10)[rltrout]

        # figure out double counted values
        r_mod = np.sqrt( (z*shape[:,2])**2 + (r*shape[:,1])**2)
        shape[np.logical_and(rltrout, r_mod > rout),2] = 1e30

        return shape
