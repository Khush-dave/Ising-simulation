# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 17:22:00 2021

@author: khush
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from tqdm import tqdm

plt.rcParams['animation.ffmpeg_path'] = 'C:\\ffmpeg\\bin\\ffmpeg.exe'

class IsingLattice:

    def __init__(self, temperature, initial_state, size):
        self.size = size
        self.T = temperature
        self.system = self._build_system(initial_state)

    @property
    def sqr_size(self):
        return (self.size, self.size)

    def _build_system(self, initial_state):
        """Build the system
        Build either a randomly distributed system or a homogeneous system (for
        watching the deterioration of magnetization
        Parameters
        ----------
        initial_state : str: "r" or other
            Initial state of the lattice.  currently only random ("r") initial
            state, or uniformly magnetized, is supported
        """

        if initial_state == 'r':
            system = np.random.choice([-1, 1], self.sqr_size)
        elif initial_state == 'u':
            system = np.ones(self.sqr_size)
        else:
            raise ValueError(
                "Initial State must be 'r', random, or 'u', uniform"
            )

        return system

    def _bc(self, i):
        """Apply periodic boundary condition
        Check if a lattice site coordinate falls out of bounds. If it does,
        apply periodic boundary condition
        Assumes lattice is square
        Parameters
        ----------
        i : int
            lattice site coordinate
        Return
        ------
        int
            corrected lattice site coordinate
        """
        if i >= self.size:
            return 0
        if i < 0:
            return self.size - 1
        else:
            return i

    def energy(self, N, M):
        """Calculate the energy of spin interaction at a given lattice site
        i.e. the interaction of a Spin at lattice site n,m with its 4 neighbors
        - S_n,m*(S_n+1,m + Sn-1,m + S_n,m-1, + S_n,m+1)
        Parameters
        ----------
        N : int
            lattice site coordinate
        M : int
            lattice site coordinate
        Return
        ------
        float
            energy of the site
        """
        return -self.system[N, M]*(
            K*self.system[self._bc(N - 1), M] + K*self.system[self._bc(N + 1), M]
            + L*self.system[N, self._bc(M - 1)] + L*self.system[N, self._bc(M + 1)]
            + D*self.system[self._bc(N-1),self._bc(M-1)] + D*self.system[self._bc(N+1),self._bc(M+1)] 
            + H
        )
    
    @property
    def heat_capacity(self):
        e = 0
        E = 0
        E_2 = 0

        for i in range(self.size):
            for j in range(self.size):
                e = self.energy(i, j)
                E += e
                E_2 += e**2

        U = E/(self.size**2)
        U_2 = E_2/(self.size**2)    
        return U_2 - U**2

    @property
    def magnetization(self):
        """Find the overall magnetization of the system
        """
        #return np.sum(self.system)/self.size**2
        e = 0
        E = 0

        for i in range(self.size):
            for j in range(self.size):
                e = self.system[i, j]
                if (i+j)%2 == 0:
                    E = E + e
                else:
                    E = E - e
                

        
        U= np.abs(E)/(self.size**2)    
        return U
    @property
    def susceptibility(self):
        e = 0
        E = 0
        E_2 = 0

        for i in range(self.size):
            for j in range(self.size):
                e = self.system[i, j]
                if (i+j)%2 == 0:
                    E = E + e
                    E_2 = E_2 + e**2
                else:
                    E = E - e
                    E_2 = E_2 + e**2

        U = E/(self.size**2)
        U_2 = E_2/(self.size**2)    
        return U_2 - U**2
    @property
    def capital_M2(self):
        """Find the overall magnetization of the system
        """
        return np.sum(self.system**2)
    @property
    def capital_M1(self):
        """Find the overall magnetization of the system
        """
        return np.sum(self.system)
    

def run(lattice,T, epochs,video=True):
    """Run the simulation
    """

    M_av=0
    m_av=0
    n_av=0
    z=0
    #with writer.saving(fig, "ising.mp4", 100):
    for epoch in range(epochs):
        # Randomly select a site on the lattice
        #N, M = np.random.randint(0, lattice.size, 2)
        N=np.random.randint(0,lattice.size)
        M=np.random.randint(0,lattice.size)
            # Calculate energy of a flipped spin
        E = -1*lattice.energy(N, M)
            # "Roll the dice" to see if the spin is flipped
        if E < 0.:
            lattice.system[N, M] = -1*lattice.system[N,M]
        elif np.exp(-(E - lattice.energy(N, M))/T) >= np.random.rand():
            lattice.system[N, M] = -1*lattice.system[N,M]
        
    print (T)
    '''to returen specific heat '''
    #return ((lattice.heat_capacity)/(T*T))
    '''to return the state of the system'''
    #return lattice.system
    #return np.abs(np.sum(lattice.system)/(lattice.size**2))
    ''' to return the staggered magnetization '''
    #return lattice.magnetization
    '''to return the susceptibility'''
    return lattice.susceptibility/(lattice.size*lattice.size*T)

    

def chart_plot(T,Q1,Q2,Q3):
    plt.figure()
    plt.title('Staggerred Susceptibility vs Temperature')
    plt.xlabel('Temperature (J/Kb)')
    plt.ylabel('Staggered Susceptibility')
    plt.plot(T,Q1,'--bo',color='blue',label='Size 10x10') 
    plt.plot(T,Q2,'--bo',color='green',label='Size 20x20')
    plt.plot(T,Q3,'--bo',color='red',label='Size 40x40')
    plt.legend()

K=-1
L=-1   
D=-1

H=1
#temperature = np.linspace(2.5,0.1,30)
initial_state='u'
size=10
epochs=10000
video=True
Mhash=[]
M10=[]
M20=[]
M40=[]
T=[]

for i in range(3):   
    lattice = IsingLattice(0.5, initial_state=initial_state, size=size)
    Mhash.append(run(lattice,0.5,500000,video))
    for temperature in np.arange(0.5,4,0.1):
        #lattice = IsingLattice(temperature, initial_state=initial_state, size=size)
        if i==0:
            M10.append(run(lattice,temperature, 2*epochs,video))
            T.append(temperature)
        elif i==1:
            M20.append(run(lattice,temperature ,epochs,video))
        else :
            M40.append(run(lattice, temperature,epochs,video))
    size*=2


chart_plot(T,M10,M20,M40)

