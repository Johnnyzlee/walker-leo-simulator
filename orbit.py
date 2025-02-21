#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: orbit.py
Author: Li ZENG @ HKUST ECE
License: MIT License
"""

import numpy as np
from satellite import Satellite

class Orbit:
    def __init__(self, constellation, id, radius, inclination, num_sats):
        # Initialize the orbit's properties
        self.id = id  # Orbit ID
        self.constellation = constellation
        self.radius = radius # Radius of the orbit in km
        self.inclination = inclination # Inclination of the orbit in radians
        self.num_sats = num_sats # Number of satellites in the orbit
        self.phase_offset = self.constellation.phasediff * (self.id - 1) # Phase offset for the orbit compared to orbit 1
        if self.constellation.type == "Walker Delta Constellation":
            self.right_ascension = 2 * np.pi * (id - 1) / constellation.num_orbits  # Longitude of the ascending node
        elif self.constellation.type == "Walker Star Constellation":
            self.right_ascension = 1 * np.pi * (id - 1) / constellation.num_orbits  # Longitude of the ascending node
        self.sats = [Satellite(i + 1, self) for i in range(num_sats)] 
        return

    def get_info(self):
        # Get the orbit's information
        return {
            'id': self.id,
            'radius': self.radius,
            'altitude': self.radius - 6371,
            'inclination': np.degrees(self.inclination),
            'num_sats': self.num_sats,
            'right_ascension': np.degrees(self.right_ascension)
        }
        
    def get_satellite_positions(self):
        # Get the ECEF positions of all satellites in the orbit (x, y, z)
        positions = [sat.position_ecef for sat in self.sats] 
        return np.array(positions) 
    
    def update_orbit(self, delta_t):
        # Update the positions of all satellites in the orbit based on the time increment
        for sat in self.sats:
            sat.update_satellite(delta_t)
        return