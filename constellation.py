#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: constellation.py
Author: Li ZENG @ HKUST ECE
License: MIT License
"""

from abc import ABC, abstractmethod
import numpy as np
from orbit import Orbit

class WalkerConstellation(ABC):
    def __init__(self, num_orbits, num_sats_per_orbit, radius):
        self.time = 0
        self.radius = radius
        self.altitude = radius - 6371
        self.num_orbits = num_orbits
        self.num_sats_per_orbit = num_sats_per_orbit
        return
    
    def update_constellation(self, delta_t):
        # Update the positions of all satellites in the constellation based on the time increment
        self.time += delta_t
        for orbit in self.orbits:
            orbit.update_orbit(delta_t)
        return


class StarConstellation(WalkerConstellation):
    def __init__(self, num_orbits=12, num_sats_per_orbit=30, radius=6371 + 550):
        super().__init__(num_orbits, num_sats_per_orbit, radius)
        self.type = "Walker Star Constellation"
        self.inclination = np.radians(90)
        self.phasediff = 0
        self.orbits = [Orbit(self, orbit_id, self.radius, self.inclination, self.num_sats_per_orbit) for orbit_id in range(1, num_orbits + 1)]       
        return


class DeltaConstellation(WalkerConstellation):
    def __init__(self, num_orbits=12, num_sats_per_orbit=30, radius=6371.0 + 550.0, inclination=53.0):
        if inclination <= 0 or inclination >= 90:
            raise ValueError("Inclination must be between 0 and 90 degrees.") 
        
        super().__init__(num_orbits, num_sats_per_orbit, radius)
        self.type = "Walker Delta Constellation" 
        self.inclination = np.radians(inclination)
        self.phasediff = self.calculate_phasediff()
        self.orbits = [Orbit(self, orbitd_id, self.radius, self.inclination, self.num_sats_per_orbit) for orbitd_id in range(1, num_orbits + 1)]
        return
        
    def calculate_phasediff(self):
        # Calculate the phase difference between the ascending nodes of adjacent orbits, which helps find the neighboring satellites with inter-orbit inter-satellite links
        # Elements: (all in geodetic coordinates, i.e., latitude, longitude, height)
        #   - orbital plane 0: the equator
        #   - orbital plane 1: of which the ascending node is at (0, 0, radius)
        #   - orbital plane 2: of which the ascending node is at (omega, 0, radius)
        #   - orbital plane 3: perpendicular to plane 1 with the intersection point at (0, 0, radius), of which the ascending node is at (0, 0, radius) and the inclination is inclination + 90
        # Notes:
        #   - equation of orbital plane 0: x^2 + y^2 = radius^2, z = 0
        #   - transform from orbital plane 0 to orbital plane 2: rotate by phi around the x-axis (anti-clockwise) and then rotate by omega around the z-axis (anti-clockwise)
        #   - transform from orbital plane 0 to orbital plane 3: rotate by (phi + 90) around the x-axis (anti-clockwise)
        #   - transform from orbital plane 2 to orbital plane 0: rotate by -omega around the z-axis and then rotate by -phi around the x-axis
        #   - transform from orbital plane 3 to orbital plane 0: rotate by -(phi + 90) around the x-axis
        omega = 2 * np.pi / self.num_orbits # omega represents the longitude of the ascending node of plane 2 (in radians)
        phi = self.inclination # phi represents the inclination angle (in radians)
        # Steps: 
        #   1. write the normal vector of plane 3: v3 = [0, -cos(phi), -sin(phi)]
        v3 = np.array([0, -np.cos(phi), -np.sin(phi)])
        #   2. write the normal vector of plane 2: v2 = [sin(omega) * sin(phi), -cos(omega) * sin(phi), cos(phi)]
        v2 = np.array([np.sin(omega) * np.sin(phi), -np.cos(omega) * np.sin(phi), np.cos(phi)])
        #   3. compute the cross product of v2 and v3 to obtain the intersection line's direction vector: v = v2 x v3 = [cos(phi)**2 + cos(omega) * sin(phi)**2, sin(omega) * sin(phi)**2, -sin(omega) * cos(phi) * sin(phi)]
        v = np.cross(v3, v2)
        #   4. compute phase difference between v and the direction of the ascending node of plane 2, i.e., [cos(omega), sin(omega), 0]
        u = np.array([np.cos(omega), np.sin(omega), 0])
        dot_product = np.dot(v, u)
        magnitude_v, magnitude_u = np.linalg.norm(v), np.linalg.norm(u)
        phasediff = np.arccos(dot_product / (magnitude_v * magnitude_u)) 
        phasediff = np.pi / 2 - np.abs(np.pi / 2 - phasediff) # prevent phasediff from being larger than pi / 2
        # Note: theta in orbit p corresponds to theta - phasediff in orbit p + 1 (since we have enforced the inclination to be in range (0, pi / 2))
        
        # Use phasediff to find the neighboring satellites with inter-orbit inter-satellite links: 
        #   - The satellite in plane 2 with an ISL to the ascending node of plane 1 should the satellite with the most similar phase to the phase difference
        #   - The satellite in plane N - 1 with an ISL to the ascending node of plane 1 should the satellite with the most similar phase to negative phase difference
        return phasediff