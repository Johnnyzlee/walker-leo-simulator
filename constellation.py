#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: constellation.py
Author: Li ZENG @ HKUST ECE
License: MIT License

Description:
This module defines the WalkerConstellation abstract base class and its derived classes for managing satellite constellations.
It includes methods for initializing constellations, updating satellite positions, and calculating phase differences.
"""

from abc import ABC, abstractmethod
import numpy as np
from orbit import Orbit

EARTH_RADIUS = 6371 # Define Earth's radius as a constant

class WalkerConstellation(ABC):
    """
    Abstract base class representing a satellite constellation.
    """
    def __init__(self, num_orbits, num_sats_per_orbit, radius):
        """
        Initialize the constellation's properties.

        Parameters:
        - num_orbits: Number of orbital planes.
        - num_sats_per_orbit: Number of satellites per orbital plane.
        - radius: Radius of the orbits in km.
        """
        self.time = 0  # Simulation time in seconds
        self.radius = radius  # Radius of the orbits in km
        self.altitude = radius - EARTH_RADIUS  # Altitude of the orbits in km (Earth's radius is EART_RADIUS = 6371 km)
        self.num_orbits = num_orbits  # Number of orbital planes
        self.num_sats_per_orbit = num_sats_per_orbit  # Number of satellites per orbital plane

    def update_constellation(self, delta_t):
        """
        Update the positions of all satellites in the constellation based on the time increment.

        Parameters:
        - delta_t: Time increment in seconds.
        """
        self.time += delta_t
        for orbit in self.orbits:
            orbit.update_orbit(delta_t)


class StarConstellation(WalkerConstellation):
    """
    Class representing a Walker Star Constellation.
    """
    def __init__(self, num_orbits=12, num_sats_per_orbit=30, radius=EARTH_RADIUS + 550.0):
        """
        Initialize the Walker Star Constellation.

        Parameters:
        - num_orbits: Number of orbital planes (default: 12).
        - num_sats_per_orbit: Number of satellites per orbital plane (default: 30).
        - radius: Radius of the orbits in km (default: EARTH_RADIUS + 550).
        """
        super().__init__(num_orbits, num_sats_per_orbit, radius)
        self.type = "Walker Star Constellation"
        self.inclination = np.radians(90)  # Inclination of 90 degrees for polar orbits
        self.phasediff = 0  # No phase difference for Walker Star Constellation
        self.orbits = [Orbit(self, orbit_id, self.radius, self.inclination, self.num_sats_per_orbit) for orbit_id in range(1, num_orbits + 1)]


class DeltaConstellation(WalkerConstellation):
    """
    Class representing a Walker Delta Constellation.
    """
    def __init__(self, num_orbits=12, num_sats_per_orbit=30, radius=EARTH_RADIUS + 550.0, inclination=53.0):
        """
        Initialize the Walker Delta Constellation.

        Parameters:
        - num_orbits: Number of orbital planes (default: 12).
        - num_sats_per_orbit: Number of satellites per orbital plane (default: 30).
        - radius: Radius of the orbits in km (default: EART_RADIUS + 550.0).
        - inclination: Inclination of the orbits in degrees (default: 53.0).
        """
        if inclination <= 0 or inclination >= 90:
            raise ValueError("Inclination must be between 0 and 90 degrees.")

        super().__init__(num_orbits, num_sats_per_orbit, radius)
        self.type = "Walker Delta Constellation"
        self.inclination = np.radians(inclination)  # Convert inclination to radians
        self.phasediff = self.calculate_phasediff()  # Calculate phase difference
        self.orbits = [Orbit(self, orbit_id, self.radius, self.inclination, self.num_sats_per_orbit) for orbit_id in range(1, num_orbits + 1)]

        
    def calculate_phasediff(self):
        """
        Calculate the phase difference between the ascending nodes of adjacent orbits.

        This helps in finding the neighboring satellites with inter-orbit inter-satellite links.

        Returns:
        - phasediff: Phase difference in radians.
        """
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