#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: satellite.py
Author: Li ZENG @ HKUST ECE
License: MIT License
"""

import numpy as np

'''
Satellite class has the following properties:
    - id: Satellite ID ranging from 1 to num_sats
    - orbit: Orbit object to which the satellite belongs
    - radius: Radius of the orbit in km (6371 + altitude)
    - inclination (phi): Inclination of the orbit in radians
    - right_ascension (omega): Longitude of the ascending node in radians of the orbit
    - phase: Phase of the satellite in the orbit, sat 1 has initial phase -orbit.phase_offset, etc.
    - rotation_matrix: Rotation matrix from the orbital plane to the ECEF coordinate system
    - position_orbit: Position of the satellite in the orbital plane (x, y, z) (in km)
    - position_ecef: Position of the satellite in the ECEF coordinate system (x, y, z) (in km)
        - x-axis points to the intersection of the prime meridian and the equator
        - y-axis points to the 90 degrees East longitude on the equator
        - z-axis points to the North Pole
    - position_geodetic: Position of the satellite in geodetic coordinates (latitude, longitude, altitude) 
        - latitude, longitude are in degrees, altitude is in km
        - north latitude and east longitude are positive
        - latitude ranges from -90 to 90 degrees, longitude ranges from -180 to 180 degrees
    - time: Current time in seconds
'''

class Satellite:
    def __init__(self, id, orbit):
        # Initialize the satellite's properties
        self.id = id
        self.orbit = orbit
        self.radius = orbit.radius
        self.inclination = orbit.inclination
        self.right_ascension = orbit.right_ascension
        self.phase = 2 * np.pi * (id - 1) / orbit.num_sats - orbit.phase_offset # in radians
        self.angular_velocity = np.sqrt(398600 / self.radius**3) # in rad/s
        self.rotation_matrix = self._get_rotation_matrix()
        self.position_orbit = self.calculate_position_orbit()
        self.position_ecef = self.calculate_position_ecef()
        self.position_geodetic = self.ecef_to_geodetic()
        self.time = 0
        return
        
    def get_info(self):
        return {
            'id': self.id, # Satellite ID ranging from 1 to num_sats
            'orbit_id': self.orbit.id, # Orbit ID ranging from 1 to num_orbits
            'altitude': self.radius - 6371, # Altitude of the satellite in km
            'inclination (in degree)': np.degrees(self.inclination), # Inclination of the orbit in degrees
            'right_ascension (in degree)': np.degrees(self.right_ascension), # Longitude of the ascending node in degrees
            'angular velocity (in rad/s)': self.angular_velocity, # Angular velocity of the satellite in rad/s
            'phase in orbit (in degree) (ascending node has phase 0)': np.degrees(self.phase), # Phase of the satellite in the orbit in degrees
            'position_ecef': self.position_ecef, # ECEF coordinates (x, y, z)
            'position_geodetic': self.position_geodetic, # Geodetic coordinates (latitude, longitude, altitude)
            'current_time': self.time # Current time in seconds
        }

    def calculate_position_orbit(self):
        # Calculate the satellite's position in the orbital plane
        x = self.radius * np.cos(self.phase)
        y = self.radius * np.sin(self.phase)
        z = 0
        return np.array([x, y, z])

    def calculate_position_ecef(self):
        # Calculate the satellite's position in the Earth-Centered, Earth-Fixed (ECEF) coordinate system
        position_ecef = np.dot(self.rotation_matrix, self.position_orbit) # Rotate the position from the orbital plane to the ECEF coordinate system
        return position_ecef
    
    def ecef_to_geodetic(self):
        # Convert ECEF coordinates to geodetic coordinates (latitude, longitude, altitude)
        x, y, z = self.position_ecef
        hyp = np.sqrt(x**2 + y**2 + z**2) # should be equal to self.radius
        lon = np.arctan2(y, x)
        if y < 0:
            lon -= np.pi # Adjust the range of longitude to be (-pi, pi]
        lat = np.arctan2(z, hyp)
        alt = hyp - 6371  # Earth's radius is 6371 km
        lon = np.degrees(lon)
        lat = np.degrees(lat)
        return np.array([lat, lon, alt]) 

    def update_satellite(self, delta_t):
        self.time += delta_t
        self._update_position(delta_t)
        return
    
    def _update_position(self, delta_t):
        # Update the satellite's position every delta_t seconds
        angular_velocity = self.angular_velocity
        self.phase += angular_velocity * delta_t
        self.phase %= (2 * np.pi)
        self.position_orbit = self.calculate_position_orbit()
        self.position_ecef = self.calculate_position_ecef()
        self.position_geodetic = self.ecef_to_geodetic()
        return

    def _get_rotation_matrix(self):
        '''
        Calculate the rotation matrix from the orbital plane to the ECEF coordinate system
        Note that from ECEF to the orbital plane, need to:
        - firstly rotate by the inclination angle around the x-axis (anti-clockwise)
        - then rotate by the longitude of the ascending node around the z-axis (anti-clockwise)
        Hence, from the orbital plane to ECEF, need to:
        - firstly rotate by the longitude of the ascending node (omega) around the z-axis (clockwise, hence negative)
        - then rotate by the negative inclination angle around the x-axis (clockwise, hence negative)
        '''
        Rz_omega = np.array([
            [np.cos(-self.right_ascension), -np.sin(-self.right_ascension), 0],
            [np.sin(-self.right_ascension), np.cos(-self.right_ascension), 0],
            [0, 0, 1]
        ])
        Rx_phi = np.array([
            [1, 0, 0],
            [0, np.cos(-self.inclination), -np.sin(-self.inclination)],
            [0, np.sin(-self.inclination), np.cos(-self.inclination)]
        ])
        return np.dot(Rx_phi, Rz_omega)