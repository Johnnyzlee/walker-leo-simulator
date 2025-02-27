#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: satellite.py
Author: Li ZENG @ HKUST ECE
License: MIT License

Description:
This module defines the Satellite class for managing individual satellites in an orbit.
It includes methods for initializing satellite properties, calculating positions, and updating satellite states.
"""

import numpy as np

EARTH_RADIUS = 6371  # Define Earth's radius as a constant

class Satellite:
    """
    Class representing a satellite in an orbit.
    """
    def __init__(self, id, orbit):
        """
        Initialize the satellite's properties.

        Parameters:
        - id: Satellite ID ranging from 1 to num_sats.
        - orbit: Orbit object to which the satellite belongs.
        """
        self.id = id  # Satellite ID
        self.orbit = orbit  # Orbit object
        self.radius = orbit.radius  # Radius of the orbit in km
        self.inclination = orbit.inclination  # Inclination of the orbit in radians
        self.right_ascension = orbit.right_ascension  # Longitude of the ascending node in radians
        self.phase = 2 * np.pi * (id - 1) / orbit.num_sats - orbit.phase_offset  # Phase of the satellite in the orbit in radians
        self.angular_velocity = np.sqrt(398600 / self.radius**3)  # Angular velocity in rad/s
        self.rotation_matrix = self._get_rotation_matrix()  # Rotation matrix from the orbital plane to the ECEF coordinate system
        self.position_orbit = self.calculate_position_orbit()  # Position in the orbital plane (x, y, z) in km
        self.position_ecef = self.calculate_position_ecef()  # Position in the ECEF coordinate system (x, y, z) in km
        self.position_geodetic = self.ecef_to_geodetic()  # Position in geodetic coordinates (latitude, longitude, altitude)
        self.time = 0  # Current time in seconds

    def get_info(self):
        """
        Get the satellite's information.

        Returns:
        - Dictionary containing satellite information.
        """
        return {
            'id': self.id,  # Satellite ID
            'orbit_id': self.orbit.id,  # Orbit ID
            'altitude': self.radius - EARTH_RADIUS,  # Altitude in km
            'inclination (in degree)': np.degrees(self.inclination),  # Inclination in degrees
            'right_ascension (in degree)': np.degrees(self.right_ascension),  # Longitude of the ascending node in degrees
            'angular velocity (in rad/s)': self.angular_velocity,  # Angular velocity in rad/s
            'phase in orbit (in degree)': np.degrees(self.phase),  # Phase in the orbit in degrees
            'position_ecef': self.position_ecef,  # ECEF coordinates (x, y, z)
            'position_geodetic': self.position_geodetic,  # Geodetic coordinates (latitude, longitude, altitude)
            'current_time': self.time  # Current time in seconds
        }

    def calculate_position_orbit(self):
        """
        Calculate the satellite's position in the orbital plane.

        Returns:
        - Numpy array of the position (x, y, z) in km.
        """
        x = self.radius * np.cos(self.phase)
        y = self.radius * np.sin(self.phase)
        z = 0
        return np.array([x, y, z])

    def calculate_position_ecef(self):
        """
        Calculate the satellite's position in the Earth-Centered, Earth-Fixed (ECEF) coordinate system.

        Returns:
        - Numpy array of the position (x, y, z) in km.
        """
        position_ecef = np.dot(self.rotation_matrix, self.position_orbit)  # Rotate the position from the orbital plane to the ECEF coordinate system
        return position_ecef
    
    def ecef_to_geodetic(self):
        """
        Convert ECEF coordinates to geodetic coordinates (latitude, longitude, altitude).

        Returns:
        - Numpy array of the position (latitude, longitude, altitude).
        """
        x, y, z = self.position_ecef
        hyp = np.sqrt(x**2 + y**2 + z**2)  # Should be equal to self.radius
        lon = np.arctan2(y, x)
        if y < 0:
            lon -= np.pi  # Adjust the range of longitude to be (-pi, pi]
        lat = np.arctan2(z, hyp)
        alt = hyp - EARTH_RADIUS  # Earth's radius is 6371 km
        lon = np.degrees(lon)
        lat = np.degrees(lat)
        return np.array([lat, lon, alt])

    def update_satellite(self, delta_t):
        """
        Update the satellite's position based on the time increment.

        Parameters:
        - delta_t: Time increment in seconds.
        """
        self.time += delta_t
        self._update_position(delta_t)

    def _update_position(self, delta_t):
        """
        Update the satellite's position every delta_t seconds.

        Parameters:
        - delta_t: Time increment in seconds.
        """
        self.phase += self.angular_velocity * delta_t
        self.phase %= (2 * np.pi)
        self.position_orbit = self.calculate_position_orbit()
        self.position_ecef = self.calculate_position_ecef()
        self.position_geodetic = self.ecef_to_geodetic()

    def _get_rotation_matrix(self):
        """
        Calculate the rotation matrix from the orbital plane to the ECEF coordinate system.

        Returns:
        - Numpy array representing the rotation matrix.
        """
        # Note that from ECEF to the orbital plane, need to:
        # - firstly rotate by the inclination angle around the x-axis (anti-clockwise)
        # - then rotate by the longitude of the ascending node around the z-axis (anti-clockwise)
        Rz_omega = np.array([
            [np.cos(self.right_ascension), -np.sin(self.right_ascension), 0],
            [np.sin(self.right_ascension), np.cos(self.right_ascension), 0],
            [0, 0, 1]
        ])
        Rx_phi = np.array([
            [1, 0, 0],
            [0, np.cos(self.inclination), -np.sin(self.inclination)],
            [0, np.sin(self.inclination), np.cos(self.inclination)]
        ])
        return np.dot(Rz_omega, Rx_phi)
        