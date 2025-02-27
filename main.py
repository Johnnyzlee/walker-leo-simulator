#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: main.py
Author: Li ZENG @ HKUST ECE
License: MIT License
"""

from constellation import StarConstellation, DeltaConstellation
from network import SatNetStar, SatNetDelta
import time

EARTH_RADIUS = 6371  # Define Earth's radius as a constant

def main(constellation_type, *args, **kwargs):
    """
    Initialize and run the constellation network.

    Parameters:
    - constellation_type (str): Type of the constellation ("Walker Star Constellation" or "Walker Delta Constellation").
    - **kwargs: Additional keyword arguments for the constellation.

    Required keyword arguments for both constellations:
    - duration (int): Duration of the simulation in seconds (default: 3600).
    - delta_t (int): Time step for the simulation in seconds (default: 1).

    Required keyword arguments for "Walker Star Constellation":
    - num_orbits (int): Number of orbits.
    - num_sats_per_orbit (int): Number of satellites per orbit.
    - altitude (float): Altitude of the orbits.
    - radius (float): Radius of the orbits.

    Required keyword arguments for "Walker Delta Constellation":
    - num_orbits (int): Number of orbits.
    - num_sats_per_orbit (int): Number of satellites per orbit.
    - altitude (float): Altitude of the orbits.
    - radius (float): Radius of the orbits.
    - inclination (float): Inclination of the orbits.
    """
    # Get duration and delta_t from kwargs with default values
    duration = kwargs.pop('duration', 3600)
    delta_t = kwargs.pop('delta_t', 1)

    # Required parameters for each constellation type
    required_params_star = ['num_orbits', 'num_sats_per_orbit', 'altitude']
    required_params_delta = ['num_orbits', 'num_sats_per_orbit', 'altitude', 'inclination']

    # Initialize the constellation and network based on the constellation type
    if constellation_type == "Walker Star Constellation":
        for param in required_params_star:
            if param not in kwargs:
                raise ValueError(f"Missing required parameter for Star Constellation: {param}")
        kwargs['radius'] = EARTH_RADIUS + kwargs.pop('altitude')
        constellation = StarConstellation(*args, **kwargs)
        network = SatNetStar(constellation)
    elif constellation_type == "Walker Delta Constellation":
        for param in required_params_delta:
            if param not in kwargs:
                raise ValueError(f"Missing required parameter for Delta Constellation: {param}")
        kwargs['radius'] = EARTH_RADIUS + kwargs.pop('altitude')
        constellation = DeltaConstellation(*args, **kwargs)
        network = SatNetDelta(constellation)
    else:
        raise ValueError("Unsupported constellation type")

    print(f"Initialized {constellation_type} with parameters: {args}, {kwargs}")

    for step in range(0, duration, delta_t):
        print(f"Time step: {step} / {duration / delta_t}, delta_t: {delta_t} seconds")
        start_time = time.time()
        network.update_network(delta_t)
        end_time = time.time()
        elapsed_time = end_time - start_time
        # print(f"Elapsed time: {elapsed_time} seconds")

    return


if __name__ == "__main__":
    main("Walker Delta Constellation",  # Type of the constellation
         duration=3600,                 # Duration of the simulation in seconds
         delta_t=1,                     # Time step for the simulation in seconds 
         num_orbits=16,                 # Number of orbits
         num_sats_per_orbit=30,         # Number of satellites per orbit
         altitude=550,                  # Altitude of the orbits in km
         inclination=53)                # Inclination of the orbits in degrees
    
    # main("Walker Star Constellation",   # Type of the constellation
    #      duration=3600,                 # Duration of the simulation in seconds
    #      delta_t=1,                     # Time step for the simulation in seconds
    #      num_orbits=10,                 # Number of orbits
    #      num_sats_per_orbit=48,         # Number of satellites per orbit
    #      altitude=550)                  # Altitude of the orbits in km