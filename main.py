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

def main(constellation_type, duration=3600, delta_t=1, *args, **kwargs):
    """
    Initialize and run the constellation network.

    Parameters:
    - constellation_type (str): Type of the constellation ("Walker Star Constellation" or "Walker Delta Constellation").
    - duration (int): Duration of the simulation in seconds.
    - delta_t (int): Time step for the simulation in seconds.
    - *args: Additional positional arguments for the constellation.
    - **kwargs: Additional keyword arguments for the constellation.

    Required keyword arguments for "Walker Star Constellation":
    - num_orbits (int): Number of orbits.
    - num_sats_per_orbit (int): Number of satellites per orbit.
    - radius (float): Radius of the orbits.

    Required keyword arguments for "Walker Delta Constellation":
    - num_orbits (int): Number of orbits.
    - num_sats_per_orbit (int): Number of satellites per orbit.
    - radius (float): Radius of the orbits.
    - inclination (float): Inclination of the orbits.
    """
    # Required parameters for each constellation type
    required_params_star = ['num_orbits', 'num_sats_per_orbit', 'radius']
    required_params_delta = ['num_orbits', 'num_sats_per_orbit', 'radius', 'inclination']

    # Initialize the constellation and network based on the constellation type
    if constellation_type == "Walker Star Constellation":
        for param in required_params_star:
            if param not in kwargs:
                raise ValueError(f"Missing required parameter for Star Constellation: {param}")
        constellation = StarConstellation(*args, **kwargs)
        network = SatNetStar(constellation)
    elif constellation_type == "Walker Delta Constellation":
        for param in required_params_delta:
            if param not in kwargs:
                raise ValueError(f"Missing required parameter for Delta Constellation: {param}")
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
    main("Walker Delta Constellation", # Type of the constellation
         3600 * 1, # Duration of the simulation in seconds
         1, # Time step for the simulation in seconds 
         num_orbits=36, num_sats_per_orbit=48, radius=6371 + 550, inclination=53)