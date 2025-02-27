#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: network.py
Author: Li ZENG @ HKUST ECE
License: MIT License

Description:
This module defines the satellite network classes for Walker Star and Walker Delta constellations.
It includes methods for building the network graph, updating satellite positions, and calculating shortest paths.
"""

from constellation import StarConstellation, DeltaConstellation
from abc import ABC, abstractmethod
import numpy as np
import networkx as nx

EARTH_RADIUS = 6371  # Define Earth's radius as a constant

class SatNet(ABC):
    """
    Abstract base class for satellite networks.
    Provides common methods and properties for satellite networks.
    """
    __LIGHT_SPEED = 299792.458  # Speed of light in km/s

    @property
    def LIGHT_SPEED(self):
        """Speed of light in km/s (read-only)"""
        return self.__LIGHT_SPEED
    
    def __init__(self, constellation):
        """
        Initialize the satellite network.

        Parameters:
        - constellation: WalkerConstellation object representing the satellite constellation.
        """
        self.constellation = constellation
        self.time = 0  # Initial time in seconds

    @abstractmethod
    def _build_graph(self):
        """Abstract method to build the network graph."""
        pass
    
    @abstractmethod
    def update_graph(self):
        """Abstract method to update the network graph."""
        pass
    
    def get_distance(self, vertex_key1, vertex_key2):
        """
        Calculate the distance between two satellites.

        Parameters:
        - vertex_key1: Tuple (orbit_id, sat_id) of the first satellite.
        - vertex_key2: Tuple (orbit_id, sat_id) of the second satellite.

        Returns:
        - Distance between the two satellites in km.
        """
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")
        pos1 = sat1.position_ecef
        pos2 = sat2.position_ecef
        return np.linalg.norm(pos1 - pos2)

    def update_network(self, delta_t):
        """
        Update the network by advancing the simulation time and updating satellite positions.

        Parameters:
        - delta_t: Time increment in seconds.
        """
        self.time += delta_t
        self.constellation.update_constellation(delta_t)
        self.update_graph()
        return
    
    def get_shortest_path(self, source, target, weight='weight'):
        """
        Find the shortest path between two satellites in the network.

        Parameters:
        - source: Tuple (orbit_id, sat_id) of the source satellite.
        - target: Tuple (orbit_id, sat_id) of the target satellite.
        - weight: Edge attribute to use as weight (default: 'weight').

        Returns:
        - path: List of vertex keys representing the shortest path.
        - latency: End-to-end propagation latency in seconds.
        """
        try:
            if source not in self.graph or target not in self.graph:
                raise ValueError("Source or target node not found in the network")

            # Find the shortest path using Dijkstra's algorithm
            path = nx.shortest_path(self.graph, source, target, weight='weight') 
            
            # Calculate total distance along the path
            total_distance = 0
            for i in range(len(path)-1):
                total_distance += self.graph[path[i]][path[i+1]]['weight']
            
            # Calculate latency (distance / speed of light)
            latency = total_distance / self.LIGHT_SPEED

            return path, latency

        except nx.NetworkXNoPath:
            raise ValueError(f"No path exists between satellites {source} and {target}")
        except Exception as e:
            raise Exception(f"Error finding shortest path: {str(e)}")
        
    def get_single_source_paths(self, source=(1, 1)):
        """
        Find shortest paths and latencies from a source satellite to all other satellites.

        Parameters:
        - source: Tuple (orbit_id, sat_id) of the source satellite.

        Returns:
        - paths: Dictionary {target: path} containing paths to all other satellites.
        - latencies: Dictionary {target: latency} containing latencies to all other satellites.
        """
        try:
            if source not in self.graph:
                raise ValueError("Source node not found in the network")

            # Get distances and paths to all nodes using single-source Dijkstra
            distances, paths = nx.single_source_dijkstra(
                self.graph, 
                source, 
                weight='weight'
            )

            # Convert distances to latencies
            latencies = {
                target: distance / self.LIGHT_SPEED 
                for target, distance in distances.items()
            }

            return paths, latencies

        except Exception as e:
            raise Exception(f"Error computing single-source paths: {str(e)}")
    

class SatNetStar(SatNet):
    """
    Satellite network class for Walker Star Constellations.
    """
    def __init__(self, constellation: StarConstellation, gamma_deg=80.0):
        """
        Initialize the Walker Star satellite network.

        Parameters:
        - constellation: StarConstellation object representing the satellite constellation.
        - gamma_deg: Maximum latitude that supports the inter-plane links (default: 80.0 degrees).
        """
        if constellation.type != "Walker Star Constellation":
            raise ValueError("SatNetStar requires a Walker Star Constellation")
        
        super().__init__(constellation)
        self.gamma_deg = gamma_deg  # Maximum latitude that supports the inter-plane links
        self.graph = self._build_graph()
        return
    
    def _build_graph(self):
        """
        Build the network graph for the Walker Star Constellation.

        Returns:
        - graph: NetworkX graph representing the satellite network.
        """
        graph = nx.Graph()
        
        # Add nodes for each satellite in the constellation
        for orbit in self.constellation.orbits:
            for sat in orbit.sats:
                # Add node with key as (orbit_id, satellite_id) and store the satellite object
                graph.add_node((orbit.id, sat.id), sat=sat)
        self.graph = graph
        
        
        # Add edges for laser inter-satellite links
        # Intra-plane LISL: add edges between neighboring satellites in the same orbit
        for orbit in self.constellation.orbits:
            if orbit.num_sats < 2:
                continue # No intra-plane links for single-satellite orbits
            for sat_id in range(1, orbit.num_sats + 1):
                sat1 = (orbit.id, sat_id)
                sat2 = (orbit.id, (sat_id % orbit.num_sats) + 1)
                if self._check_isl_feasibility(sat1, sat2):
                    distance = self.get_distance(sat1, sat2)
                    self.graph.add_edge(sat1, sat2, weight=distance)
                    
        # Inter-plane LISL: add edges between neighboring satellites in adjacent orbits
        if self.constellation.num_orbits < 2:
            return self.graph # No inter-plane links for single-orbit constellations
        for orbit_id in range(1, self.constellation.num_orbits + 1):
            for sat_id in range(1, self.constellation.num_sats_per_orbit + 1):
                sat1 = (orbit_id, sat_id)
                sat2 = ((orbit_id % self.constellation.num_orbits) + 1, sat_id)
                if self._check_isl_feasibility(sat1, sat2):
                    distance = self.get_distance(sat1, sat2)
                    self.graph.add_edge(sat1, sat2, weight=distance)

        return graph
    
    def _check_isl_feasibility(self, vertex_key1, vertex_key2):
        """
        Check if a laser inter-satellite link (LISL) is feasible between two satellites.

        Parameters:
        - vertex_key1: Tuple (orbit_id, sat_id) of the first satellite.
        - vertex_key2: Tuple (orbit_id, sat_id) of the second satellite.

        Returns:
        - True if the LISL is feasible, False otherwise.
        """
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")
        lat1 = sat1.position_geodetic[0]
        lat2 = sat2.position_geodetic[0]
        
        # Check if the LISL is blocked by the Earth
        mid_point_ecef = (sat1.position_ecef + sat2.position_ecef) / 2
        earth_radius = EARTH_RADIUS
        alt_mid_point = np.linalg.norm(mid_point_ecef) - earth_radius
        if alt_mid_point < 0:
            # The LISL is blocked by the Earth
            return False
        
        # Check if the LISL subjects to the rules of the constellation
        if vertex_key1[0] == vertex_key2[0]:
            # Intra-plane LISL, assert that the satellites' ID are neighboring
            if abs((vertex_key1[1] % self.constellation.num_sats_per_orbit) - (vertex_key2[1] % self.constellation.num_sats_per_orbit)) == 1:
                return True
        else:
            # Inter-plane LISL, assert that:
            #   - the satellites are in adjacent orbits
            #   - the satellites are in the same latitude band, i.e., abs(lat1 - lat2) <= 1
            #   - the satellites are not in the polar regions, i.e., abs(lat1 + lat2) / 2 <= gamma_deg
            if abs((vertex_key1[0] % self.constellation.num_orbits) - (vertex_key2[0] % self.constellation.num_orbits)) == 1:
                if abs(lat1 - lat2) <= 1 and abs(lat1 + lat2) / 2 <= self.gamma_deg:
                    # Check if the latitudes are close and in range [-gamma_deg, gamma_deg]
                    return True
                
        return False
    
    
class SatNetDelta(SatNet):
    """
    Satellite network class for Walker Delta Constellations.
    """
    def __init__(self, constellation: DeltaConstellation):
        """
        Initialize the Walker Delta satellite network.

        Parameters:
        - constellation: DeltaConstellation object representing the satellite constellation.
        """
        if constellation.type != "Walker Delta Constellation":
            raise ValueError("SatNetDelta requires a Walker Delta Constellation")
        
        super().__init__(constellation)
        self.graph = self._build_graph()
        return
    
    def update_graph(self):
        """Update the network graph for the Walker Delta Constellation."""
        self.graph = self._build_graph()
        return
        
    def _build_graph(self):
        """
        Build the network graph for the Walker Delta Constellation.

        Returns:
        - graph: NetworkX graph representing the satellite network.
        """
        graph = nx.Graph()
        # Add nodes for each satellite in the constellation
        for orbit in self.constellation.orbits:
            for sat in orbit.sats:
                # Add node with key as (orbit_id, satellite_id) and store the satellite object
                graph.add_node((orbit.id, sat.id), sat=sat)
        self.graph = graph
        
        # Add edges for laser inter-satellite links
        # Intra-plane LISL: add edges between neighboring satellites in the same orbit
        for orbit in self.constellation.orbits:
            if orbit.num_sats < 2:
                continue # No intra-plane links for single-satellite orbits
            for sat_id in range(1, orbit.num_sats + 1):
                sat1 = (orbit.id, sat_id)
                sat2 = (orbit.id, (sat_id % orbit.num_sats) + 1)
                if self._check_isl_feasibility(sat1, sat2):
                    distance = self.get_distance(sat1, sat2)
                    self.graph.add_edge(sat1, sat2, weight=distance)
                    
        # Inter-plane LISL: add edges between neighboring satellites in adjacent orbits
        if self.constellation.num_orbits < 2:
            return self.graph # No inter-plane links for single-orbit constellations
        for orbit_id in range(1, self.constellation.num_orbits + 1):
            if orbit_id != self.constellation.num_orbits: 
                # For all orbits except the last one
                for sat_id in range(1, self.constellation.num_sats_per_orbit + 1):
                    sat1 = (orbit_id, sat_id)
                    sat2 = ((orbit_id % self.constellation.num_orbits) + 1, sat_id)
                    if self._check_isl_feasibility(sat1, sat2):
                        distance = self.get_distance(sat1, sat2)
                        self.graph.add_edge(sat1, sat2, weight=distance)
            else: 
                # For the last orbit, add links to the first orbit, the offset is calculated based on the phasediff
                for sat_id in range(1, self.constellation.num_sats_per_orbit + 1):
                    sat_id_offset = (self.constellation.num_orbits * np.degrees(self.constellation.phasediff)) / (360 / self.constellation.num_sats_per_orbit)
                    sat_id_offset = round(sat_id_offset) # Round to the nearest integer
                    sat1 = (orbit_id, (sat_id + sat_id_offset - 1) % self.constellation.num_sats_per_orbit + 1)
                    sat2 = (1, sat_id)
                    if self._check_isl_feasibility(sat1, sat2):
                        distance = self.get_distance(sat1, sat2)
                        self.graph.add_edge(sat1, sat2, weight=distance)

        return self.graph
    
    def _check_isl_feasibility(self, vertex_key1, vertex_key2):
        """
        Check if a laser inter-satellite link (LISL) is feasible between two satellites.

        Parameters:
        - vertex_key1: Tuple (orbit_id, sat_id) of the first satellite.
        - vertex_key2: Tuple (orbit_id, sat_id) of the second satellite.

        Returns:
        - True if the LISL is feasible, False otherwise.
        """
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")

        mid_point_ecef = (sat1.position_ecef + sat2.position_ecef) / 2
        earth_radius = EARTH_RADIUS
        alt_mid_point = np.linalg.norm(mid_point_ecef) - earth_radius
        if alt_mid_point < 0:
            # The LISL is blocked by the Earth
            return False
        
        return True
