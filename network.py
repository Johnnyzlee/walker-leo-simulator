#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: network.py
Author: Li ZENG @ HKUST ECE
License: MIT License
"""

from constellation import StarConstellation, DeltaConstellation
from abc import ABC, abstractmethod
import numpy as np
import networkx as nx

class SatNet(ABC):
    def __init__(self, constellation):
        self.constellation = constellation
        self.time = 0  # initial time in seconds

    @abstractmethod
    def _build_graph(self):
        pass
    
    @abstractmethod
    def update_graph(self):
        pass
    
    def get_distance(self, vertex_key1, vertex_key2):
        # Note: distance of any satellite pair can be calculated, but it does not mean that there are direct links between them
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")
        pos1 = sat1.position_ecef
        pos2 = sat2.position_ecef
        return np.linalg.norm(pos1 - pos2)

    def update_network(self, delta_t):
        self.time += delta_t
        self.constellation.update_constellation(delta_t)
        self.update_graph()
        return


class SatNetStar(SatNet):
    def __init__(self, constellation: StarConstellation, gamma_deg=80.0):
        if constellation.type != "Walker Star Constellation":
            raise ValueError("SatNetStar requires a Walker Star Constellation")
        
        super().__init__(constellation)
        self.gamma_deg = gamma_deg # maximum latitude that supports the inter-plane links
        self.graph = self._build_graph()
        return
    
    def _build_graph(self):
        graph = nx.Graph()
        
        # Add nodes for each satellite in the constellation
        for orbit in self.constellation.orbits:
            for sat in orbit.sats:
                # Add node with key as (orbit_id, satellite_id) and store the satellite object
                graph.add_node((orbit.id, sat.id), sat=sat)
        self.graph = graph
        
        # Add edges for laser inter-satellite links
        # Intra-plance LISL: add edges between neighboring satellites in the same orbit
        for orbit in self.constellation.orbits:
            for sat_id in range(1, orbit.num_sats + 1):
                if self._check_isl_feasibility((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1)):
                    distance = self.get_distance((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1))
                    graph.add_edge((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1), weight=distance)
                    
        # Inter-plane LISL: add edges between neighboring satellites in adjacent orbits
        for orbit_id in range(1, self.constellation.num_orbits + 1):
            for sat_id in range(1, self.constellation.num_sats_per_orbit + 1):
                if self._check_isl_feasibility((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id)):
                    distance = self.get_distance((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id))
                    graph.add_edge((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id), weight=distance)
        return graph
    
    def _check_isl_feasibility(self, vertex_key1, vertex_key2):
        # check if LISL are in polar regions and connects satellites
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")
        lat1 = sat1.position_geodetic[0]
        lat2 = sat2.position_geodetic[0]
        
        # Check if the LISL is blocked by the Earth
        mid_point_ecef = (sat1.position_ecef + sat2.position_ecef) / 2
        earth_radius = 6371
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
    def __init__(self, constellation: DeltaConstellation):
        if constellation.type != "Walker Delta Constellation":
            raise ValueError("SatNetDelta requires a Walker Delta Constellation")
        
        super().__init__(constellation)
        self.graph = self._build_graph()
        return
    
    def update_graph(self):
        self.graph = self._build_graph()
        return
        
    def _build_graph(self):
        graph = nx.Graph()
        # Add nodes for each satellite in the constellation
        for orbit in self.constellation.orbits:
            for sat in orbit.sats:
                # Add node with key as (orbit_id, satellite_id) and store the satellite object
                graph.add_node((orbit.id, sat.id), sat=sat)
        self.graph = graph
        
        # Add edges for laser inter-satellite links
        # Intra-plance LISL: add edges between neighboring satellites in the same orbit
        for orbit in self.constellation.orbits:
            for sat_id in range(1, orbit.num_sats + 1):
                if self._check_isl_feasibility((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1)):
                    distance = self.get_distance((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1))
                    graph.add_edge((orbit.id, sat_id), (orbit.id, (sat_id % orbit.num_sats) + 1), weight=distance)
                    
        # Inter-plane LISL: add edges between neighboring satellites in adjacent orbits
        for orbit_id in range(1, self.constellation.num_orbits + 1):
            for sat_id in range(1, self.constellation.num_sats_per_orbit + 1):
                if self._check_isl_feasibility((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id)):
                    distance = self.get_distance((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id))
                    graph.add_edge((orbit_id, sat_id), ((orbit_id % self.constellation.num_orbits) + 1, sat_id), weight=distance)
                
        return graph
    
    def _check_isl_feasibility(self, vertex_key1, vertex_key2):
        # Check if the LISL is blocked by the Earth
        try:
            sat1 = self.graph.nodes[vertex_key1]['sat']
            sat2 = self.graph.nodes[vertex_key2]['sat']
        except KeyError:
            raise ValueError("One or both vertex keys do not exist in the network")
        
        mid_point_ecef = (sat1.position_ecef + sat2.position_ecef) / 2
        earth_radius = 6371
        alt_mid_point = np.linalg.norm(mid_point_ecef) - earth_radius
        if alt_mid_point < 0:
            # The LISL is blocked by the Earth
            return False
        return True