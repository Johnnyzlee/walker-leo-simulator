#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: visualization.py
Author: Li ZENG @ HKUST ECE
License: MIT License

Description:
This module defines the Visualization class for visualizing satellite constellations.
It includes methods for computing satellite positions and links at a given timestamp and visualizing the network as a 3D graph.
"""

import plotly.graph_objects as go
import numpy as np
from network import SatNetStar, SatNetDelta


class Visualization:
    """
    Class for visualizing satellite constellations.
    """
    def __init__(self, network):
        """
        Initialize the Visualization class.

        Parameters:
        - network: SatNet object representing the satellite network.
        """
        self.network = network

    def compute_positions_and_links(self, timestamp):
        """
        Compute satellite positions and links at a given timestamp.

        Parameters:
        - timestamp: Time in seconds.

        Returns:
        - positions: Dictionary {orbit_id: {sat_id: position_ecef}} of satellite positions.
        - intra_links: List of tuples [(sat1, sat2, weight)] representing intra-orbit links and their weights.
        - inter_links: List of tuples [(sat1, sat2, weight)] representing inter-orbit links and their weights.
        """
        self.network.update_network(timestamp - self.network.time)
        positions = {}
        intra_links = []
        inter_links = []

        for orbit in self.network.constellation.orbits:
            positions[orbit.id] = {sat.id: sat.position_ecef for sat in orbit.sats}

        for edge in self.network.graph.edges(data=True):
            sat1, sat2, data = edge
            weight = data['weight']
            if sat1[0] == sat2[0]:  # Same orbit
                intra_links.append((sat1, sat2, weight))
            else:  # Different orbits
                inter_links.append((sat1, sat2, weight))

        return positions, intra_links, inter_links

    def visualize(self, timestamp, output_file='visualization/network_visualization.html'):
        """
        Visualize the satellite network at a given timestamp.

        Parameters:
        - timestamp: Time in seconds.
        - output_file: Output HTML file for the visualization.
        """
        positions, intra_links, inter_links = self.compute_positions_and_links(timestamp)

        fig = go.Figure()

        # Add satellite nodes
        for orbit_id, sats in positions.items():
            x, y, z = [], [], []
            for sat_id, pos in sats.items():
                x.append(pos[0])
                y.append(pos[1])
                z.append(pos[2])
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z,
                mode='markers',
                marker=dict(size=4),
                name=f'Orbit {orbit_id}'
            ))

        # Add intra-orbit links (edges) with weights
        for link in intra_links:
            sat1_pos = positions[link[0][0]][link[0][1]]
            sat2_pos = positions[link[1][0]][link[1][1]]
            weight = link[2]
            fig.add_trace(go.Scatter3d(
                x=[sat1_pos[0], sat2_pos[0]],
                y=[sat1_pos[1], sat2_pos[1]],
                z=[sat1_pos[2], sat2_pos[2]],
                mode='lines',
                line=dict(color='black', width=2, dash='solid'),
                name=f'Intra-Link {link[0]}-{link[1]} ({weight:.2f} km)'
            ))

        # Add inter-orbit links (edges) with weights
        for link in inter_links:
            sat1_pos = positions[link[0][0]][link[0][1]]
            sat2_pos = positions[link[1][0]][link[1][1]]
            weight = link[2]
            fig.add_trace(go.Scatter3d(
                x=[sat1_pos[0], sat2_pos[0]],
                y=[sat1_pos[1], sat2_pos[1]],
                z=[sat1_pos[2], sat2_pos[2]],
                mode='lines',
                line=dict(color='grey', width=2, dash='solid'),
                name=f'Inter-Link {link[0]}-{link[1]} ({weight:.2f} km)'
            ))

        # Add Earth sphere
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = 6371 * np.outer(np.cos(u), np.sin(v))
        y = 6371 * np.outer(np.sin(u), np.sin(v))
        z = 6371 * np.outer(np.ones(np.size(u)), np.cos(v))
        fig.add_trace(go.Surface(
            x=x, y=y, z=z,
            colorscale='Blues',
            opacity=0.1,
            showscale=False
        ))

        # Add equator
        equator_x = 6371 * np.cos(u)
        equator_y = 6371 * np.sin(u)
        equator_z = np.zeros_like(u)
        fig.add_trace(go.Scatter3d(
            x=equator_x, y=equator_y, z=equator_z,
            mode='lines',
            line=dict(color='blue', width=1),
            name='Equator'
        ))

        # Add prime meridian
        prime_meridian_x = 6371 * np.cos(u)
        prime_meridian_y = np.zeros_like(u)
        prime_meridian_z = 6371 * np.sin(u)
        fig.add_trace(go.Scatter3d(
            x=prime_meridian_x, y=prime_meridian_y, z=prime_meridian_z,
            mode='lines',
            line=dict(color='blue', width=1),
            name='Prime Meridian'
        ))

        # Add North and South Poles
        fig.add_trace(go.Scatter3d(
            x=[0], y=[0], z=[6371],
            mode='markers+text',
            marker=dict(size=5, color='red'),
            text=['North Pole'],
            textposition='top center',
            name='North Pole'
        ))
        fig.add_trace(go.Scatter3d(
            x=[0], y=[0], z=[-6371],
            mode='markers+text',
            marker=dict(size=5, color='red'),
            text=['South Pole'],
            textposition='bottom center',
            name='South Pole'
        ))

        fig.update_layout(
            scene=dict(
            xaxis_title='X (km)',
            yaxis_title='Y (km)',
            zaxis_title='Z (km)',
            xaxis_visible=False,
            yaxis_visible=False,
            zaxis_visible=False,
            aspectmode='data',
            bgcolor='rgba(0,0,0,0)'
            ),
            title=f'Satellite Network at t={timestamp} seconds',
            annotations=[
            dict(
                text=(
                f"Number of Orbits: {len(self.network.constellation.orbits)}<br>"
                f"Number of Satellites per Orbit: {len(self.network.constellation.orbits[0].sats)}<br>"
                f"Altitude: {self.network.constellation.orbits[0].radius - 6371} km<br>"
                f"Inclination: {np.degrees(self.network.constellation.orbits[0].inclination)} degrees"
                ),
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0,
                y=1,
                xanchor='left',
                yanchor='top',
                align="left",
                font=dict(
                size=12,
                color="black"
                ),
                bgcolor="white",
                bordercolor="black",
                borderwidth=1
            )
            ]
        )

        fig.write_html(output_file)
        print(f'Visualization saved to {output_file}')
        
if __name__ == "__main__":
    from constellation import DeltaConstellation
    from network import SatNetDelta
    # Example usage for Walker Delta Constellation
    constellation = DeltaConstellation(num_orbits=30, num_sats_per_orbit=48, radius=6371 + 550, inclination=53)
    network = SatNetDelta(constellation)
    viz = Visualization(network)
    # viz.visualize(timestamp=12)


    import datetime
    current_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    viz.visualize(timestamp=12, output_file=f'visualization/network_visualization_{current_time}.html')
    