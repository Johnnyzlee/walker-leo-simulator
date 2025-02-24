# Walker-LEO Satellite Constellation Simulator

A simple Python-based simulation framework for modeling and analyzing Walker Star and Walker Delta based Low Earth Orbit (LEO) satellite constellation networks.

## Features

- Simulates both Walker Star and Walker Delta constellation patterns
- Models satellite network connectivity and dynamics
- Configurable simulation parameters:
  - Number of orbits
  - Satellites per orbit
  - Orbital altitude
  - Inclination (for Delta constellations)
  - Simulation duration
  - Time step resolution
- Real-time simulation progress tracking

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Johnnyzlee/walker-leo-simulator.git
cd walker-leo-simulator
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the simulation with default parameters:
```python
python main.py
```

Customize simulation parameters:
```bash
# Walker Star Constellation
python main.py "Walker Star Constellation" duration=3600 delta_t=1 num_orbits=36 num_sats_per_orbit=48 altitude=550

# Walker Delta Constellation
python main.py "Walker Delta Constellation" duration=3600 delta_t=1 num_orbits=36 num_sats_per_orbit=48 altitude=550 inclination=53.0
```

## Configuration

### Required Parameters

- `num_orbits`: Number of orbital planes
- `num_sats_per_orbit`: Number of satellites per orbital plane
- `radius`: Orbital radius (Earth radius + altitude in km)
- `altitude`: Altitude of the orbits.
- `inclination`: Orbital inclination (only for Delta constellations)

### Optional Parameters

- `duration`: Simulation duration in seconds (default: 3600)
- `delta_t`: Time step in seconds (default: 1)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the project
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a pull request

Please make sure to update tests as appropriate.
