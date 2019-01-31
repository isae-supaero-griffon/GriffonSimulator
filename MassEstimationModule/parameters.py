# from MassEstimationModule.part import *

# avionics params
electronics_mass = 1.

# combustion params
chamber_radius = 0.1
chamber_height = 0.6
chamber_pressure = 36e5
fuel_mass = 0.75
liner_thickness = 0.005
injector_mass = 0.5
nozzle_mass = 1.

# oxidiser params
oxidant_radius = 0.1
oxidant_height = 0.3
oxidant_pressure = 60.e5
oxidant_mass = 4.

# pressuriser params
pressuriser_valve_mass = 0.1
pressuriser_tank_mass = 1.
pressuriser_tank_radius = 0.1
pressuriser_tank_height = 0.2
he_gamma = 1.66
he_gas_constant = 2077.
room_temp = 300.

# plumbing params
pipes_mass = 0.1
valve1_mass = 0.4
valve2_mass = 0.4
valve3_mass = 0.6

# structural params
fixing_mass = 0.5
rocket_structure_mass = 20.

# material params
al_density = 2700.
al_yield = 290.e6
liner_density = 1070.

safety_factor = 3.

# dictionary structure
system_dict = {
    "identifier": "hybrid engine",
    "materials": {
        "aluminium": {
            "identifier": "aluminium",
            "density": al_density,
            "yield_stress": al_yield
        },
        "liner_mat": {
            "identifier": "liner_mat",
            "density": liner_density
        }
    },
    "subsystems": {
        "avionics": {
            "identifier": "avionics",
            "parts": {
                "electronics": {
                    "identifier": "electronics",
                    "type": 'StaticPart',
                    "mass": electronics_mass
                }
            }
        },
        "combustion": {
            "identifier": "combustion",
            "parts": {
                "chamber": {
                    "identifier": "chamber",
                    "type": 'CombustionChamber',
                    "radius": chamber_radius,
                    "height": chamber_height,
                    "pressure": chamber_pressure,
                    "propellant_mass": fuel_mass,
                    "material_id": "aluminium"
                },
                "liner": {
                    "identifier": "liner",
                    "type": 'ThermalLiner',
                    "chamber_ref": {
                        "subsystem_id": "combustion",
                        "part_id": "chamber"
                    },
                    "thickness": liner_thickness,
                    "material_id": "liner_mat"
                },
                "injector": {
                    "identifier": "injector",
                    "type": 'StaticPart',
                    "mass": injector_mass
                },
                "nozzle": {
                    "identifier": "nozzle",
                    "type": 'StaticPart',
                    "mass": nozzle_mass
                }
            }
        },
        "oxidiser": {
            "identifier": "oxidiser",
            "parts": {
                "oxidant_tank": {
                    "identifier": "oxidant_tank",
                    "type": 'OxidantTank',
                    "radius": oxidant_radius,
                    "height": oxidant_height,
                    "pressure": oxidant_pressure,
                    "propellant_mass": oxidant_mass,
                    "material_id": "aluminium"
                }
            }
        },
        "pressuriser": {
            "identifier": "pressuriser",
            "parts": {
                "pressuriser_tank": {
                    "identifier": "pressuriser_tank",
                    "type": 'PressuriserTank',
                    "radius": pressuriser_tank_radius,
                    "height": pressuriser_tank_height,
                    "structure_mass": pressuriser_tank_mass,
                    "oxidiser_ref": {
                        "subsystem_id": "oxidiser",
                        "part_id": "oxidant_tank"
                    },
                    "pressuriser_gas_constant": he_gas_constant,
                    "pressuriser_gamma": he_gamma,
                    "ambient_temperature": room_temp
                },
                "pressuriser_valve": {
                    "identifier": "pressuriser_valve",
                    "type": 'StaticPart',
                    "mass": pressuriser_valve_mass
                }
            }
        },
        "plumbing": {
            "identifier": "plumbing",
            "parts": {
                "pipes": {
                    "identifier": "pipes",
                    "type": 'StaticPart',
                    "mass": pipes_mass
                },
                "valve1": {
                    "identifier": "valve1",
                    "type": 'StaticPart',
                    "mass": valve1_mass
                },
                "valve2": {
                    "identifier": "valve2",
                    "type": 'StaticPart',
                    "mass": valve2_mass
                },
                "valve3": {
                    "identifier": "valve3",
                    "type": 'StaticPart',
                    "mass": valve3_mass
                },
            }
        },
        "structures": {
            "identifier": "structure",
            "parts": {
                "rocket_structure": {
                    "identifier": "rocket_structure",
                    "type": 'StaticPart',
                    "mass": rocket_structure_mass
                },
                "fixings": {
                    "identifier": "fixings",
                    "type": 'StaticPart',
                    "mass": fixing_mass
                }
            }
        },
    },
    "safety_factor": safety_factor
}
