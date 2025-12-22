import pandas as pd

__all__ = ["joules2kilotons"]


def joules2kilotons(joules: pd.Series):
    if (joules < 0).any():
        raise ValueError("Energy cannot be negative")

    JOULES_PER_KILOTON = 4.184e12

    return joules / JOULES_PER_KILOTON


def kinetic_energy(mass: pd.Series, velocity: pd.Series):
    return 0.5 * mass * velocity**2
