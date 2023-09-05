import numpy as np

def find_elements_in_proximity(midpoints, radius):
    proximities = []
    for i, midpoint in enumerate(midpoints):
        close_elements = [j for j, m in enumerate(midpoints) if np.linalg.norm(np.array(midpoint) - np.array(m)) <= radius and j != i]
        proximities.append(close_elements)
    return proximities

def filter(values, proximities):
    smoothed = []
    for i, close_elements in enumerate(proximities):
        all_elements = [i] + close_elements
        average_value = np.mean([values[id] for id in all_elements])
        smoothed.append(average_value)
    return smoothed