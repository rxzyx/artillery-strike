import math
from json import dumps
from time import time
from math import radians, sin, cos, sqrt, atan2

# Made by rxzyx on GitHub!
# A simple (very basic) artillery strike program!

def artillery_strike(artillery_pos, target_pos, wind_speed, wind_direction):   
    artillery_pos_cartesian = lat_lon_to_cartesian(artillery_pos)
    target_pos_cartesian = lat_lon_to_cartesian(target_pos)
    
    d = vincenty_distance(artillery_pos, target_pos)
    
    azimuth = vincenty_initial_bearing(artillery_pos, target_pos)
    
    wind_direction_rad = math.radians(wind_direction)
    
    crosswind = wind_speed * math.sin(wind_direction_rad - azimuth)
    headwind = wind_speed * math.cos(wind_direction_rad - azimuth)
    
    azimuth_adjusted = azimuth + math.atan2(crosswind, d)
    
    if d > range_threshold:
        raise ValueError("Target out of range")
    
    denominator = (d * math.tan(azimuth_adjusted) - 0.5 * 9.81 * d**2)
    numerator = -(d**2 + d * math.tan(azimuth_adjusted)**2)
    
    elevation = math.asin(numerator / denominator)

    if elevation < -math.pi/2 or elevation > math.pi/2:
        raise ValueError("Invalid elevation angle")
    
    range_adjusted = d - headwind * (d / math.sqrt((d * math.tan(azimuth_adjusted))**2 + d**2))
    
    initial_velocity = math.sqrt((range_adjusted * 9.81) / (math.sin(2 * elevation)))
    
    time_of_flight = range_adjusted / (initial_velocity * math.cos(elevation))

    artillery_latitude, artillery_longitude = cartesian_to_lat_lon(artillery_pos_cartesian)

    print((target_pos_cartesian[2] - artillery_pos_cartesian[2]) / math.cos(elevation))
    
    return {
        'distance': d,
        'azimuth': math.degrees(azimuth),
        'crosswind': crosswind,
        'headwind': headwind,
        'azimuth_adjusted': math.degrees(azimuth_adjusted),
        'elevation': math.degrees(elevation),
        'range_adjusted': range_adjusted,
        'initial_velocity': initial_velocity,
        'time_of_flight': time_of_flight,
        'artillery_latitude': artillery_latitude,
        'artillery_longitude': artillery_longitude
    }

def vincenty_distance(point_a, point_b):
    # Do not be confused, it's basically the haversine formula
    # but accounted for an ellipsoid (the earth), for maximum accuracy.
    a = 6378137.0
    f = 1 / 298.257223563       # WGS-84 model

    lat1, lon1 = radians(point_a[0]), radians(point_a[1])
    lat2, lon2 = radians(point_b[0]), radians(point_b[1])

    delta_lon = lon2 - lon1

    U1 = atan2((1 - f) * sin(lat1), cos(lat1))
    U2 = atan2((1 - f) * sin(lat2), cos(lat2))
    sinU1 = sin(U1)
    cosU1 = cos(U1)
    sinU2 = sin(U2)
    cosU2 = cos(U2)

    lambda_lon = delta_lon
    lambda_prev = 2 * 3.14159265358979323846

    while abs(lambda_lon - lambda_prev) > 1e-12:
        sin_lambda = sin(lambda_lon)
        cos_lambda = cos(lambda_lon)
        sin_sigma = sqrt((cosU2 * sin_lambda) ** 2 + (cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda) ** 2)
        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_lambda
        sigma = atan2(sin_sigma, cos_sigma)

        sin_alpha = cosU1 * cosU2 * sin_lambda / sin_sigma
        cos_sq_alpha = 1 - sin_alpha ** 2

        cos_2sigma_m = cos_sigma - 2 * sinU1 * sinU2 / cos_sq_alpha
        C = f / 16 * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))
        lambda_prev = lambda_lon
        lambda_lon = delta_lon + (1 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (cos_2sigma_m + C * cos_sigma * (-1 + 2 * cos_2sigma_m ** 2))
        )

    u_sq = cos_sq_alpha * (a ** 2 - f * (2 - f) * a ** 2) / (4 * a ** 2)
    A = 1 + u_sq / 16384 * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024 * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    delta_sigma = B * sin_sigma * (
        cos_2sigma_m + B / 4 * (
            cos_sigma * (-1 + 2 * cos_2sigma_m ** 2) - B / 6 * cos_2sigma_m * (
                -3 + 4 * sin_sigma ** 2
            ) * (-3 + 4 * cos_2sigma_m ** 2)
        )
    )

    distance = a * A * (sigma - delta_sigma)
    return distance


def vincenty_initial_bearing(coord1, coord2):
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    delta_lon = math.radians(lon2 - lon1)
    cos_lat2 = math.cos(math.radians(lat2))

    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - \
        math.sin(math.radians(lat1)) * cos_lat2 * math.cos(delta_lon)

    y = math.sin(delta_lon) * cos_lat2

    initial_bearing = math.atan2(y, x)

    return initial_bearing


def lat_lon_to_cartesian(coords):
    lat, lon = coords
    R = 6371000
    x = R * math.cos(math.radians(lat)) * math.cos(math.radians(lon))
    y = R * math.cos(math.radians(lat)) * math.sin(math.radians(lon))
    z = R * math.sin(math.radians(lat))
    return x, y, z


def cartesian_to_lat_lon(coords):
    x, y, z = coords
    R = 6371000
    lat = math.degrees(math.asin(z / R))
    lon = math.degrees(math.atan2(y, x))
    return lat, lon


artillery_pos = (46.63359481660844, 7.579840953049487) # The position of the artillery
target_pos = (46.89792522703249, 7.423970599509989) # The position of the target.
wind_speed = 10                 # Wind speed in meters per second
wind_direction = 45             # Wind direction in degrees
range_threshold = 100000        # Maximum range in meters

start_time = time()
result = artillery_strike(artillery_pos, target_pos, wind_speed, wind_direction)
end_time = time()

time_took = format(end_time - start_time, '.5f')
print(f"Program took around {time_took} seconds.")

print(
    dumps(
        result,
        indent=4,
        sort_keys=True
    )
)
