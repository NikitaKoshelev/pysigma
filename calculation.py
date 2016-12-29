# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division

from datetime import timedelta
from math import sqrt, sin, cos, pi, atan2

import ephem
import numpy as np
from haversine import haversine
from sgp4.ext import jday as days


from tle import TLE

_tle = TLE()


def sun_angle(dot, t):
    obs = ephem.Observer()
    obs.lat = str(dot[0])
    obs.long = str(dot[1])
    obs.date = t
    sun = ephem.Sun(obs)
    sun.compute(obs)
    sun_angle = float(sun.alt) * 57.2957795  # Convert Radians to degrees
    return sun_angle


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


# Задает и поворачивает единичный вектор
def origin_direction(position, velocity, angle_x, angle_y, angle_z, GMST_rotation_matrix, additional_rotation=0):
    position = np.dot(position, GMST_rotation_matrix)
    velocity = np.dot(velocity, GMST_rotation_matrix)
    direction = -position
    Xiss = velocity
    Yiss = direction
    Ziss = np.matrix(((Xiss[0, 1] * Yiss[0, 2] - Xiss[0, 2] * Yiss[0, 1],
                       Xiss[0, 2] * Yiss[0, 0] - Xiss[0, 0] * Yiss[0, 2],
                       Xiss[0, 0] * Yiss[0, 1] - Xiss[0, 1] * Yiss[0, 0])))

    Xiss = normalize(Xiss)
    Yiss = normalize(Yiss)
    Ziss = normalize(Ziss)
    shift_matrix = np.matrix([[Xiss[0, 0], Yiss[0, 0], Ziss[0, 0]],
                              [Xiss[0, 1], Yiss[0, 1], Ziss[0, 1]],
                              [Xiss[0, 2], Yiss[0, 2], Ziss[0, 2]]])

    angle_x = angle_x * pi / 180.0  # rotation around X axis
    angle_y = angle_y * pi / 180.0  # rotation around Y axis
    angle_z = angle_z * pi / 180.0  # rotation around Z axis
    ca = cos(angle_x)
    sa = sin(angle_x)
    cb = cos(angle_y)
    sb = sin(angle_y)
    cc = cos(angle_z)
    sc = sin(angle_z)
    additional_rotation_matrix = np.matrix([[1, 0, 0],
                                            [0, cos(additional_rotation), -sin(additional_rotation)],
                                            [0, sin(additional_rotation), cos(additional_rotation)]])

    Mx = np.matrix([[1, 0, 0],
                    [0, ca, -sa],
                    [0, sa, ca]])

    My = np.matrix([[cb, 0, sb],
                    [0, 1, 0],
                    [-sb, 0, cb]])
    Mz = np.matrix([[cc, -sc, 0],
                    [+sc, cc, 0],
                    [0, 0, 1]])
    iss_orientation_matrix = np.dot(Mx, My)
    iss_orientation_matrix = np.dot(iss_orientation_matrix, My)
    iss_orientation_matrix = np.dot(iss_orientation_matrix, Mz)
    iss_orientation_matrix = np.dot(iss_orientation_matrix, additional_rotation_matrix)
    # Перейдем в iss
    direction = np.dot(direction, shift_matrix)
    # Повернем вектор направления
    direction = np.dot(direction, iss_orientation_matrix)
    # Перейдем обратно в ecef
    direction = np.dot(direction, shift_matrix.I)
    return position, direction


def calc_intersection(x, y, z, u, v, w):
    # u,v,w - это координаты ecef
    # define WGS-84
    a = 6378137.0 / 1000.0  # semi-major axis of the WGS-84 ellipsoid in meters
    b = a  # semi-minor axis in meters
    c = 6356752.3142 / 1000.0
    sq_eq_a = u ** 2 * b ** 2 * c ** 2 + v ** 2 * a ** 2 * c ** 2 + w ** 2 * a ** 2 * c ** 2
    sq_eq_b = 2 * (x * u * b ** 2 * a ** 2 + y * v * a ** 2 * c ** 2 + z * w * a ** 2 * b ** 2)
    sq_eq_c = (x ** 2 * b ** 2 * c ** 2
               + y ** 2 * a ** 2 * c ** 2
               + z ** 2 * a ** 2 * b ** 2
               - a ** 2 * b ** 2 * c ** 2)
    sq_eq_D = sq_eq_b ** 2 - 4 * sq_eq_a * sq_eq_c
    if sq_eq_D > 0:
        t1 = (-sq_eq_b + sqrt(sq_eq_D)) / (2 * sq_eq_a)
        t2 = (-sq_eq_b - sqrt(sq_eq_D)) / (2 * sq_eq_a)
        hit_1 = x + t1 * u, y + t1 * v, z + t1 * w
        hit_2 = x + t2 * u, y + t2 * v, z + t2 * w
        d1 = distance((x, y, z), hit_1)
        d2 = distance((x, y, z), hit_2)
        if d1 < d2:
            hit = hit_1
        else:
            hit = hit_2
    else:
        hit = None
    return hit


# Считает JD
def days(year, month, day, hours, minutes, seconds):
    dwhole = 367 * year - 7 * (year + (month + 9) // 12) // 4 + 275 * month // 9 + day - 730531.5
    dfrac = (hours + minutes / 60 + seconds / 3600) / 24
    d = dwhole + dfrac
    return d


# Считает Гринвичское среднее звездное время, то есть угол поворота вращающейся СК относительно инерциальной
def GMST(year, month, day, hours, minutes, seconds):
    d = days(year, month, day, hours, minutes, seconds)
    GMST = 280.46061837 + 360.98564736629 * d
    GMST = GMST - 360 * int(GMST / 360)
    if GMST < 0:
        GMST = 360.0 + GMST
    GMST = GMST * pi / 180
    cz = cos(GMST)
    sz = sin(GMST)
    # Поворот относительно оси z - из центра Земли в сев. полюс
    Rz = np.matrix([
        [cz, -sz, 0],
        [sz, cz, 0],
        [0, 0, 1]
    ])
    return GMST, Rz


def convert_ll_to_xyz(geo_coords):
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    (lat, lon) = geo_coords
    alt = 0
    rad = np.float64(6378137.0)  # Radius of the Earth (in meters)
    f = np.float64(1.0 / 298.257223563)  # Flattening factor WGS84 Model
    cosLat = np.cos(lat)
    sinLat = np.sin(lat)
    FF = (1.0 - f) ** 2
    C = 1 / np.sqrt(cosLat ** 2 + FF * sinLat ** 2)
    S = C * FF
    x = (rad * C + alt) * cosLat * np.cos(lon)
    y = (rad * C + alt) * cosLat * np.sin(lon)
    z = (rad * S + alt) * sinLat
    return x, y, z


def convert_xyz_ll(coords):
    # WGS-95 ellipsoid constants
    (x, y, z) = coords
    a = 6378137.0 / 1000.0
    e = 8.1819190842622e-2
    e = 0.0818191909289062

    # calculations
    b = sqrt(a ** 2 * (1 - e ** 2))
    ep = sqrt((a ** 2 - b ** 2) / b ** 2)
    p = sqrt(x ** 2 + y ** 2)
    th = atan2(a * z, b * p)
    lon = atan2(y, x)
    lat = atan2((z + ep ** 2 * b * sin(th) ** 3), (p - e ** 2 * a * cos(th) ** 3))
    N = a / sqrt(1 - e ** 2 * sin(lat) ** 2)
    alt = p / cos(lat) - N
    return lat * 180.0 / pi, lon * 180.0 / pi


def distance(dot_1, dot_2):
    return sqrt((dot_2[0] - dot_1[0]) ** 2 + (dot_2[1] - dot_1[1]) ** 2 + (dot_2[2] - dot_1[2]) ** 2)


def am_i_approaching(obj, dt, tle=None):
    if not isinstance(tle, TLE):
        tle = _tle
    iss = tle.get_iss(dt=dt)
    coords_now = hit_the_earth(iss, dt)
    coors_next = hit_the_earth(iss, dt + timedelta(seconds=1))
    return haversine(obj, coords_now) > haversine(obj, coors_next)


def find_traverz(start_datetime, end_datetime, obj):
    if end_datetime - start_datetime < timedelta(seconds=2):
        return start_datetime

    answer = None
    half = start_datetime + timedelta(seconds=(end_datetime - start_datetime).seconds / 2)
    if am_i_approaching(obj, start_datetime) and not (am_i_approaching(obj, half)):
        answer = find_traverz(start_datetime, half, obj)
    elif am_i_approaching(obj, start_datetime) and am_i_approaching(obj, half):
        answer = find_traverz(half, end_datetime, obj)
    elif not am_i_approaching(obj, start_datetime) and am_i_approaching(obj, half):
        answer = find_traverz(half, end_datetime, obj)
    return answer


def mass_center_in_time(timeline, tle):
    return [hit_the_earth(tle.get_iss(dt=dt), dt) for dt in timeline]


def geolocate(position, velocity, alfa, beta, gamma, t):
    GMST_rotation, GMST_rotation_matrix = GMST(t.year, t.month, t.day, t.hour, t.minute, t.second)
    origin, direction = origin_direction(position, velocity, alfa, beta, gamma, GMST_rotation_matrix, 0)
    hit_1 = calc_intersection(origin[0, 0], origin[0, 1], origin[0, 2],
                              direction[0, 0], direction[0, 1], direction[0, 2])
    return convert_xyz_ll(hit_1)


def hit_the_earth(satellite, t):
    position, velocity = satellite.propagate(t.year, t.month, t.day, t.hour, t.minute, t.second)
    return geolocate(position, velocity, 0, 0, 0, t)

