# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division

from datetime import timedelta
from math import sqrt, sin, cos, atan2, degrees, radians

import ephem
import numpy as np
from collections import namedtuple
from haversine import haversine
from sgp4.ext import jday as days

from tle import TLE

ECEF = namedtuple('ECEF', ('x', 'y', 'z'))
LLA = namedtuple('LLA', ('lat', 'lon', 'alt'))

_tle = TLE()


def sun_angle(dot, t):
    obs = ephem.Observer()
    obs.lat = str(dot[0])
    obs.long = str(dot[1])
    obs.date = t
    sun = ephem.Sun(obs)
    sun.compute(obs)
    return degrees(sun.alt)


def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def ecef2lla(x, y, z):
    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622e-2

    # calculations:
    b = sqrt(a ** 2 * (1 - e ** 2))
    ep = sqrt((a ** 2 - b ** 2) / b ** 2)
    p = sqrt(x ** 2 + y ** 2)
    th = atan2(a * z, b * p)
    lon = atan2(y, x)
    lat = atan2((z + ep ** 2 * b * sin(th) ** 3), (p - e ** 2 * a * cos(th) ** 3))
    N = a / sqrt(1 - e ** 2 * sin(lat) ** 2)
    alt = p / cos(lat) - N
    return LLA(degrees(lat), degrees(lon), alt)


def lla2ecef(lat, lon, alt, rad=False):
    """
    """
    if not rad:
        lat, lon = radians(lat), radians(lon)

    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622e-2

    # intermediate calculation
    # (prime vertical radius of curvature)
    N = a / sqrt(1 - e ** 2 * sin(lat) ** 2)

    # results:
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = ((1 - e ** 2) * N + alt) * sin(lat)

    return ECEF(x, y, z)


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

    angle_x = radians(angle_x)
    angle_y = radians(angle_y)
    angle_z = radians(angle_z)
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
    iss_orientation_matrix = np.dot(np.dot(np.dot(np.dot(Mx, My), My), Mz), additional_rotation_matrix)
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

    a = 6378137.0  # semi-major axis of the WGS-84 ellipsoid in meters
    b = a  # semi-minor axis in meters
    c = 6356752.3142
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


# Считает Гринвичское среднее звездное время, то есть угол поворота вращающейся СК относительно инерциальной
def GMST(year, month, day, hours, minutes, seconds):
    d = days(year, month, day, hours, minutes, seconds)
    GMST = 280.46061837 + 360.98564736629 * d
    GMST = GMST - 360 * int(GMST / 360)
    if GMST < 0:
        GMST = 360.0 + GMST
    GMST = radians(GMST)
    cz = cos(GMST)
    sz = sin(GMST)
    # Поворот относительно оси z - из центра Земли в сев. полюс
    Rz = np.matrix([
        [cz, -sz, 0],
        [sz, cz, 0],
        [0, 0, 1]
    ])
    return GMST, Rz


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
    if hit_1:
        coords = ecef2lla(*hit_1)
        return coords.lat, coords.lon


def hit_the_earth(satellite, t):
    position, velocity = satellite.propagate(t.year, t.month, t.day, t.hour, t.minute, t.second)
    return geolocate(position, velocity, 0, 0, 0, t)


if __name__ == '__main__':
    print(ecef2lla(*lla2ecef(45, -200, 1000)))
