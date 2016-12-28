# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division
import requests
import forecastio
import ephem
import logging
import numpy as np
import sys
import time

from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from haversine import haversine
from datetime import datetime, timedelta, date
from math import sqrt, sin, cos, pi, atan2

if sys.version_info.major < 3:
    try:
        from codecs import open
    except:
        pass

log = logging.getLogger(__name__)

program_start_time = time.time()

# SETUP_START
sun_limit = 20
vision_limit = 220
trace_length_in_seconds = 360
cloudcover_limit = 0.5
start_date = datetime.now() + timedelta(days=6)  # Сегодня
day_num_in_plan = 7  # Планируем на неделю
nasa_tle_url = 'http://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html'
objects_lst_file = 'uragan.lst'
# SETUP_END

# SERVICE
google_api_key = 'AIzaSyCDpJOKZ5Lk_noCtfXnFfoR2HokxQuibnU'
tle_database = None
finish_date = start_date + timedelta(days=day_num_in_plan)
forecastio_api_key = '8959e4f1cee6df9bf7b50b87d1d6b8c3'


# SERVICE_END


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


def read_lst(tle_file):
    objects = {}
    with open(tle_file, 'r', encoding='cp866') as lines:
        for params in lines:
            try:
                lon, lat = params.split()[:2]
                lon = float(lon)
                lat = float(lat)
                name = params.split('\'')[1].split('\'')[0]
                objects[name] = (lat, lon)
            except ValueError:
                pass  # Тут нет ни долготы ни широты
    return objects


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
    year = float(year)
    month = float(month)
    day = float(day)
    hours = float(hours)
    minutes = float(minutes)
    seconds = float(seconds)
    dwhole = 367 * year - int(7 * (year + int((month + 9) / 12)) / 4) + int(275 * month / 9) + day - 730531.5
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


def hit_the_earth(iss, t):
    position, velocity = iss.propagate(t.year, t.month, t.day, t.hour, t.minute, t.second)
    GMST_rotation, GMST_rotation_matrix = GMST(t.year, t.month, t.day, t.hour, t.minute, t.second)
    origin, direction = origin_direction(position, velocity, 0, 0, 0, GMST_rotation_matrix, 0)
    hit_1 = calc_intersection(origin[0, 0], origin[0, 1], origin[0, 2], direction[0, 0], direction[0, 1],
                              direction[0, 2])
    return convert_xyz_ll(hit_1)


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


def clean_line(line):
    while True:
        l = len(line)
        line = line.replace('  ', ' ')
        if len(line) == l:
            break
    return line


def tle_date(tle_line_1):
    date = clean_line(tle_line_1).split(' ')[3]
    year = 2000 + int(date[0:2])
    days_since_year_start = float(date[2:])
    return datetime(year, 1, 1, 0, 0, 0, 0) + timedelta(days=days_since_year_start)


def get_future_tle(url):  # get precise future TLE from NASA
    tle_database = {}
    url_text = requests.get(url)
    text = url_text.text.splitlines()
    for i, line in enumerate(text):
        if line.strip() == "ISS":
            tle_line_1 = text[i + 1].strip()
            tle_line_2 = text[i + 2].strip()
            tle_database[tle_date(tle_line_1)] = (tle_line_1, tle_line_2)
    return tle_database


def get_nearest_tle(t):
    sorted_keys = sorted(tle_database.keys())
    answer = tle_database[sorted_keys[0]]  # Первый TLE
    for key in sorted_keys:
        if key > t:
            return answer
        else:
            answer = tle_database[key]
    return answer


def am_i_approaching(obj, t):
    tle_1, tle_2 = get_nearest_tle(t)
    iss = twoline2rv(tle_1, tle_2, wgs84)
    pos_iss_1 = hit_the_earth(iss, t)
    pos_iss_2 = hit_the_earth(iss, t + timedelta(seconds=1))
    h1 = haversine(obj, pos_iss_1)
    h2 = haversine(obj, pos_iss_2)
    if h1 > h2:
        return True
    else:
        return False


def find_traverz(start_datetime, end_datetime, obj):
    answer = None
    if end_datetime - start_datetime < timedelta(seconds=2):
        return start_datetime
    else:
        half = start_datetime + timedelta(seconds=(end_datetime - start_datetime).seconds / 2)
        if am_i_approaching(obj, start_datetime) and not (am_i_approaching(obj, half)):
            answer = find_traverz(start_datetime, half, obj)
        elif am_i_approaching(obj, start_datetime) and am_i_approaching(obj, half):
            answer = find_traverz(half, end_datetime, obj)
        elif not am_i_approaching(obj, start_datetime) and am_i_approaching(obj, half):
            answer = find_traverz(half, end_datetime, obj)
        else:
            return None
    return answer


def slice_timeline(start_date, finish_date):
    loops = []
    t = start_date
    while t <= finish_date + timedelta(minutes=92.62):
        loops.append(t)
        t = t + timedelta(minutes=92.62)
    return loops


def build_timeline(start_date, finish_date):
    timeline = []
    cursor = start_date
    while cursor < finish_date:
        cursor = cursor + timedelta(seconds=1)
        timeline.append(cursor)
    return timeline


def mass_center_in_time(timeline, tle_database):
    i = 0
    mass_center = []
    tle_refresh_times = sorted(tle_database.keys())
    iss = twoline2rv(tle_database[tle_refresh_times[i]][0], tle_database[tle_refresh_times[i]][1], wgs84)  # Первый TLE
    for t in timeline:
        try:
            if t > tle_refresh_times[i + 1]:
                i += 1
                iss = twoline2rv(tle_database[tle_refresh_times[i]][0], tle_database[tle_refresh_times[i]][1], wgs84)
        except IndexError:
            iss = twoline2rv(tle_database[tle_refresh_times[i]][0], tle_database[tle_refresh_times[i]][1], wgs84)
        mass_center.append(hit_the_earth(iss, t))
    return mass_center


def geolocate(tle_line1, tle_line2, alfa, beta, gamma, t):
    iss = twoline2rv(tle_line1, tle_line2, wgs84)
    position, velocity = iss.propagate(t.year, t.month, t.day, t.hour, t.minute, t.second)
    GMST_rotation, GMST_rotation_matrix = GMST(t.year, t.month, t.day, t.hour, t.minute, t.second)
    origin, direction = origin_direction(position, velocity, alfa, beta, gamma, GMST_rotation_matrix, 0)
    hit_1 = calc_intersection(origin[0, 0], origin[0, 1], origin[0, 2], direction[0, 0], direction[0, 1],
                              direction[0, 2])
    return convert_xyz_ll(hit_1)


def paint_session(t, obj_name, obj_lon, obj_lat, tle_database, cloudcover):
    timeline = build_timeline(t - timedelta(seconds=trace_length_in_seconds / 2),
                              t + timedelta(seconds=trace_length_in_seconds / 2))
    trace = mass_center_in_time(timeline, tle_database)
    kml_name = '%s_%s_%s.kml' % (str(t).replace(':', '-'), obj, cloudcover)

    with open(kml_name, 'w', encoding='utf-8') as kml:

        # Рисуем точку съемки, выставляем камеру
        kml.write('''<?xml version="1.0" encoding="UTF-8"?>
        <kml xmlns="http://www.opengis.net/kml/2.2">
            <Document>
                <Style id="illuminator">
                    <LineStyle>
                        <color>7f00ff00</color>
                        <width>3</width>
                        <gx:labelVisibility>1</gx:labelVisibility>
                    </LineStyle>
                </Style>
                <Style id="fov">
                    <LineStyle>
                        <color>7f0000ff</color>
                        <width>3</width>
                        <gx:labelVisibility>1</gx:labelVisibility>
                    </LineStyle>
                </Style>
                <Style id="trace">
                    <LineStyle>
                        <color>ffffffff</color>
                        <width>3</width>
                        <gx:labelVisibility>1</gx:labelVisibility>
                    </LineStyle>
                </Style>
                <Placemark>
                    <name>%s</name>
                    <Camera>
                        <longitude>%s</longitude>
                        <latitude>%s</latitude>
                        <altitude>400000</altitude>
                    </Camera>
                    <Point>
                        <coordinates>%s,%s</coordinates>
                    </Point>
                </Placemark>
        ''' % (obj_name,
               trace[-1 + trace_length_in_seconds // 2][1],
               trace[-1 + trace_length_in_seconds // 2][0],
               obj_lat,
               obj_lon))
        # Рисуем трассу
        _trace = ('%s,%s,40000' % (t[1], t[0]) for t in trace)
        kml.write('''
                <Placemark>
                    <name>%s</name>
                    <styleUrl>#trace</styleUrl>
                    <LineString>
                        <extrude>1</extrude>
                        <tessellate>1</tessellate>
                        <coordinates>
                           %s
                        </coordinates>
                    </LineString>
                </Placemark>
        ''' % ('Трасса МКС', '\n'.join(_trace)))

        # Find latest TLE
        tle_line1, tle_line2 = tle_database[sorted(tle_database.keys())[0]]
        for tle in sorted(tle_database.keys()):
            if t < tle:
                break
            else:
                tle_line1, tle_line2 = tle_database[tle]

        # Рисуем иллюминатор:
        fov = 30.0
        _pos_geo = (geolocate(tle_line1,
                              tle_line2,
                              0.5 * fov * sin(i * pi / 180.0),
                              0,
                              0.5 * fov * cos(i * pi / 180.0),
                              t) for i in range(0, 370, 10))

        kml.write('''
                <Placemark>
                    <name>%s</name>
                    <styleUrl>#illuminator</styleUrl>
                    <LineString>
                        <extrude>1</extrude>
                        <tessellate>1</tessellate>
                        <coordinates>
                            %s
                        </coordinates>
                    </LineString>
                </Placemark>
        ''' % (
            'Иллюминатор',
            '\n'.join('%s,%s,40000' % (pos_geo[1], pos_geo[0]) for pos_geo in _pos_geo)
        ))

        # Рисуем границу отклонения оси визирования НА:
        fov = 60.0
        _pos_geo = (geolocate(tle_line1,
                              tle_line2,
                              0.5 * fov * sin(i * pi / 180.0),
                              0,
                              0.5 * fov * cos(i * pi / 180.0),
                              t)
                    for i in range(0, 370, 10))

        kml.write('''
                <Placemark>
                    <name>%s</name>
                    <styleUrl>#fov</styleUrl>
                    <LineString>
                        <extrude>1</extrude>
                        <tessellate>1</tessellate>
                        <coordinates>
                            %s
                        </coordinates>
                    </LineString>
                </Placemark>
        ''' % (
            'Граница видимости',
            '\n'.join('%s,%s,40000' % (pos_geo[1], pos_geo[0]) for pos_geo in _pos_geo)
        ))

        # Закрываем файл
        kml.write('''
            </Document>
        </kml>''')


def config_logger(logger, fmt, level=logging.DEBUG):
    console = logging.StreamHandler()
    formatter = logging.Formatter(fmt)
    console.setFormatter(formatter)
    console.setLevel(level)
    logger.addHandler(console)


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)-8s [%(asctime)s] [%(module)s:%(lineno)d] %(message)s', level=logging.INFO)

    # DOWNLOAD TLE FROM NASA
    while tle_database is None:
        try:
            tle_database = get_future_tle(nasa_tle_url)
        except Exception as e:
            log.exception(e)
            log.info('NASA не отвечает, пробуем еще раз')

    # DOWNLOADED
    nasa_time = time.time()
    log.info('NASA DOWNLOAD DONE in %.3f sec.' % (nasa_time - program_start_time))
    timeline_time = time.time()
    loops = slice_timeline(start_date, finish_date)
    log.info('TIMELINE BUILT in %.3f sec.' % (timeline_time - nasa_time))
    objects = read_lst(objects_lst_file)
    objects_time = time.time()
    log.info('OBJECTS READ in %.3f sec.' % (objects_time - timeline_time))

    for obj in objects:
        for i in range(len(loops) - 1):
            traverz = find_traverz(loops[i], loops[i + 1], objects[obj])
            if not (isinstance(traverz, (datetime, date)) and 8 < traverz.hour < 21):
                continue
            tle = get_nearest_tle(traverz)
            iss = twoline2rv(tle[0], tle[1], wgs84)
            h1 = haversine(hit_the_earth(iss, traverz), objects[obj])
            if h1 > vision_limit or sun_angle(objects[obj], traverz) < sun_limit:
                continue

            traverz_timestamp = (traverz - datetime(1970, 1, 1)).total_seconds()
            api_params = dict(location='{0},{1}'.format(objects[obj][0], objects[obj][1]),
                              timestamp=traverz_timestamp,
                              key=google_api_key)
            api_response = requests.get('https://maps.googleapis.com/maps/api/timezone/json', params=api_params)
            log.info(api_response.url)
            api_response_dict = api_response.json()
            if api_response_dict['status'] == 'OK':
                rawOffset = api_response_dict['rawOffset']
                traverz_timestamp = traverz_timestamp + rawOffset
                local_traverz = datetime.utcfromtimestamp(traverz_timestamp)
                try:
                    forecast = forecastio.load_forecast(forecastio_api_key, objects[obj][0],
                                                        objects[obj][1],
                                                        time=local_traverz)  # Запрос мог и не пройти
                    cloudcover = forecast.currently().cloudCover
                except:
                    cloudcover = -1
            else:
                local_traverz = traverz
                cloudcover = -1

            if cloudcover < cloudcover_limit:
                log.info('Object: %s; Traverz: %s; Local traverz: %s; Cloud Cover %s',
                         obj, traverz, local_traverz, cloudcover)
                paint_session(traverz, obj, objects[obj][0], objects[obj][1], tle_database, cloudcover)

    log.info('DONE ALL in %.3f sec.' % (time.time() - program_start_time))
