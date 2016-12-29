# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division

from datetime import timedelta
from math import sin, pi, cos

import settings
from calculation import mass_center_in_time, geolocate
from utils import build_timeline, open

TRACE_LENGTH_IN_SECONDS = getattr(settings, 'trace_length_in_seconds')


def paint_session(t, obj_name, obj_lon, obj_lat, tle, cloud_cover):
    timeline = build_timeline(t - timedelta(seconds=TRACE_LENGTH_IN_SECONDS / 2),
                              t + timedelta(seconds=TRACE_LENGTH_IN_SECONDS / 2))
    trace = mass_center_in_time(timeline, tle)
    kml_name = '%s_%s_%s.kml' % (str(t).replace(':', '-'), obj_name, cloud_cover)

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
               trace[-1 + TRACE_LENGTH_IN_SECONDS // 2][1],
               trace[-1 + TRACE_LENGTH_IN_SECONDS // 2][0],
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
        iss_p, iss_v = tle.get_iss_position(dt=t)

        # Рисуем иллюминатор:
        fov = 30.0
        _pos_geo = (geolocate(iss_p,
                              iss_v,
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
        _pos_geo = (geolocate(iss_p,
                              iss_v,
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
