# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division, print_function, with_statement

import logging
import re
import sys
from datetime import timedelta, datetime

import forecastio
import requests

import settings

if sys.version_info.major < 3:
    try:
        from codecs import open
    except:
        pass
    
GOOGLE_API_KEY = getattr(settings, 'google_api_key')
FORECASTIO_API_KEY = getattr(settings, 'forecastio_api_key')
CLEAN_EXP = re.compile(r' +')

log = logging.getLogger(__name__)


def clean_line(line):
    return CLEAN_EXP.sub(' ', line)


def slice_timeline(start_date, finish_date):
    step = timedelta(minutes=92.62)
    return list(dt_range(start_date, finish_date, step=step, include_right=True))


def build_timeline(start_date, finish_date):
    return list(dt_range(start_date, finish_date, include_left=False, include_right=True))


def dt_range(start_dt, end_dt, step=None, include_left=True, include_right=False):
    """
    `range()` implementation for datetime
    :param start_dt:
    :param finish_dt:
    :param step:
    :return:
    """
    if not isinstance(step, timedelta):
        step = timedelta(seconds=1)

    if include_left:
        yield start_dt

    cursor = start_dt + step
    while cursor < end_dt:
        yield cursor
        cursor += step

    if include_right:
        yield cursor


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


def config_logger(logger, fmt, level=logging.DEBUG):
    console = logging.StreamHandler()
    formatter = logging.Formatter(fmt)
    console.setFormatter(formatter)
    console.setLevel(level)
    logger.addHandler(console)


def get_local_dt(lat, lon, timestamp):
    params = dict(location='{0},{1}'.format(lat, lon),
                  timestamp=timestamp,
                  key=GOOGLE_API_KEY)
    response = requests.get('https://maps.googleapis.com/maps/api/timezone/json', params=params)
    try:
        log.info(response.url)
        api_response = response.json()
        if api_response.get('status') == 'OK':
            local_timestamp = timestamp + api_response.get('rawOffset')
            local_dt = datetime.utcfromtimestamp(local_timestamp)
            return local_dt
    finally:
        response.close()


def get_cloud_cover(lat, lon, dt=None):
    cloud_cover = -1
    try:
        # Запрос мог и не пройти
        forecast = forecastio.load_forecast(FORECASTIO_API_KEY, lat, lon, time=dt)
        cloud_cover = forecast.currently().cloudCover
    except Exception as e:
        log.exception(e, exc_info=True)
    return cloud_cover


if __name__ == '__main__':
    start_dt = datetime(2016, 12, 1)
    end_dt = datetime(2017, 1, 1)
    step = timedelta(days=1)
    print(*slice_timeline(start_dt, end_dt), sep='\n')
