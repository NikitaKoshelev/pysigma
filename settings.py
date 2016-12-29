# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division
from datetime import datetime, timedelta

# SETUP_START
sun_limit = 20
vision_limit = 220
trace_length_in_seconds = 360
cloud_cover_limit = 0.5
start_date = datetime.now() + timedelta(days=6)  # Сегодня
day_num_in_plan = 7  # Планируем на неделю
NASA_TLE_URL = nasa_tle_url = 'http://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html'
objects_lst_file = 'uragan.lst'

# SERVICE
google_api_key = 'AIzaSyCDpJOKZ5Lk_noCtfXnFfoR2HokxQuibnU'
tle_database = None
finish_date = start_date + timedelta(days=day_num_in_plan)
forecastio_api_key = '8959e4f1cee6df9bf7b50b87d1d6b8c3'
