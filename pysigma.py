# -*- coding: utf-8  -*-
from __future__ import unicode_literals, division

import logging
import time
from datetime import date

from haversine import haversine

from calculation import find_traverz, hit_the_earth, sun_angle
from kml import paint_session
from settings import *
from tle import TLE
from utils import slice_timeline, read_lst, get_local_dt, get_cloud_cover

log = logging.getLogger(__name__)

if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)-8s [%(asctime)s] [%(module)s:%(lineno)d] %(message)s', level=logging.INFO)

    program_start_time = time.time()
    tle = TLE()
    nasa_time = time.time()
    log.info('NASA DOWNLOAD DONE in %.3f sec.' % (nasa_time - program_start_time))
    timeline_time = time.time()
    loops = slice_timeline(start_date, finish_date)
    log.info('TIMELINE BUILT in %.3f sec.' % (timeline_time - nasa_time))
    objects = read_lst(objects_lst_file)
    objects_time = time.time()
    log.info('OBJECTS READ in %.3f sec.' % (objects_time - timeline_time))

    for obj_name, obj in objects.items():
        for i in range(len(loops) - 1):
            traverz = find_traverz(loops[i], loops[i + 1], obj)
            if not isinstance(traverz, (datetime, date)):
                continue

            iss = tle.get_iss(dt=traverz)

            exprs = [
                haversine(hit_the_earth(iss, traverz), obj) > vision_limit,
                sun_angle(obj, traverz) < sun_limit
            ]
            if any(exprs):
                continue

            traverz_timestamp = (traverz - datetime(1970, 1, 1)).total_seconds()
            local_traverz = get_local_dt(obj[0], obj[1], traverz_timestamp)

            if local_traverz:
                cloud_cover = get_cloud_cover(obj[0], obj[1], local_traverz)
            else:
                local_traverz = traverz
                cloud_cover = -1

            if 8 < local_traverz.hour < 21 and cloud_cover < cloud_cover_limit:
                log.info('Object: %s; Traverz: %s; Local traverz: %s; Cloud Cover %s',
                         obj_name, traverz, local_traverz, cloud_cover)
                paint_session(traverz, obj_name, obj[0], obj[1], tle, cloud_cover)

    log.info('DONE ALL in %.3f sec.' % (time.time() - program_start_time))
