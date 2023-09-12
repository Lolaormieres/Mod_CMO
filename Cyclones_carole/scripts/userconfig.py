#!/usr/bin/env python
# -*- coding: utf-8 -*-

#epsilon = 1e-12
#mask_outside=1e17
# usermodules = [{'name':'module1', 'abspath':'/path/to/module1.py'},
#                {'name':'module2', 'abspath':'/path/to/module2.py'}
#                ]
#use_footprints_as_builder = True
#usercolormaps = {'my_cmap':'/path/to/my_cmap.cmap'}
#implemented_formats = ['FA',]
noninteractive_backend = "Agg"
usercolormaps={'quantraf_kt':'/home/labadie/.epygram/palettes/quantraf_kt.cmap',
                   'quantraf_kmh':'/home/labadie/.epygram/palettes/quantraf_kmh.cmap',
                   'quantraf_ms':'/home/labadie/.epygram/palettes/quantraf_ms.cmap',
#                   'quantraf_tempete':'/home/labadie/.epygram/palettes/quantraf_tempete.cmap',
                   'quantraf_cyclone':'/home/labadie/.epygram/palettes/quantraf_tempete.cmap',
#                   'rr6hproba': '/home/labadie/.epygram/palettes/rr6hproba.cmap',
                   'rr6hquantclima': '/home/labadie/.epygram/palettes/rr6hquantclima.cmap',
                   'tsurf': '/home/labadie/.epygram/palettes/quantraf_ms.cmap',
#                   'neiquantclima': '/home/labadie/.epygram/palettes/rr6hquantclima.cmap',
                   'rrfortquantclima': '/home/labadie/.epygram/palettes/rrfortquantclima.cmap',
                   'rrfortquantantilope': '/home/labadie/.epygram/palettes/rrfortquantantilope.cmap',
                   'testrr': '/home/labadie/.epygram/palettes/testrr.json',
                   'testrr_modere': '/home/labadie/.epygram/palettes/testrr_modere.json',
                   'quantraf_tempete':'/home/labadie/.epygram/palettes/quantraf_tempete3.json',
                   'quantraf':'/home/labadie/.epygram/palettes/quantraf.json',
                   'quantrafNEC':'/home/labadie/.epygram/palettes/quantrafNEC.json',
                   'rr6hproba': '/home/labadie/.epygram/palettes/rr6hproba.json',
                   'neiquant': '/home/labadie/.epygram/palettes/neiquant.json',
                   'rrArnaud': '/home/labadie/.epygram/palettes/rrArnaud.json',
                   'rrArnaud_zoom': '/home/labadie/.epygram/palettes/rrArnaud_zoom.json',
                   }
usercolormaps_scaling={'quantraf_kt':[30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0,150.0],
                           'quantraf_tempete':[60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,300.0],
                           'quantraf_cyclone':[60.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0],
                           'quantraf_kmh':[40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,180.0,220.0],
                           'quantraf_ms':[10.0,13.0,16.0,19.0,22.0,25.0,28.0,31.0,33.0,36.0,39.0,42.0,50.0,62.0],
#                           'rr6h':[0., 0.5, 1., 2., 4., 8., 10., 15., 20., 25., 30., 40.],
                           'rr6h':[0., 1., 5., 10., 20., 30., 40., 50., 60., 80., 100., 120.],
                           'rr6hquantclima':[0., 1., 5., 10., 20., 30., 40., 50., 60., 80., 100., 120.],
# rr6hquantclima                           'tsurf':[200., 250., 255., 260., 265., 270., 275., 280., 285., 290., 295., 300.],
                           'tsurf':[200., 250., 255., 260., 265., 270., 275., 280., 285., 290., 295., 300.,350.,400.],
                           'rrfortquantantilope':[0., 0.1, 0.5, 1.0, 2.0, 5.0, 10., 20., 30., 50., 75., 100., 150., 200., 300., 500.],
                           'rrfortquantclima':[0., 1., 5., 10., 20., 30., 40., 50., 60., 80., 100., 120.,150.,200.,250.,300.,350.,400.,450.,500.],
                           'neiquantclima':[0.,0.5, 1., 2., 3., 4., 5., 10., 20.,  30.,  40., 50.],
#                           'neiquantclima':[0.,0.1,0.5, 1., 2., 3., 4., 5., 10., 20.,  30.,  40.],
#                           'rr6h':[0., 1., 2., 4., 10., 25., 30., 40., 50., 60., 80., 100.],
                           'rr6hproba':[5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.],
#moyen                            'rr6h':[0., 0.5, 2., 4., 10., 25., 50., 100., 150., 200., 250., 300.],
# fort                           'rr6h':[0., 5., 10., 15., 20., 25., 30., 40., 50., 60., 80., 100.],
                            'rrArnaud':["0", "0.1", "1", "3", "5", "7", "10", "15", "20", "30", "50", "70", "100", "150", "200", "250", "300", "350", "500"],
                            'rrArnaud_zoom':["0", "3", "5", "7", "10", "15", "20", "30", "50", "70", "100", "150", "200", "500"],
                           }
