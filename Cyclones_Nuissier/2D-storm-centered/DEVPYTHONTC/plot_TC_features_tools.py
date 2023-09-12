
def locate_tc_rain_2D(rain,zquantrr,suffix,lat,lon,iicen,ijcen,iie,ije,NX,NY,zoom,ztime,tcname,list_ech,echstep):

    import numpy
    import epygram
    from Magics.macro import mmap, mcoast, mtext, output, mgrib, mcont, mlegend, mtext, minput, msymb, mgraph, plot, page

    ######################################
    # Parameters for the 2D moving carto #
    ######################################

    zdphi = numpy.math.pi / 180.

    area = mmap(subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=lon-(zoom/111.),
        subpage_lower_left_latitude=lat-(zoom/111.),
        subpage_upper_right_longitude=lon+(zoom/111.),
        subpage_upper_right_latitude=lat+(zoom/111.),
         )

    #full_screen = mmap(subpage_x_length         = 26.,
    #                   subpage_y_length         = 16.,
    #                   subpage_x_position       = 0.,
    #                   subpage_y_position       = 0.,
    #                   page_x_length            = 30.,
    #                   page_y_length            = 20.,
    #                   page_id_line             = 'off'
    #                  )
    full_screen = page(layout='positional',
                       subpage_x_length         = 34.,
                       subpage_y_length         = 24.,
                       subpage_x_position       = 0.,
                       subpage_y_position       = 0.,
                       page_x_length            = 30.,
                       page_y_length            = 20.,
                       page_id_line             = 'off'
                      )
    coast = mcoast(map_coastline_land_shade = "off",
           map_coastline_land_shade_colour = "RGB(191,153,107)",
           map_coastline_sea_shade = "off",
           map_coastline_sea_shade_colour = "RGB(108,166,205)",
           map_coastline_thickness = 10,
           map_grid_line_style = "dash",
           map_grid_thickness = 3,
           map_grid_colour = "grey",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = 1,
           map_grid_longitude_increment = 1,
           map_coastline_colour = "RGBA(0,0,0,0.6)")

    output = output(
                output_formats = ['png'],
                output_name = suffix+'-RAIN-CENTERED-'+tcname+'-'+str(list_ech).zfill(4),
                output_width = 1400,
                output_name_first_page_number = "off"
        )


    fout = epygram.formats.resource('precsurf'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')


    for jx in range(max(iicen-NX,1),min(iicen+NX+1,iie)):
        for jy in range(max(ijcen-NY,1),min(ijcen+NY+1,ije)):
            rain.data[jy,jx] = zquantrr[jy-ijcen+NY,jx-iicen+NX]

    fout.writefield(rain)        

    actopr = mgrib(grib_input_file_name='precsurf'+str(list_ech).zfill(4)+'.grib')

    surfrain = mcont(
                    contour                        = 'off',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'off',
                    #contour_level_list             = [1.,2.,5.,10.,20.,30.,50.,75.,100.,125.,150.,200.,250.,300.,350.],
                    contour_level_list = [1.,5.,15.,30.,50.,75.,100.,150.,200.,250.,300.,350.,400.,450.,500.], # Pour RR24
                    contour_level_selection_type   = 'level_list',
                    contour_shade                  = 'on',
                    contour_shade_colour_list      = ["RGBA(0.5,1.,0.,0.6)","RGBA(0.,0.8,0.,0.6)","RGBA(0.,0.54,0.,0.6)","RGBA(0.06,0.3,0.54,0.6)","RGBA(0.12,0.56,1.,0.6)","RGBA(0.,0.93,0.93,0.6)","RGBA(1.,1.,0.,0.6)","RGBA(1.,0.5,0.,0.6)","RGBA(1.,0.,0.,0.6)","RGBA(0.93,0.63,0.68,0.6)","RGBA(0.54,0.27,0.07,0.6)","RGBA(0.84,0.77,1.,0.6)","RGBA(0.57,0.17,0.93,0.6)","RGBA(0.,0.,0.,0.6)"],
                    contour_shade_colour_method    = 'list',
                    contour_shade_technique        = 'cell_shading',
                    contour_shade_method           = 'area_fill',
                    legend                         = 'on',
                    )

    legend = mlegend(
                     legend_display_type       = "continuous",
                     legend_automatic_position = "right",
                     legend_title              = "on",
                     legend_title_text         = "Surface precipitation (mm)",
                     legend_text_font_size     = "1.",
                     legend_text_colour        = "black",
                     legend_label_frequency    = 1)

    title = mtext(
                 text_lines = [str(tcname)+': Max RR'+str(echstep)+' '+str(suffix)+'='+str("%1.f" % numpy.amax(zquantrr))+' Valid time: '+str(ztime)+'Z'],
                 text_justification = "centre",
                 text_font_size     = 1.2,
                 text_colour        = "charcoal")

    lat100=[]
    lon100=[]
    lat200=[]
    lon200=[]

    for jtheta in range(361):
        ztheta = (jtheta) * zdphi
        lat100.append(lat+(100.*numpy.sin(ztheta)/111.))
        lon100.append(lon+(100*numpy.cos(ztheta)/111.))
        lat200.append(lat+(200.*numpy.sin(ztheta)/111.))
        lon200.append(lon+(200*numpy.cos(ztheta)/111.))

    centrepos = minput(
                          input_x_values    =    [lon],
                          input_y_values    =    [lat]
                                )
    symb_centre = msymb(
                  symbol_type    =    "both",
                  symbol_text_list    =    [""],
                  symbol_marker_index    =    18,
                  symbol_colour    =    "red",
                  symbol_height    =    0.3,
                  symbol_outline = 'on',
                  symbol_outline_colour = "black",
                  symbol_outline_thickness = 4,
                  symbol_text_font_size = 0.,
                  symbol_text_font_colour = "red",
                  symbol_text_position = "bottom",
                  symbol_text_font_style = "normal",
                  legend='off',
                  )

    trackcercle100 = minput(
                       input_x_values    =    lon100,
                       input_y_values    =    lat100
                       )
    trackcercle200 = minput(
                     input_x_values    =    lon200,
                       input_y_values    =    lat200
                      )
    tracklinecercle = mgraph(
                          graph_line_colour    =    "RGB(1,0,0)",
                          graph_line_thickness    =    4,
                          graph_line_style    = "dash"
                      )

    AllTrace=[]  


    AllTrace.append(trackcercle100)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle200)
    AllTrace.append(tracklinecercle)
    AllTrace.append(centrepos)
    AllTrace.append(symb_centre)

    plot(output, area, actopr, surfrain, AllTrace, coast, legend, title)



    return

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def locate_tc_wind_2D(wind,zquantw,uin,vin,zquantuinflow,zquantvinflow,suffix,lat,lon,iicen,ijcen,iie,ije,NX,NY,zoom,ztime,tcname,list_ech,echstep,verifBT):

    import numpy
    import epygram
    from Magics.macro import mmap, mcoast, mtext, output, mgrib, mcont, mlegend, mtext, minput, msymb, mgraph, plot, page, mwind

    ######################################
    # Parameters for the 2D moving carto #
    ######################################

    zdphi = numpy.math.pi / 180.

    area = mmap(subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=lon-(zoom/111.),
        subpage_lower_left_latitude=lat-(zoom/111.),
        subpage_upper_right_longitude=lon+(zoom/111.),
        subpage_upper_right_latitude=lat+(zoom/111.),
         )

    #full_screen = mmap(subpage_x_length         = 26.,
    #                   subpage_y_length         = 16.,
    #                   subpage_x_position       = 0.,
    #                   subpage_y_position       = 0.,
    #                   page_x_length            = 30.,
    #                   page_y_length            = 20.,
    #                   page_id_line             = 'off'
    #                  )
    full_screen = page(layout='positional',
                       subpage_x_length         = 34.,
                       subpage_y_length         = 24.,
                       subpage_x_position       = 0.,
                       subpage_y_position       = 0.,
                       page_x_length            = 30.,
                       page_y_length            = 20.,
                       page_id_line             = 'off'
                      )

    coast = mcoast(map_coastline_land_shade = "off",
           map_coastline_land_shade_colour = "RGB(191,153,107)",
           map_coastline_sea_shade = "off",
           map_coastline_sea_shade_colour = "RGB(108,166,205)",
           map_coastline_thickness = 10,
           map_grid_line_style = "dash",
           map_grid_thickness = 3,
           map_grid_colour = "grey",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = 2,
           map_grid_longitude_increment = 2,
           map_coastline_colour = "RGBA(0,0,0,0.6)")

    output = output(
                output_formats = ['png'],
                output_name = suffix+'-WIND-CENTERED-'+tcname+'-'+str(list_ech).zfill(4),
                output_width = 1400,
                output_name_first_page_number = "off"
        )


    fout1 = epygram.formats.resource('windsurf'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')
    fout2 = epygram.formats.resource('windrad'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')


    for jx in range(max(iicen-NX,1),min(iicen+NX+1,iie)):
        for jy in range(max(ijcen-NY,1),min(ijcen+NY+1,ije)):
            uin.data[jy,jx] = zquantuinflow[jy-ijcen+NY,jx-iicen+NX]
            vin.data[jy,jx] = zquantvinflow[jy-ijcen+NY,jx-iicen+NX]
            wind.data[jy,jx] = zquantw[jy-ijcen+NY,jx-iicen+NX]

    fout1.writefield(wind)
    fout2.writefield(uin)
    fout2.writefield(vin)

    surfwind = mgrib(grib_input_file_name='windsurf'+str(list_ech).zfill(4)+'.grib')
    vent10ur = mgrib(grib_input_file_name='windrad'+str(list_ech).zfill(4)+'.grib')

    vent = mcont(
                    contour                        = 'off',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'off',
                    contour_level_list             = [0.,10.,20.,30.,34.,40.,45.,50.,55.,60.,64.,74.,83.,90.,96.,105.,113.,125.,137.,150.,300.],
                    contour_level_selection_type   = 'level_list',
                    contour_shade                  = 'on',
                    contour_shade_colour_list      = ['rgba(255,255,255,0.3)','rgba(51,254,253,0.3)','rgba(0,151,152,0.3)','rgba(0,0,253,0.3)','rgba(0,152,0,0.3)',
                                                            'rgba(0,219,0,0.3)','rgba(254,254,1,0.3)','rgba(254,127,0,0.3)','rgba(253,0,0,0.3)','rgba(189,0,0,0.3)',
                                                            'rgba(127,0,254,0.5)','rgba(152,51,254,0.5)','rgba(177,101,254,0.5)','rgba(203,152,253,0.5)','rgba(227,203,253,0.5)',
                                                            'rgba(253,203,228,0.5)','rgba(254,152,202,0.5)','rgba(253,102,176,0.5)','rgba(229,52,139,0.5)','rgba(204,0,103,0.5)'],
                    contour_shade_colour_method    = 'list',
                    contour_shade_technique        = 'cell_shading',
                    contour_shade_method           = 'area_fill',
                    legend                         = 'on',
                    )

    my_windur = mwind(
                wind_field_type                          = "arrows",
                legend                                   = "on",
                wind_legend_text                         = "",
                wind_arrow_legend_text                   = "",
                #wind_arrow_legend_text                   = 'kts',
                wind_arrow_unit_velocity                 = 20.,
                wind_arrow_min_speed                     = 10.,
                wind_arrow_thickness                     = 6,
                wind_arrow_head_shape                    = 0,
                wind_arrow_head_ratio                    = 1.,
                wind_thinning_factor                     = 2.,
                wind_advanced_method                     = "on",
                wind_advanced_colour_parameter           = "speed",
                wind_advanced_colour_selection_type      = "list",
                wind_advanced_colour_table_colour_method = "list",
                #wind_advanced_colour_list                = ['#f8f9fe','#33fffe','#009899','#0000fe','#009900',
                #                                            '#00dc00','#ffff01','#ff7f00','#fe0000','#be0000',
                #                                            '#7f00ff','#9933ff','#b265ff','#cc99fe','#e4ccfe',
                #                                            '#fecce5','#ff99cb','#fe66b1','#e6348c','#cd0067'],
                wind_advanced_colour_list                = ["RGBA(58,137,35,0.6)"],
                #wind_advanced_colour_level_list          = [0.,10.,20.,30.,34.,40.,45.,50.,55.,60.,64.,74.,83.,90.,96.,105.,113.,125.,137.,150.,300.],
                wind_advanced_colour_level_list          = [10.,300.],
                )


    legend = mlegend(
                     legend_display_type       = "continuous",
                     legend_automatic_position = "right",
                     legend_title              = "on",
                     legend_title_text         = "Surface wind (knots)",
                     legend_text_font_size     = "1.",
                     legend_text_colour        = "black",
                     legend_label_frequency    = 1)

    title = mtext(
                 text_lines = [str(tcname)+': Max wind '+str(suffix)+'='+str("%1.f" % numpy.amax(zquantw))+' Valid time: '+str(ztime)+'Z'],
                 text_justification = "centre",
                 text_font_size     = 1.2,
                 text_colour        = "charcoal")

    lat100=[]
    lon100=[]
    lat200=[]
    lon200=[]
    lat300=[]
    lon300=[]
    lat400=[]
    lon400=[]

    for jtheta in range(361):
        ztheta = (jtheta) * zdphi
        lat100.append(lat+(100.*numpy.sin(ztheta)/111.))
        lon100.append(lon+(100*numpy.cos(ztheta)/111.))
        lat200.append(lat+(200.*numpy.sin(ztheta)/111.))
        lon200.append(lon+(200*numpy.cos(ztheta)/111.))
        lat300.append(lat+(300.*numpy.sin(ztheta)/111.))
        lon300.append(lon+(300*numpy.cos(ztheta)/111.))
        lat400.append(lat+(400.*numpy.sin(ztheta)/111.))
        lon400.append(lon+(400*numpy.cos(ztheta)/111.))


    centrepos = minput(
                          input_x_values    =    [lon],
                          input_y_values    =    [lat]
                                )
    symb_centre = msymb(
                  symbol_type    =    "both",
                  symbol_text_list    =    [""],
                  symbol_marker_index    =    20,
                  symbol_colour    =    "red",
                  symbol_height    =    0.8,
                  symbol_outline = 'on',
                  symbol_outline_colour = "black",
                  symbol_outline_thickness = 4,
                  symbol_text_font_size = 0.,
                  symbol_text_font_colour = "red",
                  symbol_text_position = "bottom",
                  symbol_text_font_style = "normal",
                  legend='off',
                  )

    trackcercle100 = minput(
                       input_x_values    =    lon100,
                       input_y_values    =    lat100
                       )
    trackcercle200 = minput(
                     input_x_values    =    lon200,
                       input_y_values    =    lat200
                      )
    trackcercle300 = minput(
                     input_x_values    =    lon300,
                       input_y_values    =    lat300
                      )
    trackcercle400 = minput(
                     input_x_values    =    lon400,
                       input_y_values    =    lat400
                      )
    tracklinecercle = mgraph(
                          graph_line_colour    =    "RGB(0,0,0)",
                          graph_line_thickness    =    4,
                          graph_line_style    = "dash"
                      )

    AllTrace=[]
    AllCur=[]

    AllTrace.append(trackcercle100)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle200)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle300)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle400)
    AllTrace.append(tracklinecercle)
    AllTrace.append(centrepos)
    AllTrace.append(symb_centre)

    if(verifBT=='true'):
      fobs=open('BT_HISTORY_'+tcname,'r')
      data_obs=numpy.loadtxt('BT_HISTORY_'+tcname,dtype='str')[:,:]


      nbr_obs = 0
      while fobs.readline():
          nbr_obs += 1
      print('Il y a '+str(nbr_obs)+' positions analysées pour le cyclone '+tcname)

      for jobs in range(nbr_obs-1):
          if(str(data_obs[jobs,1]+' '+data_obs[jobs,2])==str(ztime)):
            latcur=float(data_obs[jobs,3])  
            loncur=float(data_obs[jobs,4])

            obs_cur = minput(
                      input_x_values    =    [loncur],
                      input_y_values    =    [latcur]
                      )

            symb_cur=msymb(
                   symbol_type="marker",
                   symbol_marker_index = 18,
                   symbol_colour = "black",
                   symbol_outline="on",
                   symbol_outline_thickness = 3,
                   symbol_height = 0.6,
                   symbol_text_font_colour = "black",
                   symbol_connect_line ="on",
                   symbol_text_font_size = 0.
                   )

            AllCur.append(obs_cur)
            AllCur.append(symb_cur)

          dateobs=data_obs[jobs,1]
          heureobs=data_obs[jobs,2]
          latobs=float(data_obs[jobs,3])
          lonobs=float(data_obs[jobs,4])
          windobs=float(data_obs[jobs,5])  

          obs_pos = minput(
                      input_x_values    =    [lonobs],
                      input_y_values    =    [latobs]
                      )

          if (windobs>=150. and windobs<300.):
              colorobs="#cd0067"
          if (windobs>=137. and windobs<150.):
              colorobs="#e6348c"
          if (windobs>=125. and windobs<137.):
              colorobs="#fe66b1"
          if (windobs>=113. and windobs<125.):
              colorobs="#ff99cb"
          if (windobs>=105. and windobs<113.):
              colorobs="#fecce5"
          if (windobs>=96. and windobs<105. ):
              colorobs="#e4ccfe"
          if (windobs>=90. and windobs<96. ):
              colorobs="#cc99fe"
          if (windobs>=83. and windobs<90. ):
              colorobs="#b265ff"
          if (windobs>=74. and windobs<83. ):
              colorobs="#9933ff"
          if (windobs>=64. and windobs<74. ):
              colorobs="#7f00ff"
          if (windobs>=60. and windobs<64. ):
              colorobs="#be0000"
          if (windobs>=55. and windobs<60. ):
              colorobs="#fe0000"
          if (windobs>=50. and windobs<55. ):
              colorobs="#ff7f00"
          if (windobs>=45. and windobs<50. ):
              colorobs="#ffff01"
          if (windobs>=40. and windobs<45. ):
              colorobs="#00dc00"
          if (windobs>=34. and windobs<40. ):
              colorobs="#009900"
          if (windobs>=30. and windobs<34. ):
              colorobs="#0000fe"
          if (windobs>=20. and windobs<30. ):
              colorobs="#009899"
          if (windobs>=10. and windobs<20. ):
              colorobs="#33fffe"
          if (windobs<10.):
              colorobs="#f8f9fe"

          symb_obs=msymb(
                   symbol_type="marker",
                   symbol_marker_index = 18,
                   symbol_colour = colorobs,
                   symbol_outline="on",
                   symbol_outline_thickness = 5,
                   symbol_height = 0.6,
                   symbol_text_font_colour = "black",
                   symbol_connect_line ="on",
                   symbol_text_font_size = 0.
                   )
     
          AllTrace.append(obs_pos)
          AllTrace.append(symb_obs)


    plot(output, area, surfwind, vent, vent10ur, my_windur, AllTrace, AllCur, coast, legend, title)
    #plot(output, area, surfwind, vent, AllTrace, AllCur, coast, legend, title)



    return
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def locate_tc_isp_2D(isp,zquantisp,suffix,lat,lon,iicen,ijcen,iie,ije,NX,NY,zoom,ztime,tcname,list_ech,echstep):

    import numpy
    import epygram
    from Magics.macro import mmap, mcoast, mtext, output, mgrib, mcont, mlegend, mtext, minput, msymb, mgraph, plot, page, mwind

    ######################################
    # Parameters for the 2D moving carto #
    ######################################

    zdphi = numpy.math.pi / 180.

    area = mmap(subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=lon-(zoom/111.),
        subpage_lower_left_latitude=lat-(zoom/111.),
        subpage_upper_right_longitude=lon+(zoom/111.),
        subpage_upper_right_latitude=lat+(zoom/111.),
         )

    #full_screen = mmap(subpage_x_length         = 26.,
    #                   subpage_y_length         = 16.,
    #                   subpage_x_position       = 0.,
    #                   subpage_y_position       = 0.,
    #                   page_x_length            = 30.,
    #                   page_y_length            = 20.,
    #                   page_id_line             = 'off'
    #                  )
    full_screen = page(layout='positional',
                       subpage_x_length         = 34.,
                       subpage_y_length         = 24.,
                       subpage_x_position       = 0.,
                       subpage_y_position       = 0.,
                       page_x_length            = 30.,
                       page_y_length            = 20.,
                       page_id_line             = 'off'
                      )
    coast = mcoast(map_coastline_land_shade = "off",
           map_coastline_land_shade_colour = "RGB(191,153,107)",
           map_coastline_sea_shade = "off",
           map_coastline_sea_shade_colour = "RGB(108,166,205)",
           map_coastline_thickness = 10,
           map_grid_line_style = "dash",
           map_grid_thickness = 3,
           map_grid_colour = "grey",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = 1,
           map_grid_longitude_increment = 1,
           map_coastline_colour = "RGBA(0,0,0,0.6)")

    output = output(
                output_formats = ['png'],
                output_name = suffix+'-ISP-CENTERED-'+tcname+'-'+str(list_ech).zfill(4),
                output_width = 1400,
                output_name_first_page_number = "off"
        )


    fout = epygram.formats.resource('isp108'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')


    for jx in range(max(iicen-NX,1),min(iicen+NX+1,iie)):
        for jy in range(max(ijcen-NY,1),min(ijcen+NY+1,ije)):
            isp.data[jy,jx] = zquantisp[jy-ijcen+NY,jx-iicen+NX]

    fout.writefield(isp)

    isp108 = mgrib(grib_input_file_name='isp108'+str(list_ech).zfill(4)+'.grib')

    radiance = mcont(
                    contour                        = 'off',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'off',
                    contour_level_list             = [-100,-80,-59,-57,-55,-53,-49,-47,-45,-43,-41,-39,-37,-35,-33,-31,-29,-27,-25,-23,-21,-19,-17,-15,-13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15,17,19,21,23,25,27,40],
                    contour_level_selection_type   = 'level_list',
                    contour_shade                  = 'on',
                    contour_shade_colour_list      = ["rgb(255,20,147)","rgb(1.,0,0)","rgb(1.,0.30,0)","rgb(1.,0.4,0)",
"rgb(1.,0.50,0)","rgb(1.,0.50,0)","rgb(1.,0.65,0)",
"rgb(1.,0.65,0)","rgb(1.,0.75,0)","rgb(1.,0.85,0)",
"rgb(1.,0.85,0)","rgb(1.,0.92,0)","rgb(0.80,0.80,0.80)",
"rgb(0.78,0.78,0.78)","rgb(0.77,0.77,0.77)","rgb(0.72,0.83,0.73)",
"rgb(0.72,0.83,0.73)","rgb(0.70,0.80,0.70)","rgb(0.72,0.81,0.64)",
"rgb(0.72,0.81,0.64)","rgb(0.74,0.83,0.61)","rgb(0.72,0.78,0.56)",
"rgb(0.72,0.78,0.56)","rgb(0.69,0.75,0.58)","rgb(0.67,0.72,0.57)",
"rgb(0.65,0.70,0.54)","rgb(0.63,0.67,0.51)","rgb(0.58,0.66,0.56)",
"rgb(0.56,0.67,0.61)","rgb(0.51,0.65,0.66)","rgb(0.45,0.64,0.69)",
"rgb(0.40,0.63,0.74)","rgb(0.32,0.63,0.78)","rgb(0.20,0.60,0.82)",
"rgb(0.18,0.57,0.78)","rgb(0.16,0.54,0.74)","rgb(0.13,0.51,0.71)",
"rgb(0.11,0.48,0.67)","rgb(0.09,0.46,0.63)","rgb(0.07,0.42,0.59)",
"rgb(0.05,0.38,0.55)","rgb(0.04,0.36,0.51)","rgb(0.03,0.34,0.47)",
"rgb(0.02,0.31,0.43)","rgb(0.01,0.29,0.39)","rgb(0.00,0.25,0.33)",
"rgb(0.00,0.23,0.28)"],
                    contour_shade_colour_method    = 'list',
                    contour_shade_technique        = 'cell_shading',
                    contour_shade_method           = 'area_fill',
                    legend                         = 'on',
                    )

    legend = mlegend(
                     legend_display_type       = "continuous",
                     legend_automatic_position = "right",
                     legend_title              = "on",
                     legend_title_text         = "Brightness Temperature IR10.8 (°C)",
                     legend_text_font_size     = "1.",
                     legend_text_colour        = "black",
                     legend_label_frequency    = 1)

    title = mtext(
                 text_lines = [str(tcname)+': Max BT10.8 '+str(suffix)+'='+str("%1.f" % numpy.amax(zquantisp))+' Valid time: '+str(ztime)+'Z'],
                 text_justification = "centre",
                 text_font_size     = 1.2,
                 text_colour        = "charcoal")

    lat100=[]
    lon100=[]
    lat200=[]
    lon200=[]

    for jtheta in range(361):
        ztheta = (jtheta) * zdphi
        lat100.append(lat+(100.*numpy.sin(ztheta)/111.))
        lon100.append(lon+(100*numpy.cos(ztheta)/111.))
        lat200.append(lat+(200.*numpy.sin(ztheta)/111.))
        lon200.append(lon+(200*numpy.cos(ztheta)/111.))

    centrepos = minput(
                          input_x_values    =    [lon],
                          input_y_values    =    [lat]
                                )
    symb_centre = msymb(
                  symbol_type    =    "both",
                  symbol_text_list    =    [""],
                  symbol_marker_index    =    18,
                  symbol_colour    =    "black",
                  symbol_height    =    0.3,
                  symbol_outline = 'on',
                  symbol_outline_colour = "black",
                  symbol_outline_thickness = 4,
                  symbol_text_font_size = 0.,
                  symbol_text_font_colour = "red",
                  symbol_text_position = "bottom",
                  symbol_text_font_style = "normal",
                  legend='off',
                  )

    trackcercle100 = minput(
                       input_x_values    =    lon100,
                       input_y_values    =    lat100
                       )
    trackcercle200 = minput(
                     input_x_values    =    lon200,
                       input_y_values    =    lat200
                      )
    tracklinecercle = mgraph(
                          graph_line_colour    =    "RGB(0,0,0)",
                          graph_line_thickness    =    4,
                          graph_line_style    = "dash"
                      )

    AllTrace=[]

    AllTrace.append(trackcercle100)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle200)
    AllTrace.append(tracklinecercle)
    AllTrace.append(centrepos)
    AllTrace.append(symb_centre)

    plot(output, area, isp108, radiance, AllTrace, coast, legend, title)



    return
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def locate_tc_innerproc_2D(dbz12km,u10,v10,cape,zprobacb,zquantuinflow,zquantvinflow,zquantcape,suffix,lat,lon,iicen,ijcen,iie,ije,NX,NY,zoom,ztime,tcname,list_ech,echstep,zshear,mb):

    import numpy
    import epygram
    from Magics.macro import mmap, mcoast, mtext, output, mgrib, mcont, mlegend, mtext, minput, msymb, mgraph, plot, page, mwind

    ######################################
    # Parameters for the 2D moving carto #
    ######################################

    zdphi = numpy.math.pi / 180.

    area = mmap(subpage_map_projection="cylindrical",
        subpage_lower_left_longitude=lon-(zoom/111.),
        subpage_lower_left_latitude=lat-(zoom/111.),
        subpage_upper_right_longitude=lon+(zoom/111.),
        subpage_upper_right_latitude=lat+(zoom/111.),
         )

    #full_screen = mmap(subpage_x_length         = 26.,
    #                   subpage_y_length         = 16.,
    #                   subpage_x_position       = 0.,
    #                   subpage_y_position       = 0.,
    #                   page_x_length            = 30.,
    #                   page_y_length            = 20.,
    #                   page_id_line             = 'off'
    #                  )
    full_screen = page(layout='positional',
                       subpage_x_length         = 34.,
                       subpage_y_length         = 24.,
                       subpage_x_position       = 0.,
                       subpage_y_position       = 0.,
                       page_x_length            = 30.,
                       page_y_length            = 20.,
                       page_id_line             = 'off'
                      )
    coast = mcoast(map_coastline_land_shade = "off",
           map_coastline_land_shade_colour = "RGB(191,153,107)",
           map_coastline_sea_shade = "off",
           map_coastline_sea_shade_colour = "RGB(108,166,205)",
           map_coastline_thickness = 10,
           map_grid_line_style = "dash",
           map_grid_thickness = 3,
           map_grid_colour = "grey",
           map_label = "on",
           map_label_height = 0.9,
           map_grid_latitude_increment = 1,
           map_grid_longitude_increment = 1,
           map_coastline_colour = "RGBA(254,163,71,1.)")

    output = output(
                output_formats = ['png'],
                output_name = suffix+'-CB-CENTERED-'+tcname+'-'+str(list_ech).zfill(4),
                output_width = 1400,
                output_name_first_page_number = "off"
        )


    fout1 = epygram.formats.resource('CB1'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')
    fout2 = epygram.formats.resource('windrad10'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')
    fout3 = epygram.formats.resource('MUCAPE'+str(list_ech).zfill(4)+'.grib', 'w', 'GRIB')


    for jx in range(max(iicen-NX,1),min(iicen+NX+1,iie)):
        for jy in range(max(ijcen-NY,1),min(ijcen+NY+1,ije)):
            dbz12km.data[jy,jx] = zprobacb[jy-ijcen+NY,jx-iicen+NX]
            u10.data[jy,jx] = zquantuinflow[jy-ijcen+NY,jx-iicen+NX]
            v10.data[jy,jx] = zquantvinflow[jy-ijcen+NY,jx-iicen+NX]
            cape.data[jy,jx] = zquantcape[jy-ijcen+NY,jx-iicen+NX]

    fout1.writefield(dbz12km)
    fout2.writefield(u10)
    fout2.writefield(v10)
    fout3.writefield(cape)


    cb1 = mgrib(grib_input_file_name='CB1'+str(list_ech).zfill(4)+'.grib')
    vent10ur = mgrib(grib_input_file_name='windrad10'+str(list_ech).zfill(4)+'.grib')
    cape = mgrib(grib_input_file_name='MUCAPE'+str(list_ech).zfill(4)+'.grib')


    convburst = mcont(
                    contour                        = 'off',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'off',
                    contour_level_list             = [5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.],
                    contour_level_selection_type   = 'level_list',
                    contour_shade                  = 'on',
                    #contour_shade_colour_list      = ["RGBA(0.,0.54,0.,0.5)","RGBA(0.,0.80,0.,0.5)","RGBA(0.50,1.,0.,0.5)","RGBA(1.,1.,0.,0.5)","RGBA(1.,0.84,0.,0.5)","RGBA(0.80,0.52,0.,0.5)","RGBA(1.,0.50,0.,0.5)","RGBA(0.80,0.,0.,0.5)","RGBA(0.54,0.,0.,0.5)","RGBA(0.54,0.,0.54,0.5)"],
                    contour_shade_colour_method     = 'calculate',
                    contour_shade_max_level_colour  = "RGBA(0.,0.,0.,1.)",
                    contour_shade_min_level_colour  = "RGBA(1.,1.,1.,1.)",

                    #contour_shade_colour_method    = 'list',
                    contour_shade_technique        = 'cell_shading',
                    contour_shade_method           = 'area_fill',
                    legend                         = 'on',
                    )

    my_windur = mwind(
                wind_field_type                          = "arrows",
                legend                                   = "on",
                wind_legend_text                         = "",
                wind_arrow_legend_text                   = "",
                #wind_arrow_legend_text                   = 'kts',
                wind_arrow_unit_velocity                 = 30.,
                wind_arrow_min_speed                     = 15.,
                wind_arrow_thickness                     = 3,
                wind_arrow_head_shape                    = 3,
                wind_arrow_head_ratio                    = 1.,
                wind_thinning_factor                     = 7.,
                wind_advanced_method                     = "on",
                wind_advanced_colour_parameter           = "speed",
                wind_advanced_colour_selection_type      = "list",
                wind_advanced_colour_table_colour_method = "list",
                #wind_advanced_colour_list                = ['#f8f9fe','#33fffe','#009899','#0000fe','#009900',
                #                                            '#00dc00','#ffff01','#ff7f00','#fe0000','#be0000',
                #                                            '#7f00ff','#9933ff','#b265ff','#cc99fe','#e4ccfe',
                #                                            '#fecce5','#ff99cb','#fe66b1','#e6348c','#cd0067'],
                wind_advanced_colour_list                = ["RGBA(255,0,255,0.6)"],
                #wind_advanced_colour_level_list          = [0.,10.,20.,30.,34.,40.,45.,50.,55.,60.,64.,74.,83.,90.,96.,105.,113.,125.,137.,150.,300.],
                wind_advanced_colour_level_list          = [15.,300.],
                )

    instalev1 = mcont(
                    contour                        = 'on',
                    contour_line_style             = 'dash',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'on',
                    contour_label_height           = 0.8,
                    contour_label_blanking         = 'on',
                    contour_level_list             = [-10.],
                    contour_shade                   = 'off',
                    contour_level_selection_type   = 'level_list',
                    contour_line_colour            = "RGBA(0,0,1,0.5)",
                    contour_line_thickness         = 10,
                    legend                         = 'on',
                    )
    instalev2 = mcont(
                    contour                        = 'on',
                    contour_line_style             = 'dash',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'on',
                    contour_label_height           = 0.8,
                    contour_label_blanking         = 'on',
                    contour_level_list             = [-12.],
                    contour_shade                   = 'off',
                    contour_level_selection_type   = 'level_list',
                    contour_line_colour            = "RGBA(0,1,0,0.5)",
                    contour_line_thickness         = 10,
                    legend                         = 'on',
                    )
    instalev3 = mcont(
                    contour                        = 'on',
                    contour_line_style             = 'dash',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'on',
                    contour_label_height           = 0.8,
                    contour_label_blanking         = 'on',
                    contour_level_list             = [-14.],
                    contour_shade                   = 'off',
                    contour_level_selection_type   = 'level_list',
                    contour_line_colour            = "RGBA(1,0,0,0.5)",
                    contour_line_thickness         = 10,
                    legend                         = 'on',
                    )

    my_instability = mcont(
                    contour                        = 'on',
                    contour_highlight              = 'off',
                    contour_hilo                   = 'off',
                    contour_label                  = 'off',
                    contour_level_list             = [-16.,-14.,-12.,-10.,-8.,-6.,-4.],
                    contour_shade                   = 'off',
                    contour_level_selection_type   = 'level_list',
                    contour_shade_colour_method     = 'calculate',
                    contour_shade_max_level_colour  = "RGBA(1.,0.,0.,0.8)",
                    contour_shade_min_level_colour  = "RGBA(0.,0.,1.,0.8)",
                    contour_line_thickness         = 6,
                    contour_shade_colour_direction  = 'clockwise',
                    contour_shade_technique        = 'cell_shading',
                    contour_shade_method           = 'area_fill',
                    legend                         = 'on',
                    )

    legend = mlegend(
                     legend_display_type       = "continuous",
                     legend_automatic_position = "right",
                     legend_title              = "on",
                     legend_title_text         = "Probability of CB ocurrence (%)",
                     legend_text_font_size     = "1.",
                     legend_text_colour        = "black",
                     legend_label_frequency    = 1)

    title = mtext(
                 text_lines = [str(tcname)+': CB + inflow + instability - Valid time: '+str(ztime)+'Z'],
                 text_justification = "centre",
                 text_font_size     = 1.2,
                 text_colour        = "charcoal")

    lat100=[]
    lon100=[]
    lat200=[]
    lon200=[]

    for jtheta in range(361):
        ztheta = (jtheta) * zdphi
        lat100.append(lat+(100.*numpy.sin(ztheta)/111.))
        lon100.append(lon+(100*numpy.cos(ztheta)/111.))
        lat200.append(lat+(200.*numpy.sin(ztheta)/111.))
        lon200.append(lon+(200*numpy.cos(ztheta)/111.))

    centrepos = minput(
                          input_x_values    =    [lon],
                          input_y_values    =    [lat]
                                )
    symb_centre = msymb(
                  symbol_type    =    "both",
                  symbol_text_list    =    [""],
                  symbol_marker_index    =    18,
                  symbol_colour    =    "black",
                  symbol_height    =    0.3,
                  symbol_outline = 'on',
                  symbol_outline_colour = "black",
                  symbol_outline_thickness = 4,
                  symbol_text_font_size = 0.,
                  symbol_text_font_colour = "red",
                  symbol_text_position = "bottom",
                  symbol_text_font_style = "normal",
                  legend='off',
                  )

    trackcercle100 = minput(
                       input_x_values    =    lon100,
                       input_y_values    =    lat100
                       )
    trackcercle200 = minput(
                     input_x_values    =    lon200,
                       input_y_values    =    lat200
                      )
    tracklinecercle = mgraph(
                          graph_line_colour    =    "RGB(0,0,0)",
                          graph_line_thickness    =    4,
                          graph_line_style    = "dash"
                      )

    AllTraceshear=[]

    for jb in range(mb):
        latshearend = lat+(200*numpy.sin(zshear[jb])/111.)
        lonshearend = lon+(200*numpy.cos(zshear[jb])/111.)

        shear_arrow = minput(
                       input_x_values    =    [lon,lonshearend],
                       input_y_values    =    [lat,latshearend],
                             )

        if (jb==0):
           colorshear="red"
        else:
           colorshear="RGBA(0,0,0,0.3)"


        shear_shape = mgraph(
                   graph_line_colour    =    colorshear,
                   graph_line_thickness    =    30,
                   graph_line_style    =   "solid"
                   )

        AllTraceshear.append(shear_arrow)
        AllTraceshear.append(shear_shape)

    AllTrace=[]

    AllTrace.append(trackcercle100)
    AllTrace.append(tracklinecercle)
    AllTrace.append(trackcercle200)
    AllTrace.append(tracklinecercle)
    AllTrace.append(centrepos)
    AllTrace.append(symb_centre)

    #plot(output, area, isp108, radiance, AllTrace, coast, legend, title)
    plot(output, area, cb1, convburst, cape, instalev1, cape, instalev2, cape, instalev3, vent10ur, my_windur, AllTrace, AllTraceshear, legend, title, coast)



    return
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




