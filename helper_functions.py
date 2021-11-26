def gdalds2cube(gdalds, timestamp=None):

    '''
    :param gdalds: Input can be either a geotiff file or a gdal dataset object
    :param timestamp: Datetime (as a string formatted %Y%m%d%H%M)
    :return: iris cube
    '''

    import sys
    from osgeo import gdal
    import numpy as np
    import numpy.ma as ma
    import datetime as dt
    import cf_units
    import iris


    if isinstance(gdalds, str):
        ds = gdal.Open(gdalds)
    elif isinstance(gdalds, gdal.Dataset):
        ds = gdalds
    else:
        sys.exit('Didn\'t recognise the input to gdalds2cube')

    ulX, xDist, rtnX, ulY, rtnY, yDist = ds.GetGeoTransform()
    xDist = np.round(xDist, 6) # Round because sometimes gdal spits out too precise coordinates
    yDist = np.round(yDist, 6)
    yDist = yDist if yDist < 0.0 else yDist * -1 # The y cell length should be negative because we use the UPPER Y to calculate LOWER Y
    ncellx = ds.RasterXSize
    ncelly = ds.RasterYSize
    lrX    = (ulX + (ncellx * xDist)) - (xDist * 0.5)
    lrY    = (ulY + (ncelly * yDist)) - (yDist * 0.5)
    ulX    = ulX + ( xDist * 0.5 ) # Assuming the GT coord refers to the ul of the pixel
    ulY    = ulY + ( yDist * 0.5) # Assuming the GT coord refers to the ul of the pixel
    # proj   = ds.GetProjection()

    latcoord = iris.coords.DimCoord(np.arange(ulY, lrY + yDist, yDist), standard_name='latitude', units=cf_units.Unit('degrees'), coord_system=iris.coord_systems.GeogCS(6371229.0))
    loncoord = iris.coords.DimCoord(np.arange(ulX, lrX + xDist, xDist), standard_name='longitude', units=cf_units.Unit('degrees'), coord_system=iris.coord_systems.GeogCS(6371229.0))
    latcoord.guess_bounds()
    loncoord.guess_bounds()

    if ds.GetLayerCount() == 0:
        band = ds.GetRasterBand(1)
        array = band.ReadAsArray()
    else:
        print('More than 1 raster band. What do we do next?')
        return

    nodatavalue = band.GetNoDataValue()
    if not nodatavalue == None:
        array = ma.masked_equal(array, nodatavalue)

    if timestamp:
        u = cf_units.Unit('hours since 1970-01-01 00:00:00', calendar=cf_units.CALENDAR_STANDARD)
        timecoord = iris.coords.DimCoord(u.date2num(dt.datetime.strptime(timestamp, '%Y%m%d%H%M')), standard_name='time', units=u)
        array = array[np.newaxis, ...]
        cube = iris.cube.Cube(array, dim_coords_and_dims=[(timecoord, 0), (latcoord, 1), (loncoord, 2)])
    else:
        cube = iris.cube.Cube(array, dim_coords_and_dims=[(latcoord, 0), (loncoord, 1)])

    return cube
