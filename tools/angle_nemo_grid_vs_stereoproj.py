#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

'''
 Compute the rotation angle between NEMO C-grid cells and polar stereographic projection
 The cosine or sine of this angle is used later on to correct the components of velocity
 or displacement vectors
'''

from sys import argv, exit
import numpy as np

from math import acos, asin, atan2, degrees

from netCDF4 import Dataset

from cartopy.crs import PlateCarree, NorthPolarStereo

from climporn  import dump_2d_field
from sitrack import chck4f

idebug=0

# Parameters of Polar Stereographic projection:
rlon0 = -45.
rlat0 =  70.


def __argument_parsing__():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='compute the rotation angle between NEMO C-grid cells and polar stereographic projection')
    #
    # Required:
    rqrdNam = parser.add_argument_group('required arguments')
    rqrdNam.add_argument('-i', '--fin' , required=True,            help='metrics (meshmask) file of the NEMO horizontal domain')
    #
    # Optional:
    parser.add_argument('-xt', '--nlont', default='glamt',         help='name of longitude (at center of grid cells) in input file (default="glamt")')
    parser.add_argument('-yt', '--nlatt', default='gphit',         help='name of latitude  (at center of grid cells) in input file (default="gphit")')
    parser.add_argument('-o',  '--fout' , default='angle_mesh.nc', help='output file (default="angle_mesh.nc")')
    #
    args = parser.parse_args()
    #
    return args.fin, args.nlont, args.nlatt, args.fout


def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]

def cross(vA, vB):
    return vA[0]*vB[1] - vA[1]*vB[0]

def ang(lineA, lineB):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
    cross_prod = cross(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cosA = dot_prod/magA/magB
    sinA = cross_prod/magA/magB
    #
    # Get angle in radians and then convert to degrees
    #anglec = acos(dot_prod/magB/magA)
    #angles = asin(cross_prod/magB/magA)
    #
    # Basically doing angle <- angle mod 360
    #zang1 = degrees(anglec)%360
    #zang2 = degrees(angles)

    zang3 = degrees( atan2(cross_prod, dot_prod) )

    return -zang3, cosA, -sinA




if __name__ == '__main__':

    print('')

    cf_in, cv_lon_t, cv_lat_t, cf_out = __argument_parsing__()

    print(' *** Input file: '+cf_in)
    print('                => name for longitude and latitude: "'+cv_lon_t+'", "'+cv_lat_t+'"\n')


    cv_lon_u, cv_lat_u = str.replace( cv_lon_t, 't', 'u'), str.replace( cv_lat_t, 't', 'u')
    cv_lon_v, cv_lat_v = str.replace( cv_lon_t, 't', 'v'), str.replace( cv_lat_t, 't', 'v')
    cv_lon_f, cv_lat_f = str.replace( cv_lon_t, 't', 'f'), str.replace( cv_lat_t, 't', 'f')


    print(cv_lon_f,cv_lat_f)

    cv_dx_u, cv_dx_v = 'e1u','e1v'
    cv_dy_u, cv_dy_v = 'e2u','e2v'

    list_v_read = [ cv_lon_t, cv_lat_t, cv_lon_u, cv_lat_u, cv_lon_v, cv_lat_v, cv_lat_f, cv_dx_u, cv_dx_v, cv_dy_u, cv_dy_v ]

    chck4f(cf_in)


    with Dataset(cf_in) as id_in:

        list_var = list(id_in.variables.keys())
        #print(list_var)
        for cvt in list_v_read:
            if not cvt in list_var:
                print('ERROR: variable `'+cvt+'` not present in input file '+cf_in)
                exit(0)

        shpLon, shpLat = id_in.variables[cv_lon_t].shape, id_in.variables[cv_lat_t].shape
        l_2d_coordinates = ( len(shpLon)==2 and len(shpLat)==2 )
        l_3d_coordinates = ( len(shpLon)==3 and len(shpLat)==3 )

        if l_2d_coordinates:
            if shpLon != shpLat :
                print('ERROR: lon and lat variables are 2D and have different shapes!'); exit(0)
                print(' *** Coordinates are 2D, irregular grid!')
                (Ny,Nx) = shpLon
        elif ( len(shpLon)==1 and len(shpLat)==1 ):
            print(' *** Coordinates are 1D, regular grid!')
        elif l_3d_coordinates:
            print(' *** Coordinates are 3D, irregular grid!')
            (_,Ny,Nx) = shpLon
        else:
            print('ERROR: could not figure out the shape of coordinates in inpute file...'); exit(0)

            
        del shpLon, shpLat

        print('       ==> domain shape: Ny, Nx =', Ny,Nx,'\n')

        xlat_t,xlon_t = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xlat_t[:,:]   = id_in.variables[cv_lat_t][:,:]
            xlon_t[:,:]   = id_in.variables[cv_lon_t][:,:]
        if l_3d_coordinates:
            xlat_t[:,:]   = id_in.variables[cv_lat_t][0,:,:]
            xlon_t[:,:]   = id_in.variables[cv_lon_t][0,:,:]
            
        xlat_f,xlon_f = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xlat_f[:,:]   = id_in.variables[cv_lat_f][:,:]
            xlon_f[:,:]   = id_in.variables[cv_lon_f][:,:]
        if l_3d_coordinates:
            xlat_f[:,:]   = id_in.variables[cv_lat_f][0,:,:]
            xlon_f[:,:]   = id_in.variables[cv_lon_f][0,:,:]

        xlat_u,xlon_u = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xlat_u[:,:]   = id_in.variables[cv_lat_u][:,:]
            xlon_u[:,:]   = id_in.variables[cv_lon_u][:,:]
        if l_3d_coordinates:
            xlat_u[:,:]   = id_in.variables[cv_lat_u][0,:,:]
            xlon_u[:,:]   = id_in.variables[cv_lon_u][0,:,:]

        xlat_v,xlon_v = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xlat_v[:,:]   = id_in.variables[cv_lat_v][:,:]
            xlon_v[:,:]   = id_in.variables[cv_lon_v][:,:]
        if l_3d_coordinates:
            xlat_v[:,:]   = id_in.variables[cv_lat_v][0,:,:]
            xlon_v[:,:]   = id_in.variables[cv_lon_v][0,:,:]

        xe1u,xe2u = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xe1u[:,:] = id_in.variables[cv_dx_u][:,:]
            xe2u[:,:] = id_in.variables[cv_dy_u][:,:]
        if l_3d_coordinates:
            xe1u[:,:] = id_in.variables[cv_dx_u][0,:,:]
            xe2u[:,:] = id_in.variables[cv_dy_u][0,:,:]

        xe1v,xe2v = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
        if l_2d_coordinates:
            xe1v[:,:] = id_in.variables[cv_dx_v][:,:]
            xe2v[:,:] = id_in.variables[cv_dy_v][:,:]
        if l_3d_coordinates:
            xe1v[:,:] = id_in.variables[cv_dx_v][0,:,:]
            xe2v[:,:] = id_in.variables[cv_dy_v][0,:,:]



    ### closing `cf_in`...


    xX_t,xY_t = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)
    xX_f,xY_f = np.zeros((Ny,Nx), dtype=np.double),np.zeros((Ny,Nx), dtype=np.double)

    crs_src = PlateCarree() ;                                                      # geographic coordinates (lat,lon)
    crs_trg = NorthPolarStereo(central_longitude=rlon0, true_scale_latitude=rlat0) ; # that's (lon,lat) to (x,y)

    zX,zY,_ =  crs_trg.transform_points( crs_src, xlon_t, xlat_t ).T / 1000. ; # [km]
    xY_t[:,:] = zY.T
    xX_t[:,:] = zX.T


    zX,zY,_ =  crs_trg.transform_points( crs_src, xlon_f, xlat_f ).T / 1000.  ; # [km]
    xY_f[:,:] = zY.T
    xX_f[:,:] = zX.T

    del zY, zX

    if idebug>0:
        dump_2d_field( 'Xt.nc', xX_t, name='X_t' )
        dump_2d_field( 'Yt.nc', xY_t, name='Y_t' )
        dump_2d_field( 'Xf.nc', xX_f, name='X_f' )
        dump_2d_field( 'Yf.nc', xY_f, name='Y_f' )



    xang_u, xcos_u, xsin_u = np.zeros((Ny,Nx), dtype=np.double), np.zeros((Ny,Nx), dtype=np.double), np.zeros((Ny,Nx), dtype=np.double)
    xang_v, xcos_v, xsin_v = np.zeros((Ny,Nx), dtype=np.double), np.zeros((Ny,Nx), dtype=np.double), np.zeros((Ny,Nx), dtype=np.double)

    

    for jj in range(Ny):
        for ji in range(Nx-1):
            zLAB  = [ [xX_t[jj,ji],xY_t[jj,ji]], [xX_t[jj,ji+1],xY_t[jj,ji+1]] ]         ; # point A -> point B
            zLAC  = [ [xX_t[jj,ji],xY_t[jj,ji]], [xX_t[jj,ji]+xe1u[jj,ji],xY_t[jj,ji]] ] ; # point A -> point A + (dx,0)
            xang_u[jj,ji], xcos_u[jj,ji], xsin_u[jj,ji] = ang(zLAB, zLAC)

    #for jj in range(Ny-1):
    #    for ji in range(Nx):
    #        zLAB  = [ [xX_t[jj,ji],xY_t[jj,ji]], [xX_t[jj+1,ji],xY_t[jj+1,ji]] ]         ; # point A -> point B
    #        zLAC  = [ [xX_t[jj,ji],xY_t[jj,ji]], [xX_t[jj,ji]+xe2v[jj,ji],xY_t[jj,ji]] ] ; # point A -> point A + (0,dy)
    #        xang_v[jj,ji], xcos_v[jj,ji], xsin_v[jj,ji] = ang(zLAB, zLAC)
            
    for jj in range(Ny):
        for ji in range(1,Nx):
            zLAB  = [ [xX_f[jj,ji-1],xY_f[jj,ji-1]], [xX_f[jj,ji],xY_f[jj,ji]] ]                 ; # point A -> point B
            zLAC  = [ [xX_f[jj,ji-1],xY_f[jj,ji-1]], [xX_f[jj,ji-1]+xe1v[jj,ji],xY_f[jj,ji-1]] ] ; # point A -> point A + (dx,0)
            xang_v[jj,ji], xcos_v[jj,ji], xsin_v[jj,ji] = ang(zLAB, zLAC)

    if idebug>0:
        dump_2d_field( 'angle_u.nc', xang_u, name='angle_u' )
        dump_2d_field( 'angle_v.nc', xang_v, name='angle_v' )


        
        
    # Writing output file:
    id_out = Dataset(cf_out, 'w', format='NETCDF4')

    id_out.createDimension('y', Ny)
    id_out.createDimension('x', Nx)

    id_lat_t  = id_out.createVariable( 'gphit' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_t[:,:] = xlat_t[:,:] ; id_lat_t.units = 'degrees_north'
    id_lon_t  = id_out.createVariable( 'glamt' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_t[:,:] = xlon_t[:,:] ; id_lon_t.units = 'degrees_east'
    #
    id_lat_u  = id_out.createVariable( 'gphiu' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_u[:,:] = xlat_u[:,:] ; id_lat_u.units = 'degrees_north'
    id_lon_u  = id_out.createVariable( 'glamu' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_u[:,:] = xlon_u[:,:] ; id_lon_u.units = 'degrees_east'
    #
    id_lat_v  = id_out.createVariable( 'gphiv' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_v[:,:] = xlat_v[:,:] ; id_lat_v.units = 'degrees_north'
    id_lon_v  = id_out.createVariable( 'glamv' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_v[:,:] = xlon_v[:,:] ; id_lon_v.units = 'degrees_east'
    #
    id_lat_f  = id_out.createVariable( 'gphif' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lat_f[:,:] = xlat_f[:,:] ; id_lat_f.units = 'degrees_north'
    id_lon_f  = id_out.createVariable( 'glamf' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_lon_f[:,:] = xlon_f[:,:] ; id_lon_f.units = 'degrees_east'
    #
    #
    id_ang_u  = id_out.createVariable( 'angle_u' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_ang_u[:,:] = xang_u[:,:] ; id_ang_u.units = 'degrees'
    #
    id_ang_v  = id_out.createVariable( 'angle_v' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_ang_v[:,:] = xang_v[:,:] ; id_ang_v.units = 'degrees'
    #
    id_cos_u  = id_out.createVariable( 'cosa_u' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_cos_u[:,:] = xcos_u[:,:] ; id_cos_u.units = ''
    id_sin_u  = id_out.createVariable( 'sina_u' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_sin_u[:,:] = xsin_u[:,:] ; id_sin_u.units = '-'
    #
    id_cos_v  = id_out.createVariable( 'cosa_v' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_cos_v[:,:] = xcos_v[:,:] ; id_cos_v.units = ''
    id_sin_v  = id_out.createVariable( 'sina_v' ,'f8',('y','x',), zlib=True, complevel=7 )
    id_sin_v[:,:] = xsin_v[:,:] ; id_sin_v.units = ''
    #
    #
    #
    id_out.about = '`angle_nemo_grid_vs_stereoproj.py` of `sitrack`, based on file "'+cf_in+'"'
    id_out.close()

    print('\n *** '+cf_out+' written!\n')
