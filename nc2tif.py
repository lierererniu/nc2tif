"""
my environments:
    python 3.8.5
    modules:
        netCDF4 1.5.4
        xarray 0.16.1
        numpy 1.19.1
        Gdal 3.1.3
ע�����
	1���ļ�·����ò�Ҫ������·��
	2������ļ�������
        	3. �����ļ���ֻ�ܴ��.nc�ļ�����Ҫ�������ļ�
"""
from osgeo import gdal, osr, gdal_array
import xarray as xr
import numpy as np
import sys
import netCDF4 as nc
from netCDF4 import num2date
import os
def GetValueBynetCDF4(in_filename):
    '''
    ͨ��netCDF4��.nc�ļ�
    Open the .nc file through netCDF4
    '''
    nc_data = nc.Dataset(in_filename)
    return nc_data

def get_time(index,nc_data,var,filename):
    '''
    ��ȡʱ���ǩֵ
    Get the time label value
    '''
    time = nc_data.variables['time'][:]
    time_units = nc_data.variables['time'].units
    timed = num2date(time, time_units, 'standard')  # 1961-01-10 00:00:00
    nowFile = var + '_' + filename[:-14] + "_" + str(timed[index])[0:10].replace('-', '_')
    return nowFile

def GetnetCDF4InfobyName(in_filename,var_name):
    '''
    Function to read the original file's projection
    '''
    #Open netCDF file
    src_ds = gdal.Open(in_filename)
    #warnings.filterwarnings('ignore')
    if src_ds is None:
        print('Open failed')
        sys.exit()

    if len(src_ds.GetSubDatasets()) > 1:
        #������ڶ������
        subdataset = 'NETCDF:"' + in_filename +'":' + var_name
        src_ds_sd = gdal.Open(subdataset)
        #warnings.filterwarnings('ignore')
        #��ȡ����
        NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
        xsize = src_ds_sd.RasterXSize
        ysize = src_ds_sd.RasterYSize
        GeoT = src_ds_sd.GetGeoTransform()#���ص�������
        Projection = osr.SpatialReference()
        # ��ȡ��������ϵͳ��Ϣ������ѡȡ��Ҫ�ĵ�������ϵͳ

        Projection.ImportFromEPSG(3857)  # �������������ϵΪ"WGS 84"��AUTHORITY["EPSG","4326"
        #Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())

        #�ر��������ݼ�
        src_ds_sd = None
        src_ds = None

        #��xrray��ȡ����
        xr_ensemble = xr.open_dataset(in_filename)
        data = xr_ensemble[var_name]
        data = np.ma.masked_array(data,mask=data==NDV,fill_value=NDV)


        return NDV,xsize,ysize,GeoT,Projection,data

    #Create GeoTiff image
def create_geotiff(suffix,Array,NDV,xsize,ysize,GeoT,Projection,nc_data,filename):
        '''
        �������д����µ�tiffͼ��
        Create a new tiff image from the array
        '''
        DataType = gdal_array.NumericTypeCodeToGDALTypeCode(Array.dtype) ###

        if type(DataType) != np.int:
            if DataType.startswith('gdal.GDT_') == False:
                DataType = eval('gdal.GDT_' + DataType)


        zsize = Array.shape[0]

        #create a driver
        driver = gdal.GetDriverByName('GTiff')
        #Set nans to the original No Data Value
        Array[np.isnan(Array)] = NDV

        #Set up the dataset with zsize
        #write each slice of the array along the zsize
        Array_shape_length = len(Array.shape)
        if Array_shape_length == 3:
            for i,image in enumerate(Array,1):
                now_file = get_time(i - 1,nc_data,var_name,filename)
                # NewFileName = outpath + var_name + netCDF4�ļ� + 1999_12_31 + .tif
                NewFileName = suffix + '\\' + now_file + '.tif'
                print('Outpath:',NewFileName,'\n')
                DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
                DataSet.SetGeoTransform(GeoT)
                DataSet.SetProjection(Projection.ExportToWkt())
                #�ж����ݵ�ά��
                #if Array_shape_length == 3:
                DataSet.GetRasterBand(1).WriteArray(image)
                #else:
                # ����Ϊn*7*11�ľ��󣬴˴���ȡ��Ϊ��0�㣬һ��n��
                #DataSet.GetRasterBand(1).WriteArray(image[0])
                DataSet.GetRasterBand(1).SetNoDataValue(NDV)
                DataSet.FlushCache()
                DataSet = None  # close the tif file
                print('The ' + now_file + '.tif' ' has finished!!!\n')
        else:
            for i,image in enumerate(Array,1):
                now_file = get_time(i - 1,nc_data,var_name,filename)
                # NewFileName = outpath + var_name + netCDF4�ļ� + 1999_12_31 + .tif
                NewFileName = suffix + '\\' + now_file + '.tif'
                print('Outpath:',NewFileName,'\n')
                DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
                DataSet.SetGeoTransform(GeoT)
                DataSet.SetProjection(Projection.ExportToWkt())
                #�ж����ݵ�ά��
                #����Ϊn*7*11�ľ��󣬴˴���ȡ��Ϊ��0�㣬һ��n��
                DataSet.GetRasterBand(1).WriteArray(image[0])
                DataSet.GetRasterBand(1).SetNoDataValue(NDV)
                DataSet.FlushCache()
                DataSet = None  # close the tif file
                print('The ' + now_file + '.tif' ' has finished!!!\n')

def GetOutPath(output_nc,filename,var):
    '''
    ��������ļ���
    Create output folder
    '''
    path = output_nc + '\\' + filename[:-14] + '_' + var  # ��������ļ���
    try:
        os.makedirs(path)
    except:
        pass
    return path

if __name__ == '__main__':
    #�����ļ��� /Input folder
    infile = r'C:\Users\Seven\Desktop\123'
    #����ļ��� /Output folder
    output_nc = r"C:\Users\Seven\Desktop\1234"
    lista = os.listdir(infile)
    for k in range(0, len(lista)):
        filename = lista[k]
        strPth = infile + '\\' + filename
        dataset = GetValueBynetCDF4(strPth)
        length = len(dataset.variables.keys())  # ��ȡ.nc�ļ��ı�ǩ
        for i in range(4, length):
            var_name = str(list(dataset.variables.keys())[i])
            output_path = GetOutPath(output_nc,filename, var_name)
            NDV,xsize,ysize,GeoT,Projection,data = GetnetCDF4InfobyName(strPth,var_name)
            create_geotiff(output_path,data,NDV,xsize,ysize,GeoT,Projection,dataset,filename)
            print('The ' + var_name +' has finished!\n')
            #warnings.filterwarnings('ignore')
        print('The ' + filename +' conversion to .tiff file is complete! \n')
    print('Done!\n')
