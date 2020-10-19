"""
my environments:
    python 3.8.5
    modules:
        netCDF4 1.5.4
        xarray 0.16.1
        numpy 1.19.1
        Gdal 3.1.3
注意事项：
	1、文件路径最好不要有中文路径
	2、输出文件夹随意
        	3. 输入文件夹只能存放.nc文件，不要有其他文件
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
    通过netCDF4打开.nc文件
    Open the .nc file through netCDF4
    '''
    nc_data = nc.Dataset(in_filename)
    return nc_data

def get_time(index,nc_data,var,filename):
    '''
    获取时间标签值
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
        #如果存在多个变量
        subdataset = 'NETCDF:"' + in_filename +'":' + var_name
        src_ds_sd = gdal.Open(subdataset)
        #warnings.filterwarnings('ignore')
        #读取变量
        NDV = src_ds_sd.GetRasterBand(1).GetNoDataValue()
        xsize = src_ds_sd.RasterXSize
        ysize = src_ds_sd.RasterYSize
        GeoT = src_ds_sd.GetGeoTransform()#返回地理坐标
        Projection = osr.SpatialReference()
        # 获取地理坐标系统信息，用于选取需要的地理坐标系统

        Projection.ImportFromEPSG(3857)  # 定义输出的坐标系为"WGS 84"，AUTHORITY["EPSG","4326"
        #Projection.ImportFromWkt(src_ds_sd.GetProjectionRef())

        #关闭整个数据集
        src_ds_sd = None
        src_ds = None

        #用xrray读取数据
        xr_ensemble = xr.open_dataset(in_filename)
        data = xr_ensemble[var_name]
        data = np.ma.masked_array(data,mask=data==NDV,fill_value=NDV)


        return NDV,xsize,ysize,GeoT,Projection,data

    #Create GeoTiff image
def create_geotiff(suffix,Array,NDV,xsize,ysize,GeoT,Projection,nc_data,filename):
        '''
        从数组中创建新的tiff图像
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
                # NewFileName = outpath + var_name + netCDF4文件 + 1999_12_31 + .tif
                NewFileName = suffix + '\\' + now_file + '.tif'
                print('Outpath:',NewFileName,'\n')
                DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
                DataSet.SetGeoTransform(GeoT)
                DataSet.SetProjection(Projection.ExportToWkt())
                #判读数据的维度
                #if Array_shape_length == 3:
                DataSet.GetRasterBand(1).WriteArray(image)
                #else:
                # 数据为n*7*11的矩阵，此代码取的为第0层，一共n层
                #DataSet.GetRasterBand(1).WriteArray(image[0])
                DataSet.GetRasterBand(1).SetNoDataValue(NDV)
                DataSet.FlushCache()
                DataSet = None  # close the tif file
                print('The ' + now_file + '.tif' ' has finished!!!\n')
        else:
            for i,image in enumerate(Array,1):
                now_file = get_time(i - 1,nc_data,var_name,filename)
                # NewFileName = outpath + var_name + netCDF4文件 + 1999_12_31 + .tif
                NewFileName = suffix + '\\' + now_file + '.tif'
                print('Outpath:',NewFileName,'\n')
                DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
                DataSet.SetGeoTransform(GeoT)
                DataSet.SetProjection(Projection.ExportToWkt())
                #判读数据的维度
                #数据为n*7*11的矩阵，此代码取的为第0层，一共n层
                DataSet.GetRasterBand(1).WriteArray(image[0])
                DataSet.GetRasterBand(1).SetNoDataValue(NDV)
                DataSet.FlushCache()
                DataSet = None  # close the tif file
                print('The ' + now_file + '.tif' ' has finished!!!\n')

def GetOutPath(output_nc,filename,var):
    '''
    创建输出文件夹
    Create output folder
    '''
    path = output_nc + '\\' + filename[:-14] + '_' + var  # 创建输出文件夹
    try:
        os.makedirs(path)
    except:
        pass
    return path

if __name__ == '__main__':
    #输入文件夹 /Input folder
    infile = r'C:\Users\Seven\Desktop\123'
    #输出文件夹 /Output folder
    output_nc = r"C:\Users\Seven\Desktop\1234"
    lista = os.listdir(infile)
    for k in range(0, len(lista)):
        filename = lista[k]
        strPth = infile + '\\' + filename
        dataset = GetValueBynetCDF4(strPth)
        length = len(dataset.variables.keys())  # 获取.nc文件的标签
        for i in range(4, length):
            var_name = str(list(dataset.variables.keys())[i])
            output_path = GetOutPath(output_nc,filename, var_name)
            NDV,xsize,ysize,GeoT,Projection,data = GetnetCDF4InfobyName(strPth,var_name)
            create_geotiff(output_path,data,NDV,xsize,ysize,GeoT,Projection,dataset,filename)
            print('The ' + var_name +' has finished!\n')
            #warnings.filterwarnings('ignore')
        print('The ' + filename +' conversion to .tiff file is complete! \n')
    print('Done!\n')
