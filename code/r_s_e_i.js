 
var roi = table.geometry().bounds()
 
Map.centerObject(roi, 5);

function remove_water(img) {
    var year = img.get('year')
    var jrc_year = ee.ImageCollection('JRC/GSW1_4/YearlyHistory')
            .filterDate('2020-01-01', '2020-12-31')
            .first()
            .clip(roi)
            .select('waterClass')
            .select('waterClass')
            .reproject('EPSG:4326',null,1000)
    // jrc_year: JRC Yearly Water Classification History, v1.3
    var Mask = jrc_year
                       // 掩膜掉大片永久水体
                       .eq(3)
    // 此时Mask中value值有1 ， 0 ， masked，把masked转化为 0
    Mask = Mask.unmask(0).not()
    return img.updateMask(Mask)
}
 
function removeCloud(image){
  var qa = image.select('BQA')
  var cloudMask = qa.bitwiseAnd(1 << 4).eq(0)
  var cloudShadowMask = qa.bitwiseAnd(1 << 8).eq(0)
  var valid = cloudMask.and(cloudShadowMask)
  return image.updateMask(valid)
}
 
var L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
           .filterBounds(roi)
           .filterDate('2020-01-01', '2020-12-31')
           .filterMetadata('CLOUD_COVER', 'less_than',50)
           .map(function(image){
                    return image.set('year', ee.Image(image).date().get('year'))                           
                  })
           .map(removeCloud)

 
var L8imgList = ee.List([])
for(var a = 2020; a < 2021; a++){
   var img = L8.filterMetadata('year', 'equals', a).median().clip(roi)
   var L8img = img.set('year', a)
   L8imgList = L8imgList.add(L8img)
 }
 
var L8imgCol = ee.ImageCollection(L8imgList)
                 .map(function(img){
                      return img.clip(roi)
                   })
                 .map(remove_water)
                 
L8imgCol = L8imgCol.map(function(img){
  
  // Wet
  var Wet = img.expression('B*(0.1509) + G*(0.1973) + R*(0.3279) + NIR*(0.3406) + SWIR1*(-0.7112) + SWIR2*(-0.4572)',{
       'B': img.select(['B2']),
       'G': img.select(['B3']),
       'R': img.select(['B4']),
       'NIR': img.select(['B5']),
       'SWIR1': img.select(['B6']),
       'SWIR2': img.select(['B7'])
     })   
  img = img.addBands(Wet.rename('WET'))
  
  
  // NDVI
  var ndvi = img.normalizedDifference(['B5', 'B4']);
  img = img.addBands(ndvi.rename('NDVI'))
  
  
  // lst 直接采用MODIS产品
  var lst = ee.ImageCollection('MODIS/006/MOD11A1').map(function(img){
                return img.clip(roi)
           })
           .filterDate('2020-01-01', '2020-12-31')
  
  var year = img.get('year')
  lst=lst.filterDate(ee.String(year).cat('-01-01'),ee.String(year).cat('-12-31')).select(['LST_Day_1km', 'LST_Night_1km']);
      
  // reproject主要是为了确保分辨率为1000
  var img_mean=lst.mean().reproject('EPSG:4326',null,1000);
  //print(img_mean.projection().nominalScale())
  
  img_mean = img_mean.expression('((Day + Night) / 2)',{
      'Day': img_mean.select(['LST_Day_1km']),
      'Night': img_mean.select(['LST_Night_1km']),
       })
  img = img.addBands(img_mean.rename('LST'))
  
  
  // ndbsi = ( ibi + si ) / 2
  var ibi = img.expression('(2 * SWIR1 / (SWIR1 + NIR) - (NIR / (NIR + RED) + GREEN / (GREEN + SWIR1))) / (2 * SWIR1 / (SWIR1 + NIR) + (NIR / (NIR + RED) + GREEN / (GREEN + SWIR1)))', {
      'SWIR1': img.select('B6'),
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'GREEN': img.select('B3')
    })
  var si = img.expression('((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
      'SWIR1': img.select('B6'),
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'BLUE': img.select('B2')
    }) 
  var ndbsi = (ibi.add(si)).divide(2)
  return img.addBands(ndbsi.rename('NDBSI'))
})
 
 
var bandNames = ['NDVI', "NDBSI", "WET", "LST"]
L8imgCol = L8imgCol.select(bandNames)
 
// 归一化
var img_normalize = function(img){
      var minMax = img.reduceRegion({
            reducer:ee.Reducer.minMax(),
            geometry: roi,
            scale: 1000,
            maxPixels: 10e13,
        })
      var year = img.get('year')
      var normalize  = ee.ImageCollection.fromImages(
            img.bandNames().map(function(name){
                  name = ee.String(name);
                  var band = img.select(name);
                  return band.unitScale(ee.Number(minMax.get(name.cat('_min'))), ee.Number(minMax.get(name.cat('_max'))));
                    
              })
        ).toBands().rename(img.bandNames()).set('year', year);
        return normalize;
}
var imgNorcol  = L8imgCol.map(img_normalize);
imgNorcol = imgNorcol.select(bandNames)


Export.image.toDrive({
  image: imgNorcol.select(["NDVI","NDBSI","WET","LST"]).mosaic(),
  description: 'Nile_rsei_2020',
  folder: 'R_S_E_I',
  scale: 30,
  region: roi,
  maxPixels: 1E13,
});
 