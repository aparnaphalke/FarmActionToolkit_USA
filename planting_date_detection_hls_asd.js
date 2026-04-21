/*******************************************************
 FAcT-USA | Alabama | 2024
 HLS-Based Corn Planting Date Detection Workflow

 Purpose:
   Estimate corn planting date (day of year; DOY) for a
   selected Alabama Agricultural Statistics District (ASD)
   using Harmonized Landsat-Sentinel (HLS) NDVI time series
   and a phenology-guided greenup detection approach.

 Important implementation note:
   This exact workflow was repeated separately for each
   Alabama ASD zone (110, 120, 130, 140, 150, 160), using
   zone-specific corn and cropland mask inputs. Final
   statewide products were created after running the same
   method independently for each ASD and combining outputs.

 Inputs:
   1) ASD-specific corn extent map
   2) ASD-specific cropland mask
   3) NASA HLSL30 v002
   4) NASA HLSS30 v002

 Method summary:
   - Restrict analysis to mapped corn pixels within cropland
   - Build regular 6-day NDVI composites
   - Linearly interpolate missing observations
   - Smooth NDVI trajectories
   - Detect peak spring NDVI slope as greenup timing
   - Estimate planting date as greenup DOY minus SHIFT_DAYS
   - Retain only plausible Alabama corn planting dates

 Output:
   Corn planting date raster in day of year (DOY)

 Author:
   Aparna R. Phalke

 Repository:
   FarmActionToolkit_USA
********************************************************/


// ================================
// 0) USER SETTINGS
// ================================

var year = 2024;
var SELECT_STASD = 110;

// Temporal settings
var STEP_DAYS = 6;     // regular HLS compositing interval
var SHIFT_DAYS = 12;   // planting proxy = greenup DOY - SHIFT_DAYS

// Plausible Alabama corn planting window
var CORN_DOY_MIN = 88;
var CORN_DOY_MAX = 162;

// Minimum seasonal NDVI amplitude threshold to suppress
// weak or noisy trajectories
var MIN_NDVI_AMP = 0.15;

// ASD-specific assets
var CORN_ASSET = 'users/aparnapm16/AL_ASD_RF_Corn_STASD_110_Y2024';
var CROPMASK_ASSET = 'users/aparnapm16/AL_ASD_RF_CropNonCrop_STASD_110_Y2024';


// ================================
// 1) LOAD ASD BOUNDARY
// ================================

var asdAll = ee.FeatureCollection('projects/ee-aparnapm16/assets/Alabama_ASD_2012_20m')
  .filter(ee.Filter.eq('STATE', '01'));

var asdOne = asdAll.filter(ee.Filter.eq('STASD_N', SELECT_STASD));
var region = asdOne.geometry();

Map.setOptions('HYBRID');
Map.centerObject(region, 9);
Map.addLayer(
  asdOne.style({
    color: '00ffff',
    width: 3,
    fillColor: '00000000'
  }),
  {},
  'ASD ' + SELECT_STASD,
  true
);


// ================================
// 2) LOAD CORN AND CROPLAND MASKS
// ================================

// Corn extent: assume 1 = corn, 0 = non-corn
var cornRaw = ee.Image(CORN_ASSET).clip(region);
var cornMask = cornRaw.eq(1).selfMask().rename('corn');

// Cropland mask: assume 1 = crop, 0 = non-crop
var cropRaw = ee.Image(CROPMASK_ASSET).clip(region);
var cropMask = cropRaw.eq(1).selfMask().rename('crop');

// Final analysis mask = corn pixels constrained to mapped cropland
var analysisMask = cornMask.updateMask(cropMask);

Map.addLayer(cropMask, {palette: ['00ff00']}, 'Cropland mask', false);
Map.addLayer(cornMask, {palette: ['ffff00']}, 'Corn mask', false);
Map.addLayer(analysisMask, {palette: ['ff0000']}, 'Corn ∩ Cropland mask', true);


// ================================
// 3) DEFINE ANALYSIS DATE RANGE
// ================================

var dateStart = ee.Date.fromYMD(year, 1, 1);
var dateEnd = ee.Date.fromYMD(year, 12, 31);


// ================================
// 4) HLS PREPARATION FUNCTIONS
// ================================

function getFmaskBit(img, bit) {
  return img.select('Fmask').bitwiseAnd(1 << bit).neq(0);
}

// HLS Fmask bits:
// 1 = cloud
// 2 = adjacent cloud/shadow
// 3 = shadow
// 4 = snow/ice
// 5 = water
function maskHLS(img) {
  var cloud = getFmaskBit(img, 1);
  var adj = getFmaskBit(img, 2);
  var shadow = getFmaskBit(img, 3);
  var snow = getFmaskBit(img, 4);
  var water = getFmaskBit(img, 5);

  var good = cloud.not()
    .and(adj.not())
    .and(shadow.not())
    .and(snow.not())
    .and(water.not());

  return img.updateMask(good);
}

function prepL30(img) {
  var sr = img.select(['B2', 'B3', 'B4', 'B5'], ['BLUE', 'GREEN', 'RED', 'NIR'])
    .multiply(0.0001);
  return sr.copyProperties(img, ['system:time_start']);
}

function prepS30(img) {
  var sr = img.select(['B2', 'B3', 'B4', 'B8'], ['BLUE', 'GREEN', 'RED', 'NIR'])
    .multiply(0.0001);
  return sr.copyProperties(img, ['system:time_start']);
}

function toNDVI(img) {
  return img.normalizedDifference(['NIR', 'RED'])
    .rename('NDVI')
    .toFloat()
    .clamp(-1, 1)
    .copyProperties(img, ['system:time_start']);
}


// ================================
// 5) LOAD HLS NDVI SCENES
// ================================

var l30 = ee.ImageCollection('NASA/HLS/HLSL30/v002')
  .filterBounds(region)
  .filterDate(dateStart, dateEnd)
  .map(maskHLS)
  .map(prepL30)
  .map(toNDVI);

var s30 = ee.ImageCollection('NASA/HLS/HLSS30/v002')
  .filterBounds(region)
  .filterDate(dateStart, dateEnd)
  .map(maskHLS)
  .map(prepS30)
  .map(toNDVI);

var ndviScenes = l30.merge(s30).select('NDVI');

print('HLS NDVI scenes:', ndviScenes.size());


// ================================
// 6) BUILD REGULAR 6-DAY NDVI COMPOSITES
// ================================

var totalDays = dateEnd.difference(dateStart, 'day');
var offsets = ee.List.sequence(0, totalDays, STEP_DAYS);

var regular = ee.ImageCollection.fromImages(
  offsets.map(function(d) {
    d = ee.Number(d);
    var ws = dateStart.advance(d, 'day');
    var we = ws.advance(STEP_DAYS, 'day');

    var win = ndviScenes.filterDate(ws, we);

    var ndvi = ee.Image(ee.Algorithms.If(
      win.size().gt(0),
      win.mean().rename('NDVI').toFloat().clamp(-1, 1),
      ee.Image.constant(0).rename('NDVI').toFloat().clamp(-1, 1)
        .updateMask(ee.Image.constant(0))
    ));

    ndvi = ndvi.updateMask(analysisMask).clip(region);

    return ndvi.set('system:time_start', ws.millis());
  })
).sort('system:time_start');

print('Regular composite count:', regular.size());


// ================================
// 7) LINEAR INTERPOLATION OF REGULAR COMPOSITES
// ================================

function linearInterpolateRegular(ic) {
  ic = ee.ImageCollection(ic).sort('system:time_start');
  var list = ic.toList(ic.size());
  var n = ic.size();

  var fwd = ee.List.sequence(0, n.subtract(1)).iterate(function(i, acc) {
    i = ee.Number(i);
    acc = ee.Dictionary(acc);

    var prevImg = ee.Image(acc.get('prevImg'));
    var prevT = ee.Number(acc.get('prevT'));
    var out = ee.List(acc.get('out'));

    var base = ee.Image(list.get(i));
    var img = base.select('NDVI').toFloat();
    var t = ee.Number(base.get('system:time_start'));
    var has = img.mask().reduce(ee.Reducer.anyNonZero()).gt(0);

    var newPrevImg = ee.Image(ee.Algorithms.If(has, img, prevImg));
    var newPrevT = ee.Number(ee.Algorithms.If(has, t, prevT));

    var pack = newPrevImg.rename('ndvi_prev')
      .addBands(ee.Image.constant(newPrevT).rename('t_prev'))
      .toFloat();

    return ee.Dictionary({
      prevImg: newPrevImg,
      prevT: newPrevT,
      out: out.add(pack)
    });
  }, ee.Dictionary({
    prevImg: ee.Image.constant(0).rename('NDVI').toFloat()
      .updateMask(ee.Image.constant(0)),
    prevT: ee.Number(-1),
    out: ee.List([])
  }));

  fwd = ee.List(ee.Dictionary(fwd).get('out'));

  var bwd = ee.List.sequence(0, n.subtract(1)).reverse().iterate(function(i, acc) {
    i = ee.Number(i);
    acc = ee.Dictionary(acc);

    var nextImg = ee.Image(acc.get('nextImg'));
    var nextT = ee.Number(acc.get('nextT'));
    var out = ee.List(acc.get('out'));

    var base = ee.Image(list.get(i));
    var img = base.select('NDVI').toFloat();
    var t = ee.Number(base.get('system:time_start'));
    var has = img.mask().reduce(ee.Reducer.anyNonZero()).gt(0);

    var newNextImg = ee.Image(ee.Algorithms.If(has, img, nextImg));
    var newNextT = ee.Number(ee.Algorithms.If(has, t, nextT));

    var pack = newNextImg.rename('ndvi_next')
      .addBands(ee.Image.constant(newNextT).rename('t_next'))
      .toFloat();

    return ee.Dictionary({
      nextImg: newNextImg,
      nextT: newNextT,
      out: out.add(pack)
    });
  }, ee.Dictionary({
    nextImg: ee.Image.constant(0).rename('NDVI').toFloat()
      .updateMask(ee.Image.constant(0)),
    nextT: ee.Number(-1),
    out: ee.List([])
  }));

  bwd = ee.List(ee.Dictionary(bwd).get('out')).reverse();

  var outList = ee.List.sequence(0, n.subtract(1)).map(function(i) {
    i = ee.Number(i);

    var base = ee.Image(list.get(i));
    var ts = ee.Number(base.get('system:time_start'));
    var orig = base.select('NDVI').toFloat();

    var prevPack = ee.Image(fwd.get(i));
    var nextPack = ee.Image(bwd.get(i));

    var ndviPrev = prevPack.select('ndvi_prev');
    var tPrev = prevPack.select('t_prev');
    var ndviNext = nextPack.select('ndvi_next');
    var tNext = nextPack.select('t_next');

    var hasObs = orig.mask().reduce(ee.Reducer.anyNonZero()).gt(0);
    var ok = tPrev.neq(-1).and(tNext.neq(-1)).and(tNext.neq(tPrev));

    var w = ee.Image.constant(ts).subtract(tPrev).divide(tNext.subtract(tPrev));

    var interp = ndviPrev.add(ndviNext.subtract(ndviPrev).multiply(w))
      .rename('NDVI')
      .toFloat();

    var filled = ee.Image(ee.Algorithms.If(
      hasObs,
      orig,
      ee.Image(ee.Algorithms.If(ok, interp, orig))
    ));

    return filled.rename('NDVI')
      .toFloat()
      .clamp(-1, 1)
      .set('system:time_start', ts);
  });

  return ee.ImageCollection.fromImages(outList).sort('system:time_start');
}

var regularInterp = linearInterpolateRegular(regular);


// ================================
// 8) SMOOTH NDVI TIME SERIES
// ================================

function smoothMovingAverage(ic, windowSteps) {
  ic = ee.ImageCollection(ic).sort('system:time_start');
  var list = ic.toList(ic.size());
  var n = ic.size();
  var half = ee.Number(windowSteps).divide(2).floor();

  var out = ee.List.sequence(0, n.subtract(1)).map(function(i) {
    i = ee.Number(i);
    var start = i.subtract(half).max(0);
    var end = i.add(half).min(n.subtract(1));

    var win = ee.ImageCollection.fromImages(
      ee.List.sequence(start, end).map(function(j) {
        return ee.Image(list.get(j)).select('NDVI').toFloat();
      })
    );

    var sm = win.mean().rename('NDVI').toFloat();
    var center = ee.Image(list.get(i));
    var ts = ee.Number(center.get('system:time_start'));

    return sm.set('system:time_start', ts);
  });

  return ee.ImageCollection.fromImages(out).sort('system:time_start');
}

var ndviSmooth = smoothMovingAverage(regularInterp, 5);


// ================================
// 9) NDVI AMPLITUDE FILTER
// ================================

var maxNDVI = ndviSmooth.max();
var minNDVI = ndviSmooth.min();
var ampNDVI = maxNDVI.subtract(minNDVI).rename('ampNDVI');

var goodAmpMask = ampNDVI.gte(MIN_NDVI_AMP).selfMask();

Map.addLayer(
  ampNDVI,
  {min: 0, max: 0.6, palette: ['white', 'yellow', 'green']},
  'NDVI amplitude',
  false
);

Map.addLayer(
  goodAmpMask,
  {palette: ['00ffff']},
  'Good amplitude mask',
  false
);


// ================================
// 10) ADD DOY BAND TO SMOOTHED NDVI
// ================================

var ndviCol = ndviSmooth.map(function(img) {
  var doy = ee.Image.constant(
    ee.Date(img.get('system:time_start')).getRelative('day', 'year').add(1)
  ).rename('DOY').toInt16();

  return img.select('NDVI')
    .toFloat()
    .addBands(doy)
    .set('system:time_start', img.get('system:time_start'));
});


// ================================
// 11) COMPUTE NDVI SLOPE
// ================================

var ndviList = ndviCol.toList(ndviCol.size());
var n = ndviCol.size();

var slopeCol = ee.ImageCollection.fromImages(
  ee.List.sequence(1, n.subtract(1)).map(function(i) {
    i = ee.Number(i);

    var prev = ee.Image(ndviList.get(i.subtract(1))).select('NDVI');
    var curr = ee.Image(ndviList.get(i)).select('NDVI');
    var tMillis = ee.Number(ee.Image(ndviList.get(i)).get('system:time_start'));

    var slope = curr.subtract(prev).divide(STEP_DAYS).rename('slope').toFloat();

    var doy = ee.Image.constant(
      ee.Date(tMillis).getRelative('day', 'year').add(1)
    ).rename('DOY').toInt16();

    return slope.addBands(doy).set('system:time_start', tMillis);
  })
);


// ================================
// 12) DETECT SPRING GREENUP DOY
// ================================

var greenStart = ee.Date.fromYMD(year, 2, 1);
var greenEnd = ee.Date.fromYMD(year, 7, 31);

var maxSlope = slopeCol.filterDate(greenStart, greenEnd).qualityMosaic('slope');

var greenupDOY = maxSlope.select('DOY')
  .rename('greenupDOY')
  .updateMask(analysisMask)
  .updateMask(goodAmpMask);


// ================================
// 13) CONVERT GREENUP DOY TO PLANTING DOY
// ================================

var plantingDOY = greenupDOY
  .subtract(SHIFT_DAYS)
  .rename('plantingDOY')
  .toInt16()
  .updateMask(analysisMask)
  .clip(region);

// Retain only planting dates within the plausible Alabama corn planting window
var plantingDOY_clamped = plantingDOY.updateMask(
  plantingDOY.gte(CORN_DOY_MIN).and(plantingDOY.lte(CORN_DOY_MAX))
);


// ================================
// 14) APPLY LIGHT SPATIAL SMOOTHING
// ================================

var plantingDOY_smooth = plantingDOY_clamped.focal_median({
  radius: 1,
  units: 'pixels',
  kernelType: 'square'
}).rename('plantingDOY');


// ================================
// 15) VISUALIZATION
// ================================

var plantViz = {
  min: CORN_DOY_MIN,
  max: CORN_DOY_MAX,
  palette: [
    '3b0f70',
    '2c7fb8',
    '41b6c4',
    'a1dab4',
    'ffffcc',
    'fed976',
    'fd8d3c',
    'e31a1c'
  ]
};

Map.addLayer(greenupDOY, plantViz, 'Greenup DOY', false);
Map.addLayer(plantingDOY, plantViz, 'Planting DOY raw', false);
Map.addLayer(plantingDOY_clamped, plantViz, 'Planting DOY clamped', false);
Map.addLayer(plantingDOY_smooth, plantViz, 'Planting DOY smooth', true);


// ================================
// 16) SUMMARY STATISTICS
// ================================

print(
  'Planting DOY summary:',
  plantingDOY_smooth.reduceRegion({
    reducer: ee.Reducer.minMax().combine({
      reducer2: ee.Reducer.mean(),
      sharedInputs: true
    }),
    geometry: region,
    scale: 30,
    maxPixels: 1e13,
    tileScale: 4
  })
);


// ================================
// 17) EXPORTS
// ================================

Export.image.toDrive({
  image: plantingDOY_smooth,
  description: 'AL_ASD110_CornPlantingDOY_HLS_2024',
  folder: 'GEE_Exports',
  fileNamePrefix: 'AL_ASD110_CornPlantingDOY_HLS_2024',
  region: region,
  scale: 30,
  maxPixels: 1e13,
  crs: 'EPSG:4326'
});

Export.image.toDrive({
  image: plantingDOY,
  description: 'AL_ASD110_CornPlantingDOY_HLS_2024_raw',
  folder: 'GEE_Exports',
  fileNamePrefix: 'AL_ASD110_CornPlantingDOY_HLS_2024_raw',
  region: region,
  scale: 30,
  maxPixels: 1e13,
  crs: 'EPSG:4326'
});


// ================================
// 18) LEGEND
// ================================

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 12px',
    backgroundColor: 'white'
  }
});

legend.add(ui.Label({
  value: 'Corn Planting Date (DOY)',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 6px 0'
  }
}));

var palette = [
  '3b0f70',
  '2c7fb8',
  '41b6c4',
  'a1dab4',
  'ffffcc',
  'fed976',
  'fd8d3c',
  'e31a1c'
];

var labels = [
  'Early (Mar end ~ DOY 88)',
  'Early Apr',
  'Mid Apr',
  'Late Apr',
  'Early May',
  'Mid May',
  'Late May',
  'Late (Jun ~ DOY 162)'
];

for (var i = 0; i < palette.length; i++) {
  var colorBox = ui.Label('', {
    backgroundColor: '#' + palette[i],
    padding: '8px',
    margin: '0 0 4px 0'
  });

  var description = ui.Label(labels[i], {
    margin: '0 0 4px 6px'
  });

  legend.add(ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  }));
}

Map.add(legend);
