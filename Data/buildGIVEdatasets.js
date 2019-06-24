const fs = require('fs')
const util = require('util')
const zlib = require('zlib')
const path = require('path')

const fsReadfilePromise = util.promisify(fs.readFile)
const fsWritefilePromise = util.promisify(fs.writeFile)
const fsMkdirPromise = util.promisify(fs.mkdir)
const fsStatPromise = util.promisify(fs.stat)
const zlibGunzipPromise = util.promisify(zlib.gunzip)

const references = {
  'human': 'hg38',
  'mouse': 'mm10'
}

// This is the path of the `loadEncode.js` relative to the server scripts
const dataRootPath = '/var/www/epialignment/'

const outputPath = '../html/assets'

const DATA_URL_BASE = 'https://epialign.ucsd.edu/'

const matchTissueJsonFileName = 'matchingTissues.json'
const filterJsonFileName = 'filters.json'
const publicDataName = 'publicData.json'

const outputEncodeFile = 'encodeData.json'
const outputPublicFile = 'publicData.json'
const outputExperimentFile = 'experimentDict.json'

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24

var mysqlPass = process.argv[2]
const giveLoadingPath = '/home/xcao3/Github/GIVe/GIVE-Toolbox/'
const removeDataScript = 'remove_data.sh'
const createGroupScript = 'add_trackGroup.sh'
const loadBedScript = 'add_track_bed.sh'
const loadBigWigScript = 'add_track_bigWig.sh'
const mysqlParam = ' -u cpwriter -p ' + mysqlPass + ' -g epialign'
const jsonSettingsBigWig = ' "height": 70, "upperProportion": 0.01'

var encodeFilePromise = fsReadfilePromise(path.format({
  dir: outputPath,
  base: outputEncodeFile
}), 'utf8').then(data => JSON.parse(data))

var publicFilePromise = fsReadfilePromise(path.format({
  dir: outputPath,
  base: outputPublicFile
}), 'utf8').then(data => JSON.parse(data))

var expDictPromise = fsReadfilePromise(path.format({
  dir: outputPath,
  base: outputExperimentFile
}), 'utf8').then(data => JSON.parse(data))

function parseResult (resultObj, expDict, result) {
  result = result || {}
  let filters = resultObj.filters
  resultObj.matchedTissues.forEach(tissue => {
    let label = tissue['.label']
    for (let species in tissue) {
      if (tissue.hasOwnProperty(species) &&
        !species.startsWith('.') && !species.startsWith('__')
      ) {
        if (!result.hasOwnProperty(species)) {
          result[species] = []
          let removeDataCommand = giveLoadingPath + removeDataScript + mysqlParam
          removeDataCommand += ' -r ' + references[species]
          result[species].push(removeDataCommand)
          
          let createGroupCommand = giveLoadingPath + createGroupScript + mysqlParam
          createGroupCommand += ' -r ' + references[species] +
            ' -l "EpiAlignment Preset Data" -o 50 -s 0' 
          result[species].push(createGroupCommand)
        }
        tissue[species].experiments.forEach(experiment => {
          let loadingCommandBed, loadingCommandBigWig
          let sampleDesc = experiment.biosampleId
            ? ' (Biosample ID: ' + experiment.biosampleId + ')'
            : ' (Dataset ID: ' + experiment.id + ')'
          let expId = experiment.id
          let filterDesc = filters[experiment.filterId]['.label']
          let shortLabelFilterShort = filterDesc.replace(/ChIP-Seq \((.*?)\)/, '$1')
          let filePath
          // find file in expDict
          if (expDict[expId].bigwig_file) {
            loadingCommandBigWig = giveLoadingPath + loadBigWigScript + mysqlParam
            let shortLabel = shortLabelFilterShort + ' ' + label
            let longLabel = filterDesc + ' on ' + label + sampleDesc
            filePath = dataRootPath + expDict[expId].bigwig_file
            loadingCommandBigWig += ' -r ' + references[species] +
              ' -t ' + expId.toLowerCase() + '-b -l "' + longLabel + '" -s "' +
              shortLabel + '" -v full -a 1 -o 10 -f ' + filePath + ' -m ' +
              '\'"cellType": "' + label + '", "dataType": "ChIP-Seq", "trackFeature": "' +
              shortLabelFilterShort + '", "experimentId": "' + experiment.id + '",' +
              jsonSettingsBigWig + '\''
          }
          if (expDict[expId].peak_file) {
            loadingCommandBed = giveLoadingPath + loadBedScript + mysqlParam
            let shortLabel = shortLabelFilterShort + ' Peaks ' + label
            let longLabel = filterDesc + ' Peaks on ' + label + sampleDesc
            filePath = dataRootPath + expDict[expId].peak_file
            loadingCommandBed += ' -r ' + references[species] +
              ' -t ' + expId.toLowerCase() + '-p -l "' + longLabel + '" -s "' +
              shortLabel + '" -v dense -o 10 -f ' + filePath + ' -m ' +
              '\'"cellType": "' + label + '", "dataType": "ChIP-Seq", "trackFeature": "' +
              shortLabelFilterShort + '", "experimentId": "' + experiment.id + '"\''
          }
          if (loadingCommandBed) {
            result[species].push(loadingCommandBed)
          }
          if (loadingCommandBigWig) {
            result[species].push(loadingCommandBigWig)
          }
        })
      }
    }
  })
}

var result = {}

Promise.all([
  Promise.all([encodeFilePromise, expDictPromise]
  ).then(resultArray => {
    let [encodeData, expDict] = resultArray
    parseResult(encodeData, expDict, result)
  }),
  Promise.all([publicFilePromise, expDictPromise]
  ).then(resultArray => {
    let [encodeData, expDict] = resultArray
    parseResult(encodeData, expDict, result)
  })
]).then(() => {
  for (let species in references) {
    // console.log(species)
    console.log()
    result[species].forEach(entry => console.log(entry))
  }
})
