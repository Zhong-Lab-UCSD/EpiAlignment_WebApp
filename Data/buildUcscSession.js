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
const scriptPath = 'Data'

const outputPath = '../html/assets'

const DATA_URL_BASE = 'https://epialign.ucsd.edu/'

const matchTissueJsonFileName = 'matchingTissues.json'
const filterJsonFileName = 'filters.json'
const publicDataName = 'publicData.json'

const outputEncodeFile = 'encodeData.json'
const outputPublicFile = 'publicData.json'
const outputExperimentFile = 'experimentDict.json'

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24

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
        }
        tissue[species].experiments.forEach(experiment => {
          let ucscEntry = ''
          let sampleDesc = experiment.biosampleId
            ? ' (Biosample ID: ' + experiment.biosampleId + ')'
            : ' (Dataset ID: ' + experiment.id + ')'
          let expId = experiment.id
          let filterDesc = filters[experiment.filterId]['.label']
          let title = filterDesc + ' on ' + label + sampleDesc
          let dataUrl
          // find file in expDict
          if (expDict[expId].bigwig_file) {
            dataUrl = DATA_URL_BASE + expDict[expId].bigwig_file
            ucscEntry = 'track type=bigWig name="' + filterDesc + ' ' +
              (experiment.biosampleId || experiment.id) + '" description="' +
              title + '" ' +
              'visibility=0 bigDataUrl=' + dataUrl
          } else if (expDict[expId].peak_file) {
            dataUrl = DATA_URL_BASE + expDict[expId].peak_file
            ucscEntry = 'track type=bed name="' + filterDesc + ' ' +
              (experiment.biosampleId || experiment.id) + '" description="' +
              title + '" ' +
              'visibility=2 bigDataUrl=' + dataUrl
          }
          if (ucscEntry) {
            result[species].push(ucscEntry)
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
    console.log(species)
    result[species].forEach(entry => console.log(entry))
  }
})
