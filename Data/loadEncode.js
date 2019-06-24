const axios = require('axios')
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
  'human': 'GRCh38',
  'mouse': 'mm10'
}

const encodeBasePath = 'encodeData'
const publicBasePath = 'publicData'

// This is the path of the `loadEncode.js` relative to the server scripts
const scriptPath = 'Data'

const outputPath = '../html/assets'

const matchTissueJsonFileName = 'matchingTissues.json'
const filterJsonFileName = 'filters.json'
const publicDataName = 'publicData.json'

const outputEncodeFile = 'encodeData.json'
const outputPublicFile = 'publicData.json'
const outputExperimentFile = 'experimentDict.json'

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24

const DELAY_MINIMUM = 500
const DELAY_MAXIMUM = 2000

const RETRY_TIMES = 5

/**
 * Workflow:
 *
 * *  Read ./matchedTissues.json for ENCODE Biosample IDs;
 *    matchedTissues.json should be an array of objects
 *    ```json
 *    { "human": ["<human_biosample_IDs>"], "mouse": ["<mouse_biosample_IDs>"] }
 *    ```
 * *  Read ./filter.json for an array of filters (experiment type, target);
 *    Filter object should be:
 *    ```json
 *    {
 *      ".id": "<filter_id>",
 *      ".label": "<human_readable_label>",
 *      ".isExpression": "<is_expression_value>",
 *      // ...
 *      // all properties from ENCODE search object that must match
 *      // sub-properties can be represented by using '.', e.g. 'target.label'
 *    }
 *    ```
 * *  For each biosample:
 *    *   Use ENCODE API to get the UUID of the biosample;
 *    *   Use ENCODE search API to get all experiments done on all matching
 *        biosamples;
 *    *   Filter the search results by filters;
 *    *   (TODO) Include archived experiments and replace them by the
 *        superseded experiments if exist, update biosample if necessary;
 *    *   Populate the biosample list with filtered experiment IDs;
 * *  For each biosample pair (matched tissue):
 *    *   Find the overlapping experiments of the biosample between the species
 *        and remove species-specific ones.
 * *  For each experiment:
 *    *   Use ENCODE API to get experiment information;
 *    *   Find the 'released' results;
 *    *   Find the most comprehensive (pooled, replicated) peaks and bigWig
 *        files
 *    *   Download those files to
 *        <ENCODE_File_ID>/<ENCODE_File_ID>.<Extension>,
 *        and add them to the experiment list.
 *
 * The result shall be a JSON file with all matching biosamples, each has an
 * 'experiment' property listing all relevant experiments.
 *
 * All experiment objects will include final released data file IDs and files
 * on the server.
 **/

/**
 * final result format
 * @property {object} filters - the filters that can be used to distinguish
 *    different experiments (for example, ChIP-Seq H3K4me3). Filters should be
 *    in the following format:
 *    ```json
 *    {
 *      "<unique_filter_id>": {
 *        ".id": "<unique_filter_id>",
 *        ".label": "<human_readable_label>",  // e.g. 'ChIP-Seq (H3K4me3)'
 *        // ...
 *        // all properties from ENCODE search object that must match
 *        // sub-properties can be represented by using '.',
 *        //   e.g. 'target.label'
 *      }
 *      // ...
 *    }
 *    ```
 * @property {Array<object>} matchedTissues - the matched tissues, should be:
 *    ```json
 *    [
 *      {
 *        ".term": "<biosample_term_name>",
 *        ".life_stage": "<life_stage>",
 *        ".label": "<human_readable_label>",
 *        "<species_name>": {
 *          "experiments": [
 *            {
 *              "biosampleId": "<ENCODE_biosample_ID>",
 *              "id": "<ENCODE_experiment_ID>",
 *              "filterId": "<unique_filter_id>"
 *            },
 *            // ...
 *          ]
 *        },
 *        // ...
 *      },
 *      // ...
 *    ]
 *    ```
 */
var encodeResult = {
  filters: {},
  matchedTissues: [],
  experiments: {}
}

var publicResult = {
  filters: {},
  matchedTissues: []
}

var expressionIdDict = {}

/**
 * Dictionary mapping experiment IDs to data files and limited meta info.
 * Should be:
 * ```json
 * {
 *   "<ENCODE_experiment_ID>": {
 *     "id": "<ENCODE_experiment_ID>",
 *     "experiment_type": "<filter label>",
 *     "biosample_label": "<biosample label>",
 *     "peak_id": "<ENCODE_ID for peak file>",
 *     "peak_file": "<path of peak file on server>",
 *     "bigwig_id": "<ENCODE_ID for bigwig file>",
 *     "bigwig_file": "<path of bigwig file on server>"
 *   },
 *   // ...
 * }
 * ```
 */
var experimentDict = {}

async function randomDelay () {
  let time =
    Math.floor(Math.random() * (DELAY_MAXIMUM - DELAY_MINIMUM) + DELAY_MINIMUM)
  console.log('Delaying ' + time + ' ms.')
  return new Promise((resolve, reject) => setTimeout(resolve, time))
}

async function axiosGetWrapper (url, criteria) {
  let retries = 0
  let resultObj
  while (retries < RETRY_TIMES) {
    try {
      resultObj = await axios.get(url)
      if (resultObj.data && (!criteria || criteria(resultObj.data))) {
        return resultObj
      }
      throw resultObj
    } catch (err) {
      retries++
      console.log('Error requesting ' + url)
      await randomDelay()
    }
  }
  throw resultObj
}

async function axiosRequestWrapper (...args) {
  let retries = 0
  let resultObj
  while (retries < RETRY_TIMES) {
    try {
      resultObj = await axios.request(...args)
      return resultObj
    } catch (err) {
      retries++
      console.log('Error in request.')
      await randomDelay()
    }
  }
  throw resultObj
}

/**
 * Get the ENCODE JSON Object of the biosample
 * @async
 * @param {string} type ENCODE type: ('biosamples', 'experiments')
 * @param {string} id ENCODE ID
 * @returns {Promise<object>} ENCODE object
 */
async function getEncodeObject (type, id, criteria) {
  // console.log('https://www.encodeproject.org/' + type +
  //   '/' + id + '?format=json')
  return axiosGetWrapper('https://www.encodeproject.org/' + type +
    '/' + id + '?format=json', criteria
  ).then(response => response.data
  ).catch(err => {
    console.log('Error getting ENCODE object: ' +
      'https://www.encodeproject.org/' + type +
      '/' + id + '?format=json')
    throw err
  })
}

/**
 * Get the ENCODE JSON Object of the biosample
 * @async
 * @param {string} biosampleId ENCODE biosample ID
 * @param {boolean} [notReleaseOnly] Set to `true` to include non-released
 *    datasets
 * @returns {Promise<Array<string>>} An array including all experiments
 */
async function getExperimentsUsingBiosample (biosampleObj, notReleaseOnly) {
  // first get the uuid of the biosampleObj
  let uuid = biosampleObj.uuid
  let queryUrlBase = 'https://www.encodeproject.org/search/?' +
    'format=json&type=Experiment&replicates.library.biosample.uuid=' + uuid +
    '&limit=all'
  let status = ''
  if (!notReleaseOnly) {
    // release only
    status = 'status=released'
  }
  try {
    let resultObj = await axiosGetWrapper(
      queryUrlBase + (status ? '&' + status : '')
    )
    console.log('Found ' + resultObj.data['@graph'].length +
      ' experiments for biosample ' + biosampleObj.accession + '.')
    return resultObj.data
  } catch (err) {
    if (err.response && err.response.status === 404) {
      // no result
      console.log('No results for biosample ' + biosampleObj.accession +
        '. UUID = ' + biosampleObj.uuid)
      return null
    }
    console.log('Error getting ENCODE object for biosample ' +
      biosampleObj.accession + ': ' +
      queryUrlBase + (status ? '&' + status : ''))
    console.log(err)
    throw err
  }
}

const getNestedObject = function getNestedObject (nestedObj, pathArr) {
  return pathArr.reduce(
    (obj, key) => ((obj && obj[key] !== undefined) ? obj[key] : undefined),
    nestedObj)
}

/**
 * Filter the ENCODE search results
 * @param {object} searchObj Serach object returned from ENCODE
 * @param {Array<object>} filters Filters, each filter should be an object with
 *    all required string properties. Multiple filters are combined with 'OR'.
 *    properties should be matched
 * @param {string} species species of the biosample ('human' or 'mouse')
 * @returns {Array<object>} An array of experimentEntry object:
 * ```json
 * {
 *    "id": "<Experiment_IDs>",
 *    "filter": "<the filter object that this experiment passed>",
 *    "expressionId": "<the Experiment ID for expression of the same biosample"
 * }
 * ```
 * @async
 */
async function filterSearchResults (searchObj, filters, species) {
  let result = []
  if (searchObj && Array.isArray(searchObj['@graph'])) {
    for (let expIndex = 0; expIndex < searchObj['@graph'].length; expIndex++) {
      let graphObj = searchObj['@graph'][expIndex]
      for (let filterIndex = 0; filterIndex < filters.length; filterIndex++) {
        let filter = filters[filterIndex]
        let criteriaMatched = true
        for (let key in filter) {
          if (!key.startsWith('.') && filter.hasOwnProperty(key)) {
            // key may be something like 'a.b', in which case we need to descend
            // into graphObj
            if ((filter[key] + '').toLowerCase() !==
              (getNestedObject(graphObj, key.split('.')) + '').toLowerCase()
            ) {
              criteriaMatched = false
              break
            }
          }
        }
        if (criteriaMatched) {
          // Passed filter, construct an object with Experiment IDs and filter ID
          // that are passed
          result.push({
            filterId: filter['.id'],
            id: graphObj.accession,
            isExpression: filter['.isExpression'] || false
          })
          if (filter['.isExpression'] &&
            await verifyExpressionId(graphObj.accession, references[species])
          ) {
            result.expressionId = result.expressionId || graphObj.accession
          }
          break
        }
      }
    }
  }
  if (result.expressionId) {
    result.forEach(item => (item.expressionId = result.expressionId))
  }
  return result
}

/**
 * Verify the expression ID:
 * *   Use ENCODE API to get experiment information;
 * *   Find the 'released' results;
 * *   Find the 'gene quantifications' files
 * *   Save the ID of those files (for the final downloading) and add them
 *     to the expressionIdDict list.
 * @async
 * @param {string} experimentId ENCODE ID for the experiment
 * @param {string} assembly the assembly of the file
 * @returns {object} returns an experiment object (see "Final result format")
 */
async function verifyExpressionId (experimentId, assembly) {
  let experimentObj =
    await getEncodeObject('experiments', experimentId, data => data.files)
  let availableFiles = experimentObj.files.filter(file => (
    file.assembly === assembly && file.output_category !== 'raw data' &&
    file.status === 'released' &&
    file.output_type === 'gene quantifications'
  ))

  let expressionFiles = availableFiles.filter(file => (
    file.output_type === 'gene quantifications' && (
      !file.superseded_by || !file.superseded_by.length
    )
  ))

  let experimentResult = null

  // Then just use bedFiles[0] and bwFiles[0] to populate the peak and bigwig
  // fields of the experiment

  if (expressionFiles.length) {
    experimentResult = { id: experimentId }
    experimentResult.experiment_file_ids =
      expressionFiles.map(file => file.accession)
    experimentResult.experiment_file_objs = expressionFiles

    expressionIdDict[experimentId] = experimentResult
  }

  return experimentResult
}

/**
 * For each experiment:
 * *   Use ENCODE API to get experiment information;
 * *   Find the 'released' results;
 * *   Find the most comprehensive (pooled, replicated) peaks and bigWig
 *     files
 * *   Save the ID of those files (for the final downloading) and add them
 *     to the experiment list.
 * @async
 * @param {string} experimentId ENCODE ID for the experiment
 * @param {string} assembly the assembly of the file
 * @param {object} [additionalMeta] additional meta data that needs to be added
 *    to the experiment entry.
 * @returns {object} returns an experiment object (see "Final result format")
 */
async function processExperimentFileObj (
  experimentId, assembly, additionalMeta
) {
  additionalMeta = additionalMeta || {}
  let experimentObj = await getEncodeObject(
    'experiments', experimentId, data => data.files)
  let availableFiles = experimentObj.files.filter(file => (
    file.assembly === assembly && file.output_category !== 'raw data' &&
    file.status === 'released' && (
      file.file_format.toLowerCase() === 'bed' ||
      file.file_format.toLowerCase() === 'bigwig'
    )
  ))

  // For bedfiles, the ones not deriving any other file should be used
  let bedFiles = availableFiles.filter(file => (
    file.file_format.toLowerCase() === 'bed' &&
    file.file_type.toLowerCase().includes('narrowpeak')
  ))
  bedFiles = bedFiles.filter(
    file => bedFiles.every(
      fileToCompare => (
        !fileToCompare.derived_from ||
        fileToCompare.derived_from.every(
          derivedTerm => !derivedTerm.includes(file.accession)
        )
      )
    )
  )

  // For bwFiles, anything involving p-value is not included
  // Anything deriving other files are not used
  // Anything whose `derived_from` is a subset of something else are not used
  let bwFiles = availableFiles.filter(file => (
    file.file_format.toLowerCase() === 'bigwig' &&
    !file.output_type.includes('p-value')
  ))
  bwFiles = bwFiles.filter(
    file => bwFiles.every(
      fileToCompare => (
        !fileToCompare.derived_from ||
        fileToCompare.derived_from.every(
          derivedTerm => !derivedTerm.includes(file.accession)
        )
      )
    )
  ).filter(
    file => bwFiles.every(
      fileToCompare => (
        !file.derived_from || !fileToCompare.derived_from ||
        file['@id'] === fileToCompare['@id'] ||
        file.derived_from.some(
          derivedId => fileToCompare.derived_from.indexOf(derivedId) < 0
        )
      )
    )
  )

  let experimentResult = null

  // Then just use bedFiles[0] and bwFiles[0] to populate the peak and bigwig
  // fields of the experiment

  if (bedFiles.length) {
    experimentResult = Object.assign({ id: experimentId }, additionalMeta)
    experimentResult.peak_id = bedFiles[0].accession
    experimentResult.peak_file_obj = bedFiles[0]
    if (bwFiles.length) {
      experimentResult.bigwig_id = bwFiles[0].accession
      experimentResult.bigwig_file_obj = bwFiles[0]
    }
  }

  return experimentResult
}

/**
 * Finish the "for each biosample part"
 * For each biosample:
 * *   Use ENCODE API to get the UUID of the biosample;
 * *   Use ENCODE search API to get all experiments done on all matching
 *     biosamples;
 * *   Filter the search results by filters;
 * *   (TODO) Include archived experiments and replace them by the
 *     superseded experiments if exist, update biosample if necessary;
 * *   Return filtered experiment IDs;
 * @param {string} biosampleId ENCODE biosample ID
 * @param {Array<object>} filters Filters
 * @param {object} tissueResult tissueResult object
 *    (to fill '.term', '.life_stage' and '.label' element only)
 * @param {string} species species of the biosample ('human' or 'mouse')
 * @returns {Array<object>} An array of experimentEntry object:
 * ```json
 * {
 *    "id": "<Experiment_IDs>",
 *    "filter": "<the filter object that this experiment passed>"
 * }
 * ```
 */
async function processSingleBiosample (
  biosampleId, filters, tissueResult, species
) {
  let biosampleObj = await getEncodeObject('biosamples', biosampleId)

  if (!tissueResult.hasOwnProperty('.label')) {
    // initialize the properties with this biosample
    tissueResult['.term'] = biosampleObj.biosample_ontology.term_name
    tissueResult['.life_stage'] = biosampleObj.life_stage
    tissueResult['.label'] = (
      tissueResult['.life_stage'] ? tissueResult['.life_stage'] + ' ' : ''
    ) + tissueResult['.term']
  }
  try {
    // TODO: set `notReleaseOnly` to `true` and update superseded experiments
    let searchObj = await getExperimentsUsingBiosample(biosampleObj)
    // TODO: replace superseded experiments

    return await filterSearchResults(searchObj, filters, species)
  } catch (err) {
    console.log('Error processSingleBiosample (' + biosampleId + ').')
    return []
  }
}

/**
 * For each biosample pair (matched tissue):
 * *   Find the overlapping experiments of the biosample between the species
 *     and remove species-specific ones.
 * @param {object} tissueObj tissue object (see Workflow)
 * @param {Array<object>} filters filters (see Workflow)
 * @returns {object} tissueResult object (see encodeResult.matchedTissues)
 */
async function processSingleTissue (tissueObj, filters) {
  // For each species in tissue:
  //    use `processSingleBiosample` to get a list of experiments
  let tissueResult = {}
  if (tissueObj['.label']) {
    tissueResult['.label'] = tissueObj['.label']
  }
  let speciesFilterSet = []
  for (let species in tissueObj) {
    if (!species.startsWith('__') && !species.startsWith('.') &&
      tissueObj.hasOwnProperty(species)
    ) {
      let experimentEntryArray = await Promise.all(
        tissueObj[species].map(
          biosampleId => processSingleBiosample(
            biosampleId, filters, tissueResult, species
          ))
      )
      // construct a preliminary list of experiments
      tissueResult[species] = {}
      experimentEntryArray.some(entryArray => {
        if (entryArray[0] && entryArray[0].expressionId) {
          tissueResult[species].expressionId = entryArray[0].expressionId
          return true
        }
        return false
      })

      tissueResult[species].experiments = experimentEntryArray.reduce(
        (prev, current, index) => prev.concat(
          current.map(experimentEntry => ({
            biosampleId: tissueObj[species][index],
            id: experimentEntry.id,
            filterId: experimentEntry.filterId,
            expressionId: experimentEntry.expressionId ||
              tissueResult[species].expressionId || null,
            isExpression: experimentEntry.isExpression
          }))
        ), []
      )

      tissueResult[species].experiments =
        tissueResult[species].experiments.filter(expEntry => {
          if (expEntry.isExpression) {
            return false
          } else {
            delete expEntry.isExpression
            return true
          }
        })
      console.log(tissueResult[species].experiments.length +
        ' filtered experiment(s) for "' + tissueResult['.label'] +
        '" in species ' + species)

      // build a filter set (for debug purposes)
      let currSpeciesFilterSet = new Set()
      tissueResult[species].experiments.forEach(
        experimentObj => currSpeciesFilterSet.add(experimentObj.filterId)
      )

      speciesFilterSet.push(currSpeciesFilterSet)
    }
  }

  // find the overlapping experiments (filterId that appear in all species)
  return tissueResult
}

/**
 * Filter the experiment that does not have proper data available.
 * Proper data means:
 * *  Have data mapped to the desired assembly
 *    (determined by `references[speciesName]`)
 * *  Have peak files called and available (in "released" state)
 *
 * @param {object} tissueResult The tissueResult object
 * @param {object} filterDict Filter object (dictionary) for filter labels
 * @returns
 */
async function filterInvalidDataSingleTissue (tissueResult, filterDict) {
  for (let species in tissueResult) {
    if (!species.startsWith('__') && !species.startsWith('.') &&
      tissueResult.hasOwnProperty(species)
    ) {
      // First filter out duplicates
      let expIdSet = new Set()
      tissueResult[species].experiments =
        tissueResult[species].experiments.filter(
          experimentObj => {
            if (expIdSet.has(experimentObj.id)) {
              return false
            }
            expIdSet.add(experimentObj.id)
            return true
          }
        )
      // filter this preliminary list to remove the ones without proper data
      // files
      let experimentResultArray = await Promise.all(
        tissueResult[species].experiments.map(
          experimentObj => processExperimentFileObj(
            experimentObj.id, references[species], {
              'experiment_type': filterDict[experimentObj.filterId]['.label'],
              'biosample_label': tissueResult['.label'],
              'expressionId': tissueResult[species].expressionId || null
            })
        )
      )
      tissueResult[species].experiments =
        tissueResult[species].experiments.filter(
          (experimentObj, index) => !!experimentResultArray[index]
        )
      tissueResult[species]._expResultDict = experimentResultArray
        .filter(result => !!result)
        .reduce((prev, curr) => {
          prev[curr.id] = curr
          return prev
        }, {})
    }
  }
  return tissueResult
}

function filterMismatchedSingleTissue (tissueResult) {
  // reconstruct speciesFilterSet
  let speciesFilterSet = []
  for (let species in tissueResult) {
    if (!species.startsWith('__') && !species.startsWith('.') &&
      tissueResult.hasOwnProperty(species)
    ) {
      let currSpeciesFilterSet = new Set()
      tissueResult[species].experiments.forEach(experimentObj =>
        currSpeciesFilterSet.add(experimentObj.filterId))
      speciesFilterSet.push(currSpeciesFilterSet)
    }
  }

  // find the overlapping experiments (filterId that appear in all species)
  let intersectSet = new Set(speciesFilterSet.reduce(
    (setArray, currSet) => setArray.filter(id => currSet.has(id)),
    speciesFilterSet[0] ? [...speciesFilterSet[0]] : []
  ))

  // remove species-unique experiments
  for (let species in tissueResult) {
    if (!species.startsWith('__') && !species.startsWith('.') &&
      tissueResult.hasOwnProperty(species)
    ) {
      tissueResult[species].experiments =
        tissueResult[species].experiments.filter(experiment =>
          intersectSet.has(experiment.filterId)
        )
    }
  }
  return tissueResult
}

async function downloadFileObject (fileObj, filePath) {
  filePath = filePath || encodeBasePath || '.'
  let encodeId = fileObj.accession
  let extension = fileObj.href.split('.').pop()
  let needGunzip = false
  if (extension === 'gz') {
    // needs to gunzip
    needGunzip = true
    extension = fileObj.href.split('.').slice(-2)[0]
  }
  let fileName = path.format({
    dir: path.format({ dir: filePath, base: encodeId }),
    name: encodeId,
    ext: '.' + extension
  })
  let needUpdateFile = await fsMkdirPromise(
    path.format({ dir: filePath, base: encodeId })
  ).then(() => true
  ).catch(async err => {
    if (err.code === 'EEXIST') {
      // already exists, use fs.stat to test if file is too old
      return fsStatPromise(fileName)
        .then(fileStat =>
          (Date.now() - fileStat.ctime > (90 * MILLISECONDS_IN_A_DAY))
        )
        .catch(err => {
          if (err.code === 'ENOENT') {
            return true
          } else {
            console.log(err)
            throw err
          }
        })
    } else {
      console.log(err)
      throw err
    }
  }).catch(err => {
    console.log(err)
    throw err
  })

  if (needUpdateFile) {
    console.log('Downloading experiment file from: ' + fileObj.href)
  } else {
    console.log('Skipping existing file from: ' + fileObj.href)
  }
  // Download file with axios
  return needUpdateFile
    ? axiosRequestWrapper({
      url: 'https://www.encodeproject.org' + fileObj.href,
      responseType: 'arraybuffer'
    }).catch(err => {
      console.log('Error getting ENCODE file: ' +
        'https://www.encodeproject.org' + fileObj.href)
      throw err
    }).then(async response => {
      let data = response.data
      if (needGunzip) {
        // needs to gunzip
        data = await zlibGunzipPromise(data)
      }
      await fsWritefilePromise(fileName, data)
      return fileName
    })
    : fileName
}

async function downloadExperimentFiles (experimentResult) {
  if (experimentResult.peak_file_obj) {
    experimentResult.peak_file = path.format({
      dir: scriptPath,
      base: await downloadFileObject(experimentResult.peak_file_obj)
    })
    delete experimentResult.peak_file_obj
  }
  if (experimentResult.bigwig_file_obj) {
    experimentResult.bigwig_file = path.format({
      dir: scriptPath,
      base: await downloadFileObject(experimentResult.bigwig_file_obj)
    })
    delete experimentResult.bigwig_file_obj
  }
  if (experimentResult.expressionId) {
    let expressionEntry = expressionIdDict[experimentResult.expressionId]
    if (expressionEntry.experiment_file_objs) {
      expressionEntry.expression_files = []
      for (let i = 0; i < expressionEntry.experiment_file_objs.length; i++) {
        let fileObj = expressionEntry.experiment_file_objs[i]
        let base = await downloadFileObject(fileObj)
        expressionEntry.expression_files[i] = path.format({
          dir: scriptPath,
          base: base
        })
      }
    }
    experimentResult.expression_files = expressionEntry.expression_files
  }
  return experimentResult
}

async function populateTissueExperiments (
  tissueResult, experimentResultDict
) {
  for (let species in tissueResult) {
    if (!species.startsWith('__') && !species.startsWith('.') &&
      tissueResult.hasOwnProperty(species)
    ) {
      await Promise.all(
        tissueResult[species].experiments.map(
          experimentObj => {
            experimentResultDict[experimentObj.id] =
              tissueResult[species]._expResultDict[experimentObj.id]
            return downloadExperimentFiles(
              tissueResult[species]._expResultDict[experimentObj.id]
            )
          }
        )
      )
    }
  }
  return tissueResult
}

var readTissuePromise = fsReadfilePromise(
  path.format({ dir: encodeBasePath, base: matchTissueJsonFileName }),
  'utf8'
).then(result => {
  return JSON.parse(result)
}).catch(err => {
  console.log(err)
})

var readFilterPromise = fsReadfilePromise(
  path.format({ dir: encodeBasePath, base: filterJsonFileName }),
  'utf8'
).then(result => {
  return JSON.parse(result)
}).catch(err => {
  console.log(err)
})

var readPublicPromise = fsReadfilePromise(
  path.format({ dir: publicBasePath, base: publicDataName }),
  'utf8'
).then(result => {
  return JSON.parse(result)
}).catch(err => {
  console.log(err)
})

Promise.all([readTissuePromise, readFilterPromise])
  .then(results => {
    let [tissues, filters] = results
    encodeResult.filters = filters.reduce((prev, curr) => {
      if (!curr['.isExpression']) {
        prev[curr['.id']] = curr
      }
      return prev
    }, {})
    // first get raw tissue results
    return Promise.all(
      tissues.map(tissue => processSingleTissue(tissue, filters))
    )
  }).then(tissueResults => {
    console.log('===== filterInvalidDataSingleTissue =====')
    return Promise.all(tissueResults.map(tissue => {
      return filterInvalidDataSingleTissue(tissue, encodeResult.filters)
    }))
  }).then(tissueResults => {
    console.log('===== filterMismatchedSingleTissue =====')
    return tissueResults.map(tissue => filterMismatchedSingleTissue(tissue))
  }
  ).then(tissueResults => {
    console.log('===== populateTissueExperiments =====')
    // Now all tissues should be ready
    encodeResult.matchedTissues = tissueResults
    // populate and download experiments
    return Promise.all(encodeResult.matchedTissues.map(
      tissueResult =>
        populateTissueExperiments(tissueResult, encodeResult.experiments)
    ))
  }).then(() => {
    console.log('===== Clearing _expResultDict =====')
    encodeResult.matchedTissues.forEach(tissueResult => {
      for (let species in tissueResult) {
        if (!species.startsWith('__') && !species.startsWith('.') &&
          tissueResult.hasOwnProperty(species)
        ) {
          delete tissueResult[species]._expResultDict
        }
      }
    })
    console.log('===== Separating encodeResult.experiments =====')
    experimentDict = encodeResult.experiments
    delete encodeResult.experiments
  }).then(() => {
    // load public datasets and process the experiment Dict
    console.log('===== Parsing publicData.experiments =====')
    return readPublicPromise.then(publicDataObj => {
      // adding filter label and biosample label to publicDataObj
      publicDataObj.matchedTissues.forEach(tissue => {
        for (let species in tissue) {
          if (!species.startsWith('__') && !species.startsWith('.') &&
            tissue.hasOwnProperty(species)
          ) {
            let experiments = tissue[species].experiments
            experiments.forEach(experiment => {
              let experimentDictEntry =
                publicDataObj.experiments[experiment.id]
              experimentDictEntry.experiment_type =
                publicDataObj.filters[experiment.filterId]['.label']
              experimentDictEntry.biosample_label = tissue['.label']
              experimentDictEntry.peak_file = path.join(
                scriptPath, publicBasePath, experimentDictEntry.peak_file
              )
              experimentDictEntry.bigwig_file = path.join(
                scriptPath, publicBasePath, experimentDictEntry.bigwig_file
              )
              if (experimentDictEntry.expression_files) {
                experimentDictEntry.expression_files =
                  experimentDictEntry.expression_files.map(
                    expFile => path.join(
                      scriptPath, publicBasePath, expFile
                    )
                  )
              }
            })
          }
        }
      })
      experimentDict = Object.assign(
        experimentDict, publicDataObj.experiments
      )
      delete publicDataObj.experiments
      publicResult = publicDataObj
    })
  }).then(() => {
    // write encodeResults to a JSON file
    console.log('===== Writing to JSON =====')
    return fsWritefilePromise(
      path.format({ dir: outputPath, base: outputEncodeFile }),
      JSON.stringify(encodeResult, null, 2)
    ).then(() => fsWritefilePromise(
      path.format({ dir: outputPath, base: outputPublicFile }),
      JSON.stringify(publicResult, null, 2)
    )).then(() => fsWritefilePromise(
      path.format({ dir: outputPath, base: outputExperimentFile }),
      JSON.stringify(experimentDict, null, 2)
    ))
  })
