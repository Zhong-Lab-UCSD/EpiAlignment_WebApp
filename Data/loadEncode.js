const axios = require('axios')
const fs = require('fs')
const util = require('util')
const zlib = require('zlib')

const fsReadfilePromise = util.promisify(fs.readFile)
const fsWritefilePromise = util.promisify(fs.writeFile)
const fsMkdirPromise = util.promisify(fs.mkdir)
const zlibGunzipPromise = util.promisify(zlib.gunzip)

const references = {
  'human': 'GRCh38',
  'mouse': 'mm10'
}

const outputFile = 'encodeData.json'

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
 *      ".id": "<unique_filter_id>",
 *      ".label": "<human_readable_label>",  // e.g. 'ChIP-Seq (H3K4me3)'
 *      // ...
 *      // all properties from ENCODE search object that must match
 *      // sub-properties can be represented by using '.', e.g. 'target.label'
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
 *              "biosampleID": "<ENCODE_biosample_ID>",
 *              "id": "<ENCODE_experiment_ID>",
 *              "filterID": "<unique_filter_id>"
 *            },
 *            // ...
 *          ]
 *        },
 *        // ...
 *      },
 *      // ...
 *    ]
 *    ```
 * @property {object} experiments - the detailed experiments, should be:
 *    ```json
 *    {
 *      "<ENCODE_experiment_ID>": {
 *        "id": "<ENCODE_experiment_ID>",
 *        "peak_id": "<ENCODE_ID for peak file>",
 *        "peak_file": "<path of peak file on server>",
 *        "bigwig_id": "<ENCODE_ID for bigwig file>",
 *        "bigwig_file": "<path of bigwig file on server>"
 *      },
 *      // ...
 *    }
 *    ```
 */
var finalResult = {
  filters: {},
  matchedTissues: {},
  experiments: {}
}

/**
 * Get the ENCODE JSON Object of the biosample
 * @async
 * @param {string} type ENCODE type: ('biosamples', 'experiments')
 * @param {string} id ENCODE ID
 * @returns {Promise<object>} ENCODE object
 */
async function getEncodeObject (type, id) {
  return axios.get('https://www.encodeproject.org/' + type +
    '/' + id + '?format=json')
}

/**
 * Get the ENCODE JSON Object of the biosample
 * @async
 * @param {string} biosampleID ENCODE biosample ID
 * @param {boolean} [notReleaseOnly] Set to `true` to include non-released
 *    datasets
 * @returns {Promise<Array<string>>} An array including all experiments
 */
async function getExperimentsUsingBiosample (biosampleObj, notReleaseOnly) {
  // first get the uuid of the biosampleObj
  let uuid = biosampleObj.uuid
  let queryUrlBase = 'https://www.encodeproject.org/search/?' +
    'type=Experiment&replicates.library.biosample.uuid=' + uuid
  let status = ''
  if (!notReleaseOnly) {
    // release only
    status = 'status=released'
  }
  let resultObj = await axios.get(queryUrlBase + (status ? '&' + status : ''))
  return resultObj
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
 * @returns {Array<object>} An array of experimentEntry object:
 * ```json
 * {
 *    "id": "<Experiment_IDs>",
 *    "filter": "<the filter object that this experiment passed>"
 * }
 * ```
 */
function filterSearchResults (searchObj, filters) {
  let result = []
  searchObj['@graph'].forEach(
    graphObj => filters.some(filter => {
      for (let key in filter) {
        if (!key.startsWith('.') && filter.hasOwnProperty(key)) {
          // key may be something like 'a.b', in which case we need to descend
          // into graphObj
          if ((filter[key] + '').toLowerCase() !==
            (getNestedObject(graphObj, key.split('.')) + '').toLowerCase()
          ) {
            return false
          }
        }
      }
      // Passed filter, construct an object with Experiment IDs and filter ID
      // that are passed
      result.push({
        filterId: filter['.id'],
        id: graphObj.accession
      })
      return true
    })
  )
  return result
}

/**
 * For each experiment:
 * *   Use ENCODE API to get experiment information;
 * *   Find the 'released' results;
 * *   Find the most comprehensive (pooled, replicated) peaks and bigWig
 *     files
 * *   Download those files and add them to the experiment list.
 * @param {string} experimentId ENCODE ID for the experiment
 * @param {string} assembly the assembly of the file
 * @returns {object} returns an experiment object (see "Final result format")
 */
async function processExperimentFileObj (experimentId, assembly) {
  let experimentObj = await getEncodeObject('experiments', experimentId)
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

  // For bwFiles, anything involving p-value is not included for now
  let bwFiles = availableFiles.filter(file => (
    file.file_format.toLowerCase() === 'bigwig' &&
    !file.output_type.includes('p-value')
  ))

  let experimentResult = null

  // Then just use bedFiles[0] and bwFiles[0] to populate the peak and bigwig
  // fields of the experiment

  if (bedFiles.length) {
    experimentResult = {
      id: experimentId
    }
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
 * @param {string} biosampleID ENCODE biosample ID
 * @param {object} filters Filters
 * @param {object} tissueResult tissueResult object
 *    (to fill '.term', '.life_stage' and '.label' element only)
 * @returns {Array<object>} An array of experimentEntry object:
 * ```json
 * {
 *    "id": "<Experiment_IDs>",
 *    "filter": "<the filter object that this experiment passed>"
 * }
 * ```
 */
async function processSingleBiosample (biosampleID, filters, tissueResult) {
  let biosampleObj = await getEncodeObject('biosamples', biosampleID)

  if (!tissueResult.hasOwnProperty('.label')) {
    // initialize the properties with this biosample
    tissueResult['.term'] = biosampleObj.biosample_term_name
    tissueResult['.life_stage'] = biosampleObj.life_stage
    tissueResult['.label'] = (
      tissueResult['.life_stage'] ? tissueResult['.life_stage'] + ' ' : ''
    ) + tissueResult['.term']
  }
  // TODO: set `notReleaseOnly` to `true` and update superseded experiments
  let searchObj = await getExperimentsUsingBiosample(biosampleObj)
  // TODO: replace superseded experiments

  return filterSearchResults(searchObj, filters)
}

/**
 * For each biosample pair (matched tissue):
 * *   Find the overlapping experiments of the biosample between the species
 *     and remove species-specific ones.
 * @param {object} tissueObj tissue object (see Workflow)
 * @param {object} filters filters (see Workflow)
 * @returns {object} tissueResult object (see finalResult.matchedTissues)
 */
async function processSingleTissue (tissueObj, filters) {
  // For each species in tissue:
  //    use `processSingleBiosample` to get a list of experiments
  let tissueResult = {}
  let speciesFilterSet = []
  for (let species in tissueObj) {
    if (tissueObj.hasOwnProperty(species)) {
      let experimentEntryArray = await Promise.all(
        tissueObj[species].map(
          biosampleID => processSingleBiosample(
            biosampleID, filters, tissueResult
          ))
      )
      // construct a preliminary list of experiments
      tissueResult[species] = {
        experiments: experimentEntryArray.reduce(
          (prev, current, index) => prev.concat(
            current.map(experimentEntry => ({
              biosampleID: tissueObj[species][index],
              id: experimentEntry.id,
              filterID: experimentEntry.filter
            }))
          ), []
        )
      }

      // build a filter set (for debug purposes)
      let currSpeciesFilterSet = new Set()
      tissueResult[species].experiments.forEach(
        experimentObj => currSpeciesFilterSet.add(experimentObj.filterID)
      )

      speciesFilterSet.push(currSpeciesFilterSet)
    }
  }

  // find the overlapping experiments (filterID that appear in all species)
  console.log('Raw species filter set: ')
  console.log(speciesFilterSet)
  return tissueResult
}

async function filterInvalidDataSingleTissue (tissueResult) {
  for (let species in tissueResult) {
    if (tissueResult.hasOwnProperty(species)) {
      // filter this preliminary list to remove the ones without proper data
      // files
      let experimentResultArray = await Promise.all(
        tissueResult[species].experiments.map(
          experimentObj => processExperimentFileObj(
            experimentObj.id, references[species])
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

function filterMismatchedSingleTissue (tissueResult, filters) {
  // reconstruct speciesFilterSet
  let speciesFilterSet = []
  for (let species in tissueResult) {
    if (tissueResult.hasOwnProperty(species)) {
      tissueResult[species].experiments.forEach(experimentObj =>
        speciesFilterSet.add(experimentObj.filterID))
    }
  }

  // find the overlapping experiments (filterID that appear in all species)
  console.log('Processed species filter set before matching: ')
  console.log(speciesFilterSet)

  let intersectSet = new Set(speciesFilterSet.reduce(
    (setArray, currSet) => setArray.filter(id => currSet.has(id)),
    speciesFilterSet[0] ? [...speciesFilterSet[0]] : []
  ))

  // remove species-unique experiments
  for (let species in tissueResult) {
    tissueResult[species].experiments =
      tissueResult[species].experiments.filter(experiment =>
        intersectSet.has(experiment.filterID)
      )
  }
  return tissueResult
}

async function downloadFileObject (fileObj, path) {
  path = path || '.'
  let encodeId = fileObj.accession
  await fsMkdirPromise(path + '/' + encodeId).catch(err => {
    if (err.code === 'EEXIST') {
      // already exists, don't care
      return null
    } else {
      console.log(err)
      throw err
    }
  })
  // Download file with axios
  return axios.request({
    url: 'https://www.encodeproject.org' + fileObj.href,
    responseType: 'arraybuffer'
  }).then(async response => {
    let extension = fileObj.href.split('.').pop()
    let data = response.data
    if (extension === 'gz') {
      // needs to gunzip
      data = await zlibGunzipPromise(data)
      extension = fileObj.href.split('.').slice(-2)[0]
    }
    let fileName = path + '/' + encodeId + '/' + encodeId + '.' + extension
    await fsWritefilePromise(fileName, data)
    return fileName
  })
}

async function downloadExperimentFiles (experimentResult) {
  if (experimentResult.peak_file_obj) {
    experimentResult.peak_file =
      await downloadFileObject(experimentResult.peak_file_obj)
    delete experimentResult.peak_file_obj
  }
  if (experimentResult.bigwig_file_obj) {
    experimentResult.bigwig_file =
      await downloadFileObject(experimentResult.bigwig_file_obj)
    delete experimentResult.bigwig_file_obj
  }
  return experimentResult
}

async function populateTissueExperiments (
  tissueResult, experimentResultDict
) {
  for (let species in tissueResult) {
    if (tissueResult.hasOwnProperty(species)) {
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

var readTissuePromise = fsReadfilePromise('matchedTissues.json', 'utf8')
  .then(result => {
    return JSON.parse(result)
  })
  .catch(err => {
    console.log(err)
  })

var readFilterPromise = fsReadfilePromise('filter.json', 'utf8')
  .then(result => {
    return JSON.parse(result)
  })
  .catch(err => {
    console.log(err)
  })

Promise.all([readTissuePromise, readFilterPromise])
  .then(results => {
    let [tissues, filters] = results
    finalResult.filters = filters
    // first get raw tissue results
    return Promise.all(
      tissues.map(tissue => processSingleTissue(tissue, filters))
    )
  }).then(tissueResults => {
    return Promise.all(
      tissueResults.map(tissue => filterInvalidDataSingleTissue(tissue))
    )
  }).then(tissueResults =>
    tissueResults.map(tissue => filterMismatchedSingleTissue(tissue))
  ).then(tissueResults => {
    // Now all tissues should be ready
    finalResult.matchedTissues = tissueResults
    // populate and download experiments
    return Promise.all(finalResult.matchedTissues.map(
      tissueResult =>
        populateTissueExperiments(tissueResult, finalResult.experiments)
    ))
  }).then(() => {
    // write finalResults to a JSON file
    return fsWritefilePromise(outputFile, JSON.stringify(finalResult))
  })
