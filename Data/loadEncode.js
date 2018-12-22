const axios = require('axios')
const fs = require('fs')
const util = require('util')
const zlib = require('zlib')

const fsReadfilePromise = util.promisify(fs.readFile)
const fsWritefilePromise = util.promisify(fs.writeFile)

const references = {
  "human": "GRCh38",
  "mouse": "mm10"
}

/**
 * Workflow:
 *
 * *  Read matchedTissues.json for ENCODE Biosample IDs;
 *    matchedTissues.json should be an array of objects
 *    ```json
 *    { "human": ["<human_biosample_IDs>"], "mouse": ["<mouse_biosample_IDs>"] }
 *    ```
 * *  Read filter.json for an array of filters (experiment type, target);
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
 *    *   Download those files and add them to the experiment list.
 *
 * The result shall be a JSON file with all matching biosamples, each has an
 * 'experiment' property listing all relevant experiments.
 *
 * All experiment objects will include final released data file IDs and files
 * on the server.
 **/

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
 * @property {object} matchedTissues - the matched tissues, should be:
 *    ```json
 *    {  
 *      ".term": "<biosample_term_name>",  
 *      ".life_stage": "<life_stage>",  
 *      ".label": "<human_readable_label>",  
 *      "<species_name>": {  
 *        "experiments": [  
 *          {  
 *            "biosampleID": "<ENCODE_biosample_ID>",  
 *            "experimentID": "<ENCODE_experiment_ID>",  
 *            "filterID": "<unique_filter_id>"  
 *          },  
 *          // ...  
 *        ]  
 *      },  
 *      // ...  
 *    }
 *    ```
 * @property {object} experiments - the detailed experiments, should be:
 *    ```json
 *    {  
 *      "<ENCODE_experiment_ID>": {  
 *        ".peak_id": "<ENCODE_ID for peak file>",  
 *        ".peak_file": "<path of peak file on server>",  
 *        ".bigwig_id": "<ENCODE_ID for bigwig file>",  
 *        ".bigwig_file": "<path of bigwig file on server>"  
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
    '/' + id +'?format=json')
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
    ((obj, key) => (obj && obj[key] !== undefined) ? obj[key] : undefined),
    nestedObj)
}

/**
 * Filter the ENCODE search results
 * @param {object} searchObj Serach object returned from ENCODE
 * @param {Array<object>} filters Filters, each filter should be an object with
 *    all required string properties. Multiple filters are combined with 'OR'.
 *    properties should be matched
 * @returns {Array<object>} An array of object:
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
          if ((filter[key] + '').toLowerCase()
            !== (getNestedObject(graphObj, key.split('.')) + '').toLowerCase()
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
 * @returns {Array<object>} An array of object:
 * ```json
 * {
 *    "id": "<Experiment_IDs>",
 *    "filter": "<the filter object that this experiment passed>"
 * }
 * ```
 */
async function processSingleBiosample (biosampleID, filters) {
  let biosampleObj = await getEncodeObject('biosamples', biosampleID)

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
 */
async function processSingleTissue (tissueObj, filters) {
  // For each species in tissue:
  //    use `processSingleBiosample` to get a list of experiments
  let tissueResult = {}
  let speciesFilterSet = []
  for (let species in tissueObj) {
    if (tissueObj.hasOwnProperty(species)) {
      let experimentArray = await Promise.all(
        tissueObj[species].map(
          biosampleID => processSingleBiosample(biosampleID, filters))
      )
      // construct a preliminary list of experiments
      let currSpeciesFilterSet = new Set()
      tissueResult[species] = experimentArray.reduce((prev, current, index) => {
        prev.concat(current.map(experimentObj => {
          currSpeciesFilterSet.add(experimentObj.filter)
          return {
            biosampleID: tissueObj[species][index],
            experimentID: experimentObj.id,
            filterID: experimentObj.filter
          }
        }))
      }, [])
      speciesFilterSet.push(currSpeciesFilterSet)
    }
  }

  // find the overlapping experiments (filterID that appear in all species)  
  console.log('Raw species filter set: ')
  console.log(speciesFilterSet)
  return tissueResult
}

function filterSingleTissue (tissueResult, filters) {
  // reconstruct speciesFilterSet
  let speciesFilterSet = []
  for (let species in tissueResult) {
    if (tissueResult.hasOwnProperty(species)) {
      tissueResult[species].forEach(experimentObj => 
        currSpeciesFilterSet.add(experimentObj.filterID))
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
    tissueResult[species] = tissueResult[species].filter(experiment =>
      intersectSet.has(experiment.filterID)
    )
  }
  return tissueResult
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
async function processExperiment (experimentId, assembly) {
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

  // then just use bedFiles[0] and bwFiles[0] to populate the peak and bigwig
  // fields of the experiment
  await Promise.all(fsWritefilePromise()
}

Promise.all([readTissuePromise, readFilterPromise])
  .then(results => {
    let [tissues, filters] = results
    let tissuetissues.forEach(tissue => {
      let tissueObj
    })
  })
