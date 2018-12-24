/** API connecting to remote part:
 * {
 *    // General settings
 *    alignMode: String, // choices from ['promoter', 'enhancer']
 *    epiName: String,
 *
 *    // Query/search shared settings
 *    genomeAssembly: Array<String>(2), // choices from ['hg38', 'mm10']
 *    speciesText: Array<String>(2),
 *
 *    // Query specific settings
 *    speciesPeak1: Array<File>(),
 *    speciesInput1: File,
 *    promoterUp: Integer,
 *    promoterDown: Integer,
 *
 *    // Search specific settings
 *    speciesPeak2: Array<File>(),
 *    speciesInput2: File,
 *    searchRegionMode: String,
 *    // choices from ['genomeregion', 'genetype', 'genecluster', 'homoregion']
 *    genetypeSelect: String,
 *    // choices from ['protein-coding gene', 'lincRNA', 'pseudogene']
 *    clusters: String, // or Array<String>, TBD
 *    enhancerUp: Integer,
 *    enhancerDown: Integer,
 *
 *    // Parameters
 *    seqweight: Float,
 *    epiweight: Float,
 *    paras: Float,
 *    paramu: Float,
 *    parak: Float,
 *    paramu: Float,
 *    piA: Float,
 *    piC: Float,
 *    piG: Float,
 *    piT: Float,
 *    pi1: Float,
 *
 *    // Notification related
 *    emailNote: Boolean,
 *    mail: String
 * }
 */

import postAjax from './promisedAjax.js'

const PROMOTER_UP_DEFAULT = 500
const PROMOTER_DOWN_DEFAULT = 1000
const ENHANCER_FLANK_DEFAULT = 20000
const SEQ_WEIGHT_DEFAULT = 0.9
const EPI_WEIGHT_DEFAULT = 0.1
const PARA_S_DEFAULT = 0.3
const PARA_MU_DEFAULT = 0.3
const PARA_K_DEFAULT = 0.3
const PARA_PI_A_DEFAULT = 0.25
const PARA_PI_C_DEFAULT = 0.25
const PARA_PI_G_DEFAULT = 0.25
const PARA_PI_T_DEFAULT = 0.25
const PARA_PI_1_DEFAULT = 0.1

const FORM_SUBMIT_TARGET = '/backend/form_upload'

var app = new Vue({
  el: '#epialign_app',
  data: {
    showParam: false,
    hasError: false,
    submitted: false,
    formError: {
      modeNotSelected: false

    },
    showMoreParamText: 'Show more parameters...',
    submitStatus: null,

    // parameters involved in preset dataset loading
    showPreset: false,
    presetLoaded: false,
    experimentDict: {},
    speciesSupported: [
      {
        'name': 'human',
        'assembly': 'hg38',
        'shortHand': 'H'
      },
      {
        'name': 'mouse',
        'assembly': 'mm10',
        'shortHand': 'M'
      }
    ],
    // selected datasets
    selectedEntry: null,
    selectedExperimentIds: null,
    // ENCODE datasets
    encodeFilters: [],
    encodeSamples: [],
    // Other public datasets
    publicFilters: [],
    publicSamples: [],

    formParams: {
      alignMode: null,
      epiName: 'H3K4me3',
      genomeAssembly: [
        'hg38',
        'mm10'
      ],
      speciesPeak: [null, null],
      speciesText: [null, null],
      speciesInput: [null, null],
      promoterUp: PROMOTER_UP_DEFAULT,
      promoterDown: PROMOTER_DOWN_DEFAULT,
      searchRegionMode: 'genomeregion',
      genetypeSelect: null,
      clusters: null,
      enhancerUp: ENHANCER_FLANK_DEFAULT,
      enhancerDown: ENHANCER_FLANK_DEFAULT,

      seqweight: SEQ_WEIGHT_DEFAULT,
      epiweight: EPI_WEIGHT_DEFAULT,

      paras: PARA_S_DEFAULT,
      paramu: PARA_MU_DEFAULT,
      parak: PARA_K_DEFAULT,

      piA: PARA_PI_A_DEFAULT,
      piC: PARA_PI_C_DEFAULT,
      piG: PARA_PI_G_DEFAULT,
      piT: PARA_PI_T_DEFAULT,
      pi1: PARA_PI_1_DEFAULT,

      emailNote: false,
      mail: null
    }
  },
  created: function () {
    // load preset data sets from ENCODE by reading the spec JSON file
    Promise.all([
      postAjax('assets/encodeData.json', null, 'json', 'GET'),
      postAjax('assets/publicData.json', null, 'json', 'GET'),
      postAjax('assets/experimentDict.json', null, 'json', 'GET')
    ]).then(resultArray => {
      let [encodeData, publicData, experiments] = resultArray
      this.experimentDict = experiments
      // first build encodeDatasets
      this.buildDataSetFilterAndSample(
        this.encodeFilters, this.encodeSamples, encodeData
      )
      // then build publicDatasets
      this.buildDataSetFilterAndSample(
        this.publicFilters, this.publicSamples, publicData
      )
      this.presetLoaded = true
    })
  },
  methods: {
    buildDataSetFilterAndSample: function (
      filterDomArray, sampleDomArray, dataObj
    ) {
      let filterDict = this.buildFilter(filterDomArray, dataObj.filters)
      this.buildSample(
        sampleDomArray, dataObj.matchedTissues, filterDict, filterDomArray
      )
    },
    buildFilter: function (filterDomArray, filterObj) {
      let filterRevMap = new Map()
      for (let key in filterObj) {
        if (filterObj.hasOwnProperty(key)) {
          if (!filterRevMap.has(filterObj[key]['.id'])) {
            filterRevMap.set(filterObj[key]['.id'], filterRevMap.size)
            filterDomArray.push({
              id: filterObj[key]['.id'],
              label: filterObj[key]['.label']
            })
          }
        }
      }
      return filterRevMap
    },
    buildSample: function (
      sampleDomArray, sampleObj, filterRevMap, filterDomArray
    ) {
      let sampleSet = new Set()
      sampleObj.forEach(sample => {
        if (!sampleSet.has(sample['.label'])) {
          sampleSet.add(sample['.label'])
          let sampleEntry = {
            label: sample['.label']
          }
          sampleEntry.matchedFilters = Array(filterDomArray.length)
          sampleEntry.matchedFilters.fill(null)
          sampleEntry.filterDict = {}

          this.speciesSupported.forEach(species => {
            let experiments = sample[species.name].experiments
            experiments.forEach(experiment => {
              if (!filterRevMap.has(experiment.filterId)) {
                console.log('Warning: filter ID does not exist! ' +
                  experiment.filterId)
              } else {
                let filterIndex = filterRevMap.get(experiment.filterId)
                if (!sampleEntry.matchedFilters[filterIndex]) {
                  sampleEntry.matchedFilters[filterIndex] = {}
                  this.speciesSupported.forEach(spc => {
                    sampleEntry.matchedFilters[filterIndex][spc.name] = []
                  })
                }
                sampleEntry.matchedFilters[filterIndex][species.name]
                  .push(experiment)
              }
            })
          })

          sampleDomArray.push(sampleEntry)
        }
      })
    },
    selectEntry: function (entry) {
      this.selectedEntry = entry
      if (entry) {
        this.selectedExperimentIds = {}
        this.speciesSupported.forEach(species => {
          this.selectedExperimentIds[species.name] = entry[species.name][0].id
        })
      } else {
        this.selectedExperimentIds = null
      }
    },
    selectExperiment: function (experimentId, speciesName) {
      let tempSelectedIds = this.selectedExperimentIds
      tempSelectedIds[speciesName] = experimentId
      this.selectedExperimentIds = null
      this.selectedExperimentIds = tempSelectedIds
    },
    toggleShowParam: function () {
      this.showParam = !this.showParam
      this.showMoreParamText = this.showParam
        ? 'Show less parameters.' : 'Show more parameters...'
    },
    checkAlignMode: function () {
      if (this.checkModeConflict()) {
        this.searchRegionMode = 'genomeregion'
      }
    },
    checkModeConflict: function () {
      return (this.promoterSelected &&
        this.searchRegionMode === 'homoregion'
      ) || (this.enhancerSelected && (
        this.searchRegionMode === 'genetype' ||
        this.searchRegionMode === 'genecluster'
      ))
    },
    validateForm: function () {
      let hasError = false
      for (let key in this.formError) {
        if (this.formError.hasOwnProperty(key)) {
          this.formError[key] = false
        }
      }
      if (!this.formParams.alignMode) {
        this.formError.modeNotSelected = true
        hasError = true
      }

      this.hasError = hasError
      return !this.hasError
    },
    submitForm: function () {
      // TODO: validate formParams
      if (!this.validateForm()) {
        this.submitStatus = '<i class="material-icons">clear</i> ' +
          'Error encountered. Please review your ' +
          'submission before continuing.'
        return
      }
      this.submitted = false
      this.submitStatus = 'Submitting data to server. Please wait ...'
      let formData = new window.FormData(this.$refs.mainForm)
      if (this.selectedExperimentIds) {
        // populate encodeData[]
        let encodeData = []
        this.formParams.genomeAssembly.forEach(assembly => {
          this.speciesSupported.some(species => {
            if (species.assembly !== assembly) {
              return false
            }
            encodeData.push(this.selectedExperimentIds[species.name])
          })
        })
        if (encodeData.length === 2) {
          encodeData.forEach(dataId => {
            formData.append('encodeData[]', dataId)
          })
        }
      }
      postAjax(FORM_SUBMIT_TARGET, formData, 'json', 'POST')
        .then(response => {
          let runid = response.runid
          this.submitted = true
          this.submitStatus = '<i class="material-icons">' +
            'check_circle' +
            '</i> Data submitted to server. Redirecting to the ' +
            'result page ...'
          // Testing code
          // if (window.confirm('Please click "Ok" to go to the result ' +
          //   'page, click "Cancel" to remain at this page.')
          // ) {
          //   window.setTimeout(() => {
          //     window.location.href = '/result_page/' + runid
          //   }, 500)
          // } else {

          // }
          window.setTimeout(() => {
            window.location.href = '/result_page/' + runid
          }, 3500)
        })
        .catch(err => {
          this.hasError = true
          this.submitStatus = '<i class="material-icons">' +
            'clear' +
            '</i> Error encountered. Error code: ' +
            err.status
        })
    },
    selectEncodeData: function () {
      this.showPreset = true
    },
    closeEncodeDialog: function () {
      this.selectEntry(null)
      this.showPreset = false
    },
    confirmEncodeSelection: function () {
      this.showPreset = false
    }
  },
  computed: {
    promoterSelected: function () {
      return this.formParams.alignMode === 'promoter'
    },
    enhancerSelected: function () {
      return this.formParams.alignMode === 'enhancer'
    },
    genomeRegionLabel: function () {
      return 'Define searching regions with a BED file' +
        (this.enhancerSelected ? '.' : ' / a gene list.')
    },
    genomePlaceholder: function () {
      return 'Paste BED data ' +
        (this.enhancerSelected ? 'here.' : 'or gene names here.')
    },
    genomeRegionSelected: function () {
      return this.formParams.searchRegionMode === 'genomeregion'
    },
    geneTypeSelected: function () {
      return this.formParams.searchRegionMode === 'genetype'
    },
    geneClusterSelected: function () {
      return this.formParams.searchRegionMode === 'genecluster'
    },
    homologRegionSelected: function () {
      return this.formParams.searchRegionMode === 'homoregion'
    }
  }
})
