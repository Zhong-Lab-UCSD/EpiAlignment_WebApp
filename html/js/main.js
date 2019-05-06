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
 *    mail: String
 * }
 */

import postAjax from './promisedAjax.js'

const PROMOTER_UP_DEFAULT = 1000
const PROMOTER_DOWN_DEFAULT = 500
const ENHANCER_FLANK_DEFAULT = 20000
const EPI_WEIGHT_DEFAULT = 0.1
const PARA_S_DEFAULT = 0.3
const PARA_MU_DEFAULT = 0.3
const PARA_K_DEFAULT = 0.5
const PARA_PI_A_DEFAULT = 0.25
const PARA_PI_C_DEFAULT = 0.25
const PARA_PI_G_DEFAULT = 0.25
const PARA_PI_T_DEFAULT = 0.25
const PARA_PI_1_DEFAULT = 0.1

const FORM_SUBMIT_TARGET = '/backend/form_upload'
const CLUSTER_QUERY_TAGET = '/backend/get_cluster'

const CLUSTER_QUERY_DEBOUNCE = 300
const CLUSTER_QUERY_DEBOUNCE_WHEN_OPEN = 50

const NUM_UP_DOWN_STREAM = 2
const NUM_PARAMS = 8

var app = new Vue({
  el: '#epialign_app',
  data: {
    showParam: false,
    hasError: false,
    submitted: false,
    formError: {
      modeNotSelected: false,
      peakError: [false, false],
      inputError: [false, false],
      promoterLengthError: false,
      referenceError: false,
      searchRegionModeError: false,
      clusterError: false,
      enhancerLengthError: false,
      weightError: false,
      paramError: false
    },
    showMoreParamText: 'Show more parameters...',
    submitStatus: null,
    showModeHelp: false,
    showParamHelp: false,

    bedHint1On: false,
    bedHint2On: false,

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
    tempSelectedEntry: null,
    tempSelectedExperimentIds: null,
    selectedEntryDescHtml: null,
    selectedEntry: null,
    selectedExperimentIds: null,
    // ENCODE datasets
    encodeFilters: [],
    encodeSamples: [],
    // Other public datasets
    publicFilters: [],
    publicSamples: [],

    peakFiles: [[], []],
    inputFiles: [[], []],

    // Cluster related
    clusterText: null,
    postedClusterText: null,
    clusterCandidates: [],
    showClusterCandidate: false,
    selectedCluster: null,
    showClusterMessage: false,
    clusterMessage: '',

    formParams: {
      alignMode: 'enhancer',
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
      searchRegionMode: 'homoregion',
      genetypeSelect: null,
      enhancerUp: ENHANCER_FLANK_DEFAULT,
      enhancerDown: ENHANCER_FLANK_DEFAULT,

      epiweight: EPI_WEIGHT_DEFAULT,

      paras: PARA_S_DEFAULT,
      paramu: PARA_MU_DEFAULT,
      parak: PARA_K_DEFAULT,

      piA: PARA_PI_A_DEFAULT,
      piC: PARA_PI_C_DEFAULT,
      piG: PARA_PI_G_DEFAULT,
      piT: PARA_PI_T_DEFAULT,
      pi1: PARA_PI_1_DEFAULT,

      mail: null
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
      return 'Define target regions with a BED file' +
        (this.enhancerSelected ? '.' : ' / a gene list.')
    },
    targetPlaceholder: function () {
      return 'Paste BED6 data ' +
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
    },
    maxSampleFilterLength: function () {
      return Math.max(this.encodeFilters.length, this.publicFilters.length)
    },

    expectedRunTime: function () {
      return this.formParams.alignMode === 'enhancer'
        ? '2 minutes' : '10 seconds'
    },

    presetButtonText: function () {
      if (!this.selectedExperimentIds) {
        if (this.peakFiles.some(peakFile => peakFile.length)) {
          return 'Custom peak file uploaded. ' +
            'Click here to use public datasets instead.'
        } else {
          return 'Click here to select a paired public dataset.'
        }
      }
      return 'Experiment ' + this.selectedEntryDescHtml +
        ' selected. Click here to change.'
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
  mounted: function () {
    document.addEventListener('click', () => this.closeClusterPanel())
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

    tempSelectEntry: function (entry, sampleDesc, filter) {
      this.tempSelectedEntry = entry
      if (entry) {
        this.tempSelectedEntryDescHtml = '<strong>' + filter.label +
          '</strong> on <strong>' + sampleDesc + '</strong>'
        this.tempSelectedEntryDesc = filter.label + ' on ' + sampleDesc
        this.tempSelectedExperimentIds = {}
        this.speciesSupported.forEach(species => {
          this.tempSelectedExperimentIds[species.name] = entry[species.name][0].id
        })
      } else {
        this.tempSelectedExperimentIds = null
        this.tempSelectedEntryDescHtml = null
      }
    },

    selectExperiment: function (experimentId, speciesName) {
      let tempSelectedIds = this.tempSelectedExperimentIds
      tempSelectedIds[speciesName] = experimentId
      this.tempSelectedExperimentIds = null
      this.tempSelectedExperimentIds = tempSelectedIds
    },

    toggleShowParam: function () {
      this.showParam = !this.showParam
      this.showMoreParamText = this.showParam
        ? 'Show less parameters.' : 'Show more parameters...'
    },

    checkAlignMode: function () {
      // clear all input regions (peaks can be retained)
      this.clearSample()
      if (this.checkModeConflict() || (
        !this.modeNotSelected && !this.formParams.searchRegionMode
      )) {
        this.formParams.searchRegionMode =
          this.promoterSelected ? 'genecluster' : 'homoregion'
      }
    },

    checkModeConflict: function () {
      return (this.promoterSelected &&
        this.formParams.searchRegionMode === 'homoregion'
      ) || (this.enhancerSelected && (
        this.formParams.searchRegionMode === 'genetype' ||
        this.formParams.searchRegionMode === 'genecluster'
      ))
    },

    validateForm: function () {
      let hasError = !this.$refs.mainForm.checkValidity()

      for (let key in this.formError) {
        if (this.formError.hasOwnProperty(key)) {
          if (Array.isArray(this.formError[key])) {
            this.formError[key].fill(false)
          } else {
            this.formError[key] = false
          }
        }
      }
      if (!this.formParams.alignMode) {
        this.formError.modeNotSelected = true
      }

      // peakError
      if (!this.selectedExperimentIds) {
        this.peakFiles.forEach((peakFile, index) => {
          if (!peakFile.length) {
            this.formError.peakError.splice(index, 1, true)
          }
        })
      }

      // searchRegionModeError
      if (!this.formParams.searchRegionMode) {
        this.formError.searchRegionModeError = true
      }

      // inputError
      this.formParams.speciesText.forEach((text, index) => {
        if ((!text || !text.trim().length) && // no input text
          !this.inputFiles[index].length && // no input file
          (
            // first input cannot be empty in non-gene-cluster sub-mode
            (index === 0 && !this.geneClusterSelected) ||
            // second input cannot be empty in genome-region sub-mode
            (index === 1 && this.genomeRegionSelected)
          )
        ) {
          this.formError.inputError.splice(index, 1, true)
        }
      })

      // promoterLengthError
      for (let i = 0; i < NUM_UP_DOWN_STREAM; i++) {
        if (!this.$refs['promoter[' + i + ']'].checkValidity()) {
          this.formError.promoterLengthError = true
          break
        }
      }

      // referenceError
      if (this.homologRegionSelected &&
        this.formParams.genomeAssembly[0] ===
        this.formParams.genomeAssembly[1]
      ) {
        this.formError.referenceError = true
      }

      // clusterError
      if (this.geneClusterSelected && !this.selectedCluster &&
        (!this.clusterText || !this.clusterText.match(/Cluster_[1-9][0-9]*/))
      ) {
        this.formError.clusterError = true
      }

      // enhancerLengthError
      if (this.homologRegionSelected) {
        for (let i = 0; i < NUM_UP_DOWN_STREAM; i++) {
          if (!this.$refs['enhancer[' + i + ']'].checkValidity()) {
            this.formError.enhancerLengthError = true
            break
          }
        }
      }

      // weightError
      if (!this.$refs['weight'].checkValidity()) {
        this.formError.weightError = true
      }

      // paramError
      for (let i = 0; i < NUM_PARAMS; i++) {
        if (!this.$refs['param[' + i + ']'].checkValidity()) {
          this.formError.paramError = true
          break
        }
      }

      for (let key in this.formError) {
        if (this.formError.hasOwnProperty(key) && (
          Array.isArray(this.formError[key])
            ? this.formError[key].some(value => value)
            : this.formError[key]
        )) {
          hasError = true
          break
        }
      }
      this.hasError = hasError
      return !this.hasError
    },

    submitForm: function () {
      // TODO: validate formParams
      if (!this.validateForm()) {
        this.submitStatus = '<i class="material-icons iconAtLeft">clear</i> ' +
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
            return true
          })
        })
        if (encodeData.length === this.formParams.genomeAssembly.length) {
          encodeData.forEach(dataId => {
            formData.append('encodeData[]', dataId)
          })
          formData.append('publicDataDesc', this.selectedEntryDesc)
        }
      }

      if (this.selectedCluster) {
        formData.append('clusters', this.selectedCluster.id)
      } else if (this.clusterText &&
        this.clusterText.match(/Cluster_[1-9][0-9]*/)
      ) {
        formData.append('clusters', this.clusterText)
      }

      // Add timezone data
      formData.append('timeZone', moment.tz.guess())

      postAjax(FORM_SUBMIT_TARGET, formData, 'json', 'POST')
        .then(response => {
          let runid = response.runid
          this.submitted = true
          this.submitStatus = '<i class="material-icons iconAtLeft">' +
            'check_circle' +
            '</i> Data submitted to server. Redirecting to the ' +
            'result page ...'
          window.setTimeout(() => {
            window.location.href = '/result_page/' + runid
          }, 2000)
        })
        .catch(err => {
          this.hasError = true
          this.submitStatus = '<i class="material-icons iconAtLeft">' +
            'clear' +
            '</i> Error encountered. Error code: ' +
            err.status
        })
    },

    selectEncodeData: function () {
      this.showPreset = true
    },
    clearEncodeData: function () {
      this.selectedEntry = null
      this.selectedExperimentIds = null
      this.selectedEntryDescHtml = null
    },
    closeEncodeDialog: function () {
      this.tempSelectEntry(null)
      this.showPreset = false
    },
    confirmEncodeSelection: function () {
      this.showPreset = false
      this.selectedEntry = this.tempSelectedEntry
      this.selectedExperimentIds = this.tempSelectedExperimentIds
      this.selectedEntryDescHtml = this.tempSelectedEntryDescHtml
      this.selectedEntryDesc = this.tempSelectedEntryDesc
      this.peakFiles.forEach(peakFile => peakFile.splice(0))
      // Manually clear this.$ref.peakFile1 and this.$ref.peakFile2
      this.$refs.peakFile1.value = ''
      this.$refs.peakFile2.value = ''
    },

    peakFileChanged: function (event, index) {
      this.peakFiles[index].splice(0, this.peakFiles[index].length,
        ...event.target.files)
      if (event.target.files.length) {
        this.clearEncodeData()
      }
    },

    inputFileChanged: function (event, index) {
      this.inputFiles[index].splice(0, this.inputFiles[index].length,
        ...event.target.files)
    },

    clusterInputChanged: function () {
      if (this.clusterTimeOut) {
        window.clearTimeout(this.clusterTimeOut)
      }
      this.clusterTimeOut = window.setTimeout(
        () => this.postClusterText(),
        this.showClusterCandidate
          ? CLUSTER_QUERY_DEBOUNCE_WHEN_OPEN
          : CLUSTER_QUERY_DEBOUNCE
      )
    },

    postClusterText: function () {
      this.clusterTimeOut = null
      if (this.postedClusterText !== this.clusterText) {
        // post cluster text to remote server
        this.postedClusterText = this.clusterText
        if (this.postedClusterText &&
          !this.postedClusterText.startsWith('Cluster_')
        ) {
          this.selectedCluster = null
          postAjax(
            CLUSTER_QUERY_TAGET + '/' + this.postedClusterText,
            null, 'json', 'GET'
          ).then(response => {
            this.clusterCandidates.splice(0)
            response.fullMatchList.forEach(cluster =>
              this.clusterCandidates.push(cluster)
            )
            response.partialMatchList.forEach(cluster =>
              this.clusterCandidates.push(cluster)
            )
            this.showClusterMessage = !response.partialMatchList.length && (
              response.maxExceeded || !response.fullMatchList.length
            )
            this.clusterMessage = response.maxExceeded
              ? 'Continue typing to get ' +
              (response.fullMatchList.length ? 'more ' : '') + 'candidates.'
              : 'No paralogues matched in selected species.'
            this.showClusterCandidate = true
          }).catch(err => {
            console.log(err)
          })
        }
      }
    },

    selectCluster: function (cluster) {
      this.selectedCluster = cluster
    },
    applyCluster: function (clusterId) {
      this.clusterText = clusterId
      this.showClusterCandidate = false
      this.postedClusterText = ''
    },
    closeClusterPanel: function () {
      if (this.showClusterCandidate) {
        this.showClusterCandidate = false
        this.selectedCluster = null
      }
    },

    getGeneSymbolAlias: function (gene, query) {
      if (gene.symbol.toLowerCase().includes(query.toLowerCase())) {
        return this.getHighlightedString(gene.symbol, query)
      }
      let aliasResult = null
      gene.aliases.some(alias => {
        if (alias.toLowerCase() === query.toLowerCase()) {
          aliasResult = '<strong>' + alias + '</strong>'
        }
      })
      if (!aliasResult) {
        gene.aliases.some(alias => {
          if (alias.toLowerCase().includes(query.toLowerCase())) {
            aliasResult = this.getHighlightedString(alias, query)
          }
        })
      }
      return gene.symbol + (aliasResult ? ' (' + aliasResult + ')' : '')
    },
    getHighlightedString: function (haystack, needle) {
      haystack = haystack || ''
      let index = haystack.toLowerCase().indexOf(needle.toLowerCase())
      if (index < 0) {
        return haystack
      }
      return haystack.slice(0, index) + '<strong>' +
        haystack.slice(index, index + needle.length) +
        '</strong>' + haystack.slice(index + needle.length)
    },

    getDatasetHref: function (experimentId) {
      if (experimentId.startsWith('ENC')) {
        // ENOCDE
        return 'https://www.encodeproject.org/experiments/' + experimentId
      } else if (experimentId.startsWith('GS')) {
        // GEO
        return 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + experimentId
      }
    },
    getDatasetLinkTitle: function (experimentId) {
      if (experimentId.startsWith('ENC')) {
        // ENOCDE
        return 'See experiment details / download data on ENCODE website.'
      } else if (experimentId.startsWith('GS')) {
        // GEO
        return 'See experiment details / download data on the GEO.'
      }
    },

    useSample: function () {
      if (!this.formParams.alignMode) {
        console.log('Please select an align mode first!')
      }
      if (this.formParams.alignMode === 'promoter') {
        this.selectedEntryDescHtml = '<strong>ChIP-Seq (H3K4me3)</strong>' +
          ' on <strong>adult B-lymphocytes</strong>'
        this.selectedEntryDesc = 'ChIP-Seq (H3K4me3) on adult B-lymphocytes'
        this.selectedExperimentIds = {
          human: 'ENCSR057BWO',
          mouse: 'ENCSR000CGK'
        }
        this.formParams.genomeAssembly.splice(
          0, this.formParams.genomeAssembly.length, 'hg38', 'mm10')
        this.formParams.speciesText.splice(0, 1, 'SLCO4A1')
        this.clusterText = 'Cluster_6635'
        this.formParams.searchRegionMode = 'genecluster'
      } else {
        // enhancer mode
        this.selectedEntryDescHtml = '<strong>ChIP-Seq (H3K4me3)</strong>' +
          ' on <strong>round spermatids</strong>'
        this.selectedEntryDesc = 'ChIP-Seq (H3K4me3) on round spermatids'
        this.selectedExperimentIds = {
          human: 'GSM1673960',
          mouse: 'GSM1674016'
        }
        this.formParams.genomeAssembly.splice(
          0, this.formParams.genomeAssembly.length, 'mm10', 'hg38')
        this.formParams.speciesText.splice(
          0, 1, 'chr12\t8207583\t8209349\tRegion_example\t0\t+\n')
        this.formParams.searchRegionMode = 'homoregion'
        this.formParams.enhancerUp = 20000
        this.formParams.enhancerDown = 20000
      }
    },

    clearSample: function () {
      this.formParams.speciesText.splice(0, 2, '', '')
      this.$refs.speciesInputFile1.value = ''
      this.$refs.speciesInputFile2.value = ''
      this.clusterText = ''
    },

    toggleModeHelp: function () {
      this.showModeHelp = !this.showModeHelp
      if (!this.formParams.alignMode && this.showModeHelp) {
        // no align mode is chosen, turn off showModeHelp after 2 seconds.
        window.setTimeout(() => (this.showModeHelp = false), 2000)
      }
    }
  }
})
