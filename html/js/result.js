import postAjax from './promisedAjax.js'

const setTimeoutPromise = function (timeout, resolvedValue) {
  return new Promise((resolve, reject) =>
    window.setTimeout(resolve, timeout, resolvedValue)
  )
}

const INQUIRY_TARGET_PREFIX = '/backend/results/'
const INQUIRY_IMAGE_PREFIX = '/backend/result_image/'

const UCSC_TARGET = 'https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=xycao&'
const GIVE_TARGET = 'browser.html?'

const ucscRefToSession = {
  mm10: 'mm10_epialign',
  hg38: 'hg38_epialign'
}

const RUN_INFO_DICT = 'assets/runInfoDict.json'
const DEFAULT_ERROR_MSG = 'Error happened during epialignment.\n' +
  'If you wish to report the issue, please include your URL in the address ' +
  'bar in your report.'

const STATUS_RUNNING = -1
const POLLING_INTERVAL_ENHANCER = 15000 // 15 seconds
const POLLING_INTERVAL_PROMOTER = 3000 // 3 seconds

const MINIMUM_HEATMAP_GAP = 10

const ALIGN_SCORE_95 = 166.63 // 95% alignment score of random sequence

const expandedRunInfoKeyList = [
  [
    'runid',
    'submitTime',
    'completeTime'
  ], [
    'alignMode',
    'subMode',
    'publicDataDesc'
  ], [
    'queryGenomeAssembly',
    'queryInput',
    'queryPeak'
  ], [
    'targetGenomeAssembly',
    'targetPeak',
    'targetInput',
    'clusters'
  ], [
    'promoterUp',
    'promoterDown'
  ], [
    'enhancerUp',
    'enhancerDown'
  ], [
    'seqWeight',
    'epiWeight'
  ], [
    'paraS',
    'paraMu',
    'paraK',
    'piA',
    'piC',
    'piG',
    'piT',
    'pi1'
  ]
]

const collapsedRunInfoKeyList = [
  [
    'runid',
    'alignMode'
  ]
]

var app = new Vue({
  el: '#result_app',
  data: {
    expandedRunInfoList: [],
    collapsedRunInfoList: [],
    runid: null,
    loading: true,
    lastRefreshedTime: '',
    highlightLoading: false,
    showRunDetails: false,
    queryGenomeAssembly: null,
    targetGenomeAssembly: null,
    runInfoDictPromise: null,
    email: null,
    downloadLink: '',
    error: false,
    errorMessage: '',
    showWarning: false,
    imageBlobUrl: null,
    alignMode: null,
    pollingTime: null,
    promoterColumnList: [],
    epiHeatmapList: [],
    hasSequence: true,
    geneIdentifier1: null, // ['transID1', 'region_name1']
    geneIdentifier2: null, // ['transID2', 'region_name2']
    showHeatmap: false,

    heatmapDesc: {
      min: Number.MAX_VALUE,
      max: Number.MIN_VALUE,
      span: MINIMUM_HEATMAP_GAP,
      showPercentile: false,
      percentileValue: ALIGN_SCORE_95
    },
    percentileEnabled: false,

    showResultHint: true,
    doNotShowHintAgain: false,
    rowsPerPageItems: [
      10,
      20,
      50,
      { 'text': '$vuetify.dataIterator.rowsPerPageAll', 'value': -1 }
    ],
    paginationObj: {
      rowsPerPage: 20
    },
    headers: [
      {
        text: '#',
        value: 'index',
        align: 'left',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      },
      {
        text: '',
        value: null,
        align: 'left',
        sortable: false,
        class: ['dataTableCell', 'dataTableCellNoGap'],
        width: 'auto'
      },
      {
        text: 'Query coordinate',
        value: 'region1',
        align: 'left',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      },
      {
        text: 'Target coordinate',
        value: 'region2',
        align: 'left',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      },
      {
        text: 'Epi Score',
        value: 'scoreE',
        align: 'right',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      },
      {
        text: 'Seq Score',
        value: 'scoreS',
        align: 'right',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      },
      {
        text: '',
        html: '<i class="material-icons">swap_horiz</i>',
        value: 'shifted',
        align: 'left',
        sortable: true,
        class: ['dataTableCell', 'dataTableCellNoGap'],
        width: 'auto',
        tooltip: 'Whether the EpiAlignment hit is different from the sequence-only hit.'
      }
    ],
    dataEntries: []
  },
  created: function () {
    const urlParams = new window.URLSearchParams(window.location.search)
    this.runid = urlParams.get('runid')
    this.pollResult(this.runid).then(
      data => this.showData(data)
    ).catch(err => {
      // TODO: err.status should say something
      let response = err.response || err
      this.setRunParameters(response)

      this.loading = false
      this.error = true
      this.errorCode = response.status
      this.errorMessage =
        (response.errMessage || response || DEFAULT_ERROR_MSG).trim()
    })
    this.runInfoDictPromise = postAjax(RUN_INFO_DICT, null, 'json', 'GET')
  },

  mounted () {
    if (localStorage.getItem('doNotShowResultHint')) {
      this.showResultHint = false
    }
  },
  methods: {
    closeHint: function () {
      if (this.doNotShowHintAgain) {
        localStorage.setItem('doNotShowResultHint', 'true')
      }
      this.showResultHint = false
    },
    closeWarning: function () {
      this.showWarning = false
    },
    foldRunInfoDetails: function () {
      this.showRunDetails = false
    },
    showRunInfoDetails: function () {
      this.showRunDetails = true
    },
    pollResult: function (runid) {
      // TODO: validate formParams
      if (!runid) {
        return Promise.reject(new Error('✖ Error: runid does not exist! ' +
          'Please try submit your data again.'))
      }
      return postAjax(INQUIRY_TARGET_PREFIX + runid, null, 'json', 'GET')
        .then(response => {
          this.email = response.email || this.email

          if (response.status === STATUS_RUNNING) {
            this.setRunParameters(response)
            if (!this.pollingTime) {
              this.pollingTime = response.alignMode === 'promoter'
                ? POLLING_INTERVAL_PROMOTER : POLLING_INTERVAL_ENHANCER
            }
            this.lastRefreshedTime = moment().format('MMM D, Y, HH:mm:ss')
            this.highlightLoading = true
            window.setTimeout(() => (this.highlightLoading = false), 1500)
            return setTimeoutPromise(this.pollingTime, runid)
              .then(runid => this.pollResult(runid))
          } else if (response.status > 0) {
            // error, reject the promise and let .catch to handle
            return Promise.reject(response)
          } else {
            this.loading = false
            this.setRunParameters(response, true)
            return response.data
          }
        })
    },
    populateInfoList: function (keyList, infoList, response) {
      this.runInfoDictPromise.then(runInfoDict => {
        infoList.splice(0)
        keyList.forEach(keyRow => {
          let infoRow = []
          keyRow.forEach(key => {
            if (response.hasOwnProperty(key)) {
              let value = response[key]
              if (runInfoDict.hasOwnProperty(key) &&
                runInfoDict[key].values
              ) {
                runInfoDict[key].values.some(valueEntry => {
                  if (valueEntry.key === value) {
                    value = valueEntry.name
                    return true
                  }
                  return false
                })
              }
              infoRow.push({
                key: runInfoDict[key].html || runInfoDict[key].name,
                value: value,
                parameter: runInfoDict[key].parameter
              })
            }
          })
          infoList.push(infoRow)
        })
      })
    },
    setRunParameters: function (response, refresh) {
      document.title = 'EpiAlignment - Result (submitted at: ' +
        response.completeTime + ')'

      if (response.errMessage) {
        this.errorMessage = response.errMessage.trim()
        this.showWarning = true
      }

      this.alignMode = response.alignMode

      this.queryGenomeAssembly = response.queryGenomeAssembly
      this.targetGenomeAssembly = response.targetGenomeAssembly
      if (refresh || !this.expandedRunInfoList.length) {
        // Add run info to expandedRunInfoList and collapsedRunInfoList
        this.populateInfoList(
          expandedRunInfoKeyList, this.expandedRunInfoList, response
        )
        this.populateInfoList(
          collapsedRunInfoKeyList, this.collapsedRunInfoList, response
        )
      }

      if (response.publicDataDesc) {
        // used public data
        // use `response.queryPeak`, `response.targetPeak` to construct
        // datasets to be used in GIVE browser
        this.displayDataSets = [
          [
            'knownGene',
            response.queryPeak.toLowerCase() + '-p',
            response.queryPeak.toLowerCase() + '-b'
          ], [
            'knownGene',
            response.targetPeak.toLowerCase() + '-p',
            response.targetPeak.toLowerCase() + '-b'
          ]
        ]
      } else {
        this.displayDataSets = ['knownGene']
      }
      this.downloadLink = window.location.protocol + '//' +
        window.location.host + '/download/' + this.runid + '.txt'
    },
    showData: function (dataEntries) {
      if (Array.isArray(dataEntries)) {
        this.dataEntries = dataEntries
        if (dataEntries[0]) {
          this.verifyNameSource(dataEntries[0])
        }
        if (this.alignMode === 'promoter') {
          this.processPromoterData(dataEntries)
        } else {
          this.processEnhancerData(dataEntries)
        }
      }
    },
    toggleRow: function (row) {
      if (this.alignMode === 'enhancer') {
        row.expanded = !row.expanded
        Vue.set(row.item, 'image', null)
        Vue.set(row.item, 'imageError', false)
        if (row.expanded) {
          this.getEnhancerImage(row.item)
        }
      }
    },
    getEnhancerImage: function (rowItem) {
      let imageLink = INQUIRY_IMAGE_PREFIX + this.runid + '/' +
        rowItem.index + '.png'
      return postAjax(
        imageLink, null, 'blob', 'GET'
      ).then(imageBlob => {
        if (this.imageBlobUrl) {
          window.URL.revokeObjectURL(this.imageBlobUrl)
        }
        this.imageBlobUrl = window.URL.createObjectURL(imageBlob)
        Vue.set(rowItem, 'image', this.imageBlobUrl)
      }).catch(async error => {
        Vue.set(rowItem, 'imageError',
          await (new window.Response(error.response || error)).text() ||
            'Cannot load image.')
      })
    },
    verifyNameSource: function (dataEntry) {
      let headerName1Index = 2
      let headerName2Index = 4
      if (!this.geneIdentifier1) {
        this.geneIdentifier1 =
          dataEntry.transID1 === '.' ? 'region_name1' : 'transID1'
        this.headers.splice(headerName1Index, 0, {
          text: 'Query name',
          value: this.geneIdentifier1,
          align: 'left',
          sortable: true,
          class: 'dataTableCell',
          width: 'auto'
        })
        if (dataEntry.ensID1 !== '.') {
          // add corresponding row to headers
          this.headers.splice(headerName1Index + 1, 0, {
            text: 'Query Ensembl ID',
            value: 'ensID1',
            align: 'left',
            sortable: true,
            class: 'dataTableCell',
            width: 'auto'
          })
          headerName2Index++
        }
      }
      if (!this.geneIdentifier2) {
        this.geneIdentifier2 =
          dataEntry.transID2 === '.' ? 'region_name2' : 'transID2'
        this.headers.splice(headerName2Index, 0, {
          text: 'Target name',
          value: this.geneIdentifier2,
          align: 'left',
          sortable: true,
          class: 'dataTableCell',
          width: 'auto'
        })
        if (dataEntry.ensID2 !== '.') {
          // add corresponding row to headers
          this.headers.splice(headerName2Index + 1, 0, {
            text: 'Target Ensembl ID',
            value: 'ensID2',
            align: 'left',
            sortable: true,
            class: 'dataTableCell',
            width: 'auto'
          })
        }
      }
      this.hasSequence = dataEntry.scoreS !== '.'
    },
    processPromoterData: function (data) {
      let columnMap = new Map()
      let rowMap = new Map()
      data.forEach(entry =>
        this.parseSinglePromoterEntry(entry, columnMap, rowMap))
      this.processHeatmap()
      this.markRowMaxFlags()
    },
    parseSinglePromoterEntry: function (entry, columnMap, rowMap) {
      this.initPromoterEntryIfNotExists(entry, columnMap, rowMap)
      // then populate the data
      let newDataEntry = {
        epiScore: parseFloat(entry.scoreE)
      }
      if (this.hasSequence) {
        newDataEntry.seqScore = parseFloat(entry.scoreS)
      }
      this.epiHeatmapList[rowMap.get(entry[this.geneIdentifier1])]
        .values[columnMap.get(entry[this.geneIdentifier2])] = newDataEntry
    },
    initPromoterEntryIfNotExists: function (entry, columnMap, rowMap) {
      if (!columnMap.has(entry[this.geneIdentifier2])) {
        columnMap.set(
          entry[this.geneIdentifier2], this.promoterColumnList.length)
        this.promoterColumnList.push(entry[this.geneIdentifier2])
      }
      if (!rowMap.has(entry[this.geneIdentifier1])) {
        rowMap.set(
          entry[this.geneIdentifier1], this.epiHeatmapList.length)
        this.epiHeatmapList.push({
          geneId: entry[this.geneIdentifier1],
          values: []
        })
      }
    },
    processHeatmap: function () {
      let propertyList = ['epiScore']
      if (this.hasSequence) {
        propertyList.push('seqScore')
      }
      this.findMinMaxTableValue(this.epiHeatmapList, propertyList)
      this.epiHeatmapList.forEach(entry => entry.values.forEach(value => (
        value.epiScoreNorm =
          (value.epiScore - this.heatmapDesc.min) / this.heatmapDesc.span
      )))
      if (this.hasSequence) {
        this.epiHeatmapList.forEach(entry => entry.values.forEach(value => {
          value.seqScoreNorm =
            (value.seqScore - this.heatmapDesc.min) / this.heatmapDesc.span
          if (this.percentileEnabled) {
            value.insignificant =
              (this.heatmapDesc.percentileValue > value.seqScore)
          }
        }))
      }
    },
    markRowMaxFlags: function () {
      this.epiHeatmapList.forEach(entry => {
        let epiMaxIndex =
          this.markSingleRowMaxFlag(entry, 'epiScore', 'epiMax')
        if (this.hasSequence) {
          let seqMaxIndex =
            this.markSingleRowMaxFlag(entry, 'seqScore', 'seqMax')
          entry.labelMax = (epiMaxIndex !== seqMaxIndex) && (
            entry.values[epiMaxIndex].epiScore >
            entry.values[epiMaxIndex].seqScore
          ) && (
            entry.values[seqMaxIndex].seqScore >
              entry.values[seqMaxIndex].epiScore
          )
        }
      })
    },
    markSingleRowMaxFlag: function (row, property, propertyToMark) {
      let maxIndex = -1
      let maxValue = Number.MIN_VALUE
      row.values.forEach((value, index) => {
        if (value[property] > maxValue) {
          maxIndex = index
          maxValue = value[property]
        }
        delete value[propertyToMark]
      })
      row.values[maxIndex][propertyToMark] = true
      return maxIndex
    },
    findMinMaxTableValue: function (table, entryPropertyList) {
      if (this.percentileEnabled) {
        this.heatmapDesc.min = ALIGN_SCORE_95 - MINIMUM_HEATMAP_GAP
        this.heatmapDesc.max = ALIGN_SCORE_95
      }
      table.forEach(entry => entry.values.forEach(
        value => entryPropertyList.forEach(prop => {
          if (this.heatmapDesc.min > value[prop]) {
            this.heatmapDesc.min = value[prop]
          }
          if (this.heatmapDesc.max < value[prop]) {
            this.heatmapDesc.max = value[prop]
            if (this.percentileEnabled) {
              this.heatmapDesc.showPercentile = true
            }
          }
        })
      ))
      if (this.heatmapDesc.max - this.heatmapDesc.min < MINIMUM_HEATMAP_GAP) {
        this.heatmapDesc.span = MINIMUM_HEATMAP_GAP
        this.heatmapDesc.min = this.heatmapDesc.max - MINIMUM_HEATMAP_GAP
      } else {
        this.heatmapDesc.span = this.heatmapDesc.max - this.heatmapDesc.min
      }
    },
    processEnhancerData: function (data) {
      // reformat table header
      // remove target name
      this.headers.splice(4, 1)
      // remove epi score and seq score
      this.headers.splice(5, 2)
      this.headers.splice(5, 0, {
        text: 'EpiAlign hit coordinate',
        value: 'scoreE',
        align: 'left',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      }, {
        text: 'Sequence-only hit coordinate',
        value: 'scoreS',
        align: 'left',
        sortable: true,
        class: 'dataTableCell',
        width: 'auto'
      })
    },
    getUcscLink: function (assembly, region) {
      region = window.encodeURIComponent(region.replace(/\([+-]\)/, ''))
      return UCSC_TARGET + 'hgS_otherUserSessionName=' +
        ucscRefToSession[assembly] +
        '&position=' + region
    },
    getGiveLink: function (refs, regions, highlight) {
      let trackPart = ''
      if (this.displayDataSets) {
        trackPart = '&track=' +
          window.encodeURIComponent(JSON.stringify(this.displayDataSets))
      }
      return GIVE_TARGET +
        'ref=' + window.encodeURIComponent(JSON.stringify(refs)) +
        '&coordinate=' + window.encodeURIComponent(JSON.stringify(regions)) +
        '&highlight=' +
        window.encodeURIComponent(JSON.stringify(highlight)) +
        trackPart
    }
  }
})
