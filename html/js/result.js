import postAjax from './promisedAjax.js'

const setTimeoutPromise = function (timeout, resolvedValue) {
  return new Promise((resolve, reject) =>
    window.setTimeout(resolve, timeout, resolvedValue)
  )
}

const INQUIRY_TARGET_PREFIX = '/backend/results/'
const INQUIRY_IMAGE_PREFIX = '/backend/result_image/'

const STATUS_RUNNING = -1
const POLLING_INTERVAL_ENHANCER = 20000 // 20 seconds
const POLLING_INTERVAL_PROMOTER = 5000 // 5 seconds

const MINIMUM_HEATMAP_GAP = 10

var app = new Vue({
  el: '#result_app',
  data: {
    runid: null,
    loading: true,
    hasEmail: false,
    emailAddress: null,
    downloadLink: '',
    error: false,
    errorMessage: '',
    imageBlobUrl: null,
    alignMode: null,
    pollingTime: POLLING_INTERVAL_PROMOTER,
    promoterColumnList: [],
    epiHeatmapList: [],
    hasSequence: true,
    geneIdentifier1: null, // ['transID1', 'region_name1']
    geneIdentifier2: null, // ['transID2', 'region_name2']
    showHeatmap: false,
    headers: [
      {
        text: '#',
        value: 'id',
        align: 'left',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Name #1',
        value: 'region_name1',
        align: 'left',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Coordinate #1',
        value: 'region1',
        align: 'left',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Name #2',
        value: 'region_name2',
        align: 'left',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Coordinate #2',
        value: 'region2',
        align: 'left',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Epi Score',
        value: 'scoreE',
        align: 'right',
        sortable: true,
        width: 'auto'
      },
      {
        text: 'Seq Score',
        value: 'scoreS',
        align: 'right',
        sortable: true,
        width: 'auto'
      }
    ],
    dataEntries: []
  },
  created: function () {
    const urlParams = new window.URLSearchParams(window.location.search)
    this.runid = urlParams.get('runid')
    this.pollResult(this.runid).then(data => this.showData(data)
    ).catch(err => {
      // TODO: err.status should say something
      this.loading = false
      this.error = true
      this.errorMessage = err.message
    })
  },
  methods: {
    pollResult: function (runid) {
      // TODO: validate formParams
      if (!runid) {
        return Promise.reject(new Error('✖ Error: runid does not exist! ' +
          'Please try submit your data again.'))
      }
      return postAjax(INQUIRY_TARGET_PREFIX + runid, null, 'json', 'GET')
        .then(response => {
          if (response.status === STATUS_RUNNING) {
            return setTimeoutPromise(this.pollingTime, runid)
              .then(runid => this.pollResult(runid))
          } else if (response.status > 0) {
            // error, reject the promise and let .catch to handle
            return Promise.reject(response)
          } else {
            this.setRunParameters(response)
            return response.data
          }
        })
    },
    setRunParameters: function (response) {
      document.title = 'EpiAlignment - Result (runID: ' + this.runid + ')'

      this.alignMode = response.alignMode
      this.downloadLink = window.location.protocol + '//' +
        window.location.host + '/download/' + this.runid + '.txt'
      this.loading = false
    },
    showData: function (dataEntries) {
      if (Array.isArray(dataEntries)) {
        this.dataEntries = dataEntries
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
      })
    },
    verifyPromoterNameSource: function (dataEntry) {
      if (!this.geneIdentifier1) {
        this.geneIdentifier1 =
          dataEntry.transID1 === '.' ? 'region_name1' : 'transID1'
      }
      if (!this.geneIdentifier2) {
        this.geneIdentifier2 =
          dataEntry.transID2 === '.' ? 'region_name2' : 'transID2'
      }
      this.hasSequence = dataEntry.scoreS !== '.'
    },
    processPromoterData: function (data) {
      if (data[0]) {
        this.verifyPromoterNameSource(data[0])
      }
      let columnMap = new Map()
      let rowMap = new Map()
      data.forEach(entry =>
        this.parseSinglePromoterEntry(entry, columnMap, rowMap))
      this.processHeatmap()
      if (this.hasSequence) {
        this.markRowMaxFlags()
      }
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
      let normFactors =
        this.findMinMaxTableValue(this.epiHeatmapList, propertyList)
      this.epiHeatmapList.forEach(entry => entry.values.forEach(value => (
        value.epiScoreNorm =
          (value.epiScore - normFactors.min) / normFactors.span
      )))
      if (this.hasSequence) {
        this.epiHeatmapList.forEach(entry => entry.values.forEach(value => (
          value.seqScoreNorm =
            (value.seqScore - normFactors.min) / normFactors.span
        )))
      }
    },
    markRowMaxFlags: function () {
      this.epiHeatmapList.forEach(entry => {
        this.markSingleRowMaxFlag(entry, 'epiScore', 'epiMax')
        this.markSingleRowMaxFlag(entry, 'seqScore', 'seqMax')
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
    },
    findMinMaxTableValue: function (table, entryPropertyList) {
      let result = {
        min: Number.MAX_VALUE,
        max: Number.MIN_VALUE
      }
      table.forEach(entry => entry.values.forEach(
        value => entryPropertyList.forEach(prop => {
          if (result.min > value[prop]) {
            result.min = value[prop]
          }
          if (result.max < value[prop]) {
            result.max = value[prop]
          }
        })
      ))
      if (result.max - result.min < MINIMUM_HEATMAP_GAP) {
        result.span = MINIMUM_HEATMAP_GAP
      } else {
        result.span = result.max - result.min
      }
      delete result.max
      return result
    },
    processEnhancerData: function (data) {
    }
  }
})
