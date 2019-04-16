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
        tooltip: 'Whether the EpiAlign hit is different from the SeqOnly hit.'
      }
    ],
    dataEntries: [],
    expList: null,
    seqBg: null,
    figProp: {
      crossSize: 4,
      vertMargin: 8,
      lineHeight: 36,
      legendLeftX: 10.5,
      legendRightX: 240.5,
      legendBoxWidth: 30,
      width: 500,
      horizMargin: 20,
      textSize: 14,
      tickSize: 5,
      axisLabelSize: 12,
      maxNumOfTicks: 10,
      epiContribLabelSize: 50
    },
    expFigProp: {
      geneHeight: 30,
      refHeight: 20,
      axisHeight: 45,
      gap: 5,
      textAreaWidth: 80,
      textMargin: 5,
      textSize: 14,
      tickSize: 5,
      axisLabelSize: 12,
      maxNumOfTicks: 10,
      width: 350,
      horizMargin: 10
    },
    params: {
      kappa: null,
      pi1: null,
      epiWeight: null
    },
    palette: [
      0x4477AA, 0xEE6677, 0x228833,
      0xCCBB44, 0x66CCEE, 0xAA3377
    ]
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
      this.seqBg = response.seqBg
      this.params.kappa = parseFloat(response.paraK)
      this.params.pi1 = parseFloat(response.pi1)
      this.params.epiWeight = parseFloat(response.epiWeight)

      if (response.expression) {
        this.expList = response.expression
      }
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
        text: 'SeqOnly hit coordinate',
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
    },
    getXSeqScore: function (item, xValue) {
      let xMin = this.getSeqScoreMin(item)
      let xMax = this.getSeqScoreMax(item)
      let xScale = (this.figProp.width - 2 * this.figProp.horizMargin) /
        (xMax - xMin)
      return (parseFloat(xValue) - xMin) * xScale + this.figProp.horizMargin
    },
    getSeqScoreMax: function (item) {
      let xMax
      xMax = Math.max(parseFloat(this.seqBg.backgroundQ75),
        parseFloat(this.seqBg.orthoQ75))
      if (typeof item.scoreS2 === 'number') {
        xMax = Math.max(item.scoreS2, xMax)
      }
      if (typeof item.scoreS === 'number') {
        xMax = Math.max(item.scoreS, xMax)
      }
      return xMax
    },
    getSeqScoreMin: function (item) {
      let xMin = Math.min(parseFloat(this.seqBg.backgroundQ25),
        parseFloat(this.seqBg.orthoQ25))
      if (typeof item.scoreS2 === 'number') {
        xMin = Math.min(item.scoreS2, xMin)
      }
      if (typeof item.scoreS === 'number') {
        xMin = Math.min(item.scoreS, xMin)
      }
      return xMin
    },
    getTickList: function (xMin, xMax, maxNumOfTicks) {
      let length = xMax - xMin
      let span = length / maxNumOfTicks
      if (Math.ceil(span) > 0) {
        // round up to closest [1,2,5] * 10^x
        let spanExp = parseInt(Math.log(span) / Math.LN10)
        let spanHeader = span / Math.pow(10, spanExp)
        if (spanHeader > 5) {
          spanExp++
          spanHeader = 1
        } else if (spanHeader > 2) {
          spanHeader = 5
        } else if (spanHeader > 1) {
          spanHeader = 2
        }
        span = spanHeader * Math.pow(10, spanExp)
      }
      span = Math.ceil(span)
      if (span <= 0) {
        span = 1
      }

      let currValue = Math.ceil(xMin / span) * span
      let result = [currValue]
      currValue += span
      while (currValue < xMax) {
        result.push(currValue)
        currValue += span
      }
      return result
    },
    getLineY: function (lineNum) {
      return this.figProp.lineHeight * lineNum
    },
    getOneNum: function (item) {
      let oneNum = parseInt(item.oneNum)
      if (oneNum > parseInt(item.queryLength)) {
        oneNum = parseInt(item.queryLength)
      } else if (oneNum < 0) {
        oneNum = 0
      }
      return oneNum
    },
    getEpiContribMin: function (item) {
      return Math.log(1 - Math.exp(-this.params.kappa)) *
        1000 * this.params.epiWeight
    },
    getEpiContribMax: function (item) {
      let oneNum = this.getOneNum(item)
      let zeroNum = parseInt(item.queryLength) - oneNum
      let kappa = this.params.kappa
      let pi1 = this.params.pi1
      let pi0 = 1 - pi1
      return (oneNum * Math.log(
        (Math.exp(-kappa) + (1 - Math.exp(-kappa) * pi1)) / pi1
      ) + zeroNum * Math.log(
        (Math.exp(-kappa) + (1 - Math.exp(-kappa) * pi0)) / pi0
      )) / item.queryLength * 1000 * this.params.epiWeight
    },
    getXEpiScore: function (item, epiScore) {
      let xLength = this.figProp.width -
        2 * (this.figProp.horizMargin + this.figProp.epiContribLabelSize)
      let xMin = this.getEpiContribMin(item)
      let xMax = this.getEpiContribMax(item)
      let xScale = xLength / (xMax - xMin)
      return (epiScore - xMin) * xScale +
        (this.figProp.horizMargin + this.figProp.epiContribLabelSize)
    },
    getEpiContribEpiHit (item) {
      if (item.shifted === 'Y') {
        return item.scoreE - item.scoreS2
      }
      return item.scoreE - item.scoreS
    },
    getEpiContribSeqHit (item) {
      if (item.shifted === 'Y') {
        return item.scoreE2 - item.scoreS
      }
      return item.scoreE - item.scoreS
    },
    getXExpression: function (expList, expValue) {
      let xLength = this.expFigProp.width -
        (this.expFigProp.horizMargin + this.expFigProp.textAreaWidth)
      let xMin = 0
      let xMax = this.getExpressionMax(expList)
      let xScale = xLength / (xMax - xMin)
      return (expValue - xMin) * xScale + this.expFigProp.textAreaWidth
    },
    getExpressionMax: function (expList) {
      return expList.reduce((prevInRefMax, expInRef) => {
        for (let key in expInRef.expValues) {
          let geneExpMax =
            Math.max(...expInRef.expValues[key].map(expVal => expVal.FPKM))
          if (prevInRefMax < geneExpMax) {
            prevInRefMax = geneExpMax
          }
        }
        return prevInRefMax
      }, 0)
    },
    getExpressionHeight: function (expList) {
      return expList.reduce((prev, expInRef) => {
        return prev + this.getExpInRefHeight(expInRef)
      }, 0)
    },
    getExpInRefHeight: function (expInRef) {
      return Object.keys(expInRef.expValues).length *
        (this.expFigProp.geneHeight + this.expFigProp.gap) +
        this.expFigProp.refHeight
    },
    getExpInRefY: function (expList, refIndex) {
      let currY = 0
      for (let i = 0; i < refIndex; i++) {
        currY += this.getExpInRefHeight(expList[i])
      }
      return currY
    },
    capitalize: function (str) {
      if (typeof str !== 'string') {
        return ''
      }
      return str.charAt(0).toUpperCase() + str.slice(1)
    },
    getColorString: function (expList, refIndex, expIndex) {
      let colorIndex = 0
      for (let i = 0; i < refIndex; i++) {
        colorIndex += expList[i].expIds.length
      }
      return this.palette[(colorIndex + expIndex) % this.palette.length]
        .toString(16)
    },
    buildLinkFromExpressionId: function (id) {
      if (id.startsWith('ENC')) {
        // ENCODE id
        if (id.startsWith('ENCFF')) {
          return 'https://www.encodeproject.org/files/' + id + '/'
        } else if (id.startsWith('ENCSR')) {
          return 'https://www.encodeproject.org/experiments/' + id + '/'
        } else {
          return 'https://www.encodeproject.org/search/?searchTerm=' + id
        }
      } else if (id.startsWith('GSM')) {
        // GEO id
        return 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + id
      }
      return null
    }
  }
})
