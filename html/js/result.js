import postAjax from './promisedAjax.js'

const setTimeoutPromise = function (timeout, resolvedValue) {
  return new Promise((resolve, reject) =>
    window.setTimeout(resolve, timeout, resolvedValue)
  )
}

const INQUIRY_TARGET_PREFIX = '/backend/results/'
const INQUIRY_IMAGE_PREFIX = '/backend/result_image/'

const STATUS_RUNNING = -1
const POLLING_INTERVAL = 15000 // 15 seconds

const imageTag = [
  { requestId: 's', itemDataLink: 'imageSData' },
  { requestId: 'e', itemDataLink: 'imageEData' },
]

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
    imageBlobUrls: {},
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
    this.pollResult(this.runid).then(data => {
      // TODO: load the submit time from server and update title accordingly
      document.title = 'EpiAlignment - Result (runID: ' + this.runid + ')'

      this.loading = false
      this.downloadLink = window.location.protocol + '//' +
        window.location.host + '/download/' + this.runid + '.txt'
      this.showData(data)
    }).catch(err => {
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
            return setTimeoutPromise(POLLING_INTERVAL, runid)
              .then(runid => this.pollResult(runid))
          } else if (response.status > 0) {
            // error, reject the promise and let .catch to handle
            return Promise.reject(response)
          } else {
            return response.data
          }
        })
    },
    showData: function (dataEntries) {
      if (Array.isArray(dataEntries)) {
        this.dataEntries = dataEntries
      }
    },
    toggleRow: function (row) {
      row.expanded = !row.expanded
      imageTag.map(tagObj => {
        Vue.set(row.item, tagObj.itemDataLink, null)
      })
      if (row.expanded) {
        return Promise.all(imageTag.map(tagObj => {
          let imageLink = INQUIRY_IMAGE_PREFIX + this.runid + '/' +
            row.item.index + '/' + tagObj.requestId + '.png'
          return postAjax(
            imageLink, null, 'blob', 'GET'
          ).then(imageBlob => {
            if (this.imageBlobUrls[tagObj.requestId]) {
              window.URL.revokeObjectURL(this.imageBlobUrls[tagObj.requestId])
            }
            this.imageBlobUrls[tagObj.requestId] =
              window.URL.createObjectURL(imageBlob)
            Vue.set(row.item, tagObj.itemDataLink,
              this.imageBlobUrls[tagObj.requestId])
          })
        }))
      }
    }
  }
})
