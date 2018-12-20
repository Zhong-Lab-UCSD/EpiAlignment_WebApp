const setTimeoutPromise = function (timeout, resolvedValue) {
  return new Promise((resolve, reject) => 
    window.setTimeout(resolve, timeout, resolvedValue)
  )
}
import postAjax from './promisedAjax.js'

const INQUIRY_TARGET_PREFIX = '/backend/results/'

const STATUS_RUNNING = -1
const POLLING_INTERVAL = 15000    // 15 seconds

var app = new Vue({
  el: '#result_app',
  data: {
    showData: null
  },
  created: function () {
    const urlParams = new URLSearchParams(window.location.search)
    let runid = urlParams.get('runid')
    this.pollResult(runid).then(data => {
      // TODO: process data
      this.showData = JSON.stringify(data)
    })
  },
  methods: {
    pollResult: function (runid) {
      // TODO: validate formParams
      if (!runid) {
        this.errorMsg = '✖ Error: runid does not exist! ' +
          'Please try submit your data again.'
        return Promise.reject()
      }
      return postAjax(INQUIRY_TARGET_PREFIX + runid, null, 'json', 'GET')
        .then(response => {
          if (response.status === STATUS_RUNNING) {
            // TODO: display waiting message

            return setTimeoutPromise(POLLING_INTERVAL)
              .then(this.pollResult(runid))
          } else if (response.status > 0) {
            // error, reject the promise and let .catch to handle
            return Promise.reject(response)
          } else {
            return response.data
          }
        })
        .catch(err => {
          // err.status should say something
        })
    }
  },
})