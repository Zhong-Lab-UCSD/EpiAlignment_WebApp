export default function postAjax (
  target, params, responseType, method, additionalHeaders
) {
  // this is a wrapper for Ajax calls throughout
  return new Promise((resolve, reject) => {
    method = method || 'POST'
    var xhr = new window.XMLHttpRequest()
    xhr.responseType = responseType || ''
    xhr.onload = function () {
      var responses = this.response
      if (this.status >= 200 && this.status < 400) {
        if (this.responseType.toLowerCase() === 'json' &&
          (navigator.appName === 'Microsoft Internet Explorer' ||
            !!(navigator.userAgent.match(/Trident/) ||
              navigator.userAgent.match(/rv 11/)))) {
          // IE detected (should be IE 11), fix the json return issue
          let errorMsg = 'You are currently using IE 11 to visit this ' +
            'site. Some part of the site may behave differently and if ' +
            'you encounter any problems, please use the info on ' +
            '\'About us\' page to contact us.'
          console.error(errorMsg)
          responses = JSON.parse(responses)
        }
        resolve(responses)
      } else {
        let err = new Error('Connection error (' + this.status +
          ')' +
          (this.response
            ? ': ' + (typeof this.response === 'object'
              ? this.response.message || '' : this.response)
            : ''))
        err.status = this.status
        err.response = this.response
        reject(err)
      }
    }
    xhr.onerror = () => {
      let err = new Error('Connection error (' + this.status +
        ')' + (this.response ? ': ' + this.response : ''))
      err.status = this.status
      reject(err)
    }
    xhr.open(method, target)
    if (params instanceof window.FormData) {
      xhr.send(params)
    } else if (params) {
      xhr.setRequestHeader('Content-Type', 'application/json')
      xhr.send(JSON.stringify(params))
    } else {
      xhr.send()
    }
  })
}
