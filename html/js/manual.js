const SLIDE_TRANSITION_DURATION = 300

/**
 * Swipe detection
 * Adapted from https://stackoverflow.com/a/27115070
 * @param {string} el element ID
 * @param {function} func callback function
 */
function detectSwipe (el, func) {
  let swipeDetector = {
    sX: 0,
    sY: 0,
    eX: 0,
    eY: 0
  }
  let minDist = 30 // min swipe distance
  let direction = '' // [udlr] for up, down, left, right
  let elem = document.getElementById(el)
  elem.addEventListener('touchstart', function (e) {
    var t = e.touches[0]
    swipeDetector.sX = t.screenX
    swipeDetector.sY = t.screenY
  }, false)
  elem.addEventListener('touchmove', function (e) {
    var t = e.touches[0]
    swipeDetector.eX = t.screenX
    swipeDetector.eY = t.screenY
  }, false)
  elem.addEventListener('touchend', function (e) {
    // find the axis
    if (Math.abs(swipeDetector.eY - swipeDetector.sY) >
      Math.abs(swipeDetector.eX - swipeDetector.sX)
    ) {
      // vertical
      if (Math.abs(swipeDetector.eY - swipeDetector.sY) > minDist) {
        direction = swipeDetector.eY > swipeDetector.sY ? 'd' : 'u'
      }
    } else {
      // horizontal
      if (Math.abs(swipeDetector.eX - swipeDetector.sX) > minDist) {
        direction = swipeDetector.eX > swipeDetector.sX ? 'r' : 'l'
      }
    }
    if (direction) {
      if (typeof func === 'function') {
        func(el, direction)
      }
    }
    direction = ''
    swipeDetector.sX = 0; swipeDetector.sY = 0; swipeDetector.eX = 0; swipeDetector.eY = 0
  }, false)
}

var app = new Vue({
  el: '#manual_app',
  data: {
    showNav: true
  },
  created () {
    if (window.matchMedia('screen and (max-width: 800px)').matches) {
      this.showNav = false
    }
  },
  mounted () {
    detectSwipe('manual_app', (el, direction) => {
      if (direction === 'r') {
        this.toggleNav(true)
      } else if (direction === 'l') {
        this.toggleNav(false)
      }
    })

    if (window.location.hash) {
      // jump to hashed location
      this.jumpToHash(window.location.hash)
    }
  },
  methods: {
    clickOnNav: function (event) {
      event.preventDefault()
      if (window.matchMedia('screen and (max-width: 800px)').matches) {
        this.showNav = false
      }
      window.location.hash = event.target.hash
      this.jumpToHash(window.location.hash)
    },

    jumpToHash: function (hash) {
      hash = hash.slice(1)
      if (this.$refs.hasOwnProperty(hash)) {
        window.requestAnimationFrame(() =>
          this.$refs[hash].scrollIntoView(true)
        )
      }
    },

    toggleNav: function (flag) {
      this.showNav = flag
      if (!window.matchMedia('screen and (max-width: 800px)').matches) {
        // navigation will affect the size of content, needs to keep proportion
        let timeStart = null
        let topElem = null
        let findFirstVisibleNode = (node, scrollTop) => {
          if (node.hasChildNodes()) {
            if ([...node.childNodes].some(childNode =>
              findFirstVisibleNode(childNode, scrollTop))
            ) {
              return true
            }
          }
          if (node.offsetTop > scrollTop) {
            topElem = node
            return true
          }
          return false
        }
        findFirstVisibleNode(
          this.$refs.manualContent, this.$refs.manualContent.scrollTop
        )
        let topDistance = topElem
          ? this.$refs.manualContent.scrollTop - topElem.offsetTop
          : 0
        let updateScroll = timeStamp => {
          if (!timeStart) {
            timeStart = timeStamp
          }
          let progress = timeStamp - timeStart
          if (topElem) {
            this.$refs.manualContent.scrollTop = topElem.offsetTop + topDistance
          }
          if (progress <= SLIDE_TRANSITION_DURATION) {
            window.requestAnimationFrame(updateScroll)
          }
        }
        if (topElem) {
          window.requestAnimationFrame(updateScroll)
        }
      }
    }
  }
})
