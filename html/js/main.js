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
    },
  },
  methods: {
    toggleShowParam: function () {
      this.showParam = !this.showParam
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
    submitForm: function () {
      // TODO: validate formParams
      let xhr = new window.XMLHttpRequest()
      let formData = new window.FormData(this.$refs.mainForm)
      //xhr.responseType = 'json'
      xhr.onload = function () {
        if (xhr.status >= 200 && xhr.status < 400) {
          window.alert('SUCCESS: \n' + xhr.response)
        } else {
          window.alert('ERROR: ' + xhr.status)
       }
      }
      xhr.onerror = function () {
        window.alert('ERROR: ' + xhr.status)
      }
      xhr.open('POST', FORM_SUBMIT_TARGET)
      xhr.send(formData)
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
  },
})

