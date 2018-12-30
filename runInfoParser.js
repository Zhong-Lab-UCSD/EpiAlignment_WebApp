const moment = require('moment-timezone')

const fs = require('fs')
const util = require('util')

const writeFilePromise = util.promisify(fs.writeFile)

const runInfoDictPath = 'html/assets/runInfoDict.json'
const MAX_INPUT_LENGTH = 35

const runInfoArray = [
  {
    key: 'runid',
    name: 'Run ID',
    description: 'Unique ID for the run.'
  },
  {
    key: 'alignMode',
    name: 'Align mode',
    description: 'The EpiAlignment mode for the run.',
    values: [
      {
        key: 'enhancer',
        name: 'Enhancer mode'
      },
      {
        key: 'promoter',
        name: 'Promoter mode'
      }
    ]
  },
  {
    key: 'subMode',
    name: 'Target type',
    description: 'Types of target regions.',
    values: [
      {
        key: 'genecluster',
        name: 'Gene cluster',
        depends: {
          key: 'alignMode',
          value: 'promoter'
        }
      },
      {
        key: 'homoregion',
        name: 'Homologous region',
        depends: {
          key: 'alignMode',
          value: 'enhancer'
        }
      },
      {
        key: 'genomeregion',
        name: 'Genome regions or gene names'
      }
    ],
    formFields: [
      {
        key: 'searchRegionMode'
      }
    ]
  },
  {
    key: 'publicDataDesc',
    name: 'Public data used',
    description: 'Description for selected public (ENCODE / GEO) dataset.'
  },
  {
    key: 'queryGenomeAssembly',
    name: 'Query reference',
    description: 'Reference for query region.',
    formFields: [
      {
        key: 'genomeAssembly',
        index: 0
      }
    ]
  },
  {
    key: 'queryPeak',
    name: 'Query peak file',
    description: 'Peak file name, or ENCODE / GEO ID for the query.',
    formFields: [
      {
        key: 'encodeData',
        index: 0
      },
      {
        key: 'speciesPeak1[]',
        index: 0,
        property: 'originalname',
        transformFunc: (value, runInfoResults) => (
          value.length > MAX_INPUT_LENGTH
            ? value.slice(0, MAX_INPUT_LENGTH) + ' ...'
            : value
        )
      }
    ]
  },
  {
    key: 'queryInput',
    name: 'Query input file',
    description: 'Input file name, or first line of input text, for the query.',
    formFields: [
      {
        key: 'speciesInput1',
        index: 0,
        property: 'originalname',
        transformFunc: (value, runInfoResults) => (
          value.length > MAX_INPUT_LENGTH
            ? value.slice(0, MAX_INPUT_LENGTH) + ' ...'
            : value
        )
      },
      {
        key: 'speciesText',
        index: 0,
        transformFunc: (value, runInfoResults) => (
          (value || '').split(/\r?\n/, 1)[0].slice(0, MAX_INPUT_LENGTH) +
          ' ...'
        )
      }
    ]
  },
  {
    key: 'targetGenomeAssembly',
    name: 'Target reference',
    description: 'Reference for target region.',
    formFields: [
      {
        key: 'genomeAssembly',
        index: 1
      }
    ]
  },
  {
    key: 'targetPeak',
    name: 'Target peak file',
    description: 'Peak file name, or ENCODE / GEO ID for the target.',
    formFields: [
      {
        key: 'encodeData',
        index: 1
      },
      {
        key: 'speciesPeak2[]',
        index: 0,
        property: 'originalname',
        transformFunc: (value, runInfoResults) => (
          value.length > MAX_INPUT_LENGTH
            ? value.slice(0, MAX_INPUT_LENGTH) + ' ...'
            : value
        )
      }
    ]
  },
  {
    key: 'targetInput',
    name: 'Target input file',
    description: 'Input file name, or first line of input text, for the target.',
    formFields: [
      {
        key: 'speciesInput2',
        index: 0,
        property: 'originalname',
        transformFunc: (value, runInfoResults) => (
          value.length > MAX_INPUT_LENGTH
            ? value.slice(0, MAX_INPUT_LENGTH) + ' ...'
            : value
        )
      },
      {
        key: 'speciesText',
        index: 1,
        transformFunc: (value, runInfoResults) => (
          (value || '').split(/\r?\n/, 1)[0].slice(0, MAX_INPUT_LENGTH) +
          ' ...'
        )
      }
    ],
    depends: {
      key: 'subMode',
      value: 'genomeregion'
    }
  },
  {
    key: 'clusters',
    name: 'Target gene clusters',
    description: 'The gene cluster ID for target.',
    depends: {
      key: 'subMode',
      value: 'genecluster'
    }
  },
  {
    key: 'submitTime',
    name: 'Job submitted at',
    description: 'Time when the job was submitted (ISO8601 String).',
    generatorFunc: runInfoResults => moment().toISOString(),
    displayFunc: (entry, runInfoResults) =>
      moment(entry).tz(runInfoResults.timeZone || moment.tz.guess()).format(
        'MMM D, Y, HH:mm:ss'
      )
  },
  {
    key: 'promoterUp',
    name: 'Promoter upstream length',
    description: 'Number of bps upstream TSS as promoter threshold.',
    depends: {
      key: 'alignMode',
      value: 'promoter'
    },
    displayFunc: (entry, runInfoResults) => parseInt(entry) + ' bp'
  },
  {
    key: 'promoterDown',
    name: 'Promoter downstream length',
    description: 'Number of bps downstream TSS as promoter threshold.',
    depends: {
      key: 'alignMode',
      value: 'promoter'
    },
    displayFunc: (entry, runInfoResults) => parseInt(entry) + ' bp'
  },
  {
    key: 'enhancerUp',
    name: 'Enhancer upstream extension',
    description: 'Number of upstream bps for query region extension before remapping to the target reference.',
    depends: {
      key: 'subMode',
      value: 'homoregion'
    },
    displayFunc: (entry, runInfoResults) => parseInt(entry) + ' bp'
  },
  {
    key: 'enhancerDown',
    name: 'Enhancer downstream extension',
    description: 'Number of downstream bps for query region extension before remapping to the target reference.',
    depends: {
      key: 'subMode',
      value: 'homoregion'
    },
    displayFunc: (entry, runInfoResults) => parseInt(entry) + ' bp'
  },
  {
    key: 'piA',
    name: 'pi_A',
    html: '<em>&#960;</em><sub>A</sub>',
    parameter: true,
    description: 'EpiAlignment parameter pi_A',
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'piC',
    name: 'pi_C',
    html: '<em>&#960;</em><sub>C</sub>',
    parameter: true,
    description: 'EpiAlignment parameter pi_C',
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'piG',
    name: 'pi_G',
    html: '<em>&#960;</em><sub>G</sub>',
    parameter: true,
    description: 'EpiAlignment parameter pi_G',
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'piT',
    name: 'pi_T',
    html: '<em>&#960;</em><sub>T</sub>',
    parameter: true,
    description: 'EpiAlignment parameter pi_T',
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'pi1',
    name: 'pi_1',
    html: '<em>&#960;</em><sub>1</sub>',
    parameter: true,
    description: 'EpiAlignment parameter pi_1',
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'paraS',
    name: 's',
    html: '<em>s</em>',
    parameter: true,
    description: 'EpiAlignment parameter s',
    formFields: [
      {
        key: 'paras'
      }
    ],
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'paraMu',
    name: 'mu',
    html: '<em>&#956;</em>',
    parameter: true,
    description: 'EpiAlignment parameter mu',
    formFields: [
      {
        key: 'paramu'
      }
    ],
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'paraK',
    name: 'kappa',
    html: '<em>&#954;</em>',
    parameter: true,
    description: 'EpiAlignment parameter kappa',
    formFields: [
      {
        key: 'parak'
      }
    ],
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'seqWeight',
    name: 'Sequence weight',
    description: 'Weight for sequence contribution in EpiAlignment',
    formFields: [
      {
        key: 'seqweight'
      }
    ],
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'epiWeight',
    name: 'Epigenome weight',
    description: 'Weight for epigenome contribution in EpiAlignment',
    formFields: [
      {
        key: 'epiweight'
      }
    ],
    displayFunc: (value, runInfoResults) =>
      parseFloat(Number.parseFloat(value).toFixed(3))
  },
  {
    key: 'email',
    name: 'User email',
    description: 'User email',
    formFields: [
      {
        key: 'mail'
      }
    ]
  },
  {
    key: 'completeTime',
    name: 'Job completed at',
    description: 'The time when job is completed (ISO8601 String)',
    displayFunc: (entry, runInfoResults) =>
      moment(entry).tz(runInfoResults.timeZone || moment.tz.guess()).format(
        'MMM D, Y, HH:mm:ss'
      )
  },
  {
    key: 'timeZone',
    name: 'User timezone',
    description: 'User\'s time zone'
  }
]

const runInfoDict = runInfoArray.reduce((currObj, infoEntry) => {
  currObj[infoEntry.key] = infoEntry
  return currObj
}, {})

if (!fs.existsSync(runInfoDictPath)) {
  writeFilePromise(runInfoDictPath, JSON.stringify(runInfoDict, null, 2))
}

const propertyListMap = {
  'email': [
    'runid',
    'submitTime',
    'completeTime',
    'alignMode',
    'publicDataDesc',
    'queryGenomeAssembly',
    'queryPeak',
    'queryInput',
    'targetGenomeAssembly',
    'subMode',
    'targetPeak',
    'targetInput',
    'clusters'
  ],
  '_default': [
    'runid',
    'submitTime',
    'completeTime',
    'alignMode'
  ]
}

class RunInfo {
  constructor (formBody, formFiles, runid) {
    if (typeof formBody === 'string') {
      this._resultFormatted = JSON.parse(formBody)
      this._result = {}
      if (this._resultFormatted.hasOwnProperty('status')) {
        this.status = this._resultFormatted.status
      }
    } else {
      let formData = Object.assign({}, formBody, formFiles, { runid: runid })
      this._resultFormatted = {}
      this._result = {}
      runInfoArray.forEach(infoEntry => {
        let key = infoEntry.key
        if (infoEntry.generatorFunc) {
          this._result[key] = infoEntry.generatorFunc(this._result)
        } else {
          // no generator, needs to parse from formData
          let valueToWrite = null

          // check entry dependency
          if (!infoEntry.depends || (
            this._result.hasOwnProperty(infoEntry.depends.key) &&
            this._result[infoEntry.depends.key] === infoEntry.depends.value
          )) {
            if (infoEntry.formFields) {
              infoEntry.formFields.some(field => {
                if (formData.hasOwnProperty(field.key)) {
                  let fieldValue = formData[field.key]
                  if (fieldValue && field.index !== undefined) {
                    fieldValue = fieldValue[field.index]
                  }
                  if (fieldValue && field.property !== undefined) {
                    fieldValue = fieldValue[field.property]
                  }
                  if (fieldValue && field.transformFunc) {
                    fieldValue = field.transformFunc(fieldValue)
                  }
                  if (fieldValue) {
                    valueToWrite = fieldValue
                    return true
                  }
                  return false
                }
              })
            } else if (formData[key] !== undefined && formData[key] !== '') {
              valueToWrite = formData[key]
            }
            if (infoEntry.values) {
              // needs to verify if values are part of available values, and
              // their dependency (if any ) have been met
              let valueCandidate = valueToWrite
              valueToWrite = null
              infoEntry.values.some(listedCandidate => {
                if (valueCandidate === listedCandidate.key) {
                  if (!listedCandidate.depends || (
                    this._result.hasOwnProperty(listedCandidate.depends.key) &&
                    this._result[listedCandidate.depends.key] ===
                      listedCandidate.depends.value
                  )) {
                    valueToWrite = valueCandidate
                  }
                  return true
                }
                return false
              })
            }
          }

          if (valueToWrite !== null) {
            this._result[key] = valueToWrite
          }
        }
      })
    }
  }

  addProperty (propName, value) {
    if (value && !this._resultFormatted.hasOwnProperty(propName) &&
      !this._result.hasOwnProperty(propName)
    ) {
      this._result[propName] = value
    }
  }

  toJSON () {
    let jsonResult = {}
    for (let key in this._result) {
      let dispResult = this._result[key]
      if (runInfoDict.hasOwnProperty(key) && runInfoDict[key].displayFunc) {
        dispResult = runInfoDict[key].displayFunc(dispResult, this._result)
      }
      jsonResult[key] = dispResult
    }
    return Object.assign(jsonResult, this._resultFormatted)
  }

  static getRunInfo (formBody, formFiles) {
    return new this(formBody, formFiles)
  }

  static fromJSON (jsonString) {
    return new this(jsonString)
  }

  toString (groupKey) {
    groupKey = groupKey || '_default'
    if (!propertyListMap.hasOwnProperty(groupKey)) {
      groupKey = '_default'
    }
    let result = ''
    propertyListMap[groupKey].forEach(key => {
      let entry = this.getPropertyLine(key)
      if (entry) {
        result += entry + '\n'
      }
    })
    return result
  }

  getDisplayValue (property) {
    if (this._resultFormatted.hasOwnProperty(property)) {
      return this._resultFormatted[property]
    }
    if (this._result.hasOwnProperty(property)) {
      let dispResult = this._result[property]
      if (runInfoDict.hasOwnProperty(property)) {
        if (runInfoDict[property].displayFunc) {
          dispResult = runInfoDict[property].displayFunc(
            dispResult, this._result
          )
        }
        if (runInfoDict[property].values) {
          runInfoDict[property].values.some(valueEntry => {
            if (valueEntry.key === dispResult) {
              dispResult = valueEntry.name
              return true
            }
            return false
          })
        }
      }
      return dispResult
    }
  }

  getPropertyLine (property) {
    if (
      (
        this._result.hasOwnProperty(property) ||
        this._resultFormatted.hasOwnProperty(property)
      ) &&
      runInfoDict.hasOwnProperty(property)
    ) {
      return runInfoDict[property].name + ': ' +
        this.getDisplayValue(property)
    }
    return null
  }

  set status (newStatus) {
    this._result.status = newStatus
    if (newStatus >= 0 &&
      !this._resultFormatted.hasOwnProperty('completeTime') &&
      !this._result.hasOwnProperty('completeTime')
    ) {
      this._result.completeTime = moment().toISOString()
    }
  }

  get status () {
    return this._result.status
  }
}

module.exports = RunInfo
