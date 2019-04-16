const express = require('express')
const multer = require('multer')
const fs = require('fs')
const nodemailer = require('nodemailer')
const upload = multer({ dest: 'uploads/' })
const spawn = require('child_process').spawn
// express.js
const app = express()
// python
const pythonScript = 'server_agent.py'
const pythonPlotScript = 'Plot_selected_region.py'
// util
const util = require('util')
const path = require('path')
const readFilePromise = util.promisify(fs.readFile)
const writeFilePromise = util.promisify(fs.writeFile)
const readDirPromise = util.promisify(fs.readdir)
const renamePromise = util.promisify(fs.rename)

const loadSpeciesGeneAnno = require('./annoParser').loadSpeciesGeneAnno

const ClusterProcesser = require('./clusterProcessing')

const expKeyIndexEncode = {
  'ensemblGeneId': 0,
  'pme_TPM': 9,
  'pme_FPKM': 10,
  'TPM': 5,
  'FPKM': 6
}

const expKeyIndexGeo = {
  'ensemblGeneId': 1,
  'FPKM': 2
}

const species = [
  {
    name: 'human',
    latin: 'Homo_sapiens',
    reference: 'hg38',
    encode_reference: 'GRCh38'
  }, {
    name: 'mouse',
    latin: 'Mus_musculus',
    reference: 'mm10',
    encode_reference: 'mm10'
  }
]

species._map = species.reduce((prevMap, curr) => {
  prevMap[curr.reference] = curr
  return prevMap
}, {})

var clusterProc = new ClusterProcesser(
  species,
  {
    rawFilePath: 'Annotation/AnnotationFiles',
    clusterSuffix: '_clusters'
  }
)

const RunInfo = require('./runInfoParser')
const resultFolder = 'userResults'

const fallbackHostName = 'epialign.ucsd.edu'
const mailDomainPrefex = 'mail'
const mailDomainSuffix = 'givengine.org'

var runIdDict = {}

const htmlAssetPath = 'html/assets'
const expDictFileName = 'experimentDict.json'
const expDictPath = path.format({
  dir: htmlAssetPath,
  base: expDictFileName
})

const ID_LENGTH = 10

const RUNNING_CODE = -1
const PROCESS_TERMINATED_SERVER_REBOOT = 500
const E_CANNOT_READ_RESULTS = 401

const readExpDictPromise = readFilePromise(expDictPath, 'utf8')
  .then(resultText => JSON.parse(resultText))

const annotationMapPromise = species.reduce((prev, current) => {
  prev[current.reference] =
    loadSpeciesGeneAnno(current, ['symbol', 'ensemblId'], true)
  return prev
}, {})

function getMailDomain (isBeta) {
  return mailDomainPrefex + '.' + (isBeta ? 'beta.' : '') + mailDomainSuffix
}

// Server preparation step: parse the cluster file to build a
//  gene-name/ensemblID to cluster map

function makeid () {
  let text = ''
  let possible = 'abcdefghijklmnopqrstuvwxyz0123456789'

  for (var i = 0; i < ID_LENGTH; i++) {
    text += possible.charAt(Math.floor(Math.random() * possible.length))
  }

  return text
}

const cpUpload = upload.fields([
  { name: 'speciesPeak1[]', maxCount: 3 },
  { name: 'speciesPeak2[]', maxCount: 3 },
  { name: 'speciesInput1', maxCount: 1 },
  { name: 'speciesInput2', maxCount: 1 }
])

function getRunIdTempPath (runid) {
  return path.format({
    dir: resultFolder,
    name: 'tmp_' + runid
  })
}

function getRunIdFilePath (runid, ext) {
  return path.format({
    dir: getRunIdTempPath(runid),
    name: runid,
    ext: ext
  })
}

function getRunIdInfoPath (runid) {
  return getRunIdFilePath(runid, '.runInfo.json')
}

function getRunIdResultPath (runid) {
  return getRunIdFilePath(runid, '.json')
}

function getRunIdImagePath (runid, index) {
  return path.format({
    dir: getRunIdTempPath(runid),
    name: 'Image_' + index + '_' + runid,
    ext: '.png'
  })
}

/**
 * Attach the previously built additional results to Python-returned results
 * @param {string} runid
 * @param {object} [preprocessRes] the result from the preprocessing, can be
 *    `null` (in which case nothing will be added).
 * @returns {Promise<void>} returns a promise that resolves to `void`
 *    when the updated data has been written to a file
 */
function updateExpression (runid, preprocessRes) {
  // Procedures:
  // * Determine if the run is many to many
  // * Determine if there are expression files available
  // * Find all the genes' Ensembl IDs involved in the run
  // * Build an object of all expression values and attach to the result
  if (preprocessRes) {
    let resultFilePath = getRunIdResultPath(runid)
    return readFilePromise(resultFilePath, 'utf8')
      .then(result => {
        let resDataObj = JSON.parse(result)
        if (resDataObj.hasOwnProperty('expression')) {
          console.log('Warning: existing expression property detected. ' +
            'Will be overwritten.')
        }
        resDataObj.expression = preprocessRes
        return writeFilePromise(resultFilePath, JSON.stringify(resDataObj))
      })
      .catch(ignore => { })
  }
  return null
}

/**
 * Get value from the textarea element, the uploaded file, or from cluster ID.
 * The later parameters will take precedence.
 * @param {RunInfo} runInfo text from textarea element
 * @param {object} files file object from multer
 * @param {string} textareaKey key value for textarea values in `runInfo`
 * @param {string} [filesKey] key value for file object
 * @param {string} [clusterId] the ID of the cluster
 * @param {string} [speciesName] the name of the species, must be provided
 *    to use `clusterId`
 * @returns {Promise<string>} a promise resolving to the content
 * @async
 */
async function getValueFromFileOrTextArea (
  runInfo, files, textareaKey, filesKey, clusterId, speciesName
) {
  if (clusterId) {
    await clusterProc.clusterTableReadyPromise
    let result = clusterProc.getClusterById(clusterId)
    if (result) {
      return result.getGeneListBySpecies(speciesName).map(gene => gene.symbol)
    }
  }
  if (files && files[filesKey] && files[filesKey][0]) {
    return readFilePromise(files[filesKey][0].filename, 'utf8')
  }
  return runInfo.getRawPropertyValue(textareaKey)
}

/**
 * Get the gene id list from submitted forms
 * @param {RunInfo} runInfo run info object
 * @param {object} files file object returned from multer
 */
async function getGeneIdList (runInfo, files) {
  let queryTextPromise = getValueFromFileOrTextArea(
    runInfo, files, 'queryInput', 'speciesInput1')
  let targetTextPromise = getValueFromFileOrTextArea(
    runInfo, files, 'targetInput', 'speciesInput2',
    runInfo.getRawPropertyValue('clusters'),
    species._map[runInfo.getRawPropertyValue('targetGenomeAssembly')].name
  )
  let queryContent = await queryTextPromise
  let targetContent = await targetTextPromise
  let getIdFromLine = line => {
    let tokens = line.trim().split(/\t| +/)
    if (tokens.length > 1) {
      // Maybe BED files
      return tokens[4] // name part
    }
    return tokens[0]
  }
  return {
    query: queryContent
      ? (Array.isArray(queryContent)
        ? queryContent
        : (queryContent.split('\n').map(line => getIdFromLine(line))
          .filter(entry => !!entry)
        ))
      : null,
    target: targetContent
      ? (Array.isArray(targetContent)
        ? targetContent
        : (targetContent.split('\n').map(line => getIdFromLine(line))
          .filter(entry => !!entry)
        ))
      : null
  }
}

function buildExpDicts (expFiles) {
  return expFiles.reduce((prevDict, expFile, fileIndex) => {
    expFile.split('\n').forEach(currLine => {
      let tokens = currLine.split('\t')
      let expKeyIndex = null
      if (tokens.length > 10) {
        // ENCODE format
        expKeyIndex = expKeyIndexEncode
      } else {
        expKeyIndex = expKeyIndexGeo
      }
      try {
        let key = tokens[expKeyIndex.ensemblGeneId].split('.', 1)[0]
        if (!prevDict.hasOwnProperty(key)) {
          prevDict[key] = Array(expFiles.length).fill({
            'pme_TPM': 0,
            'pme_FPKM': 0,
            'TPM': 0,
            'FPKM': 0
          })
        }
        for (let expKey in prevDict[key][fileIndex]) {
          if (expKeyIndex.hasOwnProperty(expKey)) {
            if (!prevDict[key][fileIndex][expKey] ||
              prevDict[key][fileIndex][expKey] <
              parseFloat(tokens[expKeyIndex[expKey]])
            ) {
              prevDict[key][fileIndex][expKey] =
                parseFloat(tokens[expKeyIndex[expKey]])
            }
          } else {
            prevDict[key][fileIndex][expKey] = null
          }
        }
      } catch (ignore) { }
    })
    return prevDict
  }, {})
}

function buildExpObjFromFiles (list, expDicts, annotationMap) {
  return list.reduce((prev, curr) => {
    curr = curr.toLowerCase()
    if (annotationMap.has(curr)) {
      let expKey = annotationMap.get(curr).ensemblId
      let symbolKey = annotationMap.get(curr).symbol
      if (!prev.hasOwnProperty(symbolKey)) {
        prev[symbolKey] = expDicts[expKey]
      }
    }
    return prev
  }, {})
}

async function preJobRun (runInfo, files) {
  if (!runInfo ||
    runInfo.getDisplayValue('alignMode') !== 'Many-vs-many mode'
  ) {
    return null
  }
  // Determine if there are expression files available
  // If there are, resolve the promise to the loaded expression files
  let queryPeakId = runInfo.getDisplayValue('queryPeak')
  let targetPeakId = runInfo.getDisplayValue('targetPeak')
  let expDict = await readExpDictPromise
  // read expression files
  let queryExpFilesPromise = null
  if (expDict.hasOwnProperty(queryPeakId) &&
    expDict[queryPeakId].expression_files
  ) {
    queryExpFilesPromise = Promise.all(
      expDict[queryPeakId].expression_files.map(
        fileName => readFilePromise(fileName, 'utf8')
      )
    )
  }
  let targetExpFilesPromise = null
  if (expDict.hasOwnProperty(targetPeakId) &&
    expDict[targetPeakId].expression_files
  ) {
    targetExpFilesPromise = Promise.all(
      expDict[targetPeakId].expression_files.map(
        fileName => readFilePromise(fileName, 'utf8')
      )
    )
  }

  // get gene identifiers
  let geneIdentifiers = await getGeneIdList(runInfo, files)

  let result = {
    queryExpObj: null,
    targetExpObj: null
  }
  if (queryExpFilesPromise && geneIdentifiers.query) {
    let queryExpDicts =
      buildExpDicts(await queryExpFilesPromise)
    result.queryExpObj =
      buildExpObjFromFiles(geneIdentifiers.query, queryExpDicts,
        await annotationMapPromise[runInfo.getRawPropertyValue('queryGenomeAssembly')])
  }
  if (targetExpFilesPromise && geneIdentifiers.target) {
    let targetExpDicts =
      buildExpDicts(await targetExpFilesPromise)
    result.targetExpObj =
      buildExpObjFromFiles(geneIdentifiers.target, targetExpDicts,
        await annotationMapPromise[runInfo.getRawPropertyValue('targetGenomeAssembly')])
  }
  // Find all the genes' Ensembl IDs involved in the run
  // Build an object of all expression values
  if (!result.queryExpObj && !result.targetExpObj) {
    return null
  }
  return result
}

function sendEmail (runInfo) {
  let email = runInfo.getDisplayValue('email')
  if (email) {
    let hostName = runInfo.getDisplayValue('hostName') || fallbackHostName
    let message = {
      from: 'EpiAlignment Notification <messenger@' +
        getMailDomain(hostName.includes('beta')) + '>',
      to: email,
      replyTo: 'x9cao@eng.ucsd.edu'
    }
    let header = 'Dear user,\n\n'
    let footer = '\n\nThe EpiAlignment Team from Zhong Lab @ UC San Diego\n' +
      '--\nDepartment of Bioengineering, Mail Code 0412\n' +
      'University of California, San Diego\n' +
      '9500 Gilman Dr.\nLa Jolla, CA 92122-0412\nUnited States'

    let runid = runInfo.getDisplayValue('runid')
    let code = runInfo.status
    // needs to write an email
    let transporter = nodemailer.createTransport({
      sendmail: true,
      newline: 'unix',
      path: '/usr/sbin/sendmail'
    })
    if (code) {
      // some kind of error has happened
      message.subject = 'Error when computing your EpiAlignment data'
      message.text = header +
        'Unfortunately, an error occured when processing ' +
        'your data for EpiAlignment (Error Code: ' + code +
        '). We were unable to complete your request. \n\n' +
        'In some cases this may be caused by an erroneous input format, ' +
        'in which case you may try again by providing the correctly ' +
        'formatted input at https://' + hostName + '/. \n\n' +
        'Details of your run: \n\n' +
        runInfo.toString('email') + '\n' +
        'If the error keeps happening, please let us know by replying to ' +
        'this email, or sending an email to Jia Lu<jil430@eng.ucsd.edu> ' +
        'or Xiaoyi Cao <x9cao@eng.ucsd.edu>. \n\n ' +
        'Thank you for trying EpiAlignment!' + footer
    } else {
      // Normal email
      message.subject = 'Your EpiAlignment results are ready'
      message.text = header +
        'Your EpiAlignment results are now ready and can be downloaded ' +
        'by following this link:\n\n' +
        'https://' + hostName + '/result_page/' + runid + '\n\n' +
        'If the link above does not work, please copy the entire link ' +
        'and paste it into the address bar of your web browser.\n\n' +
        'Details of your run: \n\n' +
        runInfo.toString('email') + '\n' +
        'If you have any questions, please let us know by replying to ' +
        'this email, or sending an email to Jia Lu<jil430@eng.ucsd.edu> ' +
        'or Xiaoyi Cao <x9cao@eng.ucsd.edu>. \n\n' +
        'Thank you for using EpiAlignment!' + footer
    }
    transporter.sendMail(message)
  }
}

function postJobRun (
  code, runInfo, errorMsg, runInfoWritePromise, runInfoPreprocessPromise
) {
  let runid = runInfo.getDisplayValue('runid')
  runInfo.status = code

  let runInfoPath = getRunIdInfoPath(runid)

  console.log('[' + runid + '] Process quit with code : ' + code)
  if (errorMsg) {
    console.log(errorMsg)
    runInfo.addProperty('errMessage', errorMsg)
  }
  runInfoWritePromise = (runInfoWritePromise || Promise.resolve()).then(() =>
    writeFilePromise(runInfoPath, JSON.stringify(runInfo, null, 2))
  )
  let updateResultPromise = null
  if (runInfoPreprocessPromise) {
    updateResultPromise = runInfoPreprocessPromise.then(
      preprocessRes => updateExpression(runid, preprocessRes)
    )
  }
  Promise.all([runInfoWritePromise, updateResultPromise]).then(() => {
    runIdDict[runid] = code
    sendEmail(runInfo)
  })
}

function runPythonScript (
  pyMessenger, runInfoObject, runInfoPreprocessPromise,
  runInfoWritePromise, runInfoPath
) {
  // execute python code on the server.
  let scriptExecution = spawn('python', [pythonScript])

  let stdErrData = ''
  // Handle normal output
  scriptExecution.stdout.on('data', (data) => {
    console.log('[stdout] ' + data + '')
  })

  // Handle error output
  scriptExecution.stderr.on('data', data => {
    data.toString().split('\n').forEach(errLine => {
      if (errLine && errLine.trim().startsWith('[EpiAlignment]')) {
        errLine = errLine.replace('[EpiAlignment]', '').trim()
        stdErrData += (stdErrData ? '\n' : '') + errLine
        runInfoObject.addProperty('errMessage', stdErrData, true)
        runInfoWritePromise = runInfoWritePromise.then(() => writeFilePromise(
          runInfoPath, JSON.stringify(runInfoObject, null, 2)
        ))
      } else {
        console.log('[stderr] ' + data + '')
      }
    })
  })

  scriptExecution.on('exit',
    code => postJobRun(
      code, runInfoObject, stdErrData,
      runInfoWritePromise, runInfoPreprocessPromise
    )
  )
  // console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger))
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end()
}

app.post('/form_upload', cpUpload, function (req, res) {
  // req.files is an object (String -> Array) where fieldname is the key,
  //    and the value is array of files
  //
  // e.g.
  //  req.files['avatar'][0] -> File
  //  req.files['gallery'] -> Array
  //
  // req.body will contain the text fields, if there were any

  // Make new folder for results

  let runid = makeid()
  let tmpPath = getRunIdTempPath(runid)
  while (fs.existsSync(tmpPath)) {
    runid = makeid()
    tmpPath = getRunIdTempPath(runid)
  }
  fs.mkdirSync(tmpPath)

  // Move uploaded files
  let moveFilePromises = []
  for (let fileKey in req.files) {
    req.files[fileKey].forEach(file => {
      moveFilePromises.push(
        renamePromise(file.path, path.format({
          dir: getRunIdTempPath(runid),
          base: file.filename
        }))
      )
    })
  }

  // Add entry in `runIdDict`
  runIdDict[runid] = RUNNING_CODE

  // construct an object for python input
  let pyMessenger = {
    'body': req.body,
    'files': req.files,
    'runid': runid,
    'path': resultFolder
  }
  let runInfoPath = getRunIdInfoPath(runid)
  let runInfoObject = new RunInfo(req.body, req.files, runid)

  if (req.headers.host) {
    runInfoObject.addProperty('hostName', req.headers.host)
  }

  let runInfoWritePromise = writeFilePromise(
    runInfoPath, JSON.stringify(runInfoObject, null, 2))

  Promise.all(moveFilePromises).then(() => {
    let runInfoPreprocessPromise = preJobRun(runInfoObject, req.files)
    return runPythonScript(pyMessenger, runInfoObject,
      runInfoPreprocessPromise, runInfoWritePromise, runInfoPath)
  })

  res.json({ runid: runid })
})

app.get('/results/:runid', function (req, res) {
  let runid = req.params.runid
  let runInfoPath = getRunIdInfoPath(runid)
  if (!fs.existsSync(runInfoPath)) {
    // The query id does not exist.
    console.log('Invalid runid requested: ' + runid)
    res.status(401)
    res.json({
      errMessage: 'This Run ID does not exist! This may be due to a wrong' +
        ' link or an expired alignment run.\nPlease check your link and/or' +
        ' run your alignment again.'
    })
  } else {
    readFilePromise(runInfoPath, 'utf8')
      .then(infoString => {
        let infoObj = RunInfo.fromJSON(infoString)
        if (infoObj.status === undefined) {
          if (runIdDict.hasOwnProperty(runid)) {
            infoObj.status = RUNNING_CODE
            res.json(infoObj)
          } else {
            postJobRun(PROCESS_TERMINATED_SERVER_REBOOT, infoObj,
              'Server rebooted during job run. ' +
              'Please submit your job again.')
            res.status(500)
            res.json(infoObj)
          }
        } else if (infoObj.status === 0) {
          // python exited successfully.
          // read JSON file and return it.
          let resultFilePath = getRunIdResultPath(runid)
          readFilePromise(resultFilePath, 'utf8')
            .then(result => {
              let resDataObj = JSON.parse(result)
              res.json(Object.assign(resDataObj, infoObj.toJSON()))
            })
            .catch(err => {
              res.status(401)
              console.log(err)
              infoObj.status = E_CANNOT_READ_RESULTS
              infoObj.addProperty('errMessage',
                'Error loading result file. ' +
                'Please re-submit your job again.'
              )
              res.json(infoObj)
            })
        } else {
          // python exited with error.
          res.json(Object.assign({
            errMessage: 'EpiAlignment exited with an error.'
          }, infoObj.toJSON()))
        }
      }).catch(err => {
        res.status(401)
        console.log(err)
        res.json({
          status: E_CANNOT_READ_RESULTS,
          errMessage: 'Cannot read run info file.'
        })
      })
  }
})

app.get('/result_image/:runid/:index.png', function (req, res) {
  let runid = req.params.runid
  let index = req.params.index
  let imageName = getRunIdImagePath(runid, index)
  // Check if the image exists.
  if (fs.existsSync(imageName)) {
    res.sendFile(path.format({
      dir: __dirname,
      base: imageName
    }))
  } else {
    // if the image does not exist
    let scriptExecution = spawn('python', [pythonPlotScript])

    // Handle normal output
    scriptExecution.stdout.on('data', (data) => {
      console.log(data + '')
    })

    // Handle error output
    scriptExecution.stderr.on('data', (data) => {
      console.log(data + '')
    })

    scriptExecution.on('exit', (code) => {
      try {
        fs.statSync(imageName)
        res.sendFile(path.format({
          dir: __dirname,
          base: imageName
        }))
      } catch (err) {
        res.status(404).send('No image file available.')
      }
    })
    // python input
    let pyImageMessenger = {
      'index': index,
      'runid': runid,
      'path': resultFolder
    }
    scriptExecution.stdin.write(JSON.stringify(pyImageMessenger))
    // tell the node that sending inputs to python is done.
    scriptExecution.stdin.end()
  }
})

app.get('/get_cluster/:partialName', (req, res) => {
  clusterProc.getClusters(req.params.partialName).then(result => {
    res.json(result)
  }).catch(err => {
    console.log(err)
    res.status(400).send('Cannot get clusters.')
  })
})

/**
 * ***** Clean all pending (non-completed data) *****
 *
 */

readDirPromise(resultFolder)
  .then(files => files.forEach(fileName => {
    if (fileName.startsWith('tmp_')) {
      let runid = fileName.replace('tmp_', '')
      let runInfoPath = getRunIdInfoPath(runid)
      readFilePromise(runInfoPath, 'utf8')
        .then(infoString => {
          let infoObj = RunInfo.fromJSON(infoString)
          if (infoObj.status === undefined) {
            postJobRun(PROCESS_TERMINATED_SERVER_REBOOT, infoObj,
              'Server rebooted during job run. ' +
              'Please re-submit your job again.')
          }
        }).catch(err => {
          console.log(err)
        })
    }
  }))

const server = app.listen(3000, function () {
  var host = server.address().address || 'localhost'
  var port = server.address().port

  console.log('Listening http://%s:%s', host, port)
})
