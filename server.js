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

const E_CANNOT_READ_RESULTS = 401

const ClusterProcesser = require('./clusterProcessing')

var clusterProc = new ClusterProcesser(
  [
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
  ],
  {
    rawFilePath: 'Annotation/AnnotationFiles',
    clusterSuffix: '_clusters'
  }
)

const RunInfo = require('./runInfoParser')
const resultFolder = 'userResults'

var runIdDict = {}

const RUNNING_CODE = -1
const PROCESS_TERMINATED_SERVER_REBOOT = 500

// Server preparation step: parse the cluster file to build a
//  gene-name/ensemblID to cluster map

function makeid () {
  let text = ''
  let possible = 'abcdefghijklmnopqrstuvwxyz0123456789'

  for (var i = 0; i < 10; i++) {
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

function postJobRun (code, runInfo, errorMsg, writePromise) {
  let message = null
  let runid = runInfo.getDisplayValue('runid')
  let email = runInfo.getDisplayValue('email')
  runInfo.status = code

  let runInfoPath = getRunIdInfoPath(runid)
  if (email) {
    message = {
      from: 'EpiAlignment Notification <messenger@mail.beta.givengine.org>',
      to: email,
      replyTo: 'x9cao@eng.ucsd.edu'
    }
  }
  let header = 'Dear user,\n\n'
  let footer = '\n\nThe EpiAlignment Team from Zhong Lab @ UC San Diego\n' +
    '--\nDepartment of Bioengineering, Mail Code 0412\n' +
    'University of California, San Diego\n' +
    '9500 Gilman Dr.\nLa Jolla, CA 92122-0412\nUnited States'

  console.log('[' + runid + '] Process quit with code : ' + code)
  if (code) {
    runInfo.addProperty('errMessage', errorMsg)
  }
  (writePromise || Promise.resolve()).then(() =>
    writeFilePromise(runInfoPath, JSON.stringify(runInfo, null, 2))
  )
  if (message) {
    // needs to write an email
    let transporter = nodemailer.createTransport({
      sendmail: true,
      newline: 'unix',
      path: '/usr/sbin/sendmail'
    })
    runIdDict[runid] = code
    if (code) {
      // some kind of error has happened
      message.subject = 'Error when computing your EpiAlignment data'
      message.text = header +
        'Unfortunately, an error occured when processing ' +
        'your data for EpiAlignment (Error Code: ' + code +
        '). We were unable to complete your request. \n\n' +
        'In some cases this may be caused by an erroneous input format, ' +
        'in which case you may try again by providing the correctly ' +
        'formatted input at https://beta.epialign.ucsd.edu/. \n\n' +
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
        'https://beta.epialign.ucsd.edu/result_page/' + runid + '\n\n' +
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

app.post('/form_upload', cpUpload, function (req, res) {
  // req.files is an object (String -> Array) where fieldname is the key, and the value is array of files
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

  let writePromise = writeFilePromise(
    runInfoPath, JSON.stringify(runInfoObject, null, 2))

  // execute python code on the server.
  let scriptExecution = spawn('python', [pythonScript])

  let stdErrData = ''
  // Handle normal output
  scriptExecution.stdout.on('data', (data) => {
    console.log(data + '')
  })

  // Handle error output
  scriptExecution.stderr.on('data', data => {
    if (data && data.startsWith && data.startsWith('[EpiAlignment]')) {
      data = data.replace('[EpiAlignment]', '').trim()
      stdErrData += (stdErrData ? '\n' : '') + data
    } else {
      console.log(data + '')
    }
  })

  scriptExecution.on('exit',
    code => postJobRun(code, runInfoObject, stdErrData, writePromise))
  // console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger))
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end()

  res.json({ runid: runid })
})

app.get('/results/:runid', function (req, res) {
  let runid = req.params.runid
  let runInfoPath = getRunIdInfoPath(runid)
  if (!fs.existsSync(runInfoPath)) {
    // The query id does not exist.
    console.log('Invalid runid requested: ' + runid)
    res.status(401)
    res.send('This query ID does not exist!')
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
              'Please re-submit your job again.')
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
            // TODO: find a way to insert error code and msg here
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
        res.send(404, 'No image file available.')
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
    res.send(400, 'Cannot get clusters.')
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
      console.log(runid)
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
