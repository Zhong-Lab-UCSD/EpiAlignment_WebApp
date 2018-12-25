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

var runIdDict = {}
const RUNNING_CODE = -1

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
  while (fs.existsSync('tmp_' + runid)) {
    runid = makeid()
  }
  fs.mkdirSync('tmp_' + runid)
  runIdDict[runid] = RUNNING_CODE

  // construct an object for python input
  let pyMessenger = {
    'body': req.body,
    'files': req.files,
    'runid': 'tmp_' + runid
  }

  // execute python code on the server.
  let scriptExecution = spawn('python', [pythonScript])

  // Handle normal output
  scriptExecution.stdout.on('data', (data) => {
    console.log(data + '')
  })

  // Handle error output
  scriptExecution.stderr.on('data', (data) => {
    console.log(data + '')
  })

  scriptExecution.on('exit', (code) => {
    let message = null
    if (req.body.emailNote && req.body.mail) {
      message = {
        from: 'EpiAlignment Notification <messenger@mail.beta.givengine.org>',
        to: req.body.mail,
        replyTo: 'x9cao@eng.ucsd.edu'
      }
    }
    let header = 'Dear user,\n\n'
    let footer = '\n\nThe EpiAlignment Team from Zhong Lab @ UC San Diego\n' +
      '--\nDepartment of Bioengineering, Mail Code 0412\n' +
      'University of California, San Diego\n' +
      '9500 Gilman Dr.\nLa Jolla, CA 92122-0412\nUnited States'

    console.log('Process quit with code : ' + code)
    runIdDict[runid] = code
    if (message) {
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
          'formatted input at https://beta.epialign.ucsd.edu/. \n\n' +
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
          'If you have any questions, please let us know by replying to ' +
          'this email, or sending an email to Jia Lu<jil430@eng.ucsd.edu> ' +
          'or Xiaoyi Cao <x9cao@eng.ucsd.edu>. \n\n' +
          'Thank you for using EpiAlignment!' + footer
      }
      transporter.sendMail(message)
    }
  })
  // console.log(JSON.stringify(pyMessenger))
  // python input
  scriptExecution.stdin.write(JSON.stringify(pyMessenger))
  // tell the node that sending inputs to python is done.
  scriptExecution.stdin.end()

  res.json({ runid: runid })
})

app.get('/results/:runid', function (req, res) {
  let runid = req.params.runid
  if (!runIdDict.hasOwnProperty(runid)) {
    // The query id does not exist.
    console.log('Invalid runid requested: ' + runid)
    res.status(401)
    res.send('This query ID does not exist!')
  } else if (runIdDict[runid] === RUNNING_CODE) {
    res.json({ status: RUNNING_CODE })
  } else if (runIdDict[runid] === 0) {
    // python exited successfully.
    // read JSON file and return it.
    readFilePromise('tmp_' + runid + '/' + runid + '.json', 'utf8')
      .then(result => {
        let resObj = JSON.parse(result)
        resObj.status = 0
        res.json(resObj)
      })
      .catch(err => {
        res.status(401)
        res.send(err.message)
      })
  } else {
    // python exited with error.
    res.json({
      // TODO: find a way to insert error code and msg here
      status: 1,
      message: 'EpiAlignment exited with an error.'
    })
  }
})

app.get('/result_image/:runid/:index.png', function (req, res) {
  let runid = req.params.runid
  let index = req.params.index
  let imageName = 'tmp_' + runid + '/Image_' + index + '_' + runid + '.png'
  try {
    // Check if the image exists.
    fs.statSync(imageName)
    res.sendFile(path.format({
      dir: __dirname,
      base: imageName
    }))
  } catch (err) {
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
    let pyImageMessenger = { 'index': index, 'runid': runid }
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

const server = app.listen(3000, function () {
  var host = server.address().address || 'localhost'
  var port = server.address().port

  console.log('Listening http://%s:%s', host, port)
})
